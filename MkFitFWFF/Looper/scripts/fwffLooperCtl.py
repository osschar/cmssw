#!/usr/bin/env python3
"""Drive FWFFMkFitLooper in its interactive mode: send one command to its control directory and wait for the result.

    fwffLooperCtl.py DIR status                  print DIR/status.json
    fwffLooperCtl.py DIR wait                    wait until the looper waits for a command (e.g. after loop 0)
    fwffLooperCtl.py DIR rerun                   re-run the event unchanged
    fwffLooperCtl.py DIR stop                    end the job
    fwffLooperCtl.py DIR send FILE               send a command file as it is (Python or 'rerun' / 'stop')
    fwffLooperCtl.py DIR edit [--name N] [--module LABEL] [--set PARAM=PYEXPR ...] [--json PATCH]
                              [--knob GROUP.KNOB=PYEXPR ...]
                                                 send an edit made with fwffEdit(); PYEXPR is a cms-typed Python
                                                 value, e.g. --set "buildingRoutine=cms.string('cloneEngine')";
                                                 for --knob a plain value also works (1.5, 3, False), e.g.
                                                 --knob Window.dphi_trk_fac=1.5 (a -DMKFIT_V2P2_MUTABLE build)

Options: --timeout S (default 600) for the whole wait, --poll S (default 0.05) between looks at status.json,
--log FILE to append the result as one JSON line.

Protocol (python/customizeFWFFLooper.py has the full description): the client waits until status.json says the looper
is waiting, writes the command to a temporary file in DIR and renames it to DIR/command, then waits for the first
status.json with a higher 'seq' whose state is not 'running'. The latency printed is from the rename to the moment
that status was seen; 'accepted' is the time to the 'running' status, when the looper had read and applied the
command. Needs only the Python standard library; the looper reads the command with CMSSW's Python.

Exit code: 0 for the states waiting, stopped and ended; 1 for error (the command was refused, the job goes on);
2 for a timeout.
"""

import argparse
import json
import os
import sys
import time


def read_status(d):
    try:
        with open(os.path.join(d, 'status.json')) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def wait_for(d, pred, timeout, poll):
    t_end = time.time() + timeout
    while True:
        s = read_status(d)
        if s is not None and pred(s):
            return s
        if time.time() > t_end:
            return None
        time.sleep(poll)


def edit_text(name, module, sets, patch, knobs=()):
    args = [repr(name), repr(module or '')]
    if patch:
        args.append('jsonPatch = %r' % patch)
    if knobs:
        items = []
        for s in knobs:
            k, eq, v = s.partition('=')
            if not eq or '.' not in k:
                raise SystemExit('--knob needs GROUP.KNOB=PYEXPR, got %r' % s)
            items.append('%r: %s' % (k.strip(), v.strip()))
        args.append('knobs = {%s}' % ', '.join(items))
    for s in sets:
        k, eq, v = s.partition('=')
        if not eq or not k.strip():
            raise SystemExit('--set needs PARAM=PYEXPR, got %r' % s)
        args.append('%s = %s' % (k.strip(), v.strip()))
    return ('import FWCore.ParameterSet.Config as cms\n'
            'from MkFitFWFF.Looper.customizeFWFFLooper import fwffEdit\n'
            'edit = fwffEdit(%s)\n' % ', '.join(args))


def summary(s):
    ll = s.get('lastLoop') or {}
    t = ll.get('times_ms') or {}
    out = ['state=%s loop=%s seq=%s commands=%s' % (s.get('state'), s.get('loop'), s.get('seq'), s.get('commands'))]
    c = s.get('command')
    if c:
        out.append('command id=%s kind=%s name=%s' % (c.get('id'), c.get('kind'), c.get('name')))
    if s.get('error'):
        out.append('error: ' + s['error'].strip().replace('\n', '\n       '))
    if ll:
        out.append('lastLoop %s edit=%s events=%s times_ms %s' % (
            ll.get('loop'), ll.get('edit'), ll.get('events'),
            ' '.join('%s=%.1f' % (k, v) for k, v in sorted(t.items()))))
        for l in ll.get('digest', []):
            out.append('  ' + l)
    return '\n'.join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dir')
    ap.add_argument('what', choices=['status', 'wait', 'rerun', 'stop', 'send', 'edit'])
    ap.add_argument('file', nargs='?', help='send: the command file')
    ap.add_argument('--name', default='')
    ap.add_argument('--module', default='')
    ap.add_argument('--set', action='append', default=[], dest='sets')
    ap.add_argument('--json', default='', help='edit: a ConfigJsonPatcher patch of the iteration JSON')
    ap.add_argument('--knob', action='append', default=[], dest='knobs',
                    help='edit: a V2p2 knob, GROUP.KNOB=PYEXPR (a -DMKFIT_V2P2_MUTABLE build)')
    ap.add_argument('--timeout', type=float, default=600.)
    ap.add_argument('--poll', type=float, default=0.05)
    ap.add_argument('--log', default='', help='append the result as one JSON line to this file')
    a = ap.parse_args()
    d = a.dir

    if a.what == 'status':
        s = read_status(d)
        print(json.dumps(s, indent=1) if s else 'no status.json in %s' % d)
        return 0 if s else 2

    t_start = time.time()
    s = wait_for(d, lambda s: s.get('waiting') or s.get('state') in ('stopped', 'ended'), a.timeout, a.poll)
    if s is None:
        print('timeout: the looper does not wait for a command (%s)' % (read_status(d) or {}).get('state'))
        return 2
    if a.what == 'wait':
        print('waited %.0f ms\n%s' % (1000 * (time.time() - t_start), summary(s)))
        return 0
    if not s.get('waiting'):
        print('the looper has ended:\n' + summary(s))
        return 2

    if a.what in ('rerun', 'stop'):
        text = a.what + '\n'
    elif a.what == 'send':
        if not a.file:
            raise SystemExit('send needs a FILE')
        text = open(a.file).read()
    else:
        text = edit_text(a.name or 'edit', a.module, a.sets, a.json, a.knobs)

    seq0 = s['seq']
    tmp = os.path.join(d, '.command.%d.tmp' % os.getpid())
    with open(tmp, 'w') as f:
        f.write(text)
    os.replace(tmp, os.path.join(d, 'command'))
    t_sent = time.time()

    accepted = []
    def done(st):
        if st.get('seq', 0) <= seq0:
            return False
        if st.get('state') == 'running':
            if not accepted:
                accepted.append(time.time())
            return False
        return True
    r = wait_for(d, done, a.timeout, a.poll)
    t_res = time.time()
    if r is None:
        print('timeout after %.0f s; status: %s' % (t_res - t_sent, (read_status(d) or {}).get('state')))
        return 2
    lat = 1000 * (t_res - t_sent)
    acc = 1000 * (accepted[0] - t_sent) if accepted else None
    print('%s: latency_ms=%.0f accepted_ms=%s\n%s' % (a.what, lat, '%.0f' % acc if acc is not None else '-', summary(r)))
    if a.log:
        with open(a.log, 'a') as f:
            f.write(json.dumps({'what': a.what, 'text': text, 'sent_epoch': t_sent, 'latency_ms': lat,
                                'accepted_ms': acc, 'status': r}) + '\n')
    return 1 if r.get('state') == 'error' else 0


if __name__ == '__main__':
    sys.exit(main())
