# G4Snitch

## Code overview

1. `G4Snitch.h/cc` - implementation of the SimWatcher concept. Traces G4 particles down to certain kinetic energy and stores the kinematics information into a ROOT TTree. Unstable paticles (with the exception of neutrons) are tracked until final decay.  
    Tracking limits are constexpr-declared in the `.cc` file.  
    Output ROOT filename is hardcoded to `G4Snitch.root`.

2. `G4S_DataFormat.h` - definition of output structures so that dictionary can easily build in an external, ROOT-only setup. Uses REve 64-bit float 4-vectors for position (+time) and momentum (+energy) representation.

3. A combined HGCal + G4Snitch config is in `SimG4Core/HelpfulWatchers/g4s-test/runG4Snitch.py`.  
   See development notes below about dictionary generation.

## Output format notes

Tree `T`
- Branch `p` 'particles'  
  - `std::vector<G4S_Particle>` -- kine tree in vector representation  
  - Units: GeV / cm / ns - converted in functions `g4_pos/mom_to_cms()`  
  - One entry per event - can get large.  
  - Contains also particles that pass the filter but are later NOT tracked by G4.  
    -> Check the bool m_was_tracked field.  
  - Element 0 is fictitous mother particle containing the true event primaries.  
    -> Can contain "empty" / default elements (not referenced)
- Branch `i` 'info'  
  - `G4S_Info` -- meta data, some global info about event


## Possible extensions

1. Track energy deposits
2. Track physical / logical volumes of the geometry.  
   In principle can also do this from TGeo -- but think about alignment.
3. Include more "path-marks" along the track (say, every x% energy loss)

All of the above can be in a separate tree / vector / branch,


## Possible / optional improvements

1. Re-compactify vector to remove unused slots.  
   Could also store primaries separately an glue them to the end.
2. Compactify out also non-tracked particles.

Compactification requires, of course, index remapping.

3. Vector re-capacity calculation does take into account missed / skipped  primaries. It could also taper down the extra factor based on how close too 100% of primary processing we get.


## Development notes

Disabling LTO / biglib speeds up build time significantly (x10).
```
scram b disable-biglib
cmsenv
```

### Debug scram build
USER_CXXFLAGS="-g -O0" scram b -j 16

### Hacking in dictionary for ROOT output

Apparently `SimG4Core/HelpfulWatchers` is a plugin and plgins can not have ROOT dictionaries (fails in some scram check looking for `libXyzQrt.so`, instead of `pluginXyzQrt.so`). Scram dictionary generation files are left in `src/` as inactive (and now also outdated) `classes.zh` and  `classes_def.zxml`.

To bypass this we build dictionary by hand and copy lib and dict to scram lib directory.
In top-level `src/` do:
```
bash SimG4Core/HelpfulWatchers/g4s-test/make_dict.sh
```
before build (or whenever `G4SnitchDataFormat` classes change).

This does the following:
```
name=RootG4Snitch
sdir=SimG4Core/HelpfulWatchers
installpath=../lib/${SCRAM_ARCH}/

rootcling -f ${name}.cc ${sdir}/interface/G4SnitchDataFormat.h ${sdir}/g4s-test/G4S_LinkDef.h
c++ -I${ROOTSYS}/include -I. ${name}.cc -fPIC -shared -o lib${name}.so
mkdir -p ${installpath}
cp ${name}_rdict.pcm lib${name}.so ${installpath}
```
