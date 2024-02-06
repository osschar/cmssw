name=RootG4Snitch
sdir=SimG4Core/HelpfulWatchers
installpath=../lib/${SCRAM_ARCH}/

rootcling -f ${name}.cc ${sdir}/interface/G4SnitchDataFormat.h ${sdir}/g4s-test/G4S_LinkDef.h
c++ -I${ROOTSYS}/include -I. ${name}.cc -fPIC -shared -o lib${name}.so
mkdir -p ${installpath}
cp ${name}_rdict.pcm lib${name}.so ${installpath}
