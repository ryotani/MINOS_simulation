#!/bin/tcsh

make clean
make all
make -f GNUmakefileBeamIn
make -f GNUmakefileDatai
make -f GNUmakefileDataoRings
make -f GNUmakefileSetup
