#!/bin/bash

rm -rf build
mkdir build
cd build

SCOREP_WRAPPER=off cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_C_COMPILER=$(which scorep-gcc) -DCMAKE_CXX_COMPILER=$(which scorep-g++) -DENABLE_MPI=OFF -DWITH_THREADS=OFF -DPython_EXECUTABLE=$(which python) -DCMAKE_BUILD_TYPE=Release

make -j 1 
# make test
# make install


