#!/bin/bash

rm -rf build
mkdir build
cd build

# SCOREP_WRAPPER=off cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_C_COMPILER=scorep-gcc -DCMAKE_CXX_COMPILER=scorep-g++ -DENABLE_MPI=OFF -DWITH_THREADS=OFF -DPYTHON_EXECUTABLE=$(which python) -DCMAKE_BUILD_TYPE=Release

SCOREP_WRAPPER=off cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_C_COMPILER=$(which scorep-mpicc) -DCMAKE_CXX_COMPILER=$(which scorep-mpic++) -DENABLE_MPI=ON -DWITH_THREADS=OFF -DPython_EXECUTABLE=$(which python) -DCMAKE_BUILD_TYPE=Release

make -j 4
# make test
# make install
