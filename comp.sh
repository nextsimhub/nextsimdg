#!/bin/bash

rm -rf build
mkdir build
cd build

cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
	    -DCMAKE_C_COMPILER=$(which mpicc) -DCMAKE_CXX_COMPILER=$(which mpicxx) \
	        -DENABLE_MPI=ON -DWITH_THREADS=ON -DPYTHON_EXECUTABLE=$(which python)

make -j 4
make test
make install

