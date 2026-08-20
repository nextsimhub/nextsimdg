#!/bin/bash

export DD_PATH=/home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/domain_decomp
export VENV_PATH=/home/${USER}/nextsimdg
export SPACK_ROOT=/home/${USER}/spack
export SPACK_ENV_NAME=nextsimdg-scorep-gcc
export NEXTSIMDG_PATH=/home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg

source ${VENV_PATH}/.venv/bin/activate # CHANGE ME
source ${SPACK_ROOT}/share/spack/setup-env.sh # CHANGE ME
spack env activate -p ${SPACK_ENV_NAME}

spack load gcc@11
spack load python

export LD_LIBRARY_PATH=${SPACK_ROOT}/var/spack/environments/${SPACK_ENV_NAME}/.spack-env/view/lib:$LD_LIBRARY_PATH

cd ${DD_PATH}
rm -rf build
mkdir -p build
cd ${DD_PATH}/build

cmake .. -DCMAKE_CXX_COMPILER=$(which mpic++) -DCMAKE_C_COMPILER=$(which mpicc) -DCMAKE_INSTALL_PREFIX=${DD_PATH}/build -DZoltan_INCLUDE_DIRS=${SPACK_ROOT}/var/spack/environments/${SPACK_ENV_NAME}/.spack-env/view/include

export PATH=$PATH:${DD_PATH}/build

make -j 4
make test
make install

cd ${NEXTSIMDG_PATH}

rm -rf build
mkdir build
cd build

if [[ "$1" == "--enable-mpi" ]]; then
	echo "MPI Build: -DENABLE_MPI=ON -DWITH_THREADS=ON"
	cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_C_COMPILER=$(which mpicc) -DCMAKE_CXX_COMPILER=$(which mpicxx) -DENABLE_MPI=ON -DWITH_THREADS=ON -DPython_EXECUTABLE=$(which python) -DCMAKE_BUILD_TYPE=Release
elif [[ "$1" == "--enable-openmp" ]]; then
	echo "OpenMP build: -DENABLE_MPI=OFF -DWITH_THREADS=ON"
        cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++) -DENABLE_MPI=OFF -DWITH_THREADS=ON -DPython_EXECUTABLE=$(which python) -DCMAKE_BUILD_TYPE=Release
elif [[ "$1" == "--enable-scorep-mpi" ]];then
	echo "Score-P MPI Build: -DCMAKE_C_COMPILER=$(which scorep-mpicc) -DCMAKE_CXX_COMPILER=$(which scorep-mpicxx) -DENABLE_MPI=ON -DWITH_THREADS=OFF"
	SCOREP_WRAPPER=off cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_C_COMPILER=$(which scorep-mpicc) -DCMAKE_CXX_COMPILER=$(which scorep-mpicxx) -DENABLE_MPI=ON -DWITH_THREADS=OFF -DPython_EXECUTABLE=$(which python) -DCMAKE_BUILD_TYPE=Release
else
	echo "Serial Build"
	cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++) -DENABLE_MPI=OFF -DWITH_THREADS=OFF -DPython_EXECUTABLE=$(which python) -DCMAKE_BUILD_TYPE=Release
fi

make -j 4
make test
