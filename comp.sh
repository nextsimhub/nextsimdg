#!/bin/bash

export DD_PATH=/home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/domain_decomp
export VENV_PATH=/home/${USER}/nextsimdg
export SPACK_ROOT=/home/${USER}/spack
export SPACK_ENV_NAME=nextsimdg-scorep-gcc
export NEXTSIMDG_PATH=/home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg
export BENCHMARK_PATH=/home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg/run/run_june23
echo "$0 $1 $2"

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
	echo "Score-P MPI Build: -DCMAKE_C_COMPILER=$(which scorep-mpicc) -DCMAKE_CXX_COMPILER=$(which scorep-mpic++) -DENABLE_MPI=ON -DWITH_THREADS=OFF"
	SCOREP_WRAPPER=off cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_C_COMPILER=$(which scorep-mpicc) -DCMAKE_CXX_COMPILER=$(which scorep-mpic++) -DENABLE_MPI=ON -DWITH_THREADS=OFF -DPython_EXECUTABLE=$(which python) -DCMAKE_BUILD_TYPE=Release
else
	echo "Serial Build"
	cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++) -DENABLE_MPI=OFF -DWITH_THREADS=OFF -DPython_EXECUTABLE=$(which python) -DCMAKE_BUILD_TYPE=Release
fi

make -j 4
make test

cd ${BENCHMARK_PATH}

for MPI_SIZE in 64
do
	        echo "MPI Size: ${MPI_SIZE}"
		mpiexec -n ${MPI_SIZE} ${DD_PATH}/build/decomp -g ${BENCHMARK_PATH}/init_25km_NH.nc -x xdim -y ydim
		ln -sf ${BENCHMARK_PATH}/partition_metadata_${MPI_SIZE}.nc ${BENCHMARK_PATH}/partition.nc
		for NUM_THREADS in 1
		do
			export OMP_NUM_THREADS=${NUM_THREADS}
			echo "MPI Size: ${MPI_SIZE}; Number of threads ${NUM_THREADS}, ${OMP_NUM_THREADS}"
			if [[ "$1" == "--enable-scorep-mpi" ]]; then
			       mpiexec -n ${MPI_SIZE} ${NEXTSIMDG_PATH}/build/nextsim --config-file config_june23.cfg
			elif [[ "$1" == "--enable-mpi" ]]; then
			       time mpiexec -n ${MPI_SIZE} ${NEXTSIMDG_PATH}/build/nextsim --config-file config_june23.cfg
		        elif [[ "$1" == "--enable-openmp" ]]; then
				time ${NEXTSIMDG_PATH}/build/nextsim --config-file config_june23.cfg
			else
				${NEXTSIMDG_PATH}/build/nextsim --config-file config_june23.cfg
			fi
                               			       
		done
done
