#!/bin/bash
#SBATCH -J MPI76profiling-june-23-config
#SBATCH -A ICCS-SL2-CPU
#SBATCH --output=scaling_MPI76_x2048_scorep_64_32_1.out
#SBATCH --error=scaling_MPI76_x2048_scorep_64_32_1.err

#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=36:00:00
#SBATCH --mem=50000mb
#SBATCH --ntasks=76

#SBATCH -p icelake
#SBATCH --mail-type=NONE

### Set path and filenames ###
export DD_PATH=/home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/domain_decomp # CHANGE ME for domain_decomp path
export VENV_PATH=/home/${USER}/nextsimdg # CHANGE ME for virtual environment path
export SPACK_ROOT=/home/${USER}/spack # CHANGE ME for spack directory path
export SPACK_ENV_NAME=nextsimdg-scorep-gcc # CHANGE ME for a suitable environment name
export NEXTSIMDG_PATH=/home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg # CHANGE ME for path to NextSIMDG directory
export BENCHMARK_PATH=/home/nvs31/rds/rds-iccs-DKRMHAHoC3M/nvs31/nextsimdg/run/run_benchmark # CHANGE ME for Benchmark Directory
export CONFIG_FILE_NAME=config_benchmark.cfg # CHANGE ME for name of the config .cfg file in Benchmark directory
export PARTITION_DATA_FILE_NAME=init_benchmark.nc # CHANGE ME for partition data .nc file

### Intialize Virtual Environment and Spack Environment ###
source ${VENV_PATH}/.venv/bin/activate
source ${SPACK_ROOT}/share/spack/setup-env.sh
spack env activate -p ${SPACK_ENV_NAME}

### Load gcc and python modules from spack and add them to Library Path ###
spack load gcc@11
spack load python

export LD_LIBRARY_PATH=${SPACK_ROOT}/var/spack/environments/${SPACK_ENV_NAME}/.spack-env/view/lib:$LD_LIBRARY_PATH

### Build and test domain_decomp ###
cd ${DD_PATH}
rm -rf build
mkdir -p build
cd ${DD_PATH}/build

cmake .. -DCMAKE_CXX_COMPILER=$(which mpic++) -DCMAKE_C_COMPILER=$(which mpicc) -DCMAKE_INSTALL_PREFIX=${DD_PATH}/build -DZoltan_INCLUDE_DIRS=${SPACK_ROOT}/var/spack/environments/${SPACK_ENV_NAME}/.spack-env/view/include

export PATH=$PATH:${DD_PATH}/build

make -j 4
make test
make install

### Build and test NextSIMDG ###
cd ${NEXTSIMDG_PATH}

rm -rf build
mkdir build
cd build

### Set compilers as per options ###
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

if [[ "$1" == "--enable-openmp" ]]; then
	for MPI_SIZE in 1
        do
		echo "MPI Size: ${MPI_SIZE}"
		### Run domain_decomp ###
		${DD_PATH}/build/decomp -g ${BENCHMARK_PATH}/${PARTITION_DATA_FILE_NAME} -x xdim -y ydim
		ln -sf ${BENCHMARK_PATH}/partition_metadata_${MPI_SIZE}.nc ${BENCHMARK_PATH}/partition.nc
		for NUM_THREADS in 76 64 32 16 8 4 2 1
			do
				export OMP_NUM_THREADS=${NUM_THREADS}
				echo "MPI Size: ${MPI_SIZE}; Number of threads ${NUM_THREADS}, ${OMP_NUM_THREADS}"
				### Run nextsimdg ###
				time ${NEXTSIMDG_PATH}/build/nextsim --config-file ${CONFIG_FILE_NAME} >> openmp_time.txt
                         done
        done

elif [[ "$1" == "--enable-mpi" || "$1" == "--enable-scorep-mpi" ]]; then
	for MPI_SIZE in 64 32 1
	do
			echo "MPI Size: ${MPI_SIZE}"
			### Run domain_decomp ###
		        mpiexec -n ${MPI_SIZE} ${DD_PATH}/build/decomp -g ${BENCHMARK_PATH}/${PARTITION_DATA_FILE_NAME} -x x_dim -y y_dim
		        ln -sf ${BENCHMARK_PATH}/partition_metadata_${MPI_SIZE}.nc ${BENCHMARK_PATH}/partition.nc
			for NUM_THREADS in 1
				do
					export OMP_NUM_THREADS=${NUM_THREADS}
					echo "MPI Size: ${MPI_SIZE}; Number of threads ${NUM_THREADS}, ${OMP_NUM_THREADS}"
					if [[ "$1" == "--enable-scorep-mpi" ]]; then
						### Run nextsimdg ###
						mpiexec -n ${MPI_SIZE} ${NEXTSIMDG_PATH}/build/nextsim --config-file ${CONFIG_FILE_NAME}
					elif [[ "$1" == "--enable-mpi" ]]; then
						### Run nextsimdg ###
						time mpiexec -n ${MPI_SIZE} ${NEXTSIMDG_PATH}/build/nextsim --config-file ${CONFIG_FILE_NAME}
				        fi		       
		                 done
	done
else
	for MPI_SIZE in 1
	do
                        echo "MPI Size: ${MPI_SIZE}"
			### Run domain_decomp ###
                        ${DD_PATH}/build/decomp -g ${BENCHMARK_PATH}/${PARTITION_DATA_FILE_NAME} -x xdim -y ydim
                        ln -sf ${BENCHMARK_PATH}/partition_metadata_${MPI_SIZE}.nc ${BENCHMARK_PATH}/partition.nc
                        for NUM_THREADS in 1
                                do
                                        export OMP_NUM_THREADS=${NUM_THREADS}
                                        echo "MPI Size: ${MPI_SIZE}; Number of threads ${NUM_THREADS}, ${OMP_NUM_THREADS}"
                                        ### Run nextsimdg ###
					${NEXTSIMDG_PATH}/build/nextsim --config-file ${CONFIG_FILE_NAME}
                                done
        done
fi
