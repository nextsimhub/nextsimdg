#!/bin/bash
#SBATCH -J OMP76profiling-june-23-625config
#SBATCH -A ICCS-SL2-CPU
#SBATCH --output=scaling_OMP76_625.out
#SBATCH --error=scaling_OMP76_625.err

#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --mem=40000mb

#SBATCH -p icelake
#SBATCH --mail-type=NONE


cd /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/nvs31/nextsimdg/
source /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/nvs31/nextsimdg/comp-env.sh
export LD_LIBRARY_PATH=/home/${USER}/spack/var/spack/environments/nextsimdg-scorep-gcc/.spack-env/view/lib:$LD_LIBRARY_PATH
cd /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/nvs31/domain_decomp
bash /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/nvs31/domain_decomp/comp.sh
export PATH=$PATH:/home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/domain_decomp/build
cd /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/nvs31/nextsimdg/
bash /home/nvs31/rds/rds-iccs-DKRMHAHoC3M/nvs31/nextsimdg/comp-openmp.sh
# module load armforge
# module load papi/6.0.0.1/gcc/6o7wjeeo

# cd /home/nvs31/rds/rds-iccs-DKRMHAHoC3M/nvs31/nextsimdg/run/run_june23
cd /home/nvs31/rds/rds-iccs-DKRMHAHoC3M/nvs31/nextsimdg/run/run_june23_625km

for MPI_SIZE in 1
do
	echo "MPI Size: ${MPI_SIZE}"
	# mpirun -n ${MPI_SIZE} /home/nvs31/rds/rds-iccs-DKRMHAHoC3M/${USER}/domain_decomp/build/decomp -g /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg/run/run_june23/init_25km_NH.nc -x xdim -y ydim
	/home/nvs31/rds/rds-iccs-DKRMHAHoC3M/${USER}/domain_decomp/build/decomp -g /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg/run/run_june23_625km/init_6.25km_NH_open.nc -x xdim -y ydim
	ln -sf /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg/run/run_june23_625km/partition_metadata_${MPI_SIZE}.nc /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg/run/run_june23_625km/partition.nc

	for NUM_THREADS in 76
	do
		export OMP_NUM_THREADS=${NUM_THREADS}
		echo "MPI Size: ${MPI_SIZE}; Number of threads ${NUM_THREADS}, ${OMP_NUM_THREADS}"
		time /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg/build/nextsim --config-file config_NH_benchmark_6_25km.cfg
		# time mpirun -n ${MPI_SIZE} /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg/build/nextsim --config-file config_june23.cfg
		# perf-report mpirun -n ${MPI_SIZE} /home/${USER}/nextsimdg/build/nextsim --config-file config_june23.cfg
	done
done
