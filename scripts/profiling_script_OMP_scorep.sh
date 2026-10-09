#!/bin/bash
#SBATCH -J OMP76scorep-june-23-config
#SBATCH -A ICCS-SL2-CPU
#SBATCH --output=scaling_OMP76_scorep.out
#SBATCH --error=scaling_OMP76_scorep.err

#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --mem=40000mb

#SBATCH -p icelake
#SBATCH --mail-type=NONE


cd /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/nvs31/nextsimdg/
source /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/nvs31/nextsimdg/comp-env.sh
cd /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/nvs31/domain_decomp
bash /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/nvs31/domain_decomp/comp.sh
cd /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/nvs31/nextsimdg/
bash /home/nvs31/rds/rds-iccs-DKRMHAHoC3M/nvs31/nextsimdg/comp-openmp-scorep.sh
# module load armforge
# module load papi/6.0.0.1/gcc/6o7wjeeo

cd /home/nvs31/rds/rds-iccs-DKRMHAHoC3M/nvs31/nextsimdg/run/run_june23

for MPI_SIZE in 1
do
	echo "MPI Size: ${MPI_SIZE}"
	mpirun -n ${MPI_SIZE} /home/nvs31/rds/rds-iccs-DKRMHAHoC3M/${USER}/domain_decomp/build/decomp -g /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg/run/run_june23/init_25km_NH.nc -x xdim -y ydim
	# /home/nvs31/rds/rds-iccs-DKRMHAHoC3M/${USER}/domain_decomp/build/decomp -g /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg/run/run_june23/init_25km_NH.nc -x xdim -y ydim
	ln -sf /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg/run/run_june23/partition_metadata_${MPI_SIZE}.nc /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg/run/run_june23/partition.nc

	for NUM_THREADS in 76
	do
		export OMP_NUM_THREADS=${NUM_THREADS}
		echo "MPI Size: ${MPI_SIZE}; Number of threads ${NUM_THREADS}, ${OMP_NUM_THREADS}"
		# time /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg/build/nextsim --config-file config_june23.cfg
		time mpirun -n ${MPI_SIZE} /home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg/build/nextsim --config-file config_june23.cfg
		# perf-report mpirun -n ${MPI_SIZE} /home/${USER}/nextsimdg/build/nextsim --config-file config_june23.cfg
	done
done
