#!/bin/bash
#SBATCH -J profiling-june-23-config
#SBATCH -A ICCS-SL2-CPU
#SBATCH --output=profiler.out
#SBATCH --error=profiler.err

#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=10:00:00
#SBATCH --mem=40000mb

#SBATCH -p icelake
#SBATCH --mail-type=NONE


source /home/${USER}/nextsimdg/comp-env.sh # Location of Environment setup file
bash /home/${USER}/nextsimdg/comp.sh # Build nextsimdg
cd /home/${USER}/domain_decomp
bash /home/${USER}/domain_decomp/comp.sh # Build domain_decomp
module load armforge

cd /home/${USER}/nextsimdg/run/run_june23

for MPI_SIZE in 1
do
	echo "MPI Size: ${MPI_SIZE}"
	mpirun -n ${MPI_SIZE} /home/${USER}/domain_decomp/build/decomp -g /home/${USER}/nextsimdg/run/run_june23/init_25km_NH.nc -x xdim -y ydim
	ln -sf /home/${USER}/nextsimdg/run/run_june23/partition_metadata_${MPI_SIZE}.nc /home/${USER}/nextsimdg/run/run_june23/partition.nc

	for OMP_NUM_THREADS in 1 2 4 8 16 32 64 76
	do
		echo "MPI Size: ${MPI_SIZE}; Number of threads ${OMP_NUM_THREADS}"
		perf-report mpirun -n ${MPI_SIZE} /home/${USER}/nextsimdg/build/nextsim --config-file config_june23.cfg
	done
done
