#!/bin/bash
#SBATCH -J profiling-june-23-config
#SBATCH -A ICCS-SL2-CPU
#SBATCH --output=profiler.out
#SBATCH --error=profiler.err

#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=03:00:00
#SBATCH --mem=40000mb

#SBATCH -p icelake
#SBATCH --mail-type=ALL


source /home/${USER}/nextsimdg/comp-env.sh # Location of Environment setup file
bash /home/${USER}/nextsimdg/comp.sh # Build nextsimdg
cd /home/${USER}/domain_decomp
bash /home/${USER}/domain_decomp/comp.sh # Build domain_decomp

cd /home/${USER}/nextsimdg/run/run_june23
time mpirun -n 4 /home/${USER}/nextsimdg/build/nextsim --config-file config_june23.cfg

