#!/usr/bin/env bash

# module purge
module load gcc/11
module load python/3.11.0-icl
# module load scorep
# scorep --version

source /home/${USER}/nextsimdg/.venv/bin/activate # CHANGE ME
# source /home/${USER}/fenics_branch/share/spack/setup-env.sh # CHANGE ME
source /home/${USER}/spack/share/spack/setup-env.sh # CHANGE ME
spack env activate -p nextsimdg-scorep-gcc
# spack env activate -p nextsimdg_test_env # CHANGE ME
# spack install --add scorep %gcc@11
export PATH=$PATH:/home/nvs31/rds/rds-iccs-DKRMHAHoC3M/nvs31/domain_decomp/build # CHANGE ME

