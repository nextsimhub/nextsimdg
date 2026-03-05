#!/usr/bin/env bash

# module purge
module load gcc/11
module load python/3.11.0-icl

source /home/${USER}/nextsimdg/.venv/bin/activate # CHANGE ME
source /home/${USER}/fenics_branch/share/spack/setup-env.sh # CHANGE ME
spack env activate -p nextsimdg_test_env # CHANGE ME
spack add scorep
spack install scorep
export PATH="$PATH:/home/nvs31/domain_decomp/build" # CHANGE ME
