#!/usr/bin/env bash
#SBATCH --account=ICCS-SL2-CPU
#SBATCH --partition=icelake
#SBATCH --ntasks=32
#SBATCH --time=12:00:00

set -euo pipefail

ENVFILE="/home/mak60/rds/hpc-work/sasip/activate_env.sh"
if [[ ! -r "$ENVFILE" ]]; then
    printf 'Error: environment setup file not found or unreadable: %s\n' "$ENVFILE" >&2
    exit 1
fi
source "$ENVFILE"

SCRIPT_DIR="/home/mak60/rds/hpc-work/sasip/tripolar-verification/nextsimdg/test_case_tripolar"
REPO_ROOT="/home/mak60/rds/hpc-work/sasip/tripolar-verification/nextsimdg"

NTASKS="${SLURM_NTASKS:-1}"

DECOMP_EXEC="$REPO_ROOT/domain_decomp/build/decomp"
NEXTSIM_EXEC="$REPO_ROOT/build/nextsim"
GRIDFILE="$SCRIPT_DIR/grid_tnx1v4_20240627.nc"
CONFIG_FILE="$SCRIPT_DIR/tripolar_test.cfg"
PARTITION_METADATA="partition_metadata_${NTASKS}.nc"

if ! command -v mpirun >/dev/null 2>&1; then
    printf 'Error: mpirun was not found in PATH\n' >&2
    exit 1
fi

for executable in "$DECOMP_EXEC" "$NEXTSIM_EXEC"; do
    if [[ ! -x "$executable" ]]; then
        printf 'Error: executable not found or not executable: %s\n' "$executable" >&2
        exit 1
    fi
done

for input_file in "$GRIDFILE" "$CONFIG_FILE"; do
    if [[ ! -f "$input_file" ]]; then
        printf 'Error: input file not found: %s\n' "$input_file" >&2
        exit 1
    fi
done

if [[ -e "$SCRIPT_DIR/$PARTITION_METADATA" ]]; then
    printf 'Error: refusing to use existing partition metadata: %s\n' \
        "$SCRIPT_DIR/$PARTITION_METADATA" >&2
    exit 1
fi

cd "$SCRIPT_DIR"

# generate partition data
mpirun -np "$NTASKS" "$DECOMP_EXEC" --grid "$GRIDFILE" -x x_dim -y y_dim --tripolar
mv -- "$PARTITION_METADATA" partition.nc

mpirun -np "$NTASKS" "$NEXTSIM_EXEC" "--config-file=$CONFIG_FILE"
