#!/bin/env bash
# ========================================================================= #
# Activate the Python and Spack environments used by nextSIM-DG.            #
# ========================================================================= #

# Activate the Python virtual environment
# NOTE: Modify as appropriate for the location you want to install the virtual environment
if [ ! -e "${HOME}/.virtualenvs/nextsim" ]; then
  echo "Creating a new virtual environment in ${HOME}/.virtualenvs/nextsim."
  mkdir -p "${HOME}/.virtualenvs"
  python3 -m venv "${HOME}/.virtualenvs/nextsim"
fi
source "${HOME}/.virtualenvs/nextsim/bin/activate"

# Activate the Spack environment
if ! spack env activate -p nextsim; then
  echo "Failed to activate Spack environment." >&2
  exit 1
fi

# Source the environment variables for XIOS
SOFTWARE="${HOME}/software"
if ! source "${SOFTWARE}/nextSIM-DG/.github/xios_arch/ubuntu.env"; then
  echo "Failed to source XIOS environment variables." >&2
  exit 1
fi

# Set XIOS directory
# NOTE: Modify as appropriate for the location you installed XIOS
export xios_DIR="${SOFTWARE}/XIOS"

# # Put decomp in the path
# # NOTE: Modify as appropriate for the location you installed domain_decomp
# export PATH="${SOFTWARE}/tools/domain_decomp/build":${PATH}
