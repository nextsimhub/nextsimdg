#!/bin/env bash
# ========================================================================= #
# Build script for nextSIM-DG.                                              #
# ========================================================================= #

set -e # Exit immediately if a command exits with a non-zero status

BUILD_DIR="build"

# Check if a Python virtual environment is active
if [ -z "${VIRTUAL_ENV}" ]; then
  echo "No virtual environment is active. Please activate a virtual environment \
before running this script."
  exit 1
fi

# Check if a Spack environment is active
if [ -z "${SPACK_ENV}" ]; then
  echo "No Spack environment is active. Please activate a Spack environment \
before running this script."
  exit 1
fi

# Parse command line arguments
FRESH_BUILD=false
for arg in "$@"; do
  case $arg in
  --fresh | -f)
    FRESH_BUILD=true
    shift
    ;;
  *) ;;
  esac
done

# Check if a fresh build is requested
if [ "${FRESH_BUILD}" = true ]; then
  echo "Creating a fresh build..."
  rm -rf "${BUILD_DIR}"
else
  echo "Rebuilding..."
fi

# Install Python dependencies
python3 -m pip install -r requirements.txt

# Create build directory and navigate into it
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Different path to XIOS if running in a Docker container
if [ -f /.dockerenv ]; then
  xios_DIR="/xios"
fi

# Check if cmake and make are available
command -v cmake >/dev/null 2>&1 || {
  echo >&2 "cmake is required but it's not installed. Aborting."
  exit 1
}
command -v make >/dev/null 2>&1 || {
  echo >&2 "make is required but it's not installed. Aborting."
  exit 1
}

# Provide MPI configuration
# NOTE: Modify as appropriate
MPICC=mpicc
MPICXX=mpicxx
MPIF90=mpif90

# Build the model with XIOS support in Debug mode
cmake \
  -DCMAKE_BUILD_TYPE=Debug \
  -DENABLE_XIOS=ON \
  -Dxios_DIR="${xios_DIR}" \
  -DENABLE_MPI=ON \
  -DCMAKE_C_COMPILER="${MPICC}" \
  -DCMAKE_CXX_COMPILER="${MPICXX}" \
  -DCMAKE_Fortran_COMPILER="${MPIF90}" \
  -DENABLE_OASIS=ON .. \
  -DBUILD_TESTS=ON
make -j8
