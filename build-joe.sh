#!/bin/env bash
# ========================================================================= #
# Build script for nextSIM-DG.                                              #
# ========================================================================= #

set -e # Exit immediately if a command exits with a non-zero status

# Function to display help text
show_help() {
  echo "Usage: $0 [--prod | -p] [--fresh | -f] [--help | -h]"
  echo
  echo "Options:"
  echo "  --prod  | -p    Compile in Release mode."
  echo "  --fresh | -f    Create a fresh build before compiling."
  echo "  --help  | -h    Show this help message and exit."
}

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

# Set the build directory
BUILD_DIR="build"

# Parse command line arguments
PROD=false
FRESH_BUILD=false
HELP=false
for arg in "$@"; do
  case $arg in
  --prod | -p)
    PROD=true
    shift
    ;;
  --fresh | -f)
    FRESH_BUILD=true
    shift
    ;;
  --help | -h)
    HELP=true
    shift
    ;;
  *) ;;
  esac
done

# Check for --help option
if [ "${HELP}" = true ]; then
  show_help
  exit 0
fi

# Install Python dependencies
python3 -m pip install -r requirements.txt

# Check if a fresh build is requested
if [ "${FRESH_BUILD}" = true ]; then
  echo "Creating a fresh build..."
  rm -rf "${BUILD_DIR}"
else
  echo "Rebuilding..."
fi

# Use a different build directory in release mode
if [ "${PROD}" = true ]; then
  BUILD_DIR="${BUILD_DIR}_prod"
fi

# Create build directory and navigate into it
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Different path to XIOS if running in a Docker container
if [ -f /.dockerenv ]; then
  xios_DIR="/xios"
fi

# Check if CMake and GNU Make are available
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

if [ "${PROD}" = true ]; then
  # Build the model with XIOS support in Release mode
  cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DENABLE_XIOS=ON \
    -Dxios_DIR="${xios_DIR}" \
    -DENABLE_MPI=ON \
    -DENABLE_OASIS=ON .. \
    -DBUILD_TESTS=ON \
    -DCMAKE_C_COMPILER="${MPICC}" \
    -DCMAKE_CXX_COMPILER="${MPICXX}" \
    -DCMAKE_Fortran_COMPILER="${MPIF90}"
  make -j8
else
  # Build the model with XIOS support in Debug mode
  cmake \
    -DCMAKE_BUILD_TYPE=Debug \
    -DENABLE_XIOS=ON \
    -Dxios_DIR="${xios_DIR}" \
    -DENABLE_MPI=ON \
    -DENABLE_OASIS=ON .. \
    -DBUILD_TESTS=ON \
    -DCMAKE_C_COMPILER="${MPICC}" \
    -DCMAKE_CXX_COMPILER="${MPICXX}" \
    -DCMAKE_Fortran_COMPILER="${MPIF90}"
  # -DCMAKE_CXX_FLAGS="-O0 -Wall -enable=all"
  make VERBOSE=1 -j8
fi
