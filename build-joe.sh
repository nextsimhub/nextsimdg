#!/bin/env bash
# ========================================================================= #
# Build script for nextSIM-DG.                                              #
# ========================================================================= #

set -e # Exit immediately if a command exits with a non-zero status

# Function to display help text
show_help() {
  echo "Usage: $0 [--debug | -d] [--no-xios | -nx] [--no-mpi | -nm] [--fresh | -f] [--help | -h]"
  echo
  echo "Options:"
  echo "  --debug   | -d    Compile in Debug mode."
  echo "  --no-xios | -nx   Compile without XIOS support."
  echo "  --no-mpi  | -nm   Compile without MPI support."
  echo "  --fresh   | -f    Create a fresh build before compiling."
  echo "  --help    | -h    Show this help message and exit."
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
BUILD_TYPE=Release
COMPILE_MPI=ON
COMPILE_XIOS=ON
FRESH_BUILD=false
HELP=false
for arg in "$@"; do
  case $arg in
  --debug | -d)
    BUILD_DIR="${BUILD_DIR}_debug"
    BUILD_TYPE=Debug
    shift
    ;;
  --no-xios | -nx)
    BUILD_DIR="${BUILD_DIR}_noxios"
    COMPILE_XIOS=OFF
    if [ "${COMPILE_MPI}" == "OFF" ]; then
      echo "Note: the --no-mpi flag already implies --no-xios"
    fi
    shift
    ;;
  --no-mpi | -nm)
    BUILD_DIR="${BUILD_DIR}_nompi"
    COMPILE_MPI=OFF
    if [ "${COMPILE_XIOS}" == "OFF" ]; then
      echo "Note: the --no-mpi flag already implies --no-xios"
    fi
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
mkdir -p "${BUILD_DIR}"

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

# Build the model
cmake -S. -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DENABLE_MPI="${COMPILE_MPI}" \
  -DENABLE_XIOS="${COMPILE_XIOS}" \
  -Dxios_DIR="${xios_DIR}" \
  -DBUILD_TESTS=ON \
  -DCMAKE_C_COMPILER="${MPICC}" \
  -DCMAKE_CXX_COMPILER="${MPICXX}"
#   -DCMAKE_CXX_FLAGS="-O0 -Wall -enable=all"
if [ "${BUILD_TYPE}" = Debug ]; then
  cd "${BUILD_DIR}"
  make VERBOSE=1 -j8
else
  cmake --build "${BUILD_DIR}"
fi
