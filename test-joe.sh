#!/bin/env bash
# ====================================================================== #
# Run tests for nextSIM-DG.                                              #
# ====================================================================== #

set -e # Exit immediately if a command exits with a non-zero status

# Function to display help text
show_help() {
  echo "Usage: $0 [--debug | -d] [--no-xios | -nx] [--no-mpi | -nm]"
  echo "          [--integration-only | -i] [--unit-only | -u] [--help | -h]"
  echo
  echo "Options:"
  echo "  --debug            | -d    Compile in Debug mode."
  echo "  --no-xios          | -nx   Compile without XIOS support."
  echo "  --no-mpi           | -nm   Compile without MPI support."
  echo "  --integration-only | -i    Only run the integration tests."
  echo "  --unit-only        | -u    Only run the unit tests."
  echo "  --help             | -h    Show this help message and exit."
  echo
  echo "Any other options will be passed to CTest."
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
RUN_INTEGRATION_TESTS=true
RUN_UNIT_TESTS=true
HELP=false
for arg in "$@"; do
  case $arg in
  --debug | -d)
    BUILD_DIR="${BUILD_DIR}_debug"
    shift
    ;;
  --no-xios | -nx)
    BUILD_DIR="${BUILD_DIR}_noxios"
    shift
    ;;
  --no-mpi | -nm)
    BUILD_DIR="${BUILD_DIR}_nompi"
    shift
    ;;
  --integration-only | -i)
    RUN_UNIT_TESTS=false
    shift
    ;;
  --unit-only | -u)
    RUN_INTEGRATION_TESTS=false
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

# Run unit tests
if [ "${RUN_UNIT_TESTS}" = true ]; then
  cd "${BUILD_DIR}"
  ctest -V $@
  cd -
fi

# Run integration tests
if [ "${RUN_INTEGRATION_TESTS}" = true ]; then
  cp -r test "${BUILD_DIR}"
  cd "${BUILD_DIR}/test"
  python ThermoIntegration_test.py
  ./run_test_jan2010_integration.sh python
  ./run_test_advection.sh python
  ./run_test_error_handling_test.sh python
fi
