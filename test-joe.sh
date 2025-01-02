#!/bin/env bash
# ====================================================================== #
# Run tests for nextSIM-DG.                                              #
# ====================================================================== #

set -e # Exit immediately if a command exits with a non-zero status

# Function to display help text
show_help() {
  echo "Usage: $0 [--all | -a] [--test | -t <TEST_NAME>] [--help | -h]"
  echo
  echo "Options:"
  echo "  --all  | -a    Run all MPI tests."
  echo "  --test | -t    Run a specific MPI test."
  echo "  --help | -h    Show this help message and exit."
  echo
  echo "If no options are provided, all XIOS tests will be run."
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
HELP=false
RUN_ALL=false
RUN_SINGLE_TEST=false
TEST_NAME=""
case "$#" in
1)
  case "$1" in
  --all | -a)
    RUN_ALL=true
    ;;
  --help | -h)
    HELP=true
    shift
    ;;
  *)
    echo "Invalid argument: $1"
    exit 1
    ;;
  esac
  ;;
2)
  case "$1" in
  --test | -t)
    RUN_SINGLE_TEST=true
    TEST_NAME=$2
    ;;
  *)
    echo "Invalid argument: $1"
    exit 1
    ;;
  esac
  ;;
*)
  echo "Running XIOS tests"
  ;;
esac

# Check for --help option
if [ "${HELP}" = true ]; then
  show_help
  exit 0
fi

# Run additional tests if --all option is provided
if [ "${RUN_ALL}" = true ]; then
  # Run MPI-parallel tests
  for COMPONENT in core physics; do
    cd ${BUILD_DIR}/${COMPONENT}/test
    for FILE in $(find test* -maxdepth 0 -type f); do
      echo ${FILE}
      NP=$(echo ${FILE} | sed -r "s/.*MPI([0-9]+)/\1/")
      if [ ${FILE} == ${NP} ]; then
        NP=1
      fi
      mpirun --allow-run-as-root --oversubscribe -n ${NP} ./${FILE}
    done
    cd -
  done
elif [ "${RUN_SINGLE_TEST}" = true ]; then
  # Run a single MPI-parallel test
  cd ${BUILD_DIR}/core/test
  echo test${TEST_NAME}_MPI2
  mpiexec --allow-run-as-root --oversubscribe -np 2 ./test${TEST_NAME}_MPI2
else
  # Only run XIOS tests
  cd ${BUILD_DIR}/core/test
  for FILE in $(find testXios* -maxdepth 0 -type f); do
    echo ${FILE}
    mpiexec --allow-run-as-root --oversubscribe -np 2 ./${FILE}
  done
fi
