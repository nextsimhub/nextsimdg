#!/bin/sh
# Prepares and executes nextSIM-DG integration test based on the January 2010
# Arctic simulation, but with erroneous forcing data. Tests if the return code
# is non-zero.

# Log commands
set -x

# Get the name of the Python executable
if [ $# -lt 1 ]; then
  PYTHON=python
else
  PYTHON=$1
fi

# Setup
${PYTHON} make_init25kmNH_test_data.py
${PYTHON} era5_topaz4_test_data.py --errors

# Run the integration test and check it raises the expected errors
./run_integration_test.sh

# The test passes if the model returns a non-zero exit code
if [ $? -ne 0 ]; then
  echo "Success: Error caught and model exits with a non-zero exit code"
  exit 0
else
  exit 1
fi
