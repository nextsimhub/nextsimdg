#!/bin/sh

# Get the name of the python executable
if [ $# -lt 1 ]
then
    PYTHON=python
else
    PYTHON=$1
fi


# Prepares and executes nextSIM-DG integration test based on the January 2010
# Arctic simulation, but with erroneous forcing data. Tests if the return code
# is non-zero.

$PYTHON make_init25kmNH_test_data.py
$PYTHON era5_topaz4_test_data.py --errors
echo run_integration_test.sh
./run_integration_test.sh

# The test passes if the model returns a non-zero exit code
if [ $? -ne 0 ]
then
  echo "Success: Error caught and model exits with a non-zero exit code"
  exit 0
else
  exit 1
fi
