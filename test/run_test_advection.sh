#!/bin/sh
# Prepares, executes and test the nextSIM-DG advection test

# Log commands
set -x

# Get the name of the Python executable
if [ $# -lt 1 ]; then
  PYTHON=python
else
  PYTHON=$1
fi

# Set filenames
restart_file=advection_test_init.nc
final_file=advection_test.restart.nc
diag_file=advection_test.diagnostic.nc

# Setup
${PYTHON} make_advection_test_data.py

# Run the advection test and check it corrects the expected files
./run_advection_test.sh
$PYTHON test_advection.py
test_return_value=$?

# Cleanup
rm $restart_file $out_file $diag_file
rm nextsim.*.log

exit $test_return_value
