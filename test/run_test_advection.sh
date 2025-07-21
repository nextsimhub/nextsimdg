#!/bin/sh

# Get the name of the python executable
if [ $# -lt 1 ]
then
    PYTHON=python
else
    PYTHON=$1
fi


# Prepares, executes and test the nextSIM-DG advection test

restart_file=advection_test_init.nc
final_file=advection_test.restart.nc
diag_file=advection_test.diagnostic.nc

$PYTHON make_advection_test_data.py
echo run_advection_test.sh
./run_advection_test.sh
$PYTHON test_advection.py
test_return_value=$?
rm $restart_file $out_file $diag_file
rm nextsim.*.log

exit $test_return_value
