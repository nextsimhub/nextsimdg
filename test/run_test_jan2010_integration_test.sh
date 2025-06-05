#!/bin/sh

# Get the name of the python executable
if [ $# -lt 1 ]
then
    PYTHON=python
else
    PYTHON=$1
fi


# Prepares, executes and test the nextSIM-DG integration test based on the
# January 2010 Arctic simulation

restart_file=init_25km_NH.nc
era5_file=25km_NH.ERA5_2010-01-01T000000_test_data.nc
topaz_file=25km_NH.TOPAZ4_2010-01-01T000000_test_data.nc
out_file=out.integration_test.nc

$PYTHON make_init25kmNH_test_data.py
$PYTHON era5_topaz4_test_data.py
echo run_integration_test.sh
./run_integration_test.sh
$PYTHON test_integration_test.py
test_return_value=$?
rm $restart_file $era5_file $topaz_file $out_file
rm nextsim.*.log

exit $test_return_value
