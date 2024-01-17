#! /bin/bash

# Test the output of the NextSIM configuration help

# Check that nextsim exists
if [ ! -f nextsim ]; then
    echo "nextsim binary not found!"
    pwd
    exit 1
fi

# Get the help text. This will be produced in the case of unrecognized command
# line arguments.
HELP_TEXT=`./nextsim --help 2>&1`

# Check that the --help-config option exists. That is, doesn't produce the same
# output as HELP_TEXT.
HELPCFG_TEXT=`./nextsim --help-config 2>&1`

if [ "$HELP_TEXT" == "$HELPCFG_TEXT" ]; then
    echo "Test failed: --help-config option not recognized."
    exit 1
else
    echo "Test passed: --help-config option recognized."
fi

# Check the text produced by an implicit argument and by an argument of 'all'
# are the same

HELPALL_TEXT=`./nextsim --help-config all 2>&1`

if [ "$HELPALL_TEXT" == "$HELPCFG_TEXT" ]; then
    echo "Test passed: implicit value accepted."
else
    echo "Test failed: implicit value not accepted."
    exit 1
fi

# Check that CCSMIceAlbedo is listed in the 'avail'able configurable classes

HELPAVAIL_TEXT=`./nextsim --help-config avail 2>&1`

if [[ "$HELPAVAIL_TEXT" == *"CCSMIceAlbedo"* ]]; then
    echo "Test passed: CCSMIceAlbedo has configuration help."
else
    echo "Test failed: CCSMIceAlbedo is not listed as having configuration help."
    exit 1
fi

# Check the CCSMIceAlbedo help. This does depend on the CCSMIceAlbedo class not
# being changed too much 
HELPCCSM_TEXT=`./nextsim --help-config CCSMIceAlbedo 2>&1`

if [[ "$HELPCCSM_TEXT" == *"CCSMIceAlbedo.iceAlbedo"* ]] &&
    [[ "$HELPCCSM_TEXT" == *"CCSMIceAlbedo.snowAlbedo"* ]]; then
    echo "Test passed: CCSMIceAlbedo configuration help contains the expected fields."
else
    echo "Test failed: CCSMIceAlbedo configuration help does not contain the expected fields."
    exit 1
fi