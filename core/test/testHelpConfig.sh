#! /bin/bash

# Test the output of the NextSIM configuration help

# Get the help text. This will be produced in the case of unrecognized command
# line arguments.
HELP_TEXT=`./nextsim --help 2>&1`

# Check that the --help-config option exists. That is, doesn't produce the same
# output as HELP_TEXT.
HELPCFG_TEXT=`./nextsim --help-config 2>&1`

if [ "$HELP_TEXT" == "$HELPCFG_TEXT" ]; then
    echo "Test 1 failed :("
    exit 1
else
    echo "Test 1 passed :)"
fi

# Check the text produced by an implicit argument and by an argument of 'all'
# are the same

HELPALL_TEXT=`./nextsim --help-config all 2>&1`

if [ "$HELPALL_TEXT" == "$HELPCFG_TEXT" ]; then
    echo "Test 2 passed :)"
else
    echo "Test 2 failed :("
    exit 1
fi

# Check that CCSMIceAlbedo is listed in the 'avail'able configurable classes

HELPAVAIL_TEXT=`./nextsim --help-config avail 2>&1`

if [[ "$HELPAVAIL_TEXT" == *"CCSMIceAlbedo"* ]]; then
    echo "Test 3 passed :)"
else
    echo "Test 3 failed :("
    exit 1
fi

# Check the CCSMIceAlbedo help. This does depend on the CCSMIceAlbedo class not
# being changed too much 
HELPCCSM_TEXT=`./nextsim --help-config CCSMIceAlbedo 2>&1`

if [[ "$HELPCCSM_TEXT" == *"CCSMIceAlbedo.iceAlbedo"* ]] &&
    [[ "$HELPCCSM_TEXT" == *"CCSMIceAlbedo.snowAlbedo"* ]]; then
    echo "Test 4 passed :)"
else
    echo "Test 4 failed :("
    exit 1
fi