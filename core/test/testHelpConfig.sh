#! /bin/bash

# Test the output of the NextSIM configuration help

# Get the help text. This will be produced in the case of unrecognized command
# line arguments.
HELP_TEXT=`./nextsim --help 2>&1`

# Check that the --help-config option exists. That is, doesn't produce the same
# output as HELP_TEXT.
HELPCFG_TEXT=`./nextsim --help-config 2>&1`

if [ "$HELP_TEXT" == "$HELPCFG_TEXT" ]; then
    echo "Test failed :("
    exit 1
else
    echo "Test passed :)"
fi