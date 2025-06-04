#! /bin/sh
if [ ! -f ../nextsim ]; then
    echo "Copy or link the nextsim executableinto the parent directory from the build directory"
fi
../nextsim --config-file advection_test.cfg "$@"
