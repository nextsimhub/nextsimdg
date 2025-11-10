#!/usr/bin/env bash

for EXTENSION in cpp hpp ipp; do
  find . -name "*.${EXTENSION}" -exec doxygen_prep.sh {} \;
done
cd docs && doxygen

if [ -z "$1" ]; then
  JOBS=1
else
  JOBS=$1
fi
sphinx-build -M html . _build --jobs $JOBS
