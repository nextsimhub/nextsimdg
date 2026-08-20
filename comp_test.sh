#!/bin/bash


echo "$1 $2 $0"

if [[ "$1" == "--build-openmp" ]]; then
	echo "OpenMP build"
elif [[ $1 == "--build-mpi" ]]; then
	echo "MPI build"
else
	echo "Serial build"
fi


