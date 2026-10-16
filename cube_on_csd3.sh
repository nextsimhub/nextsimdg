#!/bin/bash
module purge
module unuse /usr/local/Cluster-Config/modulefiles
module use /usr/local/Cluster-Config/modulefiles-test
module load rhel8/cpu-x86_64_v4/legacy
module load cube/4.8.2/gcc/irx73yh5

cube
