# Scaling Profiling Scripts for neXtSIM_DG

As part of development of neXtSIM_DG, scaling and profiling tests are performed to identify computational bottlenecks and to improve scalability of the code. We used [Score-P](https://www.vi-hps.org/projects/score-p/overview/overview.html) for profiling neXtSIM-DG code.

This repository provides scripts used for scaling and profiling tests for nextsimDG on IceLake CSD3 node. This includes
 - Spack environment file ```spack.yaml``` which provides specification for creating spack environment. We used [Spack](https://github.com/spack/spack) version 1.2.2 in this work which installs Score-P version 9.4. [Cube](https://www.scalasca.org/software/cube-4.x/download.html) is not listed in Spack specifications but can be downloaded using ```spack install --add cube```. However, since it takes some time to install cube, we used the one provided by ```module load``` on CSD3.
 - ```comp.sh``` for actual benchmark testing which activates the spack environment, creates ```domain_decomp``` and ```nextsimdg``` executables and perfoms benchmark testing. In this file paths need to be speacified where ```CHANGE ME``` is mentioned and ```SBATCH``` directives need to be adjusted as per user requirements. It has options to execute MPI scaling test, OpenMP scaling test and Score-P MPI profiling tests. 
  - ```SBATCH comp.sh --enable-mpi``` is used for MPI scaling test.
  - ```SBATCH comp.sh --enable-openmp``` is used for OpenMP scaling test.
  - ```SBATCH comp.sh --enable-scorep-mpi``` is used for Score-P MPI profiling test.