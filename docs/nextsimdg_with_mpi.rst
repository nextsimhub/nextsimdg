.. _halo-exchange:

neXtSIM_DG with MPI
===================

Introduction
------------

This section provides detailed instructions for neXtSIM_DG MPI cases. Specifically, it includes details on

* Installation instruction with MPI using spack
* Instructions on domain_decomp
* A step-by-step example on executing MPI cases

Installation instruction with MPI using spack
----------------------

[Spack](www.spack.io) is package manager mostly used on High Performance Computing (HPC) systems. We first find compiler on system using

.. code-block::
    spack compiler find
    spack compiler list

This should show the available ``gcc`` and ``llvm`` compilers on the system which will be used by spack. Next we use below spack specifications for creating spack environment.

.. code-block::

spack:
  packages:
    llvm:
      externals:
      - spec: llvm@21.1.8+clang~flang~lld~lldb
        prefix: /usr
        extra_attributes:
          compilers:
            c: /usr/bin/clang
            cxx: /usr/bin/clang++
    tar:
      externals:
      - spec: tar@1.30
        prefix: /usr
    c:
      require:
      - gcc@11.3.0
    cxx:
      require:
      - gcc@11.3.0
    fortran:
      require:
      - gcc@11.3.0
  specs:
  - boost@1.90+log+program_options
  - 'cmake@3.20:'
  - eigen@3.4.0
  - gmake
  - hdf5@1.14.3
  - openmpi@5.0.10
  - netcdf-c@4.9.2+mpi+parallel-netcdf
  - netcdf-cxx4@4.3.1
  - netcdf-fortran
  - scorep@9.4
  - zlib-ng@2.1.4
  - zoltan
  - py-netcdf4
  concretizer:
    unify: true
    reuse: false
  view: true

Once spack environment is created and activated, environment variables PATH and LD_LIBRARY_PATH should be adjusted to point towards domain_decomp executable and use spack installed packages respectively.

.. code-block::

    export PATH=$PATH:${DD_PATH}/build #Path to Domain Decomp (DD) executable
    export LD_LIBRARY_PATH=${SPACK_ROOT}/var/spack/environments/${SPACK_ENV_NAME}/.spack-env/view/lib:$LD_LIBRARY_PATH #Adding spack libraries to the LD_LIBRARY_PATH

During the build of Domain Decomp, it might be necessary to specify:

.. code-block::

    -DZoltan_INCLUDE_DIRS=${SPACK_ROOT}/var/spack/environments/${SPACK_ENV_NAME}/.spack-env/view/include

.. code-block::

   left_neighbour_ids = 0, 1
   left_neighbour_halos = 3, 2



Notes
-----

``Halo.hpp`` is still under active development. This documentation applies in a similar way for
corner neighbours (currently being implemented). Also, this diagram represents face-centered values
like those in ``ModelArray`` ``HField``'s for example. A similar approach is used for
``VertexField``'s and ``CGVector``'s (currently being implemented).
