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

(Spack)[www.spack.io] is package manager mostly used on High Performance Computing (HPC) systems. We first find compiler on system using

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

The spack installations then can be done with

.. code-block::
    spack env create SPACK_ENV_NAME
    spack env activate SPACK_ENV_NAME
    spack concretize -f
    spack install
    spack env deactivate

In the above step spack environment with necessary dependencies is now created.

After creating spack environment ``SPACK_ENV_NAME``, we set some environment variables:

.. code-block::
    ### Set path and filenames ###
    export DD_PATH=/home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/domain_decomp # CHANGE ME for domain_decomp path
    export VENV_PATH=/home/${USER}/nextsimdg # CHANGE ME for virtual environment path
    export SPACK_ROOT=/home/${USER}/spack # CHANGE ME for spack directory path
    export SPACK_ENV_NAME=nextsimdg-test2 # nextsimdg-scorep-gcc # CHANGE ME for a suitable environment name
    export NEXTSIMDG_PATH=/home/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg # CHANGE ME for path to NextSIMDG directory
    export BENCHMARK_PATH=/homei/${USER}/rds/rds-iccs-DKRMHAHoC3M/${USER}/nextsimdg/run/run_benchmark # CHANGE ME for Benchmark Directory
    export CONFIG_FILE_NAME=config_benchmark.cfg # CHANGE ME for name of the config .cfg file in Benchmark directory
    export PARTITION_DATA_FILE_NAME=init_benchmark.nc # CHANGE ME for partition data .nc file

We also need to create virtual python environment:

.. code-block::
    python -m venv VENV_NAME

Next, virtual environment and spack environment are intialised:

.. code-block::
    source ${VENV_PATH}/.venv/bin/activate
    source ${SPACK_ROOT}/share/spack/setup-env.sh
    spack env activate -p ${SPACK_ENV_NAME}

and update ``LD_LIBRARY_PATH`` as:

.. code-block::
    export LD_LIBRARY_PATH=${SPACK_ROOT}/var/spack/environments/${SPACK_ENV_NAME}/.spack-env/view/lib:$LD_LIBRARY_PATH #Adding spack libraries to the LD_LIBRARY_PATH

Next, domain_decomp is installed with

.. code-block::
    
    ### Build and test domain_decomp ###
    cd ${DD_PATH}
    rm -rf build
    mkdir -p build
    cd ${DD_PATH}/build
    
    cmake .. -DCMAKE_CXX_COMPILER=$(which mpic++) -DCMAKE_C_COMPILER=$(which mpicc) -DCMAKE_INSTALL_PREFIX=${DD_PATH}/build -DZoltan_INCLUDE_DIRS=${SPACK_ROOT}/var/spack/environments/${SPACK_ENV_NAME}/.spack-env/view/include
    
    export PATH=$PATH:${DD_PATH}/build #Path to Domain Decomp (DD) executable
    
    make -j $(nprocs)
    make test
    make install

Next, we create install or create nextsimdg executable.

.. code-block::
    
    ### Build and test NextSIMDG ###
    cd ${NEXTSIMDG_PATH}
    
    rm -rf build
    mkdir build
    cd build

We set the compiler options for nextsimdg. Here, ``-DWITH_THREADS=ON`` allows openmp thread parallelization and ``-DENABLE_MPI`` allows MPI parallelization.

.. code-block::
    
    cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_C_COMPILER=$(which mpicc) -DCMAKE_CXX_COMPILER=$(which mpicxx) -DENABLE_MPI=ON -DWITH_THREADS=ON -DPython_EXECUTABLE=$(which python) -DCMAKE_BUILD_TYPE=Release

    make -j $(nprocs)
    make install #This includes testing