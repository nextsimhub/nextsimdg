.. Copyright (c) 2021, Nansen Environmental and Remote Sensing Center

.. raw:: html



Installation
============

First step to install neXtSIM is to download this repository :

.. code::

    git clone https://github.com/nextsimdg/nextsimdg.git
    
or a specific version :

.. code::

    git clone -b v1.0 https://github.com/nextsimdg/nextsimdg.git


Dependencies
------------

To compile neXtSIM, you need to install first some libraries :

  - `NetCDF`_
  - `Boost`_
  - `CMake`_

**Installing dependencies on on MAC OS**

If your package manager is `Homebrew`_ :

.. code::

        brew install netcdf
        brew install boost
        brew install cmake
        
        
**Installing dependencies on Ubuntu**

You must have root privilege :

.. code::

        sudo apt-get update
        sudo apt-get install netcdf-bin libnetcdf-c++4-dev libboost-all-dev cmake
        

**Installing dependencies via conda**

Install conda via anaconda or miniconda (no root privileges required)

.. code::

        conda create --name nextsimdg
        conda activate nextsimdg
        conda install netCDF4
        conda -c conda-forge boost
        conda -c anaconda cmake

**Installing dependencies via spack**

``spack`` is a package manager which was designed for HPC software. It is language agnostic and
compiles binaries from source. To install ``spack`` please follow instruction from their
documentation (see `install instructions
<https://spack.readthedocs.io/en/latest/getting_started.html#installation>`_).

In the root directory of the repository there is a ``spack.yaml`` file which can be used to install
all required dependencies. The following command will create a ``spack`` environment (called
``nextsim`` in this example), load the ``spack`` environment and install the dependencies:

.. code::

   spack env create nextsim spack.yaml
   spack env activate nextsim
   spack install

.. note::

   This will install by building from source. This can take a significant amount of time, depending
   on packages available on your machine and how your spack install is configured. If you use a
   clean install of ``spack`` this could take 15-50 minutes on a modern CPU with 4-8 cores. It is
   hard to give an accurate measurement but the idea is to be patient. This only needs to be done
   once however, and then it can be loaded in seconds whenever you need it.

To use the installed dependencies, you will need to load them first. If you are in the same bash
session as the install they will be already loaded, but in future you can use the following command
to activate the environment.

.. code::

   spack env activate nextsim

To unload an environment you can use the handy (and very fun) alias of ``despacktivate``. For more
information on using ``spack`` environments please see `using environments
<https://spack.readthedocs.io/en/latest/environments.html#using-environments>`_.

Building the code
-----------------
After all dependencies have been installed, we can build the code:

.. code::

        cd nextsimdg
        mkdir -p build
        cd build
        cmake ..
        make

Dependencies and Build for MPI Parallelisation
----------------------------------------------

To build the code with MPI support, we need to install the respective compiler as well as parallel NetCDF support.

For example, on Debian-based Linux we need to also do:

.. code::

        sudo apt-get install libnetcdf-mpi-dev 
        sudo apt-get install openmpi-bin libopenmpi-dev 

The cmake call has to enable MPI support:

.. code::

        cmake .. -DENABLE_MPI=ON 

You might need to tell cmake which compiler to use, e.g.

.. code::

        cmake .. -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx -DENABLE_MPI=ON 

Using Dockerfiles for Dependencies and Building
-----------------------------------------------

In the Dockerfiles directory we provide two versions of Dockerfiles: with and without MPI support. You need to have _Docker installed on your computer. Then you can build and run it from the root of your NextSimDG repository, for example, with 

.. code::

        docker build . --file Dockerfiles/Dockerfile_Ubuntu -t nextsim
        docker run -it --entrypoint bash nextsim
        
The Ubuntu Dockerfiles create a full Ubuntu OS environment, copy over the content of the directory it is run from (which should be your nextsim root directory) and install and build nextsim in there for you. That way, you can easily try out your current code changes.

With the last command, you are inside your running container on the bash and can use your usual command line commands to look around and run the code.
    
    
.. _NetCDF: https://www.unidata.ucar.edu/software/netcdf/
.. _Boost: https://www.boost.org/
.. _CMake: https://cmake.org/
.. _Homebrew: https://brew.sh/
.. _Docker: https://www.docker.com/
