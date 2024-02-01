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
