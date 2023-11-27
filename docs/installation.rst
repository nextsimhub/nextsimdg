.. Copyright (c) 2021, Nansen Environmental and Remote Sensing Center

.. raw:: html



Installation
============

First step to install neXtSIM is to download this repository :

.. code::

    git clone https://github.com/nextsimhub/nextsimdg.git
    
or a specific version :

.. code::

    git clone -b v1.0 https://github.com/nextsimhub/nextsimdg.git


Dependencies
------------

To compile neXtSIM, you need to install first some libraries :

  - `NetCDF`_
  - `Boost`_
  - `CMake`_

Compilation on MAC OS
---------------------

If your package manager is `Homebrew`_ :

.. code::

        brew install netcdf
        brew install netcdf-cxx
        brew install boost
        brew install cmake
        brew install eigen
        
        cd nextsimdg
        mkdir -p build
        cd build
        cmake ..
        make
        
Compilation on a Debian-based Linux distribution (Debian, Ubuntu, etc)
----------------------------------------------------------------------

You must have root privilege :

.. code::

        sudo apt-get update
        sudo apt-get install netcdf-bin libnetcdf-c++4-dev libboost-all-dev cmake libeigen3-dev

        cd nextsimdg
        mkdir -p build
        cd build
        cmake ..
        make
        

Compilation with dependencies installation via conda
----------------------------------------------------

Install conda via anaconda or miniconda (no root privileges required)

.. code::

        conda create --name nextsimdg
        conda activate nextsimdg
        conda install -c conda-forge netcdf-cxx4
        conda install -c conda-forge boost
        conda install -c anaconda cmake
        
        cd nextsimdg
        mkdir -p build
        cd build
        cmake ..
        make
    
.. _NetCDF: https://www.unidata.ucar.edu/software/netcdf/
.. _Boost: https://www.boost.org/
.. _CMake: https://cmake.org/
.. _Homebrew: https://brew.sh/
