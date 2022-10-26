.. Copyright (c) 2021, Nansen Environmental and Remote Sensing Center

.. raw:: html



Installation
============

First step to install neXtSIM is to download `this repository`_ :

.. code::

    git clone https://github.com/nextsimdg/nextsimdg.git 
    or 
    git clone git@github.com:nextsimdg/nextsimdg.git
    
or a specific version :

.. code::

    git clone -b v0.1.0 https://github.com/nextsimdg/nextsimdg.git 
    or 
    git clone -b v0.1.0 git@github.com:nextsimdg/nextsimdg.git


Dependencies
------------

To compile neXtSIM, you need to install first some libraries :

  - `NetCDF`_
  - `Boost`_
  - `Catch2`_
  - `CMake`_
  - `Eigen`_

Compilation on MAC OS
---------------------

If your package manager is `Homebrew`_ :

.. code::

        brew install netcdf-cxx
        brew install boost
        brew install catch2
        brew install cmake
        brew install eigen
        
        cd nextsimdg
        mkdir -p build
        cd build
        cmake ..
        make
        
Compilation on Ubuntu
---------------------

You must have root priviledge :

.. code::

        sudo apt-get update
        sudo apt-get install netcdf-bin libnetcdf-c++4-dev libboost-all-dev cmake libeigen3-dev
        git clone -b v2.x https://github.com/catchorg/Catch2.git
        cd Catch2
        cmake -Bbuild -H. -DBUILD_TESTING=OFF
        sudo cmake --build build/ --target install
        cd ..

        cd nextsimdg
        mkdir -p build
        cd build
        cmake ..
        make
        

Compilation with dependencies installation via conda
----------------------------------------------------

Install conda via `anaconda`_ or `miniconda`_ (no root priviledges required)

.. code::

        conda create --name nextsimdg
        conda activate nextsimdg
        conda install libgcc
        conda install netCDF4
        conda install -c conda-forge boost
        conda install -c conda-forge eigen
        conda install -c anaconda cmake
        conda install -c conda-forge catch2
        conda install -c conda-forge netcdf-cxx4
        
        cd nextsimdg
        mkdir -p build
        cd build
        cmake ..
        make
        
Extra dependencies needed to run the simple example
---------------------------------------------------

To run the simple example described `here`_, you need to install python and netCDF4 library. If you already installed the dependencies for compilation with conda, nothing else is needed.

Otherwise, install them via `anaconda`_ or `miniconda`_ (no root priviledges required)

.. code::

        conda create --name nextsimdg
        conda activate nextsimdg
        conda install netCDF4

.. _`this repository`: https://github.com/nextsimdg/nextsimdg    
.. _NetCDF: https://www.unidata.ucar.edu/software/netcdf/
.. _Boost: https://www.boost.org/
.. _Catch2: https://github.com/catchorg/Catch2
.. _Eigen: https://eigen.tuxfamily.org/
.. _CMake: https://cmake.org/
.. _Homebrew: https://brew.sh/
.. _here: https://nextsim-dg.readthedocs.io/en/latest/getting_started.html
.. _anaconda: https://www.anaconda.com/products/individual
.. _miniconda: https://docs.conda.io/en/latest/miniconda.html
