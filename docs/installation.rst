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
  - `Catch2`_
  - `CMake`_

Compilation on MAC OS
---------------------

If your package manager is `Homebrew`_ :

.. code::

        brew install netcdf
        brew install boost
        brew install catch2
        brew install cmake
        
        cd nextsimdg
        mkdir -p build
        cd build
        cmake ../..
        make
        
Compilation on Ubuntu
---------------------

You must have root priviledge :

.. code::

        sudo apt-get update
        sudo apt-get install netcdf-bin libnetcdf-c++4-dev libboost-all-dev cmake
        git clone -b v2.x https://github.com/catchorg/Catch2.git
        cd Catch2
        cmake -Bbuild -H. -DBUILD_TESTING=OFF
        sudo cmake --build build/ --target install
        cd ..

        cd nextsimdg
        mkdir -p build
        cd build
        cmake ../..
        make
        

        
    
.. _NetCDF: https://www.unidata.ucar.edu/software/netcdf/
.. _Boost: https://www.boost.org/
.. _Catch2: https://github.com/catchorg/Catch2
.. _CMake: https://cmake.org/
.. _Homebrew: https://brew.sh/
