.. Copyright (c) 2021, Nansen Environmental and Remote Sensing Center

.. raw:: html



Installation
============

First step to install neXtSIM is to download this repository :

.. code::

    git clone https://github.com/nextsimhub/nextsimdg.git
    
You will get the main branch of the code, if you need a specific version :

.. code::

    git clone -b v1.0 https://github.com/nextsimhub/nextsimdg.git

It may be easier to use either the docker file (see below), or the ``spack`` installation instructions.


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
        brew install netcdf-cxx
        brew install boost
        brew install cmake
        brew install eigen


**Installing dependencies on Ubuntu**

Compilation on a Debian-based Linux distribution (Debian, Ubuntu, etc)
----------------------------------------------------------------------

You must have root privilege :

.. code::

        sudo apt-get update
        sudo apt-get install netcdf-bin libnetcdf-c++4-dev libboost-all-dev cmake subversion libeigen3-dev
        svn checkout http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk xios
        cd xios
        ./make_xios --arch <your_architecture>


**Installing dependencies via conda**

Install conda via anaconda or miniconda (no root privileges required)

.. code::

        conda create --name nextsimdg
        conda activate nextsimdg
        conda install -c conda-forge netcdf-cxx4
        conda install -c conda-forge eigen
        conda install -c conda-forge boost
        conda install -c anaconda cmake

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

Configuring the dynamics
------------------------
The dynamics for nextSIM are chosen at the point of configuring CMake. This is in contrast to most of the model configuration, which is done at model run time. The dynamics are set through the configuration option ``DynamicsType``. The available options for the dynamics are

* ``DG1``: First order discontinuous Galerkin dynamics on a 2D rectangular grid. Advection calculations are performed with 3 DG components.

  * This is the default option if no other option is provided to CMake.

* ``DG2``: Second order discontinuous Galerkin dynamics on a 2D rectangular grid. Advection calculations are performed with 6 DG components.

The syntax for chosing the dynamics via CMake is the standard method of providing options to CMake. For example, to compile the model with second order discontinuous Galerkin dynamics (``DG2``), the CMake command line with the dynamics argument would be

.. code::

        cmake -DDynamicsType=DG2 ..

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

Using Dockerfiles for Development or Production Runs
----------------------------------------------------

In the ``Dockerfiles`` directory we provide two ``Dockerfile``'s are provided.

- ``Dockerfile.devenv`` - This is the ``Dockerfile`` used to build the development image (``ghcr.io/nextsimhub/nextsimdg-dev-env``), that is
  used in the GitHub CI.
- ``Dockerfile.production`` - This ``Dockerfile`` is based off of the development image and it
  additionally installs ``nextsim`` so that you can run on any machine with ``docker`` installed.

Development Dockerfile
^^^^^^^^^^^^^^^^^^^^^^

A development image is provided on the nextsimhub `GitHub container registry
<https://github.com/orgs/nextsimhub/packages>`_ because it is needed for the CI.

If in future, this needs to be replaced. Please see `instructions
<https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry>`_
on the GitHub website.

To build the docker image, please use these instructions:

.. code-block:: console

    docker build --file Dockerfile.devenv . -t ghcr.io/nextsimhub/nextsimdg-dev-env:latest

.. note::
   The formatting of the image name **is important**. The format is
   ghcr.io/NAME_OF_REPOSITORY/NAME_OF_IMAGE:TAG

If you want to test or use the image locally, use the following command:

.. code-block:: console

   docker pull ghcr.io/nextsimhub/nextsimdg-dev-env:latest

Production Dockerfile
^^^^^^^^^^^^^^^^^^^^^

The production image is not stored on the nextsimhub `GitHub container registry
<https://github.com/orgs/nextsimhub/packages>`_ because it is not needed for the CI. Users of the
code may be interested in building their own. The instructions are as follows:

.. code-block:: console

   docker build --file Dockerfiles/Dockerfile.production . -t nextsim-production:latest

This will build a local image of the nextsim code. The production ``Dockerfile`` supports additional
build arguments (``--build-arg`` that can be specified at build time). For example, to build with
``MPI`` enabled, using 4 processes to compile, use the following command,

.. code-block:: console

   docker build --file Dockerfiles/Dockerfile.production --build-arg mpi=ON --build-arg jobs=4 . -t nextsim-production:latest

For a full list of options, please see ``Dockerfile.production``. By default ``MPI`` and ``xios``
options are disabled and the number of build jobs is 1.
