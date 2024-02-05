.. Copyright (c) 2021, Nansen Environmental and Remote Sensing Center

.. raw:: html

Getting Started
===============

This short guide explains how to get started with nextsimdg once you have installed it with one of the methods described in the `Installation`_ section.

Simple Example
--------------

To run the simple example you need a configuration file, where the model settings are defined, and a file of the initial state from which to start the simulation.

Configuration file
~~~~~~~~~~~~~~~~~~

The control of nextsimdg is done using configuration files. One or more of these can be specified on the command line using the ``--config-file`` (for a single file) or ``--config-files`` (for several files) options. 

These files specify the configuration of the model, including the start ``model.start``, stop ``model.stop`` and time step ``model.time_step`` values, formatted as dates, and file of the initial state  ``model.init_file``. 

An example of a configuration file ``config_simple_example.cfg`` can be found in the ``run`` directory.

Initial file
~~~~~~~~~~~~

As part of the 0.2.0 release, the model can be operated on a simple example with a rectangular 30x30 grid.  

The initial state file can be generated using the Python script ``make_init.py`` available in the ``run`` directory. Be sure to have installed python and the netCDF4 library (see `Installation`_ section) and run it with ``python make_init.py``. This generates a file for the initial state ``init_rect30x30.nc`` of the correct format, which is a netCDF file of the correct structure. The desired data can be provided by editing the python script.

An example of a configuration file ``init_rect30x30.nc`` can be found in the ``run`` directory.

Running the model
~~~~~~~~~~~~~~~~~

All is happenning in the ``run`` directory : you have the configuration file and the initial state file there, you need to make a link to the compiled model ``nextsim``: ``ln -sf ../build/nextsim`` (we assume you already compiled it with one of the methods described in the `Installation`_ section).

With the value of the ``model.init_file`` variable set to the name of the correct initial state file, add the name of the configuration file as a ``config-file`` argument to the command line and execute. 

An example of a launching script ``run_simple_example.sh`` can be found in the ``run`` directory, be sure it is executable (``chmod +x run_simple_example.sh``) before executing it.

The model will produce an output file named ``output.nc``, which contains the result of applying the model physics to the initial state over the specified number of time steps.

Running with MPI
~~~~~~~~~~~~~~~~
To run the model with MPI you have to have it built with MPI or use one of the `*_MPI` Dockerfiles. In addition, you will need a partition NetCDF file as created by the `Domain Decomposition tool <https://github.com/nextsimhub/domain_decomp>`_. This partition file has to either be added to the config file (`*.cfg`) by adding the line `partition_file = partition.nc` or provided as an argument to `nextsim`.

.. code::

    mpirun -n 1 ./nextsim --config-file config_simple_example.cfg --model.partition_file partition.nc

An example partition file for the simple example is added to the ``run`` directory. For more information on the format of the partition files and how to create them, see `the documentation for the decomposition tool <https://github.com/nextsimhub/domain_decomp>`_ . Note that you will need a specific partition file depending on the number of processes you want to run, which will, by default, be named ``partition_metadata_<num_mpi_processes>.nc``.



.. _Installation: https://nextsim-dg.readthedocs.io/en/latest/installation.html

