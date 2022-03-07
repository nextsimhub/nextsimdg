.. Copyright (c) 2021, Nansen Environmental and Remote Sensing Center

.. raw:: html

Getting Started
===============

This short guide explains how to get started with nextsimdg once you have installed it with one of the methods described in the `Installation`_ section.

Simple Example
--------------

To run the simple example you need a configuration file and a restart file.

Configuration file
~~~~~~~~~~~~~~~~~~

Control of nextsimdg is done using configuration files. One or more of these can be specified on the command line using the ``--config-file`` (for a single file) or ``--config-files`` (for several files) options. 

These files specify the configuration of the model, including the initial restart file (``model.init_file``) and the start (``model.start``), stop (``model.stop``) and time step (``model.time_step``) values, formatted as simple integers. 

An example of a configuration file (``dev1.cfg``) can be found in the ``run`` directory.

Restart file
~~~~~~~~~~~~

As part of the 0.1.0 release, the model operates on a simple fixed 10x10 grid of data.  

The restart file can be generated using the Python script ``dev_res.py`` available in run. Be sure to have installed python and the netCDF4 library (see `Installation`_ section) and run it with ``python dev_res.py``. This generates an initial restart file ``dev1.res.nc`` of the correct format, which is a netCDF file of the correct structure. The desired data can be provided by editing the python script.

An example of a configuration file (``dev1.res.nc``) can be found in the ``run`` directory.

Running the model
~~~~~~~~~~~~~~~~~

All is happenning in the ``run`` directory : you have the configuration file and the restart file there, you need to make a link to the compiled model ``nextsim``: ``ln -sf ../build/nextsim``.

With the value of the ``model.init_file`` variable set to the name of the correct initialization file, add the name of the configuration file as a ``config-file`` argument to the command line and execute. 

An example of a launching script (``dev1.sh``) can be found in the ``run`` directory, be sure it is executable (``chmod +x dev1.sh``) before executing it.

The model will produce a restart file named ``restart.nc``. The results of applying the model physics to the initial data over the specified number of time steps will be found here.

.. _Installation: https://nextsim-dg.readthedocs.io/en/latest/installation.html
