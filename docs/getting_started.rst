.. Copyright (c) 2021, Nansen Environmental and Remote Sensing Center

.. raw:: html

Getting Started
===============

This short guide explains how to get started with nextsimdg once you have installed it with one of the methods described in the installation section.

Simple Example
--------------

Control of nextsimdg is done using configuration files. One or more of these can be specified on the command line using the `--config-file` (for a single file) or `--config-files` (for several files) options. These files specify the configuration of the model, including the initial restart file (`model.init_file`) and the start (`model.start`), stop (`model.stop`) and time step (`model.time_step`) values, formatted as simple integers. The configuration of parts of the model can also be changed, but this is beyond the scope of a simple example.

As part of the 0.1.0 release, the model operates on a simple fixed 10x10 grid of data.  The restart file can be generated using the Python script `dev_res.py`. This generates an initial restart file of the correct format, which is a netCDF file of the correct structure. The desired data can be provided by editing the python script.

With the value of the `model.init_file` variable set to the name of the correct initialization file, add the name of the configuration file as a `config-file` argument to the command line and execute. The model will produce a restart file named `restart.nc`. The results of applying the model physics to the initial data over the specified number of time steps will be found here.

An example config file (`dev1.cfg`) and shell script (`dev1.sh`) to run the model can be found in the `run` directory.

First Example
-------------

In this first example, we are going to ...
