.. Copyright (C) 2023, Nansen Environmental and Remote Sensing Center

.. raw:: html

Getting Help
============

You can get help on compiling and running NextSIM-DG in several ways:

- The documentation that can be found alongside this document.
- The built-in help system.
- By emailing the NextSIM-DG team at nextsim@nersc.no

Documentation
-------------

Introduction
    Links to all the other available documentation, including the Doxygen API
    listing.

Installation
    How to install the prerequisite libraries for NextSIM-DG, where to obtain
    the code and how to compile the model.

Getting Started
    How to set up the model for its first run.

Getting Help
    This document. Where to find help with running the model beyond what is
    covered in these documentation files.

Built-in Help
-------------

``--help``, ``-h``
    The Unix standard  help arguments will print a message to the terminal listing the command line options.
    Currently this reads::

        neXtSIM_DG command line options::
          -h [ --help ]         print help message
          --help-config arg     arg = [avail | moduleName | all]
                                print help associated with one or more modules.
                                avail    list all available module names.
                                moduleName    print the help message associated with 
                                the given module.
                                all    print all help messages associated with all 
                                available modules.
          --config-file arg     specify a configuration file
          --config-files arg    specify a list of configuration files

``--help-config``
    This argument will display helpful information about the possible configuration of the model. There are three values that the argument can take:

- ``avail``
        This lists the names of all the available options.
- *moduleName*
        Given the name of a module, print all relevant configuration options,
        along with some text describing the option, the range of permitted
        values and the unconfigured default value.
- ``all``
        List the information described above for all `avail`able configuration
        options.
