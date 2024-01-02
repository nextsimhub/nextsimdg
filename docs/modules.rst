Modules
=======

Rationale
---------

One of the key features of nextSIM-DG is the run-time modularity of the code. This is designed to allow a user to change the parameterizations the model uses without having to recompile. Instead, many parameterizations and model components can be changed by changing the configuration, either in the config file or by providing overrides on the command line. The nextSIM-DG module system is designed to do this with a minimal impact on the performance of the model.

Usage
-----

To use the module system to change the model configuration from its defaults, the user needs to provide additional configuration values to the model. These can be permanently edited into the initial config file or provided as override on the command line. The modules section of the config file starts with the section title ``[Modules]``. Below this are key-value pairs, one per line with the format ``key = value``. The key is the name of the directory that contains the module. This will be a subdirectory of one of the ``modules`` subdirectories in the nextSIM-DG code tree. An example of this is the subdirectory ``core/src/modules/DiagnosticOutputModule`` which would be in the config file as the key ``DiagnosticOutputModule =``. The value following the equals sign is then one of the implementations of the module. The names of these can be obtained from the online help system accessible using the command line argument ``--config-help`` as well as the relevant file in this documentation directory.

