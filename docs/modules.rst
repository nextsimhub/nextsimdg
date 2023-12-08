Modules
=======

Rationale
---------

One of the key features of nextSIM-DG is the run-time modularity of the code. This is designed to allow a user to change the parameterizations the model uses without having to recompile. Instead, many parameterizations and model components can be changed by changing the configuration, either in the config file or by providing overrides on the command line. The nextSIM-DG module system is designed to do this with a minimal impact on the performance of the model.

Usage
-----

To use the module system to change the model configuration from its defaults, the user needs to provide additional configuration values to the model. These can be permanently edited into the initial config file or provided as override on the command line. The modules section of the config file starts with the section title ``[Modules]``. Below this are key-value pairs, one per line with the format ``key = value``. The key is
