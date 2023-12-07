Model grids
===========

Model grids in nextSIM-DG are derived classes of the `IStructure` class. Currently, two are defined: ``RectangularGrid`` and ``ParametricGrid``. 

Simple rectangular
------------------

The first gird type supports simple rectangular which are defined only by the number of grid points in the x and y dimensions. The lack of any spatial scale makes it unsuitable for use with the dynamics. It is intended to be used for developing and testing the column physics of the model where the actual size of the grid cells is unimportant.

Parametric rectangular
----------------------

The second grid type holds parameters that define the position and size of each grid cell. This means it is useable with the DG dynamics and in fact can only be used with the DG dynamics as the input and output functions for this grid assume use of the `ModelArray` types that support DG and CG data arrays.

GridIO
======

The bulk of the important code is in the ``*GridIO`` classes, which define how the information transferred from the netCDF file to the model for the relevant grid. There are two classes, ``RectGridIO`` for ``RectangularGrid`` and ``ParaGridIO`` for ``ParametricGrid``. Although these two classes are not derived from a common base class to provide flexibility in the arguments taken by the member functions, they are otherwise broadly similar. In both cases the data in an instance of ``ModelState`` is written to or read from file. This class provides the interface between IO code and the rest of the model.

* ``getModelState()`` constructs an instance of `ModelState` given a file path.
* ``dumpModelState()`` writes the data in a given instance of `ModelState` to the given file path.

``RectGridIO``
--------------

* ``getModelState()`` sets the size of the ``H``, ``U``, ``V`` and ``Z`` ``ModelArrays`` using the ``dimensionSetter()`` function. With the arrays correctly sized, the function reads the data for the land mask, the ice thickness, the ice concentration, the snow thickness and the ice temperature. The netCDF variable names are those taken from ``gridnames.hpp` and mapped to the same names within the ``ModelState`` data map. Close the netCDF file.

* ``dumpModelState()`` writes the comman restart metadata using the ``CommonRestartMetadata`` class. The dimensions are written using the sizes of ``H``and ``Z``arrays. All of the data in passed ``ModelState`` is then written into the file based on whether it is an ``H`` array
or a ``Z`` array. Close the netCDF file.

``ParaGridIO``
--------------

The largest difference in ``ParaGridIO`` is the handling of files with a time axis. These will generally be forcing data files when reading and diagnostic output files when writing. Time-dependent reading is performed using the function ``readForcingTimeStatic()``, which is called directly when required. Writing files with a time axis occurs when ``ParametericGrid::dumpModelState()`` is called with the argument ``isRestart`` set to ``false`` and uses a specific code path, the function ``writeDiagnosticTime()``.

* ``getModelState()`` reads a file without a time axis, usually a restart file. The netCDF dimensions are mapped to the ``ModelArray`` dimensions defined specifically for discontinuous Galerkin numerics. Then read the netCDF variables from the file, adding them to the map by the variable name. If a variable has a type that maps to ``ModelArray::Type::Z`` then the dimensions are re-ordered and the correct extent vector used to call ``getVar()``. Close the netCDF file.

* ``readForcingTimeStatic()`` reads a file with a time axis, such as a forcing file. The index of the latest time point earlier than the passed target time is found. Extent and offset arrays for the data extraction are created. Get all of the variables in the netCDF file and add them to the returned ``ModelState`` data map. Close the netCDF file.

* ``dumpModelState()`` writes a restart file. First, write the common restart data. Write the dimensions from their ``ModelArray`` definitions to the netCDF file. For the hardcoded set of restart field names write that field to disk, using the appropriate extent and offset vector. Close the netCDF file.

* ``writeDiagnosticTime()`` writes data to a file with a time axis. Check if an open file by that name exists. Create it if required and note the fact to pick the 'new file' code branches in the rest of the function. Either write or get the metadata and data node names. Write common metadata if required, as well as creating the dimensions, including time. Get the dimensions out of the file and create extent and offset arrays for each type. A new file requires a landmask, which is (assumed) static with time, so create dimensions for that. Loop over all fields in the supplied ``ModelState``, writing them using the extent and offset masks, and using a special condition for the landmask field. Close the netCDF file.
