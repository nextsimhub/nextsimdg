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

RectGridIO
----------

* ``getModelState()`` sets the size of the ``H``, ``U``, ``V`` and ``Z`` ``ModelArrays`` using the ``dimensionSetter()`` function. With the arrays correctly sized, the function reads the data for the land mask, the ice thickness, the ice concentration, the snow thickness and the ice temperature. The netCDF variable names are those taken from ``gridnames.hpp` and mapped to the same names within the ``ModelState`` data map. Close the netCDF file.

* ``dumpModelState()`` writes the comman restart metadata using the ``CommonRestartMetadata`` class. The dimensions are written using the sizes of ``H``and ``Z``arrays. All of the data in passed ``ModelState`` is then written into the file based on whether it is an ``H`` array
or a ``Z`` array. Close the netCDF file.
