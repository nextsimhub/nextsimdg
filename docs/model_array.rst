.. Copyright (c) 2023, Nansen Environmental and Remote Sensing Center

ModelArray and ModelArrayRef
============================

Introduction
------------

The ModelArray and ModelArrayRef classes provide the data arrays that are used
throughout NextSIM_DG. They are designed to act as multi-dimensional arrays
with a size fixed at run-time, but shared between arrays of the same type. The
ModelArray class provides the actual data storage, while the ModelArrayRef
class provides shared access for data throughout the model between different
ModelComponents.

ModelArray
----------

The ModelArray class is the data storage class within the nextSIM_DG model. It
is based on the Eigen::Array class, which provides the underlying memory and
computational functionality. The size of a ModelArray is dynamic and set at
run-time. The arrays are assigned to various Types, which define the number and
sizes of the array dimensions, which translate from the underlying
one-dimensional array to the multi-dimensional array indexing provided by the
ModelArray (and ModelArrayRef) interface.

The interface that ModelArray provides includes multi-dimensional indexing,
access to non-spatial components such as finite element components. Special
access to the last spatial dimension is provided by the function
zIndexAndLayer(). This is used when one ModelArray::Type defines an extra
spatial dimension compared to another. An example might be when one Type
defines a two-dimensional spatial array and another defines an array that
covers the same horizontal domain but adds vertical levels. The function takes
the one dimensional index of an element in the two-dimensional array plus an
index for the vertical level. The data reference returned is to the data at the
given vertical level above the horizontal position referred to by the
one-dimensional index.

The class
also provides a large set of mathematical functions which act upon the data.
These include the four basic arithmetic functions, interacting with both other
ModelArrays and double-precision scalars, including assigned arithmetic
operators such as the += operator. Also included are element-wise maximum
(a.max(b)) and minimum (a.min(b)) functions and functions to clamp all elements
to a maximum or minimum value (a.clampXXX(b)).

The class also provides information on the type and number and length of
dimensions of any given array. One important pair of functions when developing
new code using ModelArray are the size() and trueSize() functions. The first,
size(), provides the logical size of the array from the definitions of the
array dimensions (see below). The second, trueSize(), provides the current size
of the array in memory, which may be different if the defining array dimensions
have been resized, since this information cannot be efficiently propagated
automatically, and so is not.

ModelArrayDetails
^^^^^^^^^^^^^^^^^

The Types of a ModelArray are defined at compile time in the headers
ModelArrayDetails.hpp and ModelArrayTypedefs.hpp and the source file
ModelArrayDetails.cpp.

The first of these, ModelArrayDetails.hpp, defines a set
of dimensions that will be used by the ModelArray::Types. This is a enum class
of the names of the dimensions. These enum values will be useable throughout
the code that includes this particular ModelArrayDetails information. There is
also an enum class of the Types of array that can be used. These will each use
some of the dimensions in their definitions, which are defined in the
ModelArrayDetails.cpp source file.

The ModelArrayDetails.cpp source file defines a map between the
ModelArray::Dimension enum members defined in the header and a pair of the
string name of the dimension and the length. It is the length here that defines
the length of the dimension. The string name is used to name the dimension when
reading or writing netCDF files. A map from ModelArray::Type to a vector of
Dimensions named typeDimensions is defined. This defines the dimensions of each
type. There is also a map from ModelArray::Type to a string name for the Type.
This source file also defines the default ModelArray constructor, which sets
the default type for a ModelArray. A function named hasDoF() returns true when
passed a ModelArray::Type that has non-spatial dimension, such as discontinuos
Galerkin degrees of freedom. There is also the definition of the constructor
for the ModelArray size map m_sizes. This is necessary as the entries in the
map must be the ModelArray::Types defined in the header. The definition of the
default constructor for the ModelArray::DimensionMap class is also included
here. Care must be take to stay consistent with the contents of typeDimensions.
Finally, a map from Type to Dimension records the dimension that defines the
number of components for the types that have them.

The other header ModelArrayTypedefs.hpp defines a typedef for each
ModelArray::Type to allow direct definition of that type of array, rather than
using the syntax ModelArray(ModelArray::Type::TYPE).

ModelArrayRef
-------------

The ModelArrayRef class is a way of sharing data around the model which is more
informative than passing around raw pointers or references. The key is a store
(an instance of MARStore) that holds the pointers to the data. A template class
then uses this store along with a string defining the array required to get the
pointer, wrapped in a class that provides indexed array access, some arithmetic
operators and direct data access. In this way data from various components of
the model (instances of ModelComponent) can access data stored in others, while
maintaining a common interface between local and referenced data arrays.

A ModelArrayRef is created by in two parts. The source of the data is created
by registering a ModelArray as the provider for a given dataset name. This name
is wrapped in a TextTag templated utility class to allow the string to be
constexpr. The provided data may be shared as read-write (RW) or read-only
(RO). The other part is the declaration of a ModelArrayRef variable with the
given name, similarly wrapped in a TextTag. This reference may also be defined
as being read-write or read-only using the same names. The instance then needs
to be constructed with the store that contains the data pointer as an argument.

The design of the ModelArrayRef and MARStore classes allows ModelArrayRef
instances to be created both before and after the registry of the array that
provides the data. The data source may even be re-registered and the references
will afterwards refer to the new source.

The data access is such that an array registered as read-only will only be
available to read-only ModelArrayRefs using the same data field name. If a
read-write ModelArrayRef tries to refer to a dataset that has been registered
as read-only, the pointer will not be connected and a segmentation violation
will occur at run time. An array registered as a read-write data source can be
referred to by either read-write or read-only ModelArrayRefs. A read-only data
source will function correctly when referred to by a read-only ModelArrayRef.