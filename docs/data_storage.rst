Model data storage
====================

Most of the data used by the model for simulating sea ice is stored in objects of the type
``ModelArray``. This class provides n-dimensional array indexing and an interface to the DG
components as well as the element-mean values that are used by the column physics calculations. The
class provides a limit number of types, enumerated by the enum class ``ModelArray::Type``, each of which
is related to a set number of axis dimensions, enumerated by ``ModelArray::Dimension``. This mechanism
ensures that arrays over the same domain have the same extent, including arrays with and without DG
components. The set of dimensions and types is set at compile time, but the lengths of the dimensions
themselves can be set at run time.

To avoid sharing data indiscriminatly, ``ModelArray`` data can be shared by ``ModelArrayRef`` objects,
which store a pointer to the ``ModelArray`` and provide pass-though functions for indexing and the four
basic mathematical operations. The ``ModelArrayRef`` class is largely designed for the column physics
part of the model, and so only provides indexing and arithmetic on the first, index 0 component of the
referenced ``ModelArray``. The one exception is the function ``allComponents()`` which returns a reference
to the underlying data type (currently ``Eigen::Array``), which includes all available components. This
allows the arrays with full DG components to be passed to the ice dynamics sub-model.

Objects of the ``ModelArrayRef`` class interact with a ``ModelArrayReferenceStore`` which allows
registration and retrieval of the pointers by name. This class also allows write access restriction
to an object. A ``ModelArray`` must be registered by name to the ``ModelArrayReferenceStore``, and optionally
can be declared RO (read-only) or read-write. If this argument is not provided, the default is RO
(read-only). A ``ModelArrayRef`` has the name of the array it targets and an access specification as template
parameters (for example, ``"hice"``, ice thickness, accessed as ``RO``, read-only). If a read-write
``ModelArrayRef`` attempts to point to the name of an array that has only been registered as read-only, a
run-time error occurs.

Four ``ModelArrayReferenceStore`` s currently exist in the model. Two exist to share the data collected
from reading the external coupling data sources, which copy data to the relevant entry of the
``ModelArrayReferenceStore`` s described below, and will not be discussed further. The main model uses two
``ModelArrayReferenceStore`` s the first contains most of the data shared around the column physics model, allowing
the array containing the sea surface latent heat flux, from the flux calculation routines, to be used in the
thermodynamics implementation. All of the ``ModelArray`` s in this store are element-mean values and have no higher
DG components. The second is the store of data to be advected. In the first instance, this includes values such as
the ice thickness and concentration, but also any arrays from any module that needs to be advected along with the
ice. This store of advected fields can then provide ``ModelArrayRef`` s to components of the column physics model
which only allow access to the element-mean value of that field. The same store can also provide, through the
``allComponents()`` member function of ``ModelArrayRef`` access to the full array including all DG components to
the ice dynamics implementation. This reference obtained is to the underlying ``DataType`` of ``ModelArray``, which
can be passed to the dynamics and ``reinterpret_cast`` ed to a ``dgVector`` for direct use in the dynamics, without
any copying of data.
