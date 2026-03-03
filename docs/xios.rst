.. Copyright (c) 2025, Nansen Environmental and Remote Sensing Center

Using XIOS for reading and writing files
========================================

Introduction
------------

`XIOS`_ stands for *XML Input/Output Server*. It is an asynchronous I/O tool
widely used in the climate modelling community. While it is traditionally
configured using user-provided XML files, nextSIM-DG configures it directly with
API calls to reduce user requirements.

.. _XIOS: https://forge.ipsl.fr/ioserver

The integration of XIOS into nextSIM-DG is built around a static ``Xios``
handler class, which provides a C++ API for the various XIOS functions. When the
handler is instantiated, the configuration sections (see section below) are
parsed. Based on the values that are parsed, the handler object will
automatically create data structures for XIOS and associate these as
appropriate. Further details on how this works can be found in the following.

Configuration
-------------

XIOS is configured automatically from two sources, as detailed in the following.

XML configuration
^^^^^^^^^^^^^^^^^

Given that XIOS is an XML-based approach, it does require a small XML file to
configure information that is known at compile time, as well as how errors and
logs are handled. As such, you will need to provide an ``iodef.xml`` file, which
needs to exist in the directory you intend to run from. A Jinja2 template is
provided as ``core/src/iodef.xml.jinja``, which accepts ``DGCOMP`` and
``DGSTRESSCOMP`` as inputs. A helper script is provided for using the template
to generate an XML file:

.. code-block::

   python3 core/src/generate_iodef.py \
       <DGCOMP> \
       <DGSTRESSCOMP> \
       core/src/iodef.xml.jinja
       <OUTPUT_FILE_NAME>

where ``DGCOMP`` is the integer number of DG components (e.g., 6),
``DGSTRESSCOMP`` is the integer number of DG stress components (e.g., 8), and
``OUTPUT_FILE_NAME`` is the output file name, including its path, e.g.,
``build/core/test/iodef.xml`` or ``run/iodef.xml``.

NextSIM-DG configuration
^^^^^^^^^^^^^^^^^^^^^^^^

As mentioned above, the XML configuration handles information that is known at
compile time. Information that is not known at compile time is configured by
nextSIM-DG through calls to XIOS' API. This information is provided in the same
way as information is provided to the rest of nextSIM-DG - from files with
``.cfg`` extension that are passed to the ``nextsim`` executable. There are
several configuration sections that are relevant to XIOS, as detailed in the
following.

The ``xios`` section contains a single entry, which determines whether or not to
build nextSIM-DG with XIOS as the I/O driver.

.. code-block::

  [xios]
  enable = true

However, there is no need to explicitly add this to the configuration because it
will be added automatically during the ``Model`` configuration step provided
nextSIM-DG has been built with XIOS support.

The ``model`` section, which is used elsewhere in nextSIM-DG, contains several
entries relevant to XIOS. The ``start``, ``stop``, and ``time_step`` entries are
used to configure the calendar used by XIOS. For example,

.. code-block::

  [model]
  start = 1970-01-01T00:00:00Z
  stop = 1970-01-01T01:00:00Z
  time_step = P0-0T01:00:00

The filename and period for restart files is configured in the same way as when
building without XIOS. That is, the ``model`` section should include
``init_file``, ``restart_file`` and ``restart_period`` entries:

.. code-block::

  [model]
  init_file = my_init_file.nc
  restart_file = my_restart_file.nc
  restart_period = P0-0T02:00:00

Information related to fields to be written to files are currently configured
via the ``XiosOutput``, and ``XiosDiagnostic`` sections, where the first two
refer to restarts. Note that both of these sections are optional. Moreover, they
provide an interim solution. Eventually, it will not be necessary to provide
XIOS-specific configuration information.

For restart files, the filename is determined from the ``restart_file`` entry in
the ``model`` section mentioned above. In the case of diagnostics files, the
file names are configured differently, via the ``filename`` entry in the
``XiosDiagnostic`` section. For example,

.. code-block::

   [XiosDiagnostic]
   filename = my_diag_file.nc

Restart and diagnostic file names may include format strings such as
``restart%Y-%m-%dT%H:%M:%SZ.nc`` or ``diagnostic%Y-%m-%dT%H:%M:%SZ.nc`` (in
fact, these are the defaults). When writing out restarts and diagnostics, a
separate file is produced for each restart period, with filename of the format
``<FILENAME>_<START_TIME>-<END_TIME>.nc``, where ``<START_TIME>`` and
``<END_TIME>`` are the start and end of the associated period, written using the
provided format string. Restart files contain the state at the beginning of the
time window, whereas diagnostics files are averaged over the timesteps in the
window.

In addition to specifying file names, we need to specify which fields are to be
written using each I/O type. For example, we could specify that two fields
labelled ``field_A`` and ``field_B`` are to be written into restart files as
follows:

.. code-block::

  [XiosOutput]
  field_names = field_A,field_B

The ``field_names`` entry may contain a single field name or a comma-separated
list. Note that this is an interim solution. Eventually, it will not be
necessary to provide field names to configure XIOS.

As elsewhere in the model, the configuration values above are all parsed by
calling the ``Model.configure()`` member function. Note that this function will
automatically initialize the model state based on the ``input_file`` entry and
all field variables found in that file, including variables defining the grid
i.e., ``longitude`` and ``latitude`` and possibly ``coords`` and
``grid_azimuth``.

Forcing files
-------------

Forcing files are handled in a pre-defined fashion. The ``ERA5Atmosphere``
config section is used for atmospheric forcings and the ``TOPAZOcean`` config
section is used for ocean forcings. The only entry that they support is ``file``
to indicate the corresponding forcing file name. In both cases, pre-processed
files are currently expected, with the following properties:

* Pre-processed ERA5 forcing files use a 1 hour frequency and include variables
  ``dew2m``, ``lw_in``, ``sw_in``, ``pair``, ``tair``, ``wind_speed``, ``u``,
  and ``v``.
* Pre-processed TOPAZ forcing files use a 1 day frequency and include variables
  ``mld``, ``sss``, ``sst``, ``ssh``, ``u``, and ``v``.

In both cases, it is assumed that the variables are defined to have ``HField``
type based on the same grid as used in the model.

Note that this is an interim approach. Eventually, it will be possible to read
forcing data with XIOS directly from original climatology files defined on
different grids.

Order of operations for XIOS setup
----------------------------------

The required order of operations to set everything up correctly for reading and
writing files using XIOS is demonstrated in the
``core/test/XiosReadRestart_test.cpp`` and
``core/test/XiosWriteRestart_test.cpp`` worked examples and is elaborated in the
following:

1. Set configuration options as for other parts of the model. (See section above
   for options.)
2. Configure the ``Model`` by constructing it and calling its ``configure``
   member function.
3. Construct a ``ParametricGrid`` and a new ``ParaGridIO`` pointer.
4. Associate the ``ParaGridIO`` pointer with the ``ParametricGrid`` instance
   using the latter's ``setIO`` member function.
5. Get the ``Xios`` handler singleton using ``Xios::getInstance``.
6. For each field to be written to file with XIOS, call the
   ``setPrognosticFieldType`` or ``setDiagnosticFieldType`` member function of
   ``Xios`` as appropriate, providing the field name as the first argument
   and the ``ModelArray::Type::<TYPE>`` enum as the second argument (replacing
   ``<TYPE>>`` with the desired type).

There is no need to explicitly close the XIOS context definition because this
happens automatically upon reading or writing. As such, all XIOS configuration
must take place before any files are read from or written to.

XIOS concepts
-------------

XIOS has several key concepts, five of which are used in nextSIM-DG. Whilst it's
not important to understand all of them in detail, it's useful to have an idea
of how they are used.

XIOS Domain concept
^^^^^^^^^^^^^^^^^^^

The Domain type defines the horizontal domain and its MPI parallel
decomposition. The Domain definition requires the global size of the :math:`x`-
and :math:`y`-dimensions of the domain, as well as the local sizes and start
indices (corner) in each dimension for each MPI rank. This information is
provided to XIOS automatically upon configuring the ``Model``.

Different domains are used for different field types, depending on the number of
degrees of freedom (DoFs) they have in the horizontal. The ``HDomain`` has one
DoF per cell and is used for the ``HField``, ``DGField``, and ``DGSField``
types. The ``VertexDomain`` has one DoF per vertex and is used for the
``VertexField`` type. Finally, the ``CGDomain`` has the same (degree-dependant)
number of DoFs as the continuous Galerkin discretisation and is used by
``CGField``.

XIOS Axis concept
^^^^^^^^^^^^^^^^^

The Axis type is used to define an additional dimension for some field types.
For ``VertexField``, it gives rise to a vector field with as many components as
the spatial dimension, i.e., two. For ``DGField`` and ``DGSField``, the Axis
concept is used to give rise to vector fields with as many components as
``dg_comp`` and ``dgstress_comp``, respectively. These Axes are set up
automatically based on the XML configuration.

XIOS Grid concept
^^^^^^^^^^^^^^^^^

The Grid type is used to define the discretisation for each field. That is, it
associates a Domain and (in some cases) an Axis with a Field.

XIOS Field concept
^^^^^^^^^^^^^^^^^^

An instance of the Field type is based on a Grid and is associated with a
``ModelArray::Type``. Fields are created automatically upon configuring the
``Model``. For fields that are read from file, the field type is determined
automatically during the file read. For fields that are written, it's required
to call the ``setPrognosticFieldType`` or ``setDiagnosticFieldType`` member
function of the ``Xios`` singleton, providing the field name and
``ModelArray::Type``.

XIOS File concept
^^^^^^^^^^^^^^^^^

The File concept holds metadata related to input and output files, including
which Fields are associated with it and whether it is to be used for reading or
writing. Files are created automatically upon configuring the ``Model``.

Developer notes
---------------

* While XIOS and nextSIM-DG are both written in C++, the majority of XIOS' users
  work in Fortran and interact with it via its Fortran bindings. These are based
  off its C interface using the ``iso_c_binding`` Fortran module. Instead of
  ``include``-ing the XIOS source directly in nextSIM-DG, we decided to wrap the
  C interface. The header file ``core/src/include/xios_c_header.hpp`` provides
  an interface for the functions that are used.

* The nextSIM-DG XIOS integration is set up such that the filename and field
  names coincide with the fileId and fieldId of the corresponding ``File`` and
  ``Field`` objects.

* Allowing fields to be both read and written in the same run requires
  additional fields. This is because the two may have different I/O modes, e.g.,
  reading or writing at a specified frequency (``"instant"```), writing after
  applying a reduction operator (e.g., ``"average"``), or reading or writing
  once at the start of the time window (``"once"``). We use the 'base' field
  associated with the field ID for writing. For reading a field, we create a
  separate 'inherited' field, which has the same name but with ID appended by
  ``"_input"`` or ``"_forcing"``, depending on whether the field is to be read
  as a restart or as a forcing. The inherited field references the base one so
  that they share the same data, but differ in how they are handled for I/O. All
  of this happens automatically in the ``createField`` function.
