.. Copyright (c) 2025, Nansen Environmental and Remote Sensing Center

Using XIOS for reading and writing files
========================================

Introduction
------------

Section under construction.

XIOS concepts
-------------

Section under construction.

Configuration
-------------

XIOS is configured automatically from files with ``.cfg`` extension that are
passed to the ``nextsim`` executable. There are several configuration sections
that are relevant to XIOS.

The ``xios`` section contains a single entry, which determines whether or not to
build nextSIM-DG with XIOS as the I/O driver.

.. code-block::

  [xios]
  enable = true

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

In the cases of forcing and diagnostics files, the file names are configured
differently, via the ``filename`` entry in the ``XiosForcing`` and/or
``XiosDiagnostic`` sections, respectively. For example,

.. code-block::

   [XiosForcing]
   filename = my_forcing_file.nc

The fields to be read from and written to files are configured via the
``XiosInput``, ``XiosOutput``, and ``XiosForcing`` sections, where the first two
refer to restarts. For example, we could specify that two fields labelled
``field_A`` and ``field_B`` are to be written into restart files as follows:

.. code-block::

  [XiosOutput]
  field_names = field_A,field_B

The ``field_names`` entry may contain a single field name or a comma-separated
list. Note that all of the ``XiosOutput``, ``XiosInput``, ``XiosDiagnostic``,
and ``XiosForcing`` sections are optional.

As elsewhere in the model, these configuration values are parsed by calling the
``Model.configure`` member function. Since building with XIOS implies also
building with MPI, you will need to pass the MPI communicator to this member
function when calling it.

The XIOS I/O implementation supports fields of ``HField``, ``VertexField``,
``DGField``, ``DGSField``, and ``CGField``. When reading from file, the field
type will be deduced from its dimension. When writing to file, the field type
should be set using the ``setFieldType`` member function of the ``Xios`` handler
class.

Developer notes
^^^^^^^^^^^^^^^

The integration of XIOS into nextSIM-DG's is built around a static ``Xios``
handler class, which provides a C++ API for the various XIOS functions. When the
handler is instantiated, the config sections above will be parsed. Based on the
values that are parsed, the handler object will automatically create ``Field``
and ``File`` data structures for XIOS and associate these as appropriate. If the
``XiosInput`` section is used, an attempt will be made to open the file provided
by the ``filename`` entry. Note that the nextSIM-DG XIOS integration is set up
such that the filename and field names coincide with the fileId and fieldId of
the corresponding ``File`` and ``Field`` objects.

XIOS requires information on the domain decomposition for it to be able to read
and write data in parallel. This information is held by the ``ModelMetadata``
class, which may be constructed based off a partition metadata file. Upon
construction of this object, if XIOS is enabled then the XIOS handler will
automatically create a ``Domain`` data structure with the appropriate local and
global sizes and offsets.
