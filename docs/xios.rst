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

The ``model`` section, which is used elsewhere in nextSIM-DG, contains two
entries relevant to XIOS: ``start`` and ``time_step``. These are used to
configure the calendar used by XIOS and take the following default values.

.. code-block::
  [model]
  start = 1970-01-01T00:00:00Z
  time_step = P0-0T01:00:00

Files and fields to be read and written to files are configured via
``XiosInput`` and ``XiosOutput`` sections, which accept the same entries. For
example, we could set the following values, which would specify two fields
labelled ``field_A`` and ``field_B``, which are written to file
``my_output_file.nc`` every two (simulated) hours.

.. code-block::
  [XiosOutput]
  period = P0-0T02:00:00
  filename = my_output_file.nc
  fields = field_A,field_B

Note that both the ``XiosInput`` and ``XiosOutput`` sections are optional.

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
