.. _clts:

Command Line Tools
##################

Opacplot2 supports a number of Opacity/EoS table formats that are listed in the
:ref:`introduction`.

It is possible to directly convert between two table formats with a suitable
Python script. A more systematic way consist in converting the source file to
the HDF5 format which is then converted to the output format::

    input table -> HDF5 EoS/Opacity format  ->  output table

Convenience scripts exist to perform this function, installed as:
``opac_table2hdf`` and ``opac_hdf2table``.

After installation of ``opacplot2``, these two scripts will be placed in the
system's path and can be called using ``opac_table2hdf`` and similarly
``opac_hdf2table``.


``opac_table2hdf``
******************

Located in ``opacplot2.scripts.opac_table2hdf``.

Usage::

    opac_table2hdf [-h] [-t {multi,ionmix,hdf5,propaceos}] [-o OUTNAME] [--Znum ZNUM] [--Xnum XNUM] input_file

This will take a file, ``input_file`` of the type specified by the ``-t`` tag
and convert it to an HDF5 file with a name specified by the ``-o`` tag.

In general, one will have to specify the atomic numbers and material fractions
by using ``--Znum`` and ``--Xnum``, respectively. These are comma separated
lists.


``opac_hdf2table``
******************

Located in ``opacplot2.scripts.opac_hdf2table``.

Usage::

    opac_hdf2table [-h] [-t {multi,ionmix,ascii}] [-o OUTFILE] input_file

This will take a file, ``input_file`` of HDF5 format and convert it to a file of
the type specified by ``-t``, with a name specified by ``-o``.

Other Tools
***********

``opacplot2``
=============

**Has not been tested!**

Located in ``opacplot2.scripts.main``.

This script is used to automate EoS/opacity tables generation for *FLASH* from the
SESAME database.

Usage::

    opacplot2 [-h] [-d DBDIR] [-t DBTYPE] -n TABLENUM [-o OUT]

where::

    -d DBDIR, --dbdir DBDIR             Path to the database. Default: current directory.

    -t DBTYPE, --dbtype DBTYPE          Database type. Currently only supports SESAME (default).

    -n TABLENUM, --tablenum TABLENUM    Table ID.

    -o OUT, --out OUT                   Output Filename. Default: DBTYPE-eos-TABLENUM


``opacdump``
============

Located in ``opacplot2.scripts.main``.

Currently not working, fails with ``KeyError: 'rho'``.

``opac_checkhdf``
=================

Located in ``opacplot2.scripts.main``.

Currently not working, fails with ``KeyError: 'dens'``.


``opac_diff``
=============

Located in ``opacplot2.scripts.opac_diff``.

Currently not working, fails with ``KeyError: 'rho'``.
