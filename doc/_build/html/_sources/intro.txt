.. _introduction:

Introduction
############

*FLASH* uses the IONMIX4 format for its equations of state (EOS) and opacities.
However, several other file formats for EOS and opacity data are popular, such
as Propaceos and SESAME that are incompatible with *FLASH.* So, the Python
module ``opacplot2`` was created. It contains read/write/write2hdf capabilities for
various table formats, and as such it can be used to convert between the
formats as needed for *FLASH* and other MHD Codes. Listed below are the current
capabilities of ``opacplot2``:

+--------------------------+----------+----------+-----------+----------------------------------------------+
| Name                     | Reader   | Writer   | Write2hdf | Comments                                     |
+==========================+==========+==========+===========+==============================================+
| Native HDF5              |    ✔     |   ✔      |           |                                              |
+--------------------------+----------+----------+-----------+----------------------------------------------+
| IONMIX (.cn4)            |    ✔     |   ✔      |           |                                              |
+--------------------------+----------+----------+-----------+----------------------------------------------+
| MULTI (.opp, .opr, .eps) |    ✔     |   ✔      |    ✔      |                                              |
+--------------------------+----------+----------+-----------+----------------------------------------------+
| SESAME ASCII (.ses)      |    ✔     |          |     ✔     |                                              |
+--------------------------+----------+----------+-----------+----------------------------------------------+
| SESAME Binary (.sesb)    |          |          |           | `see pyeospac`_                              |
+--------------------------+----------+----------+-----------+----------------------------------------------+
| TOPS                     |    ✔     |          |           |                                              |
+--------------------------+----------+----------+-----------+----------------------------------------------+
| INFERNO                  |    ✔     |          |           |                                              |
+--------------------------+----------+----------+-----------+----------------------------------------------+
| UW EoS                   |    ✔     |          |           |                                              |
+--------------------------+----------+----------+-----------+----------------------------------------------+
| PROPACEOS ASCII (.prp)   |    ✔     |          |      ✔    | not distributed                              |
+--------------------------+----------+----------+-----------+----------------------------------------------+
| Opal                     |          |          |           | `see pystellar`_                             |
+--------------------------+----------+----------+-----------+----------------------------------------------+

.. _see pyeospac: http://github.com/luli/pyeospac
.. _see pystellar: https://github.com/alexrudy/pystellar/blob/master/pystellar/opacity.py

Read/Write/Write2hdf
********************

If a format can be *read* by ``opacplot2``, it means that ``opacplot2`` can browse a
file of that particular format and it will return a dictionary object
containing keys corresponding to particular data types (such as the planck
emmisivity, ``emp_mg``) that match up with the ``numpy`` arrays of the data.
*Write* cabailities are not as simple. Some format submodules have different
methods of writing to files than others. *Write2hdf* means that the format
submodule can open up and read a file of that format, and then write that data
to an HDF5 (.h5) file. The HDF5 format is somewhat of a middleman for most
conversions, since it can store any kind of data from a file. Then, using the
data stored in the HDF5 file, format submodules with *write* capabilities can
write to a file.


Installation
************

.. _prerequisites:

Prerequisites
=============

Python Packages
---------------

Required packages:

#. ``numpy``
#. ``pytables``
#. ``six``
#. ``setuptools``
#. ``periodictable``
#. ``matplotlib``
#. ``scipy``

On Linux
--------

You could just use your system Python distribution and install numpy, pytables,
six, periodictable with the system package manager. On Debian-based systems
(Ubuntu, etc) this can be achieved with::

    sudo apt-get install python-numpy python-tables python-six python-setuptools
    sudo pip install periodictable matplotlib


On OSX/Windows
--------------

The recommended way is to use Anaconda Python distribution (although as long as
you install the required packages this hardly matters):

#. (optional)  `Download <https://www.continuum.io/downloads>`_ and install Anaconda Python
    3.5 (or 2.7).
#.  Install the required packages::

        conda install numpy tables six setuptools matplotlib
        pip install periodictable


Install Opacplot2
=================

#.  Setup a scientific Python environment including numpy and PyTables (see :ref:`prerequisites`)
#.  Download the LULI fork of ``opacplot2`` `here <https://github.com/luli/opacplot2>`__.
#.  (optional) if you have the opg_propaceos.py script (not distributed on github),
    copy it to the ``opacplot2-master/opacplot2/`` folder, then make sure it is
    ignored by git::

        cp opg_propaceos.py opacplot2-master/opacplot2/

#.  Install ``opacplot2``::

        cd opacplot2-master/
        python setup.py install

