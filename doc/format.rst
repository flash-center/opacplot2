File Format Functions and Classes
#################################

.. _opg_ionmix:

IONMIX (.cn4)
*************

.. autoclass:: opacplot2.OpacIonmix
   :members:

.. autofunction:: opacplot2.writeIonmixFile
   
.. autoclass:: opacplot2.adapt.EosMergeGrids
   :members:

MULTI (.opp, .opr, .eps)
************************

.. autoclass:: opacplot2.OpgMulti
   :members:


.. autofunction:: opacplot2.get_related_multi_tables

SESAME ASCII (.ses)
*******************

.. autoclass:: opacplot2.OpgSesame
   :members:

.. note:: There are currently no handling functions for 500 series entries in the SESAME tables.
        
Data Prefixes
=============

From `Los Alamos National Laboratory <http://www.lanl.gov/org/padste/adtsc/theoretical/physics-chemistry-materials/_assets/docs/LAUR-92-3407.pdf>`_,
each SESAME table has five different parts of their EoS tables:

   +-----------+-----------------------------------------+------------------------+
   | Table #   |   Data Type                             | ``opacplot2`` prefix   |
   +===========+=========================================+========================+
   | Table 301 |   TotalEOS(304+305+306)                 | ``total_``             |
   +-----------+-----------------------------------------+------------------------+
   | Table 303 |   Ion EOS Plus Cold Curve (305 + 306)   | ``ioncc_``             |
   +-----------+-----------------------------------------+------------------------+
   | Table 304 |   Electron EOS                          | ``ele_``               |
   +-----------+-----------------------------------------+------------------------+
   | Table 305 |   Ion EOS (Including Zero Point)        | ``ion_``               |
   +-----------+-----------------------------------------+------------------------+
   | Table 306 |   Cold Curve (No Zero Point)            | ``cc_``                |
   +-----------+-----------------------------------------+------------------------+

If we wanted to print out all of the EoS data for the
electrons::

   >>> import opacplot2 as opp
   >>> op = opp.OpgSesame('sesame.ses', opp.OpgSesame.SINGLE)
   >>> data = op.data[13719]
   >>> for key in data.keys(): # Table ID for aluminum.
   ...     if 'cc_' == key[:3]:
   ...         print(key+':')
   ...         print(data[key])

Data Points
===========

There are several top-level data points that are not included in a specific
curve of the SESAME tables:

   +--------------------------+----------------------------+
   | ``opacplot2`` Abbr.      | Phys. Meaning              |
   +==========================+============================+
   | ``abar``                 | Mean Atomic Mass           |
   +--------------------------+----------------------------+
   | ``bulkmod``              | Bulk Modulus               |
   +--------------------------+----------------------------+
   | ``excoef``               | Exchange Coefficient       |
   +--------------------------+----------------------------+
   | ``rho0``                 | Normal Density             |
   +--------------------------+----------------------------+
   | ``zmax``                 | Mean Atomic number         |
   +--------------------------+----------------------------+

Furthermore, each table (301, ..., 305) prefixes these data points:

   +-----------------------------+-------------------------------+
   | ``opacplot2`` Abbr.         | Phys. Meaning                 |
   +=============================+===============================+
   | ``dens``                    | Density                       |
   +-----------------------------+-------------------------------+
   | ``eint``                    | Energy                        |
   +-----------------------------+-------------------------------+
   | ``ndens``                   | Number of Densities           |
   +-----------------------------+-------------------------------+
   | ``ntemp``                   | Number of Temperatures        |
   +-----------------------------+-------------------------------+
   | ``pres``                    | Pressures                     |
   +-----------------------------+-------------------------------+
   | ``temps``                   | Temperatures                  |
   +-----------------------------+-------------------------------+

From the example in the previous section, if we wanted to print out the 
density array of electrons, we would simply type::

   >>> print(data['ele_dens'])
   array([...]) # Electron density array.

PROPACEOS ASCII (.prp)
**********************

.. warning:: Handling for Propaceos EoS tables is not publicly distributed with ``opacplot2`` and so its documentation will not be presented here.

HDF5
****

.. autoclass:: opacplot2.OpgHdf5
   :members:

Data Points
===========

Depending on what kind of data has been written to an HDF5 file, each HDF5 
file may vary widely. Despite this fact, listed below are the naming conventions
``opacplot2`` uses for HDF5 data points that are relevant to IONMIX. The
abbreviations are keys unless otherwise stated to be an attribute to the 
``OpgHdf5`` object.
   
   +----------------------------+---------------------------+
   | ``opacplot2`` Abbr         | Physical Meaning          |                                                                               
   +============================+===========================+
   | ``Znum``                   | Atomic numbers            |                                                                               
   +----------------------------+---------------------------+
   | ``Xnum``                   | Element fractions         |                                                                               
   +----------------------------+---------------------------+
   | ``idens``                  | Number densities          |                                                                               
   +----------------------------+---------------------------+
   | ``temp``                   | Temperature array         |                                                                               
   +----------------------------+---------------------------+
   | ``Ng`` (attr)              | Number of energy groups   |                                                                               
   +----------------------------+---------------------------+
   | ``groups``                 | Energy group boundaries   |                                                                               
   +----------------------------+---------------------------+
   | ``opp_mg``                 | Planck absorption         |                                                                               
   +----------------------------+---------------------------+
   | ``opr_mg``                 | Rosseland opacities       |                                                                               
   +----------------------------+---------------------------+
   | ``emp_mg``                 | Planck emissivity         |                                                                               
   +----------------------------+---------------------------+

For example, if we wanted to open up an HDF5 file named ``input.h5`` and write it to an IONMIX file
named ``output.cn4``, we could write::

   >>> import opacplot2 as opp
   >>> op = OpgHdf5.open_file('input.h5')
   >>> opp.writeIonmixFile(outfile,
                    op['Znum'], op['Xnum'],
                    numDens=op['idens'][:], temps=op['temp'][:],
                    ngroups=op.Ng,
                    opac_bounds=op['groups'][:],
                    planck_absorb=op['opp_mg'][:],
                    rosseland=op['opr_mg'][:],
                    planck_emiss=op['emp_mg'][:])
