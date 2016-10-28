#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import absolute_import
#from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals

import tables
import numpy as np
import sys
import types
from six import iteritems

class OpgHdf5(dict):
    @classmethod
    def open_file(cls, filename, explicit_load=False):
        """
        Open an HDF5 file containing opacity data.

        Parameters
        ----------
        filename : str
                  Name of file to open.
        explicit_load : bool
                  Option to load the whole file to memory.
        
        Examples
        --------
        To open a file::

           >>> import opacplot2 as opp
           >>> op = opp.OpgHdf5.open_file('infile.h5')
           >>> print(op.keys())
           dict_keys(['Anum', ..., 'Znum']) # OpgHdf5 is a dictionary.
           >>> print(op['Zf_DT'])
           array([...]) # Array for average ionization.
        
        Notes
        -----
        Loading the entire file into memory can be potentially dangerous. Use
        ``explicit_load`` with caution.
        """
        self = cls()
        self.f = f = tables.open_file(filename, 'r')
        for el in f.root:
            if isinstance(el, tables.group.Group):
                # loading ion_frac group
                key = el._v_name
                self[key] = {}
                for subel in el:
                    self[key][subel.name] = subel
            else: #isinstance(tables.carray.CArray)
                self[el.name] = el
        for key in f.root._v_attrs._f_list():
            self[key] = f.root._v_attrs[key]

        if explicit_load: self.force_eval()
        self._compute_ionization()
        self.Nr =  self['dens'].shape[0]
        self.Nt =  self['temp'].shape[0]
        self.Ng =  self['groups'].shape[0] - 1
        return self

    def write2file(self, filename, **args):
        """Write to an HDF5 output file.
        
        Parameters
        ----------
        filename : str
           Name of output file.
        args : dict
           Dictionary of data to write.
        
        """
        import tables
        h5filters = tables.Filters(complib='blosc', complevel=7, shuffle=True)
        f = tables.open_file(filename, 'w', filters=h5filters)

        ATTR_LIST = ['BulkMod', 'ElemNum', 'Abar', 'Zmax']

        for key in sorted(self.keys()):
            if key in ATTR_LIST or key in ["ion_frac"]: continue
            if key in args:
                val = args[key]
            else:
                if self[key] is None: continue
                val = self[key][:]
            atom = tables.Atom.from_dtype(val.dtype)
            ds = f.createCArray(f.root, key, atom, val.shape, filters=h5filters)
            ds[:] = val
        
        # I believe we should only put ion_frac in the table if it was already
        # in the data. -JT
        if 'ion_frac' in self:
            f.createGroup(where='/', name='ion_frac', filters=h5filters)
            for  ion_frac_key,  ion_frac_val in iteritems(self['ion_frac']):
                atom = tables.Atom.from_dtype(ion_frac_val.dtype)
                ds = f.createCArray(f.root.ion_frac, ion_frac_key, atom, ion_frac_val.shape)
                ds[:] = ion_frac_val[:]

        # writing attributes
        for attr in ['BulkMod', 'ElemNum', 'Abar', 'Zmax']:
            if attr in self:
                setattr(f.root._v_attrs,attr, self[attr])
        f.close()

    def force_eval(self):
        """
        Load the whole table into memory.
        """
        for key, val in iteritems(self):
            if type(val) is tables.carray.CArray:
                self[key] = val[:]
            elif type(val) is dict:
                for key_in, val_in in iteritems(val):
                    print(key, key_in)
                    self[key][key_in] = val_in[:]
    
    def _compute_ionization(self):
        """
        Compute ionization from populations if available.
        """
        if 'ion_frac' not in self:
            self['Zfo_DT'] = None
        else:
            DT_shape = self['Zf_DT'].shape
            self['Zfo_DT'] = np.zeros(DT_shape)
            self['ion_frac_sum'] = np.zeros(DT_shape)
            for Zel, Zfrac  in iteritems(self['ion_frac']):
                IonLvls = np.arange(int(Zel[1:])+1)
                IonLvls_arr = np.tile(IonLvls, DT_shape).reshape(DT_shape +(-1,))
                self['Zfo_DT'] += (Zfrac*IonLvls_arr).sum(axis=-1)
                self['ion_frac_sum'] += np.sum(Zfrac, axis=-1)





