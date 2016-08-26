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
from .tests.suite import testsuite

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
        self.f = f = tables.openFile(filename, 'r')
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
        self._initialize_ionfrac_tests()
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
        f = tables.openFile(filename, 'w', filters=h5filters)

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

        f.createGroup(where='/', name='ion_frac', filters=h5filters)
        if 'ion_frac' in self:
            for  ion_frac_key,  ion_frac_val in iteritems(self['ion_frac']):
                atom = tables.Atom.from_dtype(ion_frac_val.dtype)
                ds = f.createCArray(f.root.ion_frac, ion_frac_key, atom, ion_frac_val.shape)
                ds[:] = ion_frac_val[:]

        # writing attributes
        for attr in ['BulkMod', 'ElemNum', 'Abar', 'Zmax']:
            if attr in self:
                setattr(f.root._v_attrs,attr, self[attr])
        f.close()

    # Would we want force_eval() to count in coverage tests, since it can
    # be dangerous to load entire file into memory? - JT
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

    test_keys=['dens', 'temp']
    
    def run_testsuite(self, mode='short'):
        for attr in dir(self):
            if 'check_' in attr:
                getattr(self, attr)(mode)
        print('')

    @testsuite(test_keys, var='Zf_DT', mode='short')
    def check_Zf_gt_0(self, mode):
        """Checking that Zf_DT>0!"""
        return self['Zf_DT'][:] > 0

    @testsuite(test_keys, var='Zf_DT', mode='short')
    def check_Zf_lt_Zmax(self, mode):
        """Checking that Zf_DT<Zmax!"""
        return self['Zf_DT'][:] <= self['Zmax']

    @testsuite(test_keys, var='Zfo_DT', mode='short')
    def check_Zfo_lt_Zmax(self, mode):
        """Checking that Zfo_DT<Zmax!"""
        return self['Zfo_DT'][:]<= self['Zmax']

    @testsuite(test_keys, var='Anum_prp', mode='short')
    def check_Abar_prp(self, mode):
        """Checking that Anum_prp is consistant with Abar!"""
        return np.isclose(self['Anum_prp'][:], self['Anum'][:], atol=0.05)


    def _initialize_ionfrac_tests(self):
        """ Initialize tests of ionfrac """
        if 'ion_frac' in self:
            for Znum in self['Znum']:
                key = 'Z{0:02}'.format(Znum)

                @testsuite(var=lambda x: x['ion_frac'][key], mode='short')
                def check_ionfrac_lt_1(self, mode):
                    """Checking that ionfrac is less then 1!"""
                    return self['ion_frac'][key][:]<=1.0
                setattr(self, 'check_ion_frac_{0}_lt_1'.format(key),
                            types.MethodType(check_ionfrac_lt_1, self))
    
                @testsuite(var=lambda x: x['ion_frac'][key], mode='short')
                def check_ionfrac_sum_1(self, mode):
                    """Checking that ionfrac sums to 1!"""
                    return np.isclose(np.sum(self['ion_frac'][key][:], axis=-1), 1.0, atol=0.01)
                setattr(self, 'check_ion_frac_{0}_sum_1'.format(key),
                            types.MethodType(check_ionfrac_sum_1, self))
    
    # Would we want to include our error message in our coverage tests? -JT
#    def _repr_arr_error(self, idx, var, msg):
#        """
#        Print an error message when array failed to pass consistency checks.
#        """
#        out = ['-'*80]
#        out.append(' Test failed for {0}/{1} points: {2}'.format(idx.sum(),idx.size, msg))
#        out.append('-'*80)
#        arr2str = np.array2string
#
#        if idx.size >= self['temp'][:].size* self['dens'][:].size:
#            out.append(' == density mask ==')
#            out.append(arr2str(np.nonzero(idx)[0]))
#            out.append(arr2str(self['dens'][np.nonzero(idx)[0]]))
#
#            out.append(' == temperature mask ==')
#            out.append(arr2str(np.nonzero(idx)[1]))
#            out.append(arr2str(self['temp'][np.nonzero(idx)[1]]))
#
#        out.append(' ==     var     ==')
#        if var.ndim == idx.ndim:
#            out.append(arr2str(var[idx]))
#        elif var.ndim == 3 and idx.ndim == 2:
#            # probably something to do with the "ionfrac sums to 1" test
#            out.append(arr2str(np.sum(var, axis=-1)[idx]))
#
#
#        out.append('-'*80)
#        return '\n'.join(out)
#





