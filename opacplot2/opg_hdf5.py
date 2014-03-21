#!/usr/bin/python
# -*- coding: utf-8 -*-
import tables
import numexpr as ne
import numpy as np


class OpgHdf5(dict):
    def __init__(self, filename, explicit_load=False, check=True):
        """
        Open an HDF5 file containing opacity data

        Parameters:
        -----------
        filename: str
                  path to the hdf5 file
        explicit_load: bool
                  load the whole file to memory
        """
        self.f = f = tables.open_file(filename, 'r')
        for el in f.root:
            self[el.name] = el
        for key in f.root._v_attrs._f_list():
            self[key] = f.root._v_attrs[key]

        if explicit_load: self.force_eval()
        if check: self.check_consistency()


    def force_eval(self):
        """
        Load the while table into memory.
        Warning that can 
        """
        for key, val in self.iteritems():
            if type(val) is tables.carray.CArray:
                self[key] = val[:]

    def check_consistency(self):
        test_list = [('{0} > 0', 'Zf_DT', '')] 
        for el in test_list:
            self._run_test(*el)



    def _run_test(self, expr, var, err_message):
        arr = self[var]
        mask = ne.evaluate(expr.format('arr'))
        if len(np.nonzero(mask)[0]):
            print self._repr_arr_error(mask, var, err_message)


    def _repr_arr_error(self, idx, var, msg):
        """
        Print an error message when array failed to pass consistency checks.
        """
        out = ['='*80]
        out.append(' Test failed for {0} :{1}'.format(var, msg))
        out.append('='*80)
        arr2str = np.array2string
        out.append(' == density mask ==')
        out.append(arr2str(idx[0]))
        out.append(arr2str(self.dens[idx[0]]))

        out.append(' == temperature mask ==')
        out.append(arr2str(idx[1]))
        out.append(arr2str(self.temp[idx[1]]))

        out.append(' ==     {0}      =='.format(var))
        out.append(arr2str(self[var][idx]))

        out.append('='*80)
        return '\n'.join(out)
