#!/usr/bin/python
# -*- coding: utf-8 -*-
import tables


class OpgHdf5(dict):
    def __init__(self, filename, explicit_load=False ):
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


    def force_eval(self):
        """
        Load the while table into memory.
        Warning that can 
        """
        for key, val in self.iteritems():
            if type(val) is tables.carray.CArray:
                self[key] = val[:]




