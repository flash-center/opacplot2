#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os.path
import unittest
import numpy as np
import opacplot2 as opp


class test_hdf5(unittest.TestCase):
    BASE_DIR = os.path.join(os.path.dirname(__file__), 'data')
    reference_file =  os.path.join(BASE_DIR, 'Al_snp_40gr.h5')
    
    def setUp(self):
        self.fh = opp.OpgHdf5.open_file(self.reference_file)
        
    def test_hdf5_read_consistency(self):
        self.assertTrue((self.fh['Zf_DT'][:]>0).sum() \
                        == self.fh['Zf_DT'][:].size, 
                        msg='Checking that Zf_DT>0!')

        self.assertTrue((self.fh['Zf_DT'][:] <= self.fh['Zmax']).sum() \
                        == self.fh['Zf_DT'][:].size, 
                        msg='Checking that Zf_DT<Zmax!')

        self.assertTrue((self.fh['Znum'][:] <= self.fh['Zmax']).sum() \
                        == self.fh['Znum'][:].size,
                        msg='Checking that Znum <= Zmax!')
                        
    def test_hdf5_write_consistency(self):
        tmp_file =  os.path.join(self.BASE_DIR, 'Al_snp_40gr_tmp.h5')
        
        try:
            self.fh.write2file(tmp_file)
            self.fh2 = opp.OpgHdf5.open_file(tmp_file)
            
            self.assertTrue(self.fh.keys() == self.fh2.keys(),
                            msg='Checking that the read/written files '
                                'have the same keys!')
            
            # If the keys are different, we can still check that the
            # common keys lead to equal values.
            keys=[key for key in self.fh.keys() if key in self.fh2.keys()]
            
            for key in keys:
                try:
                    np.testing.assert_array_equal(
                        self.fh2[key][:],
                        self.fh[key][:], 
                        err_msg='Checking that {0} for the read/written '
                                'files is the same!'.format(key))
                except IndexError:
                    # Scalars will throw an IndexError.
                    self.assertTrue(
                        self.fh2[key] == self.fh[key], 
                        msg='Checking that {0} for the read/written '
                                'files is the same!'.format(key))
                except TypeError:
                    # Zfo_DT might throw 'TypeError: 'NoneType' object is not
                    # subscriptable'.
                    if self.fh[key] is None and self.fh2[key] is None:
                        self.assertTrue(self.fh[key] == self.fh2[key],
                                        msg='Zfo_DT was set to None!')
                            
        except:
            raise
        finally:
            self.fh.f.close()
            self.fh2.f.close()
            if os.path.exists(tmp_file):
                os.remove(tmp_file)











