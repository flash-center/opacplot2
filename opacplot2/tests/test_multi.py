#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os.path
import numpy as np
import opacplot2 as opp
import unittest

class test_multi(unittest.TestCase):
    BASE_DIR = os.path.join(os.path.dirname(__file__), 'data')
    reference_name =  'He_snp'
    reference_znum = [2]
    tmp_name = 'He_snp_tmp'
    tmp_path = os.path.join(BASE_DIR, tmp_name)
    tmp_h5 = os.path.join(BASE_DIR, tmp_name+'.h5')
    
    def setUp(self):
        self.fh = opp.OpgMulti.open_file(self.BASE_DIR, 
                                         self.reference_name, 
                                         verbose=False)
    
    def test_multi_read_consistency(self):
        # Find all of the opacity data keys.
        op_keys = [key for key in self.fh.keys() if '_mg' in key]
        
        # Check that the numbers of dens/temp points corresponds correctly
        # to the shape of the opacity data.
        for key in op_keys:
            self.assertTrue(self.fh['dens'].shape[0] == self.fh[key].shape[0],
                            msg='Checking # of densities consistency '
                                'for {0}!'.format(key))
            self.assertTrue(self.fh['temp'].shape[0] == self.fh[key].shape[1],
                            msg='Checking # of temperatures consistency '
                                'for {0}!'.format(key))
            
    def test_multi_write_consistency(self):
        
        
        try:
            self.fh.write(self.tmp_path)
            
            # Open generated file to check there are no errors.
            self.fh2 = opp.OpgMulti.open_file(
                            self.BASE_DIR,  
                            self.tmp_name, 
                            verbose=False)
            
            # Check if all the keys are equal.
            for key in self.fh2.keys():
                np.testing.assert_array_equal(
                    self.fh[key],
                    self.fh2[key],
                    err_msg='Checking written file\'s consistency for '
                            '{0}!'.format(key))

        except:
            raise
        finally:
            # Clean up the written files.
            for ext in ['opp', 'opr', 'opz', 'eps']:
                real_file = "{prefix}.{ext}.gz".format(prefix=self.tmp_path,
                                                       ext=ext)
                if os.path.exists(real_file):
                    os.remove(real_file)
    
    def test_export_hdf5(self):
        try:
            self.fh.write2hdf(self.tmp_h5, Znum=self.reference_znum)

            test = opp.OpgHdf5.open_file(self.tmp_h5, explicit_load=True)
            
            # Check if the zbar data is close enough with desired rel. diff.
            # of rtol.
            np.testing.assert_allclose(
                self.fh['zbar'],
                test['Zf_DT'][:],
                err_msg='Checking the zbar consistency with HDF5!',
                rtol=1e-05)
            
            test.f.close()
            
        except:
            raise
        finally:
            if os.path.exists(self.tmp_h5):
                os.remove(self.tmp_h5)
    

