#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os.path

import nose
from numpy.testing import assert_allclose
import numpy as np

import opacplot2 as opp
import unittest

class test_imx(unittest.TestCase):
    BASE_DIR = os.path.join(os.path.dirname(__file__), 'data')
    reference_file =  os.path.join(BASE_DIR, 'imx_sample.cn4')
    fields = ['numDens', 'temps',
              'eion', 'eele', 'pion', 'pele', 'zbar',
              'ngroups', 'opac_bounds',
              'rosseland', 'planck_absorb', 'planck_emiss']
    
    tmp_file =  os.path.join(BASE_DIR, 'imx_sample_tmp.cn4')
    tmp_file_0 =  os.path.join(BASE_DIR, 'imx_sample_tmp_0.cn4')
    
    abar = 1.0
    zmax = 1.0
    fracs = (1.0,)
    
    def setUp(self):
        self.eos_data = opp.OpacIonmix(
                            self.reference_file, 
                            self.abar/opp.NA,
                            twot=True, man=True, verbose=False)
    
    def test_ionmix_read(self):
        self.assertTrue(self.eos_data.temps.size == self.eos_data.ntemp,
                        msg='Checking ntemps!')
        self.assertTrue(self.eos_data.dens.size == self.eos_data.ndens,
                        msg='Checking ndens!')
    
    def test_ionmix_extend_to_zero(self):
        eos_data_extended = opp.OpacIonmix(
                                self.reference_file, 
                                self.abar/opp.NA,
                                twot=True, man=True, verbose=False)
        eos_data_extended.extendToZero()
    
    def test_writeIonmixFile(self):
        try:
            pars = {key: getattr(self.eos_data, key) for key in self.fields}
            opp.writeIonmixFile(self.tmp_file, (self.zmax,), 
                                self.fracs, **pars)
            # Open generated file to check there are no errors.
            self.eos_data_new = opp.OpacIonmix(self.tmp_file, 
                                               self.abar/opp.NA, 
                                               twot=True, 
                                               man=True, 
                                               verbose=False)

            for key in self.fields:
                np.testing.assert_allclose(getattr(self.eos_data, key), 
                                           getattr(self.eos_data_new, key))

        except:
            raise
        finally:
            if os.path.exists(self.tmp_file):
                os.remove(self.tmp_file)
        
    def test_ionmix_write(self):
        # Verify the OpacIonmix.write function.
        try:
            self.eos_data.write(self.tmp_file_0, (self.zmax,), self.fracs, twot=True, man=True)
            # Open generated file to check there are no errors.
            self.eos_data_new = opp.OpacIonmix(self.tmp_file_0, self.abar/opp.NA, twot=True, man=True, verbose=False)

            for key in self.fields:
                np.testing.assert_allclose(getattr(self.eos_data, key), 
                                           getattr(self.eos_data_new, key))

        except:
            raise
        finally:
            if os.path.exists(self.tmp_file_0):
                os.remove(self.tmp_file_0) 
    
