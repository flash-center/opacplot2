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
#from opacplot2.opg_hdf5 import OpgHdf5



BASE_DIR = os.path.join(os.path.dirname(__file__), 'data')
reference_file =  os.path.join(BASE_DIR, 'imx_sample.cn4')


fields = ['numDens', 'temps',
            'eion', 'eele', 'pion', 'pele', 'zbar',
            'ngroups', 'opac_bounds',
            'rosseland', 'planck_absorb', 'planck_emiss']


def test_ionmix_read():
    abar = 1.0
    eos_data = opp.OpacIonmix(reference_file, 
            abar/opp.NA,
            twot=True, man=True, verbose=False)


def test_extend_to_zero():
    abar = 1.0
    eos_data = opp.OpacIonmix(reference_file, 
            abar/opp.NA,
            twot=True, man=True, verbose=False)

    eos_data.extendToZero()


def test_ionmix_read_write():
    # Verify the writeIonmixFile function
    abar = 1.0
    zmax = 1.0
    fh = opp.OpacIonmix(reference_file, 
            abar/opp.NA,
            twot=True, man=True, verbose=False)

    fracs = (1.0,)
    tmp_file =  os.path.join(BASE_DIR, 'imx_sample_tmp.cn4')
    try:
        pars = {key: getattr(fh, key) for key in fields}
        opp.writeIonmixFile(tmp_file, (zmax,), fracs, **pars)
        # Open generated file to check there are no errors.
        fh_new = opp.OpacIonmix(tmp_file, abar/opp.NA, twot=True, man=True, verbose=False)

        for key in fields:
            yield assert_allclose, getattr(fh, key), getattr(fh_new, key)

    except:
        raise
    finally:
        if os.path.exists(tmp_file):
            os.remove(tmp_file) 


def test_ionmix_read_write_0():
    # Verify the OpacIonmix.write function
    abar = 1.0
    zmax = 1.0
    fh = opp.OpacIonmix(reference_file, 
            abar/opp.NA,
            twot=True, man=True, verbose=False)

    fracs = (1.0,)
    tmp_file =  os.path.join(BASE_DIR, 'imx_sample_tmp_0.cn4')

    try:
        fh.write(tmp_file, (zmax,), fracs, twot=True, man=True)
        # Open generated file to check there are no errors.
        fh_new = opp.OpacIonmix(tmp_file, abar/opp.NA, twot=True, man=True, verbose=False)

        for key in fields:
            yield assert_allclose, getattr(fh, key), getattr(fh_new, key)

    except:
        raise
    finally:
        if os.path.exists(tmp_file):
            os.remove(tmp_file) 

