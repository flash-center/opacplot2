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

from opacplot2.opg_hdf5 import OpgHdf5



BASE_DIR = os.path.join(os.path.dirname(__file__), 'data')
reference_file =  os.path.join(BASE_DIR, 'Al_snp_40gr.h5')



def test_hdf5_read():
    fh = OpgHdf5.open_file(reference_file)
    fh.run_testsuite(mode='full')
    fh.f.close()


def test_hdf5_write():
    fh = OpgHdf5.open_file(reference_file)
    tmp_file =  os.path.join(BASE_DIR, 'Al_snp_40gr_tmp.h5')

    try:
        fh.write2file(tmp_file)
    except:
        raise
    finally:
        fh.f.close()
        if os.path.exists(tmp_file):
            os.remove(tmp_file) 


