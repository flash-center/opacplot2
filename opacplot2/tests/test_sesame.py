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



BASE_DIR = os.path.join(os.path.dirname(__file__), 'data')
reference_name =  'matr_009999.ses'


def test_opac_read():
    fh = opp.OpgSesame(os.path.join(BASE_DIR, reference_name),
            opp.OpgSesame.SINGLE)

# SEE OpgSesame
#
# def test_OpgSesame():
#     fh = opp.OpgSesame(os.path.join(BASE_DIR, reference_name),
#             opp.OpgSesame.SINGLE)
#     fh.run_testsuite(mode='full')