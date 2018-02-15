#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import os.path
import nose
import numpy as np
import opacplot2 as opp


class test_sesame(unittest.TestCase):
    BASE_DIR = os.path.join(os.path.dirname(__file__), 'data')
    reference_name =  'matr_009999.ses'
    tab_id = 9999
    def setUp(self):
        self.fh = opp.OpgSesame(
                os.path.join(self.BASE_DIR, self.reference_name),
                opp.OpgSesame.SINGLE)
        self.data = self.fh.data[self.tab_id]

    def test_sesame_read_consistency(self):
        # Check if all of the temp points are equal for the major three curves.
        np.testing.assert_array_equal(
            self.data['total_temps'],
            self.data['ioncc_temps'],
            err_msg='Checking that the temperatures are consistent!')

        np.testing.assert_array_equal(
            self.data['total_temps'],
            self.data['ele_temps'],
            err_msg='Checking that the temperatures are consistent!')
