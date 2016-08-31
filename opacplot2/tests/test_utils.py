import sys
import opacplot2 as opp
import os.path
import unittest

class test_interpDT(unittest.TestCase):
    BASE_DIR = os.path.join(os.path.dirname(__file__), 'data')
    reference_name = os.path.join(BASE_DIR, 'imx_sample.cn4')
    mpi = 1.00794
    man=True
    twot=True
    
    def setUp(self):
        # Ionmix file info.
        self.imx = opp.OpacIonmix(
            self.reference_name, self.mpi, 
            man=self.man, twot=self.twot)
        self.pion = self.imx.pion
        self.dens = self.imx.dens
        self.temps = self.imx.temps
        self.ndens = self.imx.ndens
        self.ntemp = self.imx.ntemp
        
    def test_interp_extrap_at_zero(self):
        interp = opp.utils.interpDT(
            self.pion, self.dens, self.temps, 0, 0, 
            bcdmin=opp.BC_EXTRAP_ZERO, bctmin=opp.BC_EXTRAP_ZERO)
        self.assertTrue(interp==0,
                        msg='Checking if interpDT adds a point at 0!')
                        
    def test_interp_extrap_lower_bound(self):
        interp = opp.utils.interpDT(
            self.pion, self.dens, self.temps, 0, 0,
            bcdmin=opp.BC_BOUND, bctmin=opp.BC_BOUND)
        self.assertTrue(interp <= self.pion[0][0] and interp >= 0,
                        msg='Checking if interpDT interpolates '
                            'between 0 and the lower bounds!')
                            
    def test_interp_extrap_upper_bound(self):
        interp = opp.utils.interpDT(
            self.pion, self.dens, self.temps, 
            self.dens[self.ndens-1]+1, self.temps[self.ntemp-1],
            bcdmax=opp.BC_BOUND, bctmax=opp.BC_BOUND)
        self.assertTrue(interp == self.pion[self.ndens-1][self.ntemp-1],
                        msg='Checking if interpDT interpolates '
                            'correctly at the upper bounds!')

class test_EosMergeGrids(unittest.TestCase):
    def setUp(self):    
        self.ses_file = os.path.join(os.path.dirname(__file__), 'data/matr_009999.ses')
        self.op = opp.OpgSesame(
            self.ses_file, 
            opp.OpgSesame.SINGLE)
        self.eos_data = self.op.data[9999]
        # Setting intersect just in case the defaults change & testsuite
        # doesn't.
        self.md = opp.utils.EosMergeGrids(
                    self.eos_data,
                    filter_temps=lambda x: x > 1,
                    filter_dens=lambda x: x>1,
                    intersect=['ele', 'ioncc'])
    
    def test_docstring(self):
        self.assertTrue(self.md.__doc__,
                        msg='Checking docstring!')
            
    def test_filter_temps(self):
        self.assertTrue(
            (self.md['ele_dens'] > 1).sum() \
            == self.md['ele_dens'].size,
            msg='Check if our dens filter is actually working!')
    
    def check_filter_dens(self):
        self.assertTrue(
            (self.md['ele_temps'] > 1).sum() \
            == self.md['ele_temps'].size,
            msg='Check if our temp filter is actually working!')
        
    
    