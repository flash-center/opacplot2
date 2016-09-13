import sys
import opacplot2 as opp
import opacplot2.utils
import os.path
import unittest
import numpy as np

class test_randomize_ionmix(unittest.TestCase):
    BASE_DIR = os.path.join(os.path.dirname(__file__), 'data')
    reference_file =  os.path.join(BASE_DIR, 'imx_sample.cn4')
    tmp_file =  os.path.join(BASE_DIR, 'imx_randomize_tmp.cn4')
    
    fields = ['numDens', 'temps',
              'eion', 'eele', 'pion', 'pele', 'zbar',
              'ngroups', 'opac_bounds',
              'rosseland', 'planck_absorb', 'planck_emiss']
    
    abar = 1.0
    zmax = 1.0
    fracs = (1.0,)
    
    def setUp(self):
        self.eos_data = opp.OpacIonmix(
                            self.reference_file, 
                            self.abar/opp.NA,
                            twot=True, man=True, verbose=False)
    
    def test_randomize_ionmix(self):
        try:
            opp.utils.randomize_ionmix(
                    self.reference_file,
                    self.tmp_file)
            
            self.random_eos_data = opp.OpacIonmix(
                                        self.tmp_file,
                                        self.abar/opp.NA,
                                        twot=True,
                                        man=True,
                                        verbose=False)
            for attr in self.fields:
                if attr not in ['ngroups']:
                    self.assertTrue(
                            np.any(np.not_equal(
                                        getattr(self.random_eos_data, attr),
                                        getattr(self.eos_data, attr))),
                            msg='Checking if {0} data was'
                                'randomized!'.format(attr))
            
        except:
            raise
        finally:
            if os.path.exists(self.tmp_file):
                os.remove(self.tmp_file)

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
            self.pion, self.dens, self.temps,
            bcdmin=opp.BC_EXTRAP_ZERO, bctmin=opp.BC_EXTRAP_ZERO)
        self.assertTrue(interp(0,0)==0,
                        msg='Checking if interpDT adds a point at 0!')
                        
    def test_interp_extrap_lower_bound(self):
        interp = opp.utils.interpDT(
            self.pion, self.dens, self.temps,
            bcdmin=opp.BC_BOUND, bctmin=opp.BC_BOUND)
        self.assertTrue(interp(0,0) <= self.pion[0][0] and interp(0,0) >= 0,
                        msg='Checking if interpDT interpolates '
                            'between 0 and the lower bounds!')
                            
    def test_interp_extrap_upper_bound(self):
        interp = opp.utils.interpDT(
            self.pion, self.dens, self.temps)
        self.assertTrue(interp(self.dens[self.ndens-1]+1, self.temps[self.ntemp-1])
                               == self.pion[self.ndens-1][self.ntemp-1],
                               msg='Checking if interpDT interpolates '
                                   'correctly at the upper bounds!')
    
    def test_interp_lookup_INTERP_DFDD(self):
        interp = opp.utils.interpDT(
            self.pion, self.dens, self.temps, 
            bcdmin=opp.BC_EXTRAP_ZERO, bctmin=opp.BC_EXTRAP_ZERO,
            lookup=opp.INTERP_DFDD)
        
        self.assertTrue(interp(1,1) >= 0,
                        msg='Checking that interpDT finds positive'
                            'DFDD for ion pressure!')
            
    def test_interp_lookup_INTERP_DFDT(self):
        interp = opp.utils.interpDT(
            self.pion, self.dens, self.temps, 
            bcdmin=opp.BC_EXTRAP_ZERO, bctmin=opp.BC_EXTRAP_ZERO,
            lookup=opp.INTERP_DFDT)
        
        self.assertTrue(interp(1,1) >= 0,
                        msg='Checking that interpDT finds positive'
                            'DFDT for ion pressure!')
        

class test_EosMergeGrids(unittest.TestCase):
    def setUp(self):    
        self.ses_file = os.path.join(os.path.dirname(__file__), 
                                     'data/matr_009999.ses')
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
        
    
    