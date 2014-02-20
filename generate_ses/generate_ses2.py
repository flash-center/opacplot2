# -*- coding: utf-8 -*-
import sys
import os.path
import copy

import opacplot2 as opp
from opacplot2.smooth import CheckEosConsistency
import numpy as np
import hedp
from hedp.eos import thomas_fermi_ionization

#table_id = 7593
#table_id = 3720
table_id = 3336
#mat = 'polystyrene'
mat = 'Cu'
rho0 = 8.96
filename = '{0}-ses-{1:4d}-v2.cn4'.format(mat.lower(), int(table_id))

mat_dict = hedp.matdb(mat)


SESAME_DIR = '../../sesame/'
print 'Loading SESAME db...'
eos_sesame = opp.OpgSesame(os.path.join(SESAME_DIR, "xsesame_ascii"), opp.OpgSesame.SINGLE, verbose=False)


eos_i  = eos_sesame.data[table_id]
eos_w = copy.deepcopy(eos_i)


# merge Ion and Ele grids
eos_o = opp.adapt.EosMergeGrids(eos_w,
                filter_dens=lambda x: (x>0),
                filter_temps=lambda x: (x>0.62), # remove everything bellow 0.62 eV to avoid negative values
                thresh=[])
                

dens_arr, temp_arr = np.meshgrid(eos_o['ele_dens'], eos_o['ele_temps'] )

zbar_tf = thomas_fermi_ionization(dens_arr, temp_arr, mat_dict.Z, mat_dict.A).T

CheckEosConsistency(eos_o)

# Lookig what is happenning at room temperature and solid density
rho0_idx = np.argmin(np.abs(eos_o['ele_dens']-0.92))
rho0_nearest = eos_o['ele_dens'][rho0_idx]
temps0_idx = np.argmin(np.abs(eos_o['ele_temps']-0.025))
temps0_nearest = eos_o['ele_temps'][temps0_idx]
print 'Solid density: {0:.2f} g.cm⁻³ nearest grid point {1:.2f} g.cm⁻³'.format(rho0, rho0_nearest)
print 'Room temperature: {0:.1f} meV nearest grid point {1:.1f} meV'.format(25, temps0_nearest*1e3)
print ' '*10, 'ele      ioncc'
print 'Pressure [Bar]:  {0:.3e}    {1:.3e}'.format(
                      *[eos_o[el+'_pres'][rho0_idx, temps0_idx]*1e-6 for el in ['ele', 'ioncc']])
print 'Eint [erg]    :  {0:.3e}    {1:.3e}'.format(
                      *[eos_o[el+'_eint'][rho0_idx, temps0_idx] for el in ['ele', 'ioncc']])
print 'Zbar    :  {0:.3f}'.format(zbar_tf[rho0_idx, temps0_idx])
print 'A: {0} / {1}    Z: {2} / {3}'.format(mat_dict.A, eos_o['abar'], mat_dict.Z, eos_w['zmax'])
print u'dens[0]: {0} g.cm⁻³, temps[0]: {1} eV'.format(eos_o['ele_dens'][0], eos_o['ele_temps'][0])
print eos_o['ioncc_eint'][0,0]

# writing everything to files
numDens = opp.NA * eos_o['ele_dens'] / eos_o['abar']
#print mat_dict.snop.z, mat_dict.snop.fraction
opp.writeIonmixFile(filename, mat_dict.snop.z, mat_dict.snop.fraction, 
                        numDens=numDens, temps=eos_o['ele_temps'],
                        eion=eos_o["ioncc_eint"],
                        eele=eos_o["ele_eint"],
                        pion=eos_o["ioncc_pres"],
                        pele=eos_o["ele_pres"],
                        zbar=zbar_tf)


# rereading the output
ionmix = opp.OpacIonmix(filename, eos_o['abar']/opp.NA, twot=True, man=True, verbose=False)

test_dict = {'temps': 'ele_temps', 'dens': 'ele_dens',
             'pele': 'ele_pres', 'pion': 'ioncc_pres',
             'eele': 'ele_eint', 'eion': 'ioncc_eint'}
print 'Checking: ',
for attr, key in test_dict.iteritems():
     print attr,'...',
     np.testing.assert_allclose(getattr(ionmix, attr), eos_o[key], atol=1e-5, rtol=1e-5)
np.testing.assert_allclose(ionmix.zbar, zbar_tf, atol=1e-5, rtol=1e-5)


print '\nFile successfuly written!'
print ionmix.zbar

# <codecell>


