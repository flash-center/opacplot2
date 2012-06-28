# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Visualisation routines for tabulated EoS
# ========================================
# 
# This set of routines aims to help detect issues with tabulated EoS before they are translated into a FLASH readable format. 

# <codecell>

import os.path
import copy

import opacplot2 as opp
from opacplot2 import plot_eos_grid, plot_eos_field, plot_zbar
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import periodictable

matplotlib.rcParams.update({'savefig.dpi': 120,
                        'font.size' : 12 })

SESAME_DIR = '../../sesame/'
if 'eos_sesame' not in globals():
    eos_sesame = opp.OpgSesame(os.path.join(SESAME_DIR, "xsesame_ascii"), opp.OpgSesame.SINGLE, verbose=False)
    cond_sesame = opp.OpgSesame(os.path.join(SESAME_DIR, "sescu_ascii"), opp.OpgSesame.DOUBLE, verbose=False)

# <codecell>

# Selecting the SESAME table
#table_id = 3720
#table_id = 2140
table_id = 3336
eos_data_i  = eos_sesame.data[table_id]

# list of all zbar tables for this element
cond_keys = sorted([key for key, val in cond_sesame.data.iteritems() if val['zmax'] == eos_data_i['zmax']])

# by default select the last (i.e. newest) table avalable
cond_data_i = cond_sesame.data[cond_keys[-1]]

eos_data = copy.deepcopy(eos_data_i)

# uncoment this if you want to merge ele and ion grids
#eos_data = opp.adapt.EosMergeGrids(eos_data,
#                filter_dens=lambda x: (x!=0.81)*(x>0))
#                filter_temps=lambda x: (x>0))
eos_data.update(cond_data_i)

# <codecell>

# Lookig what is happenning at room temperature and solid density

eos_unified = opp.adapt.EosMergeGrids(eos_data)
rho0 = periodictable.elements[eos_unified['zmax']].density
rho0_idx = np.argmin(np.abs(eos_unified['ele_dens']-rho0))
rho0_nearest = eos_unified['ele_dens'][rho0_idx]
temps0_idx = np.argmin(np.abs(eos_unified['ele_temps']-0.025))
temps0_nearest = eos_unified['ele_temps'][temps0_idx]
print 'Solid density: {0:.2f} g.cm⁻³ nearest grid point {1:.2f} g.cm⁻³'.format(rho0, rho0_nearest)
print 'Room temperature: {0:.1f} meV nearest grid point {1:.1f} meV'.format(25, temps0_nearest*1e3)
print 'Pressure: Ele:  {0:.3e}  Ion: {1:.3e}  Ioncc: {2:.3e}'.format(
                      *[eos_unified[el+'_pres'][rho0_idx, temps0_idx]*1e-10 for el in ['ele', 'ion', 'ioncc']]),\
       'cc: {0:.3e}'.format(eos_unified['cc_pres'][rho0_idx][0])

# <markdowncell>

# Density and Temperature grids
# ------------------------------

# <codecell>

plot_eos_grid(eos_data, 'dens')
plot_eos_grid(eos_data, 'temp')

# <codecell>

print 'Table   len(rho)   len(T)'
print '-'*25
for el in eos_data:
    if '_dens' in el:
        key = el.split('_')[0]
        print  "{0:6}   {1: 4d}      {2: 4d}".format(key, len(eos_data[key+'_dens']), len(eos_data[key+'_temps']))

# <markdowncell>

# Electron EoS - Pressure
# -----------------------

# <codecell>

plot_eos_field(eos_data, 'ele', 'pres')

# <codecell>

# Searching for negative pressure

idx = np.nonzero(eos_data['ele_pres']<0)
print 'dens:', idx[0], eos_data['ele_dens'][idx[0]]
print 'temps:', idx[1], eos_data['ele_temps'][idx[1]]
print 'Pressure:', eos_data['ele_pres'][idx]

# <codecell>

plot_eos_field(eos_data, 'ele', 'pres', grad=True)

# <markdowncell>

# Electron EoS - Internal energy
# ------------------------------

# <codecell>

plot_eos_field(eos_data, 'ele', 'eint')

# <codecell>

plot_eos_field(eos_data, 'ele', 'eint', grad=True)

# <markdowncell>

# Ion EoS - Pressure
# ------------------

# <codecell>

plot_eos_field(eos_data, 'ioncc', 'pres')
plot_eos_field(eos_data, 'ion', 'pres')

# <codecell>

plot_eos_field(eos_data, 'ion', 'pres', grad=True)

# <markdowncell>

# Ion EoS - Internal Energy
# -------------------------

# <codecell>

plot_eos_field(eos_data, 'ion', 'eint')

# <codecell>

plot_eos_field(eos_data, 'ion', 'eint', grad=True)

# <markdowncell>

# Ionisation table
# ----------------

# <codecell>

opp.plot_zbar(cond_data_i['zbar_dens'], cond_data_i['zbar_temps'], cond_data_i['zbar'], cond_data_i['zmax'], plt.figure())

# <codecell>

ntemp = eos_data["ele_ntemp"]
ndens = eos_data["ele_ndens"]
temps = eos_data["ele_temps"]
denss = eos_data["ele_dens"]

# this part is not working as it is expected to...
zbar = np.empty((ndens, ntemp))
for jd in xrange(ndens):
    for jt in xrange(ntemp):
        zbar[jd,jt] = opp.interpDT(cond_data_i["zbar"],
                                       cond_data_i["zbar_dens"],
                                       cond_data_i["zbar_temps"],
                                       denss[jd], temps[jt],
                                       bctmin=opp.BC_EXTRAP_ZERO)
            
opp.plot_zbar(denss, temps, zbar, cond_data_i['zmax'], plt.figure())

# <codecell>

# write ionmix file with 

fracs = (1.0,)
filename = 'al-ses.cn4'
numDens = opp.NA * eos_data['ele_dens'] / eos_data["abar"]
opp.writeIonmixFile(filename, (eos_data['zmax'],), fracs,  
                        numDens=numDens, temps=eos_data['ele_temps'],
                        eion=eos_data["ioncc_eint"],
                        eele=eos_data["ele_eint"],
                        pion=eos_data["ioncc_pres"],
                        pele=eos_data["ele_pres"])

opp.writeIonmixFile(filename, (eos_data['zmax'],), fracs,  
                        numDens=numDens, temps=eos_data['ele_temps'])

ionmix = opp.OpacIonmix(filename, eos_data["abar"]/opp.NA, twot=True, man=True, verbose=False)

# <codecell>


# <codecell>


