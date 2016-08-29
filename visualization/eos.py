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
from opacplot2 import utils
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
table_id = 3717
eos_data_i  = eos_sesame.data[table_id]

# list of all zbar tables for this element
cond_keys = sorted([key for key, val in cond_sesame.data.iteritems() if val['zmax'] == eos_data_i['zmax']])

# by default select the last (i.e. newest) table avalable
cond_data_i = cond_sesame.data[cond_keys[-1]]

eos_data = copy.deepcopy(eos_data_i)

# uncoment this if you want to merge ele and ion grids
eos_data = opp.adapt.EosMergeGrids(eos_data,
                filter_dens=lambda x: (x>0),
                filter_temps=lambda x: (x>0),
                thresh=['ioncc_eint'])
eos_data.update(cond_data_i)

# <codecell>

# Looking at what is happenning at room temperature and solid density


#eos_unified = opp.adapt.EosMergeGrids(eos_data)
eos_unified = eos_data
rho0 = periodictable.elements[eos_unified['zmax']].density
rho0_idx = np.argmin(np.abs(eos_unified['ele_dens']-rho0))
rho0_nearest = eos_unified['ele_dens'][rho0_idx]
temps0_idx = np.argmin(np.abs(eos_unified['ele_temps']-0.025))
temps0_nearest = eos_unified['ele_temps'][temps0_idx]
print 'Solid density: {0:.2f} g.cm⁻³ nearest grid point {1:.2f} g.cm⁻³'.format(rho0, rho0_nearest)
print 'Room temperature: {0:.1f} meV nearest grid point {1:.1f} meV'.format(25, temps0_nearest*1e3)
print 'Pressure: Ele:  {0:.3e}  Ion: {1:.3e}  Ioncc: {2:.3e}'.format(
                      *[eos_unified[el+'_pres'][rho0_idx, temps0_idx]*1e-10 for el in ['ele', 'ion', 'ioncc']])

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

idx = np.nonzero(eos_data['ioncc_eint']<0)
print 'dens:', idx[0], eos_data['ioncc_dens'][idx[0]]
print 'temps:', idx[1], eos_data['ioncc_temps'][idx[1]]
print 'Pressure:', eos_data['ioncc_pres'][idx]

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

# Ion EoS - Internal Energy
# -------------------------

# <codecell>

plot_eos_field(eos_data, 'ioncc', 'eint')

# <codecell>

plot_eos_field(eos_data, 'ioncc', 'eint', grad=True)

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
        zbar[jd,jt] = opp.utils.interpDT(cond_data_i["zbar"],
                                       cond_data_i["zbar_dens"],
                                       cond_data_i["zbar_temps"],
                                       denss[jd], temps[jt],
                                       bctmin=opp.BC_EXTRAP_ZERO)
            

# <codecell>

import pysnop

args= dict(
    rho_a= denss,
    t_a= temps,
    dsm=2.7,
    csm=3.5e5,
    dhvapm=3e10,
    tboil=2.5e3,
    ideg=1,
    z=[13],
    drek=[3.0],
    fact=[100],
    nsigma=[0],
    )
    
zbar2 = pysnop.run(**args)['zbar']
#opp.plot_zbar(denss, temps, zbar, cond_data_i['zmax'], plt.figure())
#opp.plot_zbar(denss, temps, zbar-zbar2, cond_data_i['zmax'], plt.figure())

# <markdowncell>

# Ion EoS - Pressure
# ------------------

# <codecell>

# write ionmix file with 

#RHO, TE = np.meshgrid(eos_data.origin['ioncc_dens'], eos_data.origin['ioncc_temps'])
#mask  = 1e5*(np.abs(np.log10(RHO.T) - 1)**(-30))*np.exp(-np.abs(np.log(TE.T))**2)
#eos_data.origin["ioncc_pres"] = np.where(eos_data_i["ioncc_pres"] <= 0,
#                        mask,
#                        eos_data.origin["ioncc_pres"])

mask = eos_data_i["ioncc_pres"]>0

eos_data.origin['ioncc_pres'][~mask] =  eos_data_i['ioncc_pres'][mask].min()
#eos_data.origin['ioncc_pres'] = eos_data_i['ioncc_pres'] - eos_data_i['ioncc_pres'].min()

#print eos_data['ioncc_pres'].min()

# <codecell>

ax = plt.subplot(111)
idx = 4
semilogx(eos_data['ioncc_dens'], eos_data['ioncc_pres'][:,idx]*1e-9,'.-')
ax.set_yscale('symlog', linthreshy=1e-6)
ax.set_xlim([1e-6,100])
ax.set_ylim([-3e11,1e12])
ax.grid(which='major')

# <codecell>

plot_eos_field(eos_data, 'ioncc', 'pres')

# <codecell>

plot_eos_field(eos_data, 'ioncc', 'pres', grad='rho')

# <codecell>

fracs = (1.0,)
filename = 'al-ses-3717.cn4'
numDens = opp.NA * eos_data['ele_dens'] / eos_data["abar"]
opp.writeIonmixFile(filename, (eos_data['zmax'],), fracs,  
                        numDens=numDens, temps=eos_data['ele_temps'],
                        eion=np.fmax(eos_data["ioncc_eint"],1e-10),
                        eele=eos_data["ele_eint"],
                        pion=eos_data["ioncc_pres"],
                        pele=eos_data["ele_pres"],
                        zbar=np.fmax(zbar2,1e-3))


ionmix = opp.OpacIonmix(filename, eos_data["abar"]/opp.NA, twot=True, man=True, verbose=False)

# <codecell>



