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
from opacplot2 import plot_eos_grid, plot_eos_field
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

matplotlib.rcParams.update({'savefig.dpi': 120,
                        'font.size' : 12 })

SESAME_DIR = '../../sesame/'
if 'eos_sesame' not in globals():
    eos_sesame = opp.OpgSesame(os.path.join(SESAME_DIR, "xsesame_ascii"), opp.OpgSesame.SINGLE,verbose=False)
    cond_sesame = opp.OpgSesame(os.path.join(SESAME_DIR, "sescu_ascii"), opp.OpgSesame.DOUBLE, verbose=False)
    
# Selecting the SESAME table we want to look at
eos_data  = eos_sesame.data[3140]
cond_data = cond_sesame.data[23713]

all_data = copy.deepcopy(eos_data)
#all_data = opp.filter.EosFilterGrids(eos_data,
#                         filter_dens=lambda x: (x!=0.81)*(x>0))
all_data.update(cond_data)

# <markdowncell>

# Density and Temperature grids
# ------------------------------

# <codecell>

plot_eos_grid(all_data, 'dens')
plot_eos_grid(all_data, 'temp')

# <codecell>

print 'Table   len(rho)   len(T)'
print '-'*25
for el in all_data:
    if '_dens' in el:
        key = el.split('_')[0]
        print  "{0:6}   {1: 4d}      {2: 4d}".format(key, len(all_data[key+'_dens']), len(all_data[key+'_temps']))

# <markdowncell>

# Electron EoS - Pressure
# -----------------------

# <codecell>

plot_eos_field(all_data, 'ele', 'pres')

# <codecell>

# Searching for negative pressure

idx = np.nonzero(all_data['ele_pres']<0)
print 'dens:', idx[0], all_data['ele_dens'][idx[0]]
print 'temps:', idx[1], all_data['ele_temps'][idx[1]]
print 'Pressure:', all_data['ele_pres'][idx]

# <codecell>

plot_eos_field(all_data, 'ele', 'pres', grad='rho')
plot_eos_field(all_data, 'ele', 'pres', grad='T')

# <markdowncell>

# Electron EoS - Internal energy
# ------------------------------

# <codecell>

plot_eos_field(all_data, 'ele', 'eint')

# <codecell>

plot_eos_field(all_data, 'ele', 'eint', grad='rho')
plot_eos_field(all_data, 'ele', 'eint', grad='T')

# <markdowncell>

# Ion EoS - Pressure
# ------------------

# <codecell>

plot_eos_field(all_data, 'ion', 'pres')

# <codecell>

plot_eos_field(all_data, 'ion', 'pres', grad='rho')

plot_eos_field(all_data, 'ion', 'pres', grad='T')

# <markdowncell>

# Ion EoS - Internal Energy
# -------------------------

# <codecell>

plot_eos_field(all_data, 'ion', 'eint')

# <codecell>

plot_eos_field(all_data, 'ion', 'eint', grad='rho')
plot_eos_field(all_data, 'ion', 'eint', grad='T')

# <codecell>


