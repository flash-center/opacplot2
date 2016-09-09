#!/usr/bin/python
# -*- coding: utf-8 -*-
#from opacplot2.opg_opal import OpalTable, OpalComposition
import numpy as np
from pystellar.opacity import OpacityTable

def pretty_print_composition(comp):
    return ' '.join(['{0}={1:.3f}'.format(key,val)\
                    for key, val in comp.iteritems()])

def ionmix_grid_bounds(op, logTlim=None):
    """Given bounds on logT compute logR bounds for an OPAL table
    that would be valid for all temperatures in logTlim."""
    if logTlim is None:
        logRlim, logTlim = op.bounds()
    else:
        logRlim, _ = op.bounds()
    logRholim = []
    for idx in range(2):
        logRholim.append(op.invert_points(logR=logRlim,logT=[logTlim[idx]]*2).T[0][::-1])
    return np.array(intersection(*logRholim)[::-1]), np.array(logTlim)

def intersection(A, B):
    """ Compute the intersection of two intervals I0, I1"""

    if max(0, min(A[1], B[1]) - max(A[0], B[0])):
        # intersection is not empty!
        if A[0] < B[0]:
            return [B[0], A[1]]
        else:
            return [A[0], B[1]]
    else:
        print 'Warning: intersection is an empty set!'
        return [np.nan, np.nan]

#===============================================================================
# User defined parameters
#===============================================================================
composition_request = {'X': 0.9, 'Y': 0.02, 'dXc': 0.001}
logT_bounds = [3.8,8]
logRho_bounds = [-11, -3]

#===============================================================================
# Loading OPAL table
#===============================================================================
# 1. requires: pystellar and AstroObject python modules see: https://github.com/alexrudy
# 2. add requiered opal table to pystellar/data/
# 3. edit pystellar/OPAL.yml and add a new entry for the table
# 4. make a symbolic link or copy pystellar/data/ to ./
op = OpacityTable(fkey='Gz020',filename='/home/rth/src/pystellar/OPAL.yml')

op.composition(**composition_request)

composition_match = {key: getattr(op, key) for key in ['X', 'Y', 'Z', 'dXc', 'dXo']}
print 'Requested composition :', pretty_print_composition(composition_request)
print 'Found composition     :', pretty_print_composition(composition_match)

#===============================================================================
# Grid definition
#===============================================================================
print 'Density bounds [g/cm3]:', logRho_bounds
print 'Temperature bounds [K]:', logT_bounds
rho = np.logspace(logRho_bounds[0], logRho_bounds[1], 100) # density grid
temp = np.logspace(logT_bounds[0], logT_bounds[1], 100) # temperature grid 
Nr, Nt = len(rho), len(temp)
rho_arr, temp_arr = np.meshgrid(rho, temp, indices='ij')


#===============================================================================
# Getting Rosseland mean opacities from OPAL
#===============================================================================
opr = np.ones((Nr,Nt))*np.nan
for tidx in range(Nt):
    opr[:,tidx] = op.lookup(logrho=np.log10(rho), logT=np.log10(temp[tidx])*np.ones(Nr))
opr = 10**(opr)
#    # getting density boundaries for a given temperature
#    logRholim, logTlim = ionmix_grid_bounds(op, logTlim=[np.log10(temp[tidx])]*2)
#    rho_mask = (rho>10**logRholim[1])*(rho<10**logRholim[0])
#    rho_mask_len = len(np.nonzero(rho_mask)[0])
#
#    if rho_mask_len:
#        opr[rho_mask,tidx] =  op.kappa(rho=rho[rho_mask], T=temp[tidx]*np.ones(rho_mask_len))



#===============================================================================
# Ploting values of opr
#===============================================================================
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ax = plt.subplot(111)
cs = ax.pcolormesh(rho_arr, temp_arr, opr,
        norm=LogNorm(vmin=np.nanmin(opr), vmax=np.nanmax(opr)))
plt.colorbar(cs)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\rho$ [g.cm$^{-3}$]')
ax.set_ylabel(r'$T_e$ [K]')
ax.set_title('Rosseland mean opacity OPAL: Y=0.080 X=0.700 Z=0.020 dXc=0.100 dXo=0.100')
plt.savefig('opal_opr.png', bbox_inches='tight')

#===============================================================================
# Planck mean opacities
#===============================================================================


#===============================================================================
# Calculate emissivity
#===============================================================================


#===============================================================================
# Write everything to IONMIX format
#===============================================================================
