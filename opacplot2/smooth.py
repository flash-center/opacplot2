#!/usr/bin/python
from __future__ import absolute_import
#from __future__ import division
from __future__ import print_function

import copy
import numpy as np

class CheckEosConsistency:
    def __init__(self, eos):
        self.fail = 0
        self.num_tests = 0
        self.eos = eos
        self.check_pos()
        self.check_sound_speed()
        self.check_heat_capacity()
        if not self.fail:
            print('Sucess: passed {0}/{0} tests !'.format(self.num_tests))
        else:
            print('Failure: {0}/{1} tests failed!'.format(self.fail, self.num_tests))

    def check_pos(self):
        res = True
        for spec in ['ele', 'ioncc']:
            for tab_name in ['pres', 'eint']:
                tab = self.eos['_'.join([spec, tab_name])]
                self.num_tests += 1
                if np.any(tab < 0):
                    print("{0}_{1} table has negative values".format(spec, tab_name))
                    bad_idx = np.nonzero(tab<0)
                    print('        dens             temp           {0}_{1}'.format(spec, tab_name))
                    print(np.array([self.eos[spec+'_dens'][bad_idx[0]],
                                   self.eos[spec+'_temps'][bad_idx[1]],
                                   self.eos['_'.join([spec, tab_name])][bad_idx]]).T)
                    self.fail +=1 
        return res
                    
    def check_sound_speed(self):
        self._check_deriv('ioncc_pres', "dens")
        self._check_deriv('ele_pres', "dens")
    def check_heat_capacity(self):
        self._check_deriv('ioncc_eint', "temp")
        self._check_deriv('ele_eint', "temp")
    def _check_deriv(self, tab_name, axis):
        tab = self.eos[tab_name]
        ax_idx = {'dens': 0, 'temp': 1}[axis]
        diff = np.diff(tab, axis=ax_idx)
        self.num_tests += 1
        if np.any(diff<0):
            print('d{{{0}}}/d{{{1}}} has negative values'.format(tab_name, axis))
            self.fail +=1 

def psedo_maxwell_loops(dens, temp, pres_in):
    """
    Fix  Van der Waals loops regions.
    This is not proper Maxwell constructions, but it should be close...
    """
    pres = pres_in.copy()
    p_fact = 0.5
    p_max_prev = 0
    small_offset = np.arange(len(dens))*1e-9/len(dens)
    for t_idx in range(len(temp))[::-1]:
        p_i = pres[:,t_idx] # isotherm
        end_loop_idx =  np.nonzero(p_i< 0)[0]
        if not len(end_loop_idx):
            # no Van der Waals loop found
            continue
        end_loop_idx = end_loop_idx.max()
        p_max_idx = p_i[0:end_loop_idx+1].argmax()
        p_max = p_i[p_max_idx]
        if p_max <= 0 and p_max_prev:
            p_max = 0.5*p_max_prev
            begin_loop_idx = 0
        else:
            begin_loop_idx = np.nonzero(p_i[0:end_loop_idx+1]>p_fact*p_max)[0].min()
            
        pres[:,t_idx][begin_loop_idx:end_loop_idx+1] = p_fact*p_max +\
                                    small_offset[begin_loop_idx:end_loop_idx+1]
        # setting to the loop value points at the right of the loop if they are bellow the treshold
        right_mask = p_i[end_loop_idx:]<p_fact*p_max
        pres[end_loop_idx:,t_idx][right_mask] = p_fact*p_max +\
                    small_offset[end_loop_idx:][right_mask]
        
        p_max_prev = p_max
        
    return pres


def insure_monotonicity(dens, temp, table_in, axis='dens'):
    table = copy.deepcopy(table_in)
    if axis == "dens":
        X = dens
        Y = temp
    elif axis== "temp":
        X = temp
        Y = dens
        table = table.T
    print("Assuring monotonicity", end='')
    for y_idx in range(1,len(Y)):

        for x_idx in range(1,len(X)):
            df = table[x_idx, y_idx] - table[x_idx-1, y_idx]
            if df<0.0:
                print('.', end='')
                table[x_idx, y_idx] = table[x_idx-1, y_idx] + 1e-9
    print('')
    if axis == "temp":
        table = table.T
   
    return table

def eint_offset(table):
    if np.any(table<0):
        print("Insuring eint is positive. Adding offset of {0:.3e}".format(-table.min()))
        table = table + np.abs(table.min()) + 1e-9
    return table

def interp_isochores_1d(eos, table='ele', ref_grid='ioncc'):
    #from scipy.interpolate import interp1d
    eosi = copy.deepcopy(eos)
    temp_mask = np.array([temp_el not in eos[table+'_temps']\
                        for temp_el in eos[ref_grid+'_temps']], dtype='bool')

    eos[table+'_temps'] = eos[ref_grid+'_temps']
    eos[table+'_pres'] = np.zeros((len(eos[table+'_dens']), len(eos[table+'_temps'])))

    eos[table+'_eint'] = np.zeros(eos[table+'_pres'].shape)

    for dens_idx  in range(len(eos[table+'_dens'])):
        for par in ['pres', 'eint']:
            # copy values that are on the same grid
            eos[table+'_'+par][dens_idx, ~temp_mask] =  eosi[table+'_'+par][dens_idx, :]

            # interpolate the couple of values that are not there
            #itp = interp1d(eosi[table+'_temps'], eosi[table+'_'+par][dens_idx])
            #eos[table+'_'+par][dens_idx, temp_mask] =  itp(eos[table+'_temps'][temp_mask])
            eos[table+'_'+par][dens_idx, temp_mask] =  np.interp(eos[table+'_temps'][temp_mask],
                                                        eosi[table+'_temps'],
                                                        eosi[table+'_'+par][dens_idx])
    return eos





