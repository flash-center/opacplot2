#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
from .tests.suite import testsuite
import os.path
import opacplot2.opg_sesame

class EosMergeGrids(dict):
    """This class provides filtering capabilities for the EoS temperature and
    density grids. 
             
    For instance, SESAME tables may have some additional points 
    in the ion EoS table, compared to the electron EoS table, and as 
    FLASH requires the same density and temperature grid for all species, 
    the simplest solution is to remove those extra points.
    
    Parameters
    ----------
    eos_data : dict 
        Dictionary contraining the EoS data.
    intersect : list 
        The resulting temperature [eV] and density [g/cm⁻³]
        grids will be computed as an intersection of grids of all the
        species given in this list. Default: ['ele', 'ioncc']
    filter_dens : function
        A function that takes a grid of densities
        and returns a mask of points we don't wont to keep.
        Defaut: (lamdba x: x>0.) i.e. don't remove anything.
    filter_temps : function
        A function that takes a grid of temperatures
        and returns a mask of points we don't wont to keep.
        Defaut: (lamdba x: x>0.) i.e. don't remove anything.
    thresh : list 
        Zero threshold on following keys
    
    Returns
    -------
    out : dict 
        A dictionary with the same keys as eos_data. The species specified by
        ``intersect`` will have equal temperature and density grids.
    
    Examples
    --------
    >>> eos_sesame = opp.OpgSesame("../sesame/xsesame_ascii", opp.OpgSesame.SINGLE,verbose=False)
    >>> eos_data  = eos_sesame.data[3720]  # Aluminum
    >>> eos_data_filtered = EosMergeGrids(eos_data,
            intersect=['ele', 'ioncc'],   # Merge ele and ioncc grids
            filter_temps=lamda x: x>1.) # Remove temperatures below 1eV
    """
    def __init__(self, eos_data, filter_dens=lambda x: x>=0.,
                 filter_temps=lambda x: x>=0., intersect=['ele', 'ioncc'],
                 thresh=[]):
        
        user_filter = dict(temps=filter_temps, dens=filter_dens)
        self.origin = eos_data
        self.threshold = thresh
        # Computing the intersection of 'ele' and 'ion' grids.
        i_grids = {}
        for key in ['dens', 'temps']:
            i_grids[key] = eos_data['_'.join((intersect[0],key))]
            for el in intersect[1:]:
                i_grids[key] = np.intersect1d(i_grids[key], eos_data['_'.join((el,key))])

        # Defining indexes we want to keep.
        mask = {}
        for species in ['ele', 'ion', 'total', 'cc', 'ioncc']:
            for var in ['dens', 'temps']:
               key = species + '_' + var
               # mask based on the intersection of 'ele' and 'ion' grids
               mask[key] = np.in1d(eos_data[key], i_grids[var], assume_unique=True)
               # mask from user provided parameters
               mask[key] = mask[key]*user_filter[var](eos_data[key])
        self.mask = mask
        # Initalising dictionary.
        for key in eos_data:
            self[key] = None # The actual values returned by __getitem__().
        return

    def _get_mask(self, key):
        if '_dens' in key or '_temps' in key:
            return self.mask[key]
        elif any([key.endswith(prefix) for prefix in ['pres', 'eint', 'free']]):
            species = key.split('_')[0]
            indexes = np.meshgrid(
                        np.nonzero(self.mask[species+'_dens'])[0],\
                        np.nonzero(self.mask[species + '_temps'])[0])
            return indexes[0].T, indexes[1].T

    def __getitem__(self, key):
        if key in self.origin:
            if '_dens' in key or '_temps' in key:
                mask = self.mask[key]
                return self.origin[key][mask]
            elif '_ndens' in key:
                return len(self[key.replace('_n', '_')])
            elif '_ntemp' in key:
                # there is an incoherence between '_dens' -> '_ndens'
                # and '_temps' -> '_ntemps' that should really be fixed
                return len(self[key.replace('_n', '_')+'s'])
            elif any([key.endswith(word) for word in ['pres', 'eint', 'free']]):
                data =  self.origin[key][self._get_mask(key)]
                if key in self.threshold:
                    return np.fmax(data, 0)
                else:
                    return data
            else:
                return self.origin[key]
        else:
            # Now just in case we have added some extra keys in there,
            # reproduce a normal dict's behaviour
            return dict.__getitem__(self, key)
    
    test_keys=['ele_dens', 'ele_temps']
    
    def run_testsuite(self, mode='short'):
        for attr in dir(self):
            if 'check_' in attr:
                getattr(self, attr)(mode)
        print('')
    
    @testsuite(test_keys, var='ele_dens', mode='short')
    def check_filter_temps(self, mode):
        """Check if our dens filter is actually working!"""
        return self['ele_dens'] > 1
    
    @testsuite(test_keys, var='ele_temps', mode='short')
    def check_filter_dens(self, mode):
        """Check if our temp filter is actually working!"""
        return self['ele_temps'] > 1
            
#    def _repr_arr_error(self, idx, var, msg):
#        """
#        Print an error message when array failed to pass consistency checks.
#        """
#        out = ['-'*80]
#        out.append(' Test failed for {0}/{1} points: {2}'.format(idx.sum(),idx.size, msg))
#        out.append('-'*80)
#        arr2str = np.array2string
#
#        if idx.size >= self['ele_temps'][:].size * self['ele_dens'][:].size:
#            out.append(' == density mask ==')
#            out.append(arr2str(np.nonzero(idx)[0]))
#            out.append(arr2str(self['ele_dens'][np.nonzero(idx)[0]]))
#
#            out.append(' == temperature mask ==')
#            out.append(arr2str(np.nonzero(idx)[1]))
#            out.append(arr2str(self['ele_temps'][np.nonzero(idx)[1]]))
#
#        out.append(' ==     var     ==')
#        if var.ndim == idx.ndim:
#            out.append(arr2str(var[idx]))
#        elif var.ndim == 3 and idx.ndim == 2:
#            # probably something to do with the "ionfrac sums to 1" test
#            out.append(arr2str(np.sum(var, axis=-1)[idx]))
#
#
#        out.append('-'*80)
#        return '\n'.join(out)
#