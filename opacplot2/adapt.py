#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np

class EosMergeGrids(dict):
    def __init__(self, eos_data, filter_dens=lambda x: x>=0.,
                 filter_temps=lambda x: x>=0., intersect=['ele', 'ioncc'],
                 thresh=[]):
        """This class provides filtering capabilities for the EoS temperature and
        density grids. For instance SESAME tables may have some additionnal points 
        in the ion EoS table, compared to the elecron EoS table, and as FLASH requires
        the same density and temperature grid for all species, the simplest solution is
        to remove those extra points.

        Parameters
        ----------
        - eos_data: [dict] dictionary contraining the EoS data.
        - intersect: [list] the resulting temperature [eV] and density [g/cm⁻³]
                grids will be computed as an intersection of grids of all the
                species given in this list. Default: ['ele', 'ioncc']
        - filter_dens, filter_temps: [function] a function that takes a grid
                and returns a mask of points we don't wont to keep.
                Defaut: (lamdba x: x>0.) i.e. don't remove anything.
        - thresh: zero threshold on folowing keys
        Returns
        -------
        - out: [dict] a dictionary with the same keys a eos_data.

        Exemple
        -------
        >> eos_sesame = opp.OpgSesame("../sesame/xsesame_ascii", opp.OpgSesame.SINGLE,verbose=False)
        >> eos_data  = eos_sesame.data[3720]  # Aluminum
        >> eos_data_filtered = EosMergeGrids(eos_data,
                intersect=['ele', 'ioncc'],   # merge ele and ioncc grids
                filter_temps=lamda x: x>1.) # remove temperatures below 1eV
        """
        user_filter = dict(temps=filter_temps, dens=filter_dens)
        self.origin = eos_data
        self.threshold = thresh
        # computing the intersection of 'ele' and 'ion' grids
        i_grids = {}
        for key in ['dens', 'temps']:
            i_grids[key] = eos_data['_'.join((intersect[0],key))]
            for el in intersect[1:]:
                i_grids[key] = np.intersect1d(i_grids[key], eos_data['_'.join((el,key))])

        # defining indexes we want to keep
        mask = {}
        for species in ['ele', 'ion', 'total', 'cc', 'ioncc']:
            for var in ['dens', 'temps']:
               key = species + '_' + var
               # mask based on the intersection of 'ele' and 'ion' grids
               mask[key] = np.in1d(eos_data[key], i_grids[var], assume_unique=True)
               # mask from user provided parameters
               mask[key] = mask[key]*user_filter[var](eos_data[key])
        self.mask = mask
        # initalising dictionary
        for key in eos_data:
            self[key] = None
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

