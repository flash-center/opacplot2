import re
import random

import math
import numpy as np

from .constants import BC_BOUND, BC_EXTRAP_ZERO
from .constants import INTERP_FUNC, INTERP_DFDD, INTERP_DFDT

import os.path
import opacplot2.opg_sesame

import scipy as sp
import scipy.interpolate
import scipy.misc

def randomize_ionmix(filename, outfilename):
    """Randomizes the data from an existing ionmix file and rewrites it 
    to the outfile.
    
    Parameters
    ----------
    filename : str
        Name of file to randomize.
    outfilename : str
        Name of output file.
    """
    basepttn = "\d{6}"
    exppttn = "E[+-]\d{2}"
    rand_base = lambda mobj: "{0:06}".format(random.randint(0, 999999))
    rand_exp = lambda mobj: "E{0:+03}".format(random.randint(-99, 99))

    with open(filename) as f:
        lines = f.readlines()

    newlines = []
    for line in lines:
        newline = re.sub(basepttn, rand_base, line)
        newline = re.sub(exppttn, rand_exp, newline)
        newlines.append(newline)

    with open(outfilename, 'w') as f:
        f.write("".join(newlines))
        
def avgopac(energies_in, opacs_in, trad, ebnds, 
            weight="constant", bound="error"):
    """
    
    Parameters
    ----------
    energies_in : 
        
    opacs_in : 
        
    trad :
    
    ebnds : 
        
    weight='constant' : 
        
    bound='error' : 
        
    """
    try:
        from scipy.integrate import quad
    except ImportError:
        print('Error: Scipy not installed, cannot caclulate opacity integrals!')
        raise
    except:
        raise
    # Check for errors:

    # Make sure that none of the energy group boundaries is outside of
    # the energy range for this opacity:
    energies = energies_in.copy()
    opacs = opacs_in.copy()

    if bound == "error":
        for en in ebnds:
            if en < energies[0] or en > energies[-1]:
                raise ValueError('Energy outside'
                                 'of range {0} {1}'.format(en, energies[0]))
    elif bound == "continue":
        emin = np.min(ebnds)
        emax = np.max(ebnds)

        energies = np.empty(len(energies)+2)
        energies[0] = emin
        energies[1:-1] = energies_in[:]
        energies[-1] = emax
        
        opacs = np.empty(len(opacs)+2)
        opacs[0] = opacs_in[0]
        opacs[1:-1] = opacs_in[:]
        opacs[-1] = opacs_in[-1]
                    
    else:
        raise ValueError("Illegal boundary treatment")
    
    def op(en):

        idx = energies.searchsorted(en)
        # The energy is between index idx and idx-1
        
        de = energies[idx] - energies[idx-1]
        do = opacs[idx] - opacs[idx-1]

        ans = (en-energies[idx-1]) * do/de + opacs[idx-1]
        return ans
    
    def weight_constant(en):
        return 1.0

    def weight_planck(en):
        x = en/trad
        if x > 700.0: return 1.0
        return en**3/(1-exp(-x)) * exp(-x)

    # def weight_rosseland(en):
        
        

    # Choose the weight function:
    f = lambda en: weight(en)*op(en)
    if weight == "constant":
        weight = weight_constant
    elif weight == "planck":
        weight = weight_planck
    elif weight == "rosseland":
        weight = weight_rossland
        f = lambda en: weight(en)/op(en)

    opavg = np.empty(len(ebnds)-1)

    for i in range(len(opavg)):
        numerator = quad(f, ebnds[i], ebnds[i+1], limit=500, epsrel=0.001)
        denominator = quad(weight, ebnds[i],
                           ebnds[i+1], limit=500, 
                           epsrel=1.0e-06)
        opavg[i] = numerator[0]/denominator[0]

    return opavg

def interpDT(arr, dens, temps,
             bcdmin=BC_BOUND, bctmin=BC_BOUND, 
             lookup=INTERP_FUNC):
    """
    Depending on the choice for lookup, this function returns an interpolation
    function for values in arr, the density derivative of arr, or the
    temperature derivative of arr.
    
    What ``interpDT()`` returns is dependent upon the ``lookup`` and 
    ``bcdmin``/``bctmin`` arguments:
             
        ``INTERP_FUNC`` will return an interpolation function for any
        density/temperature point within the range.

        ``INTERP_DFDD`` will return a function for the density derivative
        at any density/temperature point within the range.

        ``INTERP_DFDT`` will return a function for the temperature derivative
        at any density/temperature point within the range.
             
        ``BC_BOUND`` is the default setting. If an input density/temp is smaller 
        than the minimum value in the ``dens`` (or ``temps``) array, then it
        will automatically be set to this minimum value.

        ``BC_EXTRAP_ZERO`` will insert zero points into the ``arr`` and either
        ``dens`` or ``temps`` for ``bcdmin``, ``bctmin`` respectively.
             
    Parameters
    ----------
    arr : numpy.ndarray
        Data array.
    dens : numpy.ndarray
        Density array.
    temps : numpy.ndarray
        Temperature array.
    bcdmin : int
        Boundary conditions for density.
    bctmin : int
        Boundary conditions for temperature.
    lookup : int
        Type of function to return.
             
    Examples
    --------
    In order to find a function that interpolates between dens/temp points
    for an IONMIX file that we had previously opened as ``imx``,
    we could use::
        
        >>> import opacplot2 as opp
        >>> f = opp.utils.interpDT(imx.pion, imx.dens, imx.temps,
        ...                        bcdmin=BC_EXTRAP_ZERO,
        ...                        bctmin=BC_EXTRAP_ZERO,
        ...                        lookup=INTERP_FUNC)
        >>> print(f(0,0)) # We added points at zero.
        0
        >>> print(f(123, 456)) # Density of 123 and temperature of 456.
        1234.5678 # Resulting ion pressure at this dens/temp.
    """
   
    # Adjust for extrapolation to zero.
    if bcdmin == BC_EXTRAP_ZERO and dens[0] != 0:
        # Density arrays should be 1D.
        dens = np.insert(dens, 0, 0)
        arr = np.insert(arr, 0, 0, axis=0)
    
    if bctmin == BC_EXTRAP_ZERO and temps[0] != 0:
        # Temp arrays should be 1D.
        temps = np.insert(temps, 0, 0)
        arr = np.insert(arr, 0, 0, axis=1)
    
    # "When on a regular grid with x.size = m and y.size = n,
    # if z.ndim ==2, then z must have shape (n, m)."
    # arr will have shape (dens.size, temps.size), so we must transpose it.
    f = sp.interpolate.interp2d(dens, temps, arr.T)
    
    if lookup == INTERP_FUNC:
        def interp_func(func):
            def wrapper(d, t):
                # Deal with data being out of range.
                # Adjust for lower bounds.
                if d < dens[0]:
                    d = dens[0]
                if t < temps[0]:
                    t = temps[0]
                # Adjust for upper bounds.
                if d > dens[-1]:
                    d = dens[-1]
                if t > temps[-1]:
                    t = temps[-1]
                return func(d,t)
            return wrapper
        
        # Return the wrapper of f that takes care of data being out of range.
        return interp_func(f)
    
    if lookup == INTERP_DFDD:
        def df_dd(func):
            # Input to df_dd(func) will be a (dens, temp) point.
            def outter_wrap(dens, temp):
                # Fix temp point but let dens point range.
                def inner_wrap(x):
                    return func(x, temp)
                
                return sp.misc.derivative(
                            inner_wrap, 
                            dens, # Evaluate derivative at dens point.
                            dx=(dens*1e-12))
            
            return outter_wrap
        
        return df_dd(f)
    
    if lookup == INTERP_DFDT:
        def df_dt(func):
            # Input to df_dt(func) will be a (dens, temp) point.
            def outter_wrap(dens, temp):
                # Fix temp point but let temps point range.
                def inner_wrap(x):
                    return func(dens, x)
                    
                return sp.misc.derivative(
                            inner_wrap, 
                            temp, # Evaluate derivative at temp point.
                            dx=(temp*1e-12))
                            
            return outter_wrap
        
        return df_dt(f)
    
    raise ValueError("lookup must be INTERP_FUNC, INTERP_DFDD, or INTERP_DFDT")
    
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
    >>> eos_sesame = opp.OpgSesame("../sesame/xsesame_ascii", 
                                   opp.OpgSesame.SINGLE,verbose=False)
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
                i_grids[key] = np.intersect1d(
                                    i_grids[key], 
                                    eos_data['_'.join((el,key))])

        # Defining indexes we want to keep.
        mask = {}
        for species in ['ele', 'ion', 'total', 'cc', 'ioncc']:
            for var in ['dens', 'temps']:
               key = species + '_' + var
               # mask based on the intersection of 'ele' and 'ion' grids
               mask[key] = np.in1d(eos_data[key], 
                                   i_grids[var], 
                                   assume_unique=True)
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


