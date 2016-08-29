import re
import random

import math
import numpy as np

from .constants import BC_BOUND, BC_EXTRAP_ZERO
from .constants import INTERP_FUNC, INTERP_DFDD, INTERP_DFDT


def munge_h5filename(filename, h5filename):
    """Gets an opacplot hdf5 filename from existing names.
    
    Parameters
    ----------
    filename : str
        File name.
    h5filename : str
        HDF5 name.
    """
    if not isinstance(h5filename, type(str(encoding='utf-8'))):
        if isinstance(filename, type(str(encoding='utf-8'))): 
            h5filename = filename.rpartition(".")[0] + '.h5'
        else:
            h5filename = "opp.h5"
    return h5filename

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
        
def avgopac(energies_in, opacs_in, trad, ebnds, weight="constant", bound="error"):
    """
    
    Parameters
    ----------
    energies_in : 
        
    opacs_in : 
        
    trad, ebnds : 
        
    weight="constant" : 
        
    bound="error" : 
        
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
                raise ValueError("Energy outside of range %f %f" % (en, energies[0]))
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
        denominator = quad(weight, ebnds[i], ebnds[i+1], limit=500, epsrel=1.0e-06)
        opavg[i] = numerator[0]/denominator[0]

    return opavg

def interpDT(arr, dens, temps, d, t, 
             bcdmin=BC_BOUND, bcdmax=BC_BOUND, 
             bctmin=BC_BOUND, bctmax=BC_BOUND, 
             log = False, lookup=INTERP_FUNC):

    # Do some error checking:
    # TODO

    # Adjust for logarithmic lookup:
    if log == True:
        d     = np.log10(d)
        t     = np.log10(t)
        dens  = np.log10(dens)
        temps = np.log10(temps)
        
    # First, find the temperature/density cell we are in.
    # The opacity will be computed using densities:
    #   dens[jd-1], dens[jd]
    # and temperatures:
    #   temp[jt-1], temp[jt]

    jd = np.searchsorted(dens, d)

    if jd == 0 and bcdmin != BC_EXTRAP_ZERO: 
        d = dens[0]
        jd += 1

    if jd == len(dens): 
        jd = jd - 1
        d = dens[-1]

    d1 = 0.0 if jd == 0 else dens[jd-1]
    d2 = dens[jd]

    jt = np.searchsorted(temps, t)
    if jt == 0 and bctmin != BC_EXTRAP_ZERO: 
        t = temps[0]
        jt += 1

    if jt == len(temps): 
        jt = jt - 1
        t = temps[-1]

    t1 = 0.0 if jt == 0 else temps[jt-1]
    t2 = temps[jt]

    if jt == 0 and log == True:
        t2 = 10**t2
        t = 10**t

    f1 = 0.0 if jt == 0 or jd == 0 else arr[jd-1,jt-1]
    f2 = 0.0 if jt == 0 else arr[jd  ,jt-1]
    f3 = 0.0 if jd == 0 else arr[jd-1,jt  ]
    f4 = arr[jd  ,jt  ]

    if lookup == INTERP_FUNC:

        # Now that the surrounding temperatures/densities have been
        # identified, the interpolation coefficients can be computed.
        # c1 -> weight for dens[jd-1] and temp[jt-1]
        # c2 -> weight for dens[jd]   and temp[jt-1]
        # c3 -> weight for dens[jd-1] and temp[jt]
        # c4 -> weight for dens[jd]   and temp[jt]
    
        delta = (d-d1)/(d2-d1)
        tau   = (t-t1)/(t2-t1)
        
        c1 = (delta-1.0)*(tau-1.0)
        c2 = delta*(1-tau)
        c3 = tau*(1-delta)
        c4 = delta * tau

        # Compute the interpolated opacity:
        return c1 * f1 + c2 * f2 + c3 * f3 + c4 * f4

    if lookup == INTERP_DFDD:
        # Interpolate along d1 and d2:  
        fa = (f3-f1) * (t-t1)/(t2-t1) + f1
        fb = (f4-f2) * (t-t1)/(t2-t1) + f2        

        ans = (fb-fa)/(d2-d1)
        if log == True and jd != 0: ans = ans/(10**d * math.log(10.0))

        return ans

    if lookup == INTERP_DFDT:
        fa = (f2-f1) * (d-d1)/(d2-d1) + f1
        fb = (f4-f3) * (d-d1)/(d2-d1) + f3

        ans = (fb-fa)/(t2-t1)
        if log == True and jt != 0: ans = ans/(10**t * math.log(10.0))

        # print( d,t, fa, fb, t2, t1, ans)

        return ans
    
    raise ValueError("lookup must be INTERP_FUNC, INTERP_DFDD, or INTERP_DFDT")

