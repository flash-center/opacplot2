from constants import BC_BOUND, BC_EXTRAP, BC_EXTRAP_ZERO
from constants import INTERP_FUNC, INTERP_DFDD, INTERP_DFDT
import numpy as np
import math

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

        # print d,t, fa, fb, t2, t1, ans

        return ans
    
    raise ValueError("lookup must be INTERP_FUNC, INTERP_DFDD, or INTERP_DFDT")
