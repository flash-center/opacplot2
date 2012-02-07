from math import exp
import numpy as np
from scipy.integrate import quad

def avgopac(energies_in, opacs_in, trad, ebnds, weight="constant", bound="error"):
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

    for i in xrange(len(opavg)):
        numerator = quad(f, ebnds[i], ebnds[i+1], limit=500, epsrel=0.001)
        denominator = quad(weight, ebnds[i], ebnds[i+1], limit=500, epsrel=1.0e-06)
        opavg[i] = numerator[0]/denominator[0]

    return opavg
