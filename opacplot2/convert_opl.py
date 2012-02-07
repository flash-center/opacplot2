import numpy as np
from opl_grid import OplGrid
from opl_list import OplList
from opl_tempgrid import OplTempGrid
from avgopac import avgopac

def listToGrid(opllist, ndens, ntemps):
    
    # Create an array containing all temperatures and densities:
    rho  = np.empty(opllist.nopacs)
    tele = np.empty(opllist.nopacs)

    for n in xrange(opllist.nopacs):
        rho[n], tele[n] = opllist.getDensTemp(n)
        # print "%15.6e%15.6e" % (rho[n],tele[n])

    rho.sort()
    tele.sort()

    rho = rho[::ntemps]
    tele = tele[::ndens]
    
    # rho and tele should now contain the densities and temperatures
    return OplGrid(rho, tele, opllist.getEnergies(0), 
                   lambda jd, jt: opllist.findExact(rho[jd],tele[jt])[1])


def avgOplList(opllist, ebds, weight="constant", bound="error"):

    getEnergies = lambda n: ebds

    def getOpac(n):
        rho_n, trad_n = opllist.getDensTemp(n)

        return avgopac(opllist.getEnergies(n), 
                       opllist.getOpac(n),
                       trad_n, 
                       ebds, 
                       weight=weight,
                       bound=bound)

    return OplList(opllist.nopacs, opllist.getDensTemp, getEnergies, getOpac)

def listToTempGrid(opllist, ntemps):
    
    # Create an array containing all of the temperatures:
    rhos     = np.empty(opllist.nopacs)
    alltrads = np.empty(opllist.nopacs)
    dt       = np.empty(opllist.nopacs-1)
    trads    = np.empty(ntemps)

    for n in xrange(opllist.nopacs):
        rhos[n], alltrads[n] = opllist.getDensTemp(n)

    alltrads.sort()
    
    for n in xrange(opllist.nopacs-1): 
        dt[n] = alltrads[n+1] - alltrads[n]
    idxs = np.argsort(dt)

    trads[0] = alltrads[-1]
    for n in xrange(1,ntemps):
        i = -ntemps + n
        trads[n] = alltrads[idxs[i]]

    trads.sort()

    # trads now contains a sorted list of the temperatures. For each
    # temperature, make a sorted list of densities.
    rhos = []
    for i in xrange(ntemps):
        current_rho = []
        for n in xrange(opllist.nopacs):
            rho, trad = opllist.getDensTemp(n)
            if abs(trad-trads[i])/trads[i] <= 1.0e-12:
                current_rho.append(rho)
                
        rhos.append(np.array(current_rho))
        rhos[-1].sort()

    # Check to ensure same group structure used:
    energies = opllist.getEnergies(0)
    for n in xrange(opllist.nopacs):
        if len(energies) != len(opllist.getEnergies(n)):
            raise ValueError("Invalid group structure")

    def getOpac(jd, jt):
        en, op = opllist.findExact(rhos[jt][jd], trads[jt], rtol=1.0e-10)
        return op

    return OplTempGrid(rhos, trads, energies, getOpac)
        
