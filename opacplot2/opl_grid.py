from opl_list import OplList
import numpy as np

class OplGrid(OplList):
    """
    Locator for opacities with a temperature/density grid structure.
    """

    def __init__(self, dens, temps, energies, getOpac):
        self.dens = dens
        self.temps = temps
        self.energies = energies
        self.argGetOpac = getOpac

        def getDensTemp(n):
            jt = n % len(self.temps)
            jd = n / len(self.temps)
            return self.dens[jd], self.temps[jt]
    
        def getListOpac(n):
            jt = n % len(self.temps)
            jd = n / len(self.temps)
            return self.go(jd,jt)

        OplList.__init__(self,
                         len(self.dens)*len(self.temps), 
                         getDensTemp,
                         lambda n: self.energies,
                         getListOpac)

    def go(self, jd, jt):
        return self.argGetOpac(jd,jt)

    def getOpac(self, jd, jt):
        return self.energies, self.go(jd,jt)

    def interp(self, rho, temp, log=False):
        dens = self.dens
        temps = self.temps

        if log == True:
            rho   = np.log10(rho)
            temp  = np.log10(temp)
            dens  = np.log10(dens)
            temps = np.log10(temps)

        # First, find the temperature/density cell we are in.
        # The opacity will be computed using densities:
        #   dens[jd-1], dens[jd]
        # and temperatures:
        #   temp[jt-1], temp[jt]

        jd = np.searchsorted(dens, rho)
        if jd == 0: 
            rho = dens[0]
            jd += 1

        if jd == len(dens): 
            jd = jd - 1
            rho = dens[-1]

        jt = np.searchsorted(temps, temp)
        if jt == 0: 
            temp = temps[0]
            jt += 1

        if jt == len(temps): 
            jt = jt - 1
            temp = temps[-1]

        # Now that the surrounding temperatures/densities have been
        # identified, the interpolation coefficients can be computed.
        # c1 -> weight for dens[jd-1] and temp[jt-1]
        # c2 -> weight for dens[jd]   and temp[jt-1]
        # c3 -> weight for dens[jd-1] and temp[jt]
        # c4 -> weight for dens[jd]   and temp[jt]

        d1 = dens[jd-1]
        d2 = dens[jd]
        t1 = temps[jt-1]
        t2 = temps[jt]
                
        delta = (rho-d1)/(d2-d1)
        tau   = (temp-t1)/(t2-t1)

        c1 = (delta-1.0)*(tau-1.0)
        c2 = delta*(1-tau)
        c3 = tau*(1-delta)
        c4 = delta * tau

        # Compute the interpolated opacity:
        return \
            c1 * self.go(jd-1,jt-1) + \
            c2 * self.go(jd  ,jt-1) + \
            c3 * self.go(jd-1,jt  ) + \
            c4 * self.go(jd  ,jt  )

