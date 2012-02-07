from opl_list import OplList
import numpy as np

class OplTempGrid(OplList):
    """
    Locator for opacities with a set of fixed temperatures, but
    varying densities at each of those temperatures. Each opacity is
    assumed to have the same energy group structure.
    """

    def __init__(self, dens, temps, energies, getOpac):
        """
        temps -> numpy array of temperatures [eV]

        dens -> list of numpy arrays. len(dens) = len(temps. dens[n]
        stores a sorted list of all of the densities [g/cc]
        corresponding to temperature temps[n]
        
        energies -> [eV] energy group structure
        """
        self.dens = dens
        self.temps = temps
        self.energies = energies
        self.go = getOpac

        def map(n):
            count = 0
            for jt in xrange(len(self.temps)):
                for jd in xrange(len(self.dens[jt])):
                    if count == n:
                        return jd, jt
                    count += 1
            raise ValueError("Bad mapping")

        def getDensTemp(n):
            jd, jt = map(n)
            return self.dens[jt][jd], self.temps[jt]
    
        def getListOpac(n):
            jd, jt = map(n)
            return self.go(jd,jt)

        nopac = 0
        for jt in xrange(len(self.temps)):
            nopac += len(self.dens[jt])

        OplList.__init__(self, nopac,
                         getDensTemp,
                         lambda n: self.energies,
                         getListOpac)


    def getOpac(self, jd, jt):
        return self.energies, self.go(jd,jt)

    def interp(self, rho, temp, log=False):
        dens = self.dens[:]
        temps = self.temps

        if log == True:
            rho   = np.log10(rho)
            temp  = np.log10(temp)
            for jt in xrange(len(temps)) : dens[jt] = np.log10(dens[jt])
            temps = np.log10(temps)

        # First, find the temperatures that straddle temp:
        jt = np.searchsorted(temps, temp)
        if jt == 0: 
            temp = temps[0]
            jt += 1

        if jt == len(temps): 
            jt = jt - 1
            temp = temps[-1]
        # temp is bounded by temps[jt] and temps[jt-1].

        # find the densities that bound rho:
        rho = min(rho, dens[jt-1][-1], dens[jt][-1])
        rho = max(rho, dens[jt-1][ 0], dens[jt][ 0])

        jdm = np.searchsorted(dens[jt-1], rho)
        if jdm == 0: jdm += 1

        jdp = np.searchsorted(dens[jt], rho)
        if jdp == 0: jdp += 1

        # The density is bounded by:
        #   dens[jt-1][jdm-1], dens[jt-1][jdm] and
        #   dens[jt][jdp-1],   dens[jt][jdp]

        # First, do interpolation along temps[jt-1]:
        opa = self.go(jdm-1, jt-1)
        opb = self.go(jdm, jt-1)
        
        rhoa = dens[jt-1][jdm-1]
        rhob = dens[jt-1][jdm]

        opm = (rho-rhoa)*(opb-opa)/(rhob-rhoa) + opa

        # Second, do interpolation along temps[jt]:
        opa = self.go(jdp-1, jt)
        opb = self.go(jdp, jt)
        
        rhoa = dens[jt][jdp-1]
        rhob = dens[jt][jdp]

        opp = (rho-rhoa)*(opb-opa)/(rhob-rhoa) + opa
        
        # Interpolate in the temperature direction to find the solution:
        tempa = temps[jt-1]
        tempb = temps[jt]

        return (temp-tempa)*(opp-opm)/(tempb-tempa) + opm

    def __str__(self):
        string = ""
        for jt in xrange(len(self.temps)):
            string += "%13.6e  |  " % self.temps[jt]
            for jd in xrange(len(self.dens[jt])):
                string += "%15.6e" % self.dens[jt][jd]

            string += "\n"
        return string
            
