from histogram import histdata
import opl_grid

class OplList:

    def __init__(self, nopacs, getDensTemp, getEnergies, getOpac):
        self.nopacs = nopacs
        self.getDensTemp = getDensTemp
        self.getEnergies = getEnergies
        self.getOpacList = getOpac

    def getOpac(self, n):
        return self.getOpacList(n)

    def findExact(self, rho, temp, rtol = 1.0e-04, ttol = 1.0e-04, hist=False, verbose=False):
        """
        Find the opacity with a given temperature and density
        according to some tolerance.
        """

        for n in xrange(self.nopacs):
            rho_n, temp_n = self.getDensTemp(n)

            if( (abs(rho  - rho_n )/rho_n  <= rtol) and
                (abs(temp - temp_n)/temp_n <= ttol)) :

                en = self.getEnergies(n)
                op = self.getOpacList(n)

                if verbose: print "findExact: %13g  %13g  %13g  %13g" % (rho, rho_n, temp, temp_n)

                if hist == True and len(en)-1 == len(op):
                    return histdata(en,op)
                return en,op

        raise ValueError("Could not find opacity")
