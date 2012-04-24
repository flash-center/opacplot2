import numpy as np
from constants import KELVIN_TO_EV, GPA_TO_ERGCC, MJKG_TO_ERGG

class OpgQeos:
    """
    This class is responsible for parsing data from the legacy format
    produced by the QEOS code used by scientists at LULI.
    """        

    def __init__(self, filename, datatype, verbose=False):
        """
        filename: name of the qeos file
        datatype: either "zstar" or "eos"
        """

        self.datatype = datatype
        self.verbose = verbose
        self.fhand = open(filename, "r")
        self.count = 0

        if self.verbose == True:
            print "Loading QEOS " + datatype + " file %s\n" % filename

        self.parse()


    def getnext(self):

        data = self.fhand.read(15)
        self.count += 1
        if self.count == 4:
            self.count = 0
            self.fhand.readline()
        
        return data

    def getblock(self):
        data = np.empty((self.ndens, self.ntemps))
        
        for jt in xrange(self.ntemps):
            for jd in xrange(self.ndens):
                data[jd,jt] = float(self.getnext())

        return data

    def parse(self):

        # Read the table id:
        self.tabid = int(self.getnext()) 

        # Throw away the next number (alway 6 for some reason):
        self.getnext()
        
        # Read the number of mass density points:
        self.ndens = int(float(self.getnext()))
        
        # Read the number of temperature points:
        self.ntemps = int(float(self.getnext()))

        # Read the densities:
        self.denss = np.empty(self.ndens)
        if self.verbose: 
            print "%i Densities [g/cc]:" % self.ndens
        for i in xrange(self.ndens):
            num = self.getnext()
            if num == "           -inf":
                self.denss[i] = 0.0
            else:
                self.denss[i] = float(num)
                if self.datatype == "zstar": self.denss[i] = 10**self.denss[i]
                
            if self.verbose:
                print "%3i  %15.6e" % (i, self.denss[i])
                
        # Read the temperatures:
        self.temps = np.empty(self.ntemps)
        if self.verbose: 
            print "\n%i Temperatures [eV]:" % self.ntemps
        for i in xrange(self.ntemps):
            num = self.getnext()
            if num == "           -inf":
                self.temps[i] = 0.0
            else:
                self.temps[i] = num

                if self.datatype == "zstar": 
                    # I believe that the temperatures in this file are
                    # in log10(kilo-Kelvin). So first convert to
                    # kilo-Kelvin:
                    self.temps[i] = 10**self.temps[i]

                    # Now convert to K:
                    self.temps[i] *= 1000.0

                # Now convert from K to eV:
                self.temps[i] *= KELVIN_TO_EV
                
            if self.verbose:
                print "%3i  %15.6e" % (i, self.temps[i])


        print "\n"
                
        if self.datatype == "zstar":
            self.zbar = 10**self.getblock()
                

        if self.datatype == "eos":
            # I believe pressures are in GPa:
            self.pres = self.getblock() * GPA_TO_ERGCC

            # I believe that energies are in MJ/kg:
            self.eint = self.getblock() * MJKG_TO_ERGG
            self.efree = self.getblock()
