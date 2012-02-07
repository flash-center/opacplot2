import numpy as np
import math
from opl_grid import OplGrid

class OpacTabop(OplGrid):
    
    def __init__(self, fn, e0, verbose=False):
        self.fn = fn
        self.verbose = verbose
        self.fhand = open(self.fn, 'r')

        # Read the entire file:
        self.data = ""
        while True:
            line = self.nextline()
            if not line: break
            self.data += line

        self.data = self.data.split()
        
        # Read the table number:
        self.data.pop(0)
        self.table_num = int(self.data.pop(0))
        if self.verbose: print "Table Number =", self.table_num

        # Read zbar:
        self.data.pop(0)
        self.zbar = int(self.data.pop(0))
        if self.verbose: print "zbar =", self.table_num

        # Read abar:
        self.data.pop(0)
        self.abar = float(self.data.pop(0))
        if self.verbose: print "abar =", self.table_num

        # Read temperatures:
        self.data.pop(0)
        n = int(self.data.pop(0))
        self.temps = np.empty(n)
        if self.verbose: print "\nNumber of temperatures =", n

        for i in xrange(n): 
            self.temps[i] = 1000.0 * math.e**float(self.data.pop(0))
            if self.verbose: print "%6i  %13.6e" % (i,self.temps[i])

        # Read densities:
        self.data.pop(0)
        n = int(self.data.pop(0))
        self.dens = np.empty(n)
        if self.verbose: print "\nNumber of Densities =", n

        for i in xrange(n): 
            self.dens[i] = math.e**float(self.data.pop(0))
            if self.verbose: print "%6i  %13.6e" % (i,self.dens[i])

        # Read energies:
        self.data.pop(0)
        n = int(self.data.pop(0))
        self.energies = np.empty(n+1)
        self.energies[0] = e0/1000.0
        if self.verbose: print "\nNumber of Energy Groups =", n

        if self.verbose: print "%6i  %13.6e" % (0,self.energies[0])
        for i in xrange(n): 
            self.energies[i+1] = float(self.data.pop(0))**2/self.energies[i]
            if self.verbose: print "%6i  %13.6e" % (i+1,self.energies[i+1])

        # Convert from keV to eV:
        self.energies *= 1000.0

        # Read opacity:
        self.data.pop(0)
        self.opac = np.empty((len(self.dens),len(self.temps),(len(self.energies)-1)))

        for g in xrange(len(self.energies)-1):
            for jd in xrange(len(self.dens)):
                for jt in xrange(len(self.temps)):
                    self.opac[jd,jt,g] = math.e**float(self.data.pop(0))

        OplGrid.__init__(self,self.dens, self.temps, self.energies, 
                         lambda jd, jt: self.opac[jd,jt,:])


    def nextline(self):
        while True:
            line = self.fhand.readline()
            if line == "": return line

            if not line.isspace() and line.strip()[0] != '*':
                return line

