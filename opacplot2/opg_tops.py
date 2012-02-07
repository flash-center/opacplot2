"""
This module contains classes used to read and process TOPS opacity files.
"""

import numpy as np
from opl_list import OplList
from constants import NA

class OpacTopsSingle():
    def __init__(self, matid, tele, rho, nele, energies, opac1, opac2, opac3):
        """
        This class represents a single opacity from the TOPS file
        """
        self.matid = matid
        self.tele  = tele
        self.rho   = rho
        self.nele  = nele

        self.energies = energies
        self.opac1    = opac1
        self.opac2    = opac2
        self.opac3    = opac3

class OpacTops():
    def __init__(self, fn, verbose=False, zmax=None, abar=None):
        """
        This initializer is used to read an existing TOPS file

        The tops format contains three opacity for each
        temperature/density/frequency point. These opacities are:

        1. Total opacity
        2. Opacity without scattering
        3. Scattering only
        """

        self.verbose = verbose
        self.fn = fn
        self.fhand = open(self.fn, "r")

        # First, read through the file and determine the number of entries:
        self.nopacs = 0
        for line in self.fhand:
            if line.find("start") == -1: continue
            self.nopacs += 1

        if(self.verbose): print "Total number of opacities =", self.nopacs

        # Rewind the file:
        self.fhand = open(self.fn, "r")

        if(self.verbose):
            print "%6s  %6s  %6s  %13s  %13s  %13s  %13s" % \
                ("# num ", "npts", "matid", "tele (eV)", "rho (g/cc)", "nele (1/cc)", "zbar")

        self.opacs = []
        self.energies = []
        self.opac1 = []
        self.opac2 = []
        self.opac3 = []

        line = self.fhand.readline()
        matid, tele, rho, nele = self.__parseHeader(line)
        for line in self.fhand:
            if line.isspace(): continue

            if line.find("start") != -1:

                # We've read in one complete opacity
                if(self.verbose == True):
                    if zmax == None:
                        print "%6i  %6i  %6i  %13.6e  %13.6e  %13.6e" % \
                            (len(self.opacs), 
                             len(self.energies),
                             matid, tele, rho, nele)
                    else:
                        nion = NA * rho / abar
                        print "%6i  %6i  %6i  %13.6e  %13.6e  %13.6e  %13g" % \
                            (len(self.opacs), 
                             len(self.energies),
                             matid, tele, rho, nele, nele/nion )
                        
                
                self.opacs.append(OpacTopsSingle(matid, tele, rho, nele, 
                                                 np.array(self.energies), 
                                                 np.array(self.opac1), 
                                                 np.array(self.opac2), 
                                                 np.array(self.opac3)))

                self.energies = []
                self.opac1 = []
                self.opac2 = []
                self.opac3 = []
                        
                # Read the next header
                matid, tele, rho, nele = self.__parseHeader(line)
                continue

            e, op1, op2, op3 = line.split()
            self.energies.append(float(e)*1000.0)
            self.opac1.append(float(op1))
            self.opac2.append(float(op2))
            self.opac3.append(float(op3))
            

        if(self.verbose == True):
            if zmax == None:
                print "%6i  %6i  %6i  %13.6e  %13.6e  %13.6e" % \
                    (len(self.opacs), 
                     len(self.energies),
                     matid, tele, rho, nele)
            else:
                nion = NA * rho / abar
                print "%6i  %6i  %6i  %13.6e  %13.6e  %13.6e  %13.6e" % \
                    (len(self.opacs), 
                     len(self.energies),
                     matid, tele, rho, nele, nele/nion)

        self.opacs.append(OpacTopsSingle(matid, tele, rho, nele, 
                                         np.array(self.energies), 
                                         np.array(self.opac1), 
                                         np.array(self.opac2), 
                                         np.array(self.opac3)))

        

    def __parseHeader(self, line):
        """
        Parse the header line for each opacity and extract all of the
        data.
        """
        
        matid = int(line[32:39])
        tele  = float(line[39:51])*1000.0
        rho   = float(line[51:63])
        nele  = float(line[63:75])

        return matid, tele, rho, nele

    def oplList(self):
        return OplList(self.nopacs, 
                       lambda n: (self.opacs[n].rho, self.opacs[n].tele), 
                       lambda n: self.opacs[n].energies,
                       lambda n: self.opacs[n].opac1)

    def oplListNoScat(self):
        return OplList(self.nopacs, 
                       lambda n: (self.opacs[n].rho, self.opacs[n].tele), 
                       lambda n: self.opacs[n].energies,
                       lambda n: self.opacs[n].opac2)

    def oplListScat(self):
        return OplList(self.nopacs, 
                       lambda n: (self.opacs[n].rho, self.opacs[n].tele), 
                       lambda n: self.opacs[n].energies,
                       lambda n: self.opacs[n].opac3)
