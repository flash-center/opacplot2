from opl_list import OplList
import numpy as np

class OpacInfernoSingle():
    """
    This class represents a single opacity from an INFERNO opacity file
    """
    def __init__(self, znuc, anum, zbar, rho, tele, eta, rmean, pmean,
                 energies, opac_tot, opac_bb, opac_bf, opac_ff, opac_scat):
        
        self.znuc      = znuc
        self.anum      = anum 
        self.zbar      = zbar
        self.rho       = rho 
        self.tele      = tele
        self.eta       = eta
        self.rmean     = rmean
        self.pmean     = pmean
        self.energies  = energies
        self.opac_tot  = opac_tot
        self.opac_bb   = opac_bb
        self.opac_bf   = opac_bf
        self.opac_ff   = opac_ff
        self.opac_scat = opac_scat

    
class OpacInferno:
    """
    Class to read data from INFERNO opacity files
    """
    
    def __init__(self, fn):
        """
        Opens and parses the file, storing all data
        """

        # Open the file:
        self.fn = fn
        self.file = open(self.fn, "r")
        
        # Start looping through the file, one line at a time...
        self.opacs = []
        while self.file.readline():
            self.__readOpac()
        
    def __splitline(self, line):
        words = line.split()

        for i in xrange(len(words)):
            if(words[i][6] == '-'):
                words[i] = words[i][:6] + "E-" + words[i][7:]
               
            if(words[i][6] == '+'):
                words[i] = words[i][:6] + "E+" + words[i][7:]
                
        return [float(word) for word in words]

    def __readOpac(self):
        """
        Read a single opacity from the file
        """
        
        # The next line contains a lot of important information:
        line = self.file.readline()
        
        znuc  = int(line[0:10]) / 100000
        npts  = int(line[0:10]) % 100000
        anum  = float(line[11:20])
        zbar  = float(line[21:30])
        rho   = float(line[31:40])
        tele  = float(line[41:50])
        eta   = float(line[51:60])
        rmean = float(line[61:70])
        pmean = float(line[71:80])

        energies  = np.empty(npts)
        opac_tot  = np.empty(npts)
        opac_bb   = np.empty(npts)
        opac_bf   = np.empty(npts)
        opac_ff   = np.empty(npts)
        opac_scat = np.empty(npts)

        # The next line contains the electron number density. This
        # isn't really needed because we have rho and zbar.
        self.file.readline()

        for n in xrange(npts):
            numbers = self.__splitline(self.file.readline())

            energies[n]  = numbers[0]
            opac_tot[n]  = numbers[1]
            opac_bb[n]   = numbers[2]
            opac_bf[n]   = numbers[3]
            opac_ff[n]   = numbers[4]
            opac_scat[n] = numbers[5]

        self.opacs.append(OpacInfernoSingle(znuc  = znuc, 
                                            anum  = anum,
                                            zbar  = zbar, 
                                            rho   = rho,
                                            tele  = tele,
                                            eta   = eta,
                                            rmean = rmean,
                                            pmean = pmean,
                                            energies = energies,
                                            opac_tot = opac_tot,
                                            opac_bb = opac_bb,
                                            opac_bf = opac_bf,
                                            opac_ff = opac_ff,
                                            opac_scat = opac_scat))

        # Skip the last line and return...
        self.file.readline()
        return

    def oplTotal(self):
        return OplList(len(self.opacs),
                       lambda n: (self.opacs[n].rho, self.opacs[n].tele), 
                       lambda n: self.opacs[n].energies,
                       lambda n: self.opacs[n].opac_tot)

    def oplScat(self):
        return OplList(len(self.opacs),
                       lambda n: (self.opacs[n].rho, self.opacs[n].tele), 
                       lambda n: self.opacs[n].energies,
                       lambda n: self.opacs[n].opac_scat)

    def oplAbsorb(self):
        return OplList(len(self.opacs),
                       lambda n: (self.opacs[n].rho, self.opacs[n].tele), 
                       lambda n: self.opacs[n].energies,
                       lambda n: self.opacs[n].opac_bb + self.opacs[n].opac_bf + self.opacs[n].opac_ff)
