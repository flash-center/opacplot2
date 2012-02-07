import numpy as np
from constants import KELVIN_TO_EV, GPA_TO_ERGCC, MJKG_TO_ERGCC


class OpgSesame:
    """
    This class is responsible for loading all SESAME formatted data files
    """
    
    SINGLE = 1
    DOUBLE = 2

    CHAR_LINE_LEN = 80    
    WORDS_PER_LINE = 5

    def __init__(self, filename, precision, verbose=False):
        self.verbose = verbose
        self.fhand = open(filename)
        if(precision == self.SINGLE):
           self.entry_len = 15
        elif(precision == self.DOUBLE):
            self.entry_len = 22
        else:
            raise ValueError("precision must be SINGLE or DOUBLE")
        

        self.fdict = { 101 : self.parseComment, 
                       102 : self.parseComment, 
                       103 : self.parseComment,
                       104 : self.parseComment,
                       201 : self.parseInfo,
                       301 : self.parseEos,
                       303 : self.parseEos,
                       304 : self.parseEos,
                       305 : self.parseEos,
                       306 : self.parseEos,
                       401 : self.parseVap,
                       411 : self.parseSolid,
                       412 : self.parseLiquid,
                       431 : self.parseShear,
                       601 : self.parseZbar,
                       602 : self.parseEcond,
                       603 : self.parseTcond,
                       604 : self.parseTcond,
                       605 : self.parseTcond,
                       }

        self.data = {}

        self.parse()

    def parse(self):

        while True:
            header = self.fhand.readline()
            if not header: break # Reached EOF
            if header[:3] == " 2 ": break

            matid = int(header.split()[1])
            recid = int(header.split()[2])
            nentries = int(header.split()[3])

            if not matid in self.data: self.data[matid] = {}

            if self.verbose and (recid > 104):
                print "Material = %8i  Record = %8i  Entries = %8i" % (matid, recid, nentries)

            if not recid in self.fdict:
                raise ValueError("No handling function for record %d" % recid)

            self.fdict[recid](nentries,matid, recid)

    def parseComment(self, nentries, matid, recid):

        nlines = (nentries-1) / self.CHAR_LINE_LEN + 1
        nchar = nentries + nlines        
        return self.fhand.read(nchar)

    def parseInfo(self, nentries, matid, recid):
        words = self.readEntries(nentries)
        self.data[matid]["zmax"] = words[0]
        self.data[matid]["abar"] = words[1]
        self.data[matid]["rho0"] = words[2]
        self.data[matid]["bulkmod"] = words[3]
        self.data[matid]["excoef"] = words[4]

        if self.verbose: print "  zbar = %g  abar = %g  rho0 = %g" % (words[0],words[1],words[2])

        return self.data

    def parseEos(self, nentries, matid, recid):
        words = self.readEntries(nentries)

        prefixes = { 301 : "total_",
                     303 : "ioncc_",
                     304 : "ele_",
                     305 : "ion_",
                     306 : "cc_" }

        prefix = prefixes[recid]

        ndens = int(words[0])
        ntemp = int(words[1])
        start = 2

        self.data[matid][prefix+"ndens"] = ndens
        self.data[matid][prefix+"ntemp"] = ntemp

        self.data[matid][prefix+"dens"] = words[start:start+ndens] # In g/cc
        start = start+ndens

        self.data[matid][prefix+"temps"] = words[start:start+ntemp]*KELVIN_TO_EV
        start = start + ntemp        

        if self.verbose and prefix == "total_":
            dens = self.data[matid]["total_dens"]
            temps = self.data[matid]["total_temps"]

            print "ndens   = %13i ntemp    = %13i" % (ndens, ntemp)
            print "dens[0] = %13.6e  dens[-1] = %13.6e" % (dens[0], dens[-1])
            print "temp[0] = %13.6e  temp[-1] = %13.6e" % (temps[0], temps[-1])


        # Read pressure array (in GPa and convert to ergs/cc):
        self.data[matid][prefix+"pres"] = words[start:start+ntemp*ndens].reshape((ntemp,ndens)).transpose()*GPA_TO_ERGCC
        start = start + ntemp*ndens

        # Read specific internal energy array (in MJ/kg and convert to ergs/g):
        self.data[matid][prefix+"eint"] = words[start:start+ntemp*ndens].reshape((ntemp,ndens)).transpose()*MJKG_TO_ERGCC
        start = start + ntemp*ndens

        if start == nentries:
            # There is no Helmholtz free energy data, this is the end
            # of this record...
            return

        # Read specific Helmholtz free energy array (in MJ/kg and convert to ergs/g):
        self.data[matid][prefix+"free"] = words[start:start+ntemp*ndens].reshape((ntemp,ndens)).transpose()*MJKG_TO_ERGCC
        start = start + ntemp*ndens

        # exit()

    def parseVap(self, nentries, matid, recid):    
        words = self.readEntries(nentries)

    def parseSolid(self, nentries, matid, recid):    
        words = self.readEntries(nentries)

    def parseLiquid(self, nentries, matid, recid):    
        words = self.readEntries(nentries)

    def parseShear(self, nentries, matid, recid):    
        words = self.readEntries(nentries)

    def parseZbar(self, nentries, matid, recid):    
        words = self.readEntries(nentries)

        prefix = "zbar"

        ndens = int(words[0])
        ntemp = int(words[1])
        start = 2

        self.data[matid][prefix+"_ndens"] = ndens
        self.data[matid][prefix+"_ntemp"] = ntemp

        self.data[matid][prefix+"_dens"] = 10**words[start:start+ndens]
        start = start+ndens

        self.data[matid][prefix+"_temps"] = 10**words[start:start+ntemp]
        start = start + ntemp        

        self.data[matid][prefix] = 10**(words[start:start+ntemp*ndens].reshape((ntemp,ndens)).transpose())

    def parseEcond(self, nentries, matid, recid):    
        words = self.readEntries(nentries)

    def parseTcond(self, nentries, matid, recid):    
        words = self.readEntries(nentries) 
   
    def readEntries(self,nentries):
        nlines = (nentries-1) / self.WORDS_PER_LINE + 1

        data = np.empty(nentries)

        string = ""
        for i in xrange(nlines):
            string += self.fhand.readline()[:self.WORDS_PER_LINE*self.entry_len]

        for i in xrange(nentries):
            word = string[self.entry_len*i:self.entry_len*(i+1)]
            if word[-4] == '-': word = word[:-4] + "E" + word[-4:]
            data[i] = float(word)

        return data
