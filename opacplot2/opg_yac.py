from constants import JOULE_TO_ERG
from StringIO import StringIO
from opl_grid import OplGrid

import numpy as np
# from util import readword
# from GridUtil import CONV_JOULES_TO_ERGS, BOLTZMANN_ERGEV
# from mat_info import OpacTable

def readword(f, advance=True):
    """
    Read one word from a file-like object. Returns None when EOF is
    reached.

    f [in/out] the file-like object to read
    
    advance [in] When True, readword will advance the file after the
    word is read to the next non-white space character
    """

    def backup():
        f.seek(-1,1)

    def eof(char):
        if len(char) == 0: return True
        return False

    def strip():
        c = f.read(1)
        while c.isspace():
            c = f.read(1)
        backup()
        return eof(c)
        
    if strip(): # Strip out leading white space
        return None # end of file reached

    word = ""

    c = f.read(1)    
    while not c.isspace():
        word += c
        c = f.read(1)

    if advance == True: strip()
    return word


class OpacYac:
    """
    Class for parsing and storing UW equation of state and opacity
    files.
    """

    def __init__(self, fn, mpi, verbose=False):
        self.fn = fn
        self.data = StringIO(open(fn,'r').read())
        self.mpi = mpi
        self.verb = verbose

        if verbose: print "Reading UW EOS file \"%s\"\n" % (fn)

        self.nheader = 24
        self.parse()

    def read_mesh_block(self, block_contents, n=None):
        verbose = block_contents and self.verb

        if n == None:
            n = int(readword(self.data))
        vals = np.zeros(n)

        if verbose: print "  Number of %s = %i" % (block_contents, n)
        for i in range(n):
            vals[i] = float(readword(self.data))
            if verbose: print "  %15.6e" % (vals[i])
        if verbose: print ""
        return vals        


    def parse(self):
        
        # Skip the header:
        for i in range(self.nheader): self.data.readline()

        # Read EOS temperatures and densities:
        self.eos_temps = self.read_mesh_block("EOS Temperatures (eV)")
        self.eos_dens = self.mpi * self.read_mesh_block("EOS Densitites (g/cc)") * self.mpi

        nteos = len(self.eos_temps)
        ndeos = len(self.eos_dens)

        # Skip rho0 (I think this number specifies the solid
        # density...):
        readword(self.data)
        
        # Skip some comment lines:
        for i in range(4): self.data.readline()
        
        # Read Opacity temperatures and densities:
        self.opac_temps = self.read_mesh_block("Opacity Temperatures (eV)")
        self.opac_dens = self.read_mesh_block("Opacity Densitites (g/cc)") * self.mpi
        ntopac = len(self.opac_temps)
        ndopac = len(self.opac_dens)

        # Read the number of opacity groups:
        self.ngroups = int(readword(self.data))

        # Skip a comment line:
        self.data.readline()

        # Read the group structure:
        self.opac_bounds = self.read_mesh_block("Opacity Group Bounds (eV)", self.ngroups+1)
        
        # Skip a comment line:
        self.data.readline()

        # Read zbar:
        self.zbar = self.read_mesh_block("", nteos*ndeos).reshape(ndeos,nteos)

        # Skip a comment line:
        self.data.readline()
        
        # Read Eint:
        self.eint = self.read_mesh_block("", nteos*ndeos).reshape(ndeos,nteos) * JOULE_TO_ERG

        # Skip a comment line:
        self.data.readline()
        
        # Read dE/dT:
        self.dedt = self.read_mesh_block("", nteos*ndeos).reshape(ndeos,nteos) * JOULE_TO_ERG

        # Skip a comment line:
        self.data.readline()
        
        # Read dE/dN:
        self.dedn = self.read_mesh_block("", nteos*ndeos).reshape(ndeos,nteos) * JOULE_TO_ERG

        # Skip a comment line:
        self.data.readline()
        
        # Read Eion:
        self.eion = self.read_mesh_block("", nteos*ndeos).reshape(ndeos,nteos) * JOULE_TO_ERG

        # Skip a comment line:
        self.data.readline()
        
        # Read Eele:
        self.eele = self.read_mesh_block("", nteos*ndeos).reshape(ndeos,nteos) * JOULE_TO_ERG
        
        # Skip a comment line:
        self.data.readline()
        
        # Read dEi/dT:
        self.deidt = self.read_mesh_block("", nteos*ndeos).reshape(ndeos,nteos) * JOULE_TO_ERG

        # Skip a comment line:
        self.data.readline()
        
        # Read dEe/dT:
        self.deedt = self.read_mesh_block("", nteos*ndeos).reshape(ndeos,nteos) * JOULE_TO_ERG

        # Skip a comment line:
        self.data.readline()
        
        # Read Pion:
        self.pion = self.read_mesh_block("", nteos*ndeos).reshape(ndeos,nteos)

        # Skip a comment line:
        self.data.readline()
        
        # Read Pele:
        self.pele = self.read_mesh_block("", nteos*ndeos).reshape(ndeos,nteos)

        # Skip a comment line:
        self.data.readline()
        
        # Read Pi/dT:
        self.dpidt = self.read_mesh_block("", nteos*ndeos).reshape(ndeos,nteos)

        # Skip a comment line:
        self.data.readline()
        
        # Read Pe/dT:
        self.dpedt = self.read_mesh_block("", nteos*ndeos).reshape(ndeos,nteos)


        # Skip a comment line:
        self.data.readline()
        
        # Read Rosseland Mean Opacity:
        data = self.read_mesh_block("", ntopac*ndopac*self.ngroups)
        self.rosseland = np.empty((ndopac, ntopac, self.ngroups))

        i = 0
        for g in xrange(self.ngroups):
            for d in xrange(ndopac):
                for t in xrange(ntopac):
                    # Read opacity and convert to units of 1/cm:
                    self.rosseland[d,t,g] = data[i]
                    i += 1

        # Skip a comment line:
        self.data.readline()
        
        # Read Planck emission opacity:
        data = self.read_mesh_block("", ntopac*ndopac*self.ngroups)
        self.planck_emiss = np.empty((ndopac, ntopac, self.ngroups))

        i = 0
        for g in xrange(self.ngroups):
            for d in xrange(ndopac):
                for t in xrange(ntopac):
                    # Read opacity and convert to units of 1/cm:
                    self.planck_emiss[d,t,g] = data[i]
                    i += 1


        # Skip a comment line:
        self.data.readline()
        
        # Read Planck absorption opacity:
        data = self.read_mesh_block("", ntopac*ndopac*self.ngroups)
        self.planck_absorb = np.empty((ndopac, ntopac, self.ngroups))

        i = 0
        for g in xrange(self.ngroups):
            for d in xrange(ndopac):
                for t in xrange(ntopac):
                    # Read opacity and convert to units of 1/cm:
                    self.planck_absorb[d,t,g] = data[i]
                    i += 1
                    
    def planckAbsorb(self):
        return OplGrid(self.opac_dens, self.opac_temps, self.opac_bounds, 
                       lambda jd, jt: self.planck_absorb[jd,jt,:])
        
