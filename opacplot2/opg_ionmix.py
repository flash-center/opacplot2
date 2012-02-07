from StringIO import StringIO
import numpy as np
import re 
import math

from opl_grid import OplGrid
from constants import ERG_TO_JOULE

class OpacIonmix:
    """
    Class to read in IONMIX EOS and Opacity Files

    All energies in this file are in Joules and must be converted to ergs
    """

    joules_to_ergs = 1.0e+07


    def __init__(self, fn, mpi, twot=False, man=False, verbose=False):
        
        self.fn = fn
        self.mpi = mpi
        self.twot = twot
        self.man  = man
        self.verb = verbose
        if verbose: print "Reading IONMIX file \"%s\"\n" % (fn)

        f = open(fn,'r')

        # Read the number of temperatures/densities:
        self.ntemp = int(f.read(10))
        self.ndens = int(f.read(10))

        # Skip the next three lines:
        for i in range(3): f.readline()

        # Setup temperature/density grid:
        if self.man == False:
            # Read information about the temperature/density grid:
            self.ddens_log10 = float(f.read(12))
            self.dens0_log10 = float(f.read(12))
            self.dtemp_log10 = float(f.read(12))
            self.temp0_log10 = float(f.read(12))

            # Compute number densities:
            self.numDens = np.logspace(self.dens0_log10, 
                                       self.dens0_log10+self.ddens_log10*(self.ndens-1), 
                                       self.ndens)

            self.temps = np.logspace(self.temp0_log10, 
                                     self.temp0_log10+self.dtemp_log10*(self.ntemp-1), 
                                     self.ntemp)

            # Read number of groups:
            self.ngroups = int(f.read(12))
        else:
            self.ngroups = int(f.read(12))
            f.readline()

        # Read the rest of the file, remove all of the white space,
        # and store the string in self.data:
        self.data = StringIO(re.sub(r'\s', '', f.read()))
                        
        if self.man == True:
            # For files where temperatures/densities are manually
            # specified, read the manual values here.
            self.temps = self.get_block(self.ntemp)
            self.numDens = self.get_block(self.ndens)

        self.dens = self.numDens * self.mpi
            
        if self.verb: 
            print "  Number of temperatures: %i" % self.ntemp
            for i in range(0, self.ntemp):
                print "%6i%27.16e" % (i, self.temps[i])

            print "\n  Number of densities: %i" % self.ndens
            for i in range(0, self.ndens):
                print "%6i%21.12e%27.16e" % (i, self.dens[i], self.numDens[i])

        self.read_eos()
        self.read_opac()

    def get_block(self,n):
        arr = np.zeros(n)
        for i in range(n):
            arr[i] = float(self.data.read(12))
        return arr

    def read_eos(self):
        nt = self.ntemp
        nd = self.ndens
        ng = self.ngroups

        self.zbar  = self.get_block(nd*nt).reshape(nd,nt)

        if self.twot == False:
            # Read in e and cv, but convert from J to ergs:
            self.etot  = self.get_block(nd*nt).reshape(nd,nt) * self.joules_to_ergs
            self.cvtot = self.get_block(nd*nt).reshape(nd,nt) * self.joules_to_ergs
            self.dedn  = self.get_block(nd*nt).reshape(nd,nt)

        else: 
            # Read in pressure, specific internal energies and
            # specific heats, but convert from J to ergs:
            self.dzdt  = self.get_block(nd*nt).reshape(nd,nt)
            self.pion  = self.get_block(nd*nt).reshape(nd,nt) * self.joules_to_ergs
            self.pele  = self.get_block(nd*nt).reshape(nd,nt) * self.joules_to_ergs
            self.dpidt = self.get_block(nd*nt).reshape(nd,nt) * self.joules_to_ergs
            self.dpedt = self.get_block(nd*nt).reshape(nd,nt) * self.joules_to_ergs
            self.eion  = self.get_block(nd*nt).reshape(nd,nt) * self.joules_to_ergs
            self.eele  = self.get_block(nd*nt).reshape(nd,nt) * self.joules_to_ergs
            self.cvion = self.get_block(nd*nt).reshape(nd,nt) * self.joules_to_ergs
            self.cvele = self.get_block(nd*nt).reshape(nd,nt) * self.joules_to_ergs
            self.deidn = self.get_block(nd*nt).reshape(nd,nt) * self.joules_to_ergs
            self.deedn = self.get_block(nd*nt).reshape(nd,nt) * self.joules_to_ergs


    def read_opac(self):
        """
        Load the opacities from the file. The opacities are arranged
        in the file so that temperature varies the fastest, then
        density, and group number varies the slowest. Note, that this
        is not the ordering of the arrays once they are loaded.
        """

        nt = self.ntemp
        nd = self.ndens
        ng = self.ngroups
        
        # Read group bounds in eV and convert to ergs:
        self.opac_bounds = self.get_block(ng+1)

        if self.verb: 
            print "\n  Number of Energy Groups: %i" % self.ngroups
            for i in range(0, self.ngroups+1):
                print "%6i%15.6e" % (i, self.opac_bounds[i])


        self.rosseland     = np.empty((nd,nt,ng))
        self.planck_absorb = np.empty((nd,nt,ng))
        self.planck_emiss  = np.empty((nd,nt,ng))

        arr_ro = self.get_block(nd*nt*ng)
        arr_pa = self.get_block(nd*nt*ng)
        arr_pe = self.get_block(nd*nt*ng)

        i = 0
        for g in xrange(ng):
            for d in xrange(nd):
                for t in xrange(nt):
                    self.rosseland[d,t,g]     = arr_ro[i]
                    self.planck_absorb[d,t,g] = arr_pa[i]
                    self.planck_emiss[d,t,g]  = arr_pe[i]
                    i += 1    
        
    def oplAbsorb(self):
        return OplGrid(self.dens, self.temps, self.opac_bounds, 
                       lambda jd, jt: self.planck_absorb[jd,jt,:])

    def oplEmiss(self):
        return OplGrid(self.dens, self.temps, self.opac_bounds, 
                       lambda jd, jt: self.planck_emiss[jd,jt,:])

    def oplRosseland(self):
        return OplGrid(self.dens, self.temps, self.opac_bounds, 
                       lambda jd, jt: self.rosseland[jd,jt,:])

    def write(self, fn, zvals, fracs, twot=None, man=None):
        if twot == None: twot = self.twot
        if twot == True and self.twot == False:
            raise ValueError("Error: Cannot write two-temperature data")

        if man == None: man = self.man
        if man == False and self.man == True:
            raise ValueError("Error: Cannot write manual temp/dens points")

        # Write the header:
        f = open(fn,'w')
        f.write("%10i%10i\n" % (self.ntemp,self.ndens))
        f.write(" atomic #s of gases: ")
        for z in zvals: f.write("%10i" % z)
        f.write("\n relative fractions: ")
        for frac in fracs: f.write("%10.2E" % frac)
        f.write("\n")

        # Write temperature/density grid and number of groups:
        def convert(num):
            string_org = "%12.5E" % (num)
            negative = (string_org[0] == "-")            
            lead = "-." if negative else "0."
            string = lead + string_org[1] + string_org[3:8] + "E"

            # Deal with the exponent:
            
            # Check for zero:
            if int(string_org[1] + string_org[3:8]) == 0:
                return string + "+00"

            # Not zero:
            expo = int(string_org[9:]) + 1
            if expo < 0:
                string += "-"
            else:
                string += "+"
            string += "%02d" % abs(expo)
            return string

        def write_block(var):
            count = 0
            for n in xrange(len(var)):
                count += 1

                f.write("%s" % convert(var[n]))
                    
                if count == 4:
                    count = 0
                    f.write("\n")

            if count != 0: f.write("\n")

        def write_opac_block(var):
            count = 0
            for g in xrange(self.ngroups):
                for jd in xrange(self.ndens):
                    for jt in xrange(self.ntemp):
                        count += 1

                        f.write("%s" % convert(var[jd,jt,g]))
                
                        if count == 4:
                            count = 0
                            f.write("\n")

            if count != 0: f.write("\n")

        if man == False:    
            f.write("%s%s%s%s" % (convert(self.ddens_log10), 
                                  convert(self.dens0_log10), 
                                  convert(self.dtemp_log10), 
                                  convert(self.temp0_log10)) )
            
        f.write("%12i\n" % self.ngroups)

        if man == True:
            write_block(self.temps)
            write_block(self.numDens)

        write_block(self.zbar.flatten())

        if twot == False:
            write_block(self.etot.flatten()/self.joules_to_ergs)
            write_block(self.cvtot.flatten()/self.joules_to_ergs)
            write_block(self.enntab.flatten())

        else:
            write_block( self.dzdt.flatten())
            write_block( self.pion.flatten()/self.joules_to_ergs)
            write_block( self.pele.flatten()/self.joules_to_ergs)
            write_block(self.dpidt.flatten()/self.joules_to_ergs)
            write_block(self.dpedt.flatten()/self.joules_to_ergs)
            write_block(self.eion.flatten()/self.joules_to_ergs)
            write_block(self.eele.flatten()/self.joules_to_ergs)
            write_block(self.cvion.flatten()/self.joules_to_ergs)
            write_block(self.cvele.flatten()/self.joules_to_ergs)
            write_block(self.deidn.flatten()/self.joules_to_ergs)
            write_block(self.deedn.flatten()/self.joules_to_ergs)

        write_block(self.opac_bounds)
        write_opac_block(self.rosseland)
        write_opac_block(self.planck_absorb)
        write_opac_block(self.planck_emiss)


    def extendToZero(self):
        """
        This routine adds another temperature point at zero
        """

        nd = self.ndens
        nt = self.ntemp
        ng = self.ngroups

        arr = np.zeros((nd,nt+1))
        arr[:,1:] = self.dzdt[:,:]
        self.dzdt = arr

        arr = np.zeros((nd,nt+1))
        arr[:,1:] = self.pion[:,:]
        self.pion = arr

        arr = np.zeros((nd,nt+1))
        arr[:,1:] = self.pele[:,:]
        self.pele = arr

        arr = np.zeros((nd,nt+1))
        arr[:,1:] = self.dpidt[:,:]
        self.dpidt = arr

        arr = np.zeros((nd,nt+1))
        arr[:,1:] = self.dpedt[:,:]
        self.dpedi = arr

        arr = np.zeros((nd,nt+1))
        arr[:,1:] = self.eion[:,:]
        self.eion = arr

        arr = np.zeros((nd,nt+1))
        arr[:,1:] = self.eele[:,:]
        self.eele = arr

        arr = np.zeros((nd,nt+1))
        arr[:,1:] = self.zbar[:,:]
        self.zbar = arr

        arr = np.zeros((nd,nt+1))
        arr[:,1:] = self.cvion[:,:]
        self.cvion = arr

        arr = np.zeros((nd,nt+1))
        arr[:,1:] = self.cvele[:,:]
        self.cvele = arr

        arr = np.zeros((nd,nt+1))
        arr[:,1:] = self.deidn[:,:]
        self.deidn = arr

        arr = np.zeros((nd,nt+1))
        arr[:,1:] = self.deedn[:,:]
        self.deedn = arr

        arr = np.zeros((nd,nt+1,ng))
        arr[:,1:,:] = self.rosseland[:,:,:]
        self.rosseland = arr

        arr = np.zeros((nd,nt+1,ng))
        arr[:,1:,:] = self.planck_absorb[:,:,:]
        self.planck_absorb = arr

        arr = np.zeros((nd,nt+1,ng))
        arr[:,1:,:] = self.planck_emiss[:,:,:]
        self.planck_emiss = arr

        # Reset temperatures:

        arr = np.zeros((nt+1))
        arr[1:] = self.temps[:]
        self.temps = arr

        self.ntemp += 1


def writeIonmixFile(fn, zvals, fracs, ndens, ntemps, numDens, temps, 
                    zbar=None,  dzdt=None, pion=None, pele=None,
                    dpidt=None, dpedt=None, eion=None, eele=None,
                    cvion=None, cvele=None, deidn=None, deedn=None,
                    ngroups=None, opac_bounds=None,
                    rosseland=None, planck_absorb=None, planck_emiss=None):

    if  zbar == None:  zbar = np.zeros((ndens,ntemps))
    if  dzdt == None:  dzdt = np.zeros((ndens,ntemps))
    if  pion == None:  pion = np.zeros((ndens,ntemps))
    if  pele == None:  pele = np.zeros((ndens,ntemps))
    if dpidt == None: dpidt = np.zeros((ndens,ntemps))
    if dpedt == None: dpedt = np.zeros((ndens,ntemps))
    if  eion == None:  eion = np.zeros((ndens,ntemps))
    if  eele == None:  eele = np.zeros((ndens,ntemps))
    if cvion == None: cvion = np.zeros((ndens,ntemps))
    if cvele == None: cvele = np.zeros((ndens,ntemps))
    if deidn == None: deidn = np.zeros((ndens,ntemps))
    if deedn == None: deedn = np.zeros((ndens,ntemps))

    if ngroups       == None: ngroups = 1
    if opac_bounds   == None: opac_bounds = (0.0,1.0)
    if rosseland     == None: rosseland = np.zeros((ndens,ntemps,ngroups))
    if planck_absorb == None: planck_absorb = np.zeros((ndens,ntemps,ngroups))
    if planck_emiss  == None: planck_emiss = np.zeros((ndens,ntemps,ngroups))

    # Write the header:
    f = open(fn,'w')
    f.write("%10i%10i\n" % (ntemps,ndens))
    f.write(" atomic #s of gases: ")
    for z in zvals: f.write("%10i" % z)
    f.write("\n relative fractions: ")
    for frac in fracs: f.write("%10.2E" % frac)
    f.write("\n")

    # Write temperature/density grid and number of groups:
    def convert(num):
        string_org = "%12.5E" % (num)
        negative = (string_org[0] == "-")            
        lead = "-." if negative else "0."
        string = lead + string_org[1] + string_org[3:8] + "E"

        # Deal with the exponent:
        
        # Check for zero:
        if int(string_org[1] + string_org[3:8]) == 0:
            return string + "+00"

        # Not zero:
        expo = int(string_org[9:]) + 1
        if expo < 0:
            string += "-"
        else:
            string += "+"
        string += "%02d" % abs(expo)
        return string

    def write_block(var):
        count = 0
        for n in xrange(len(var)):
            count += 1

            f.write("%s" % convert(var[n]))
                
            if count == 4:
                count = 0
                f.write("\n")

        if count != 0: f.write("\n")

    def write_opac_block(var):
        count = 0
        for g in xrange(ngroups):
            for jd in xrange(ndens):
                for jt in xrange(ntemps):
                    count += 1

                    f.write("%s" % convert(var[jd,jt,g]))
            
                    if count == 4:
                        count = 0
                        f.write("\n")

        if count != 0: f.write("\n")
        
    f.write("%12i\n" % ngroups)

    write_block(temps)
    write_block(numDens)

    write_block(zbar.flatten())
    
    write_block(dzdt.flatten())
    write_block(pion.flatten()*ERG_TO_JOULE)
    write_block(pele.flatten()*ERG_TO_JOULE)
    write_block(dpidt.flatten()*ERG_TO_JOULE)
    write_block(dpedt.flatten()*ERG_TO_JOULE)
    write_block(eion.flatten()*ERG_TO_JOULE)
    write_block(eele.flatten()*ERG_TO_JOULE)
    write_block(cvion.flatten()*ERG_TO_JOULE)
    write_block(cvele.flatten()*ERG_TO_JOULE)
    write_block(deidn.flatten()*ERG_TO_JOULE)
    write_block(deedn.flatten()*ERG_TO_JOULE)

    write_block(opac_bounds)
    write_opac_block(rosseland)
    write_opac_block(planck_absorb)
    write_opac_block(planck_emiss)
