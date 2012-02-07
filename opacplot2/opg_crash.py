from opl_list import OplList
import numpy as np

class OpacCrash:

    ev_to_k = 11604.55
    kgm3_to_gcm3 = 1.0e-03
    m2okg_to_cm2og = 100.0**2/1000.0

    def __init__(self, fn, energies, verbose=False):
        self.fn = fn

        self.energies = energies

        data = np.loadtxt(fn)

        self.dens = np.empty(len(data))
        self.temps = np.empty(len(data))
        self.planck = np.empty((len(data),len(self.energies)-1))
        self.rosseland = np.empty((len(data),len(self.energies)-1))

        for n,row in enumerate(data):
            self.dens[n] = 10**row[0] * self.kgm3_to_gcm3
            self.temps[n] = 10**row[1] / self.ev_to_k
            self.planck[n,:] = row[2:2+len(self.energies)-1] * self.m2okg_to_cm2og
            self.rosseland[n,:] = row[2+len(self.energies)-1:] * self.m2okg_to_cm2og

        if verbose == True:
            print "\nNumber of Opacities =", len(self.dens)

            print "\nNumber of Energy Groups =", len(self.energies)
            for i,x in enumerate(self.energies):
                print "%6i  %15.6e" % (i,x)

    def oplPlanck(self):
        return OplList(len(self.dens), 
                       lambda n: (self.dens[n], self.temps[n]), 
                       lambda n: self.energies,
                       lambda n: self.planck[n])

    def oplRosseland(self):
        return OplList(len(self.dens), 
                       lambda n: (self.dens[n], self.temps[n]), 
                       lambda n: self.energies,
                       lambda n: self.rosseland[n])


    # def opl(self):

    #     def getOpac(jd,jt):
    #         row = jd + len(self.dens)*jt
    #         opac = self.data[row,:] 
    #         #return opac[2:2+len(self.energies)-1].reshape(1,1,len(self.energies)-1)
    #         return opac[2:2+len(self.energies)-1]
        
    #     return OplGrid(self.dens, self.temps, self.energies, getOpac)


# class OpacCrash:

#     ev_to_k = 11604.55

#     def __init__(self, fn, ntemp, ndens, energies, verbose=False):
#         self.fn = fn

#         self.energies = energies

#         self.data = np.loadtxt(fn)

#         # Load densities
#         self.dens = 10**(self.data[0:ndens,0])

#         # Convert from kg/m^3 to g/cm^3
#         self.dens = self.dens / 1000.0

#         print "Number of Densities =", len(self.dens)
#         for i,x in enumerate(self.dens):
#             print "%6i  %15.6e" % (i,x)

#         # Load temperatures
#         self.temps = 10**(self.data[::ndens,1])

#         # Convert from K to eV
#         self.temps = self.temps / self.ev_to_k

#         print "\nNumber of Temperatures =", len(self.temps)
#         for i,x in enumerate(self.temps):
#             print "%6i  %15.6e" % (i,x)


#         print "\nNumber of Energy Groups =", len(self.energies)
#         for i,x in enumerate(self.energies):
#             print "%6i  %15.6e" % (i,x)


#     def opl(self):

#         def getOpac(jd,jt):
#             row = jd + len(self.dens)*jt
#             opac = self.data[row,:] 
#             #return opac[2:2+len(self.energies)-1].reshape(1,1,len(self.energies)-1)
#             return opac[2:2+len(self.energies)-1]
        
#         return OplGrid(self.dens, self.temps, self.energies, getOpac)
