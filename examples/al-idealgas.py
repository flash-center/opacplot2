"""
This script generates an IONMIX file with an ideal gas structure which
includes electron entropy, for fully ionized Aluminum.

The Sackur-Tetrode equation is used to compute the electron entropy.

The file is written out as "al-idealgas-001.cn6" and is then read back
in to ensure the reader works.
"""

from math import *
import numpy as np
import opacplot2 as opp

a = 26.9815386
z = 13.0

mpi = a / opp.NA 

temps = np.logspace(-1, 5, 13)
nt = len(temps)
print "temps =", temps

nions = np.logspace(12, 27, 16)
neles = z * nions
nni   = len(nions)
print "nions =", nions


# Make arrays:
zbar = np.empty((nni, nt))
pion = np.empty((nni, nt))
pele = np.empty((nni, nt))
eion = np.empty((nni, nt))
eele = np.empty((nni, nt))
sele = np.empty((nni, nt))

for jni, nion in enumerate(nions):
    for jt, temp in enumerate(temps):
        zbar[jni, jt] = z
        pion[jni, jt] = nion * opp.KB * temp
        pele[jni, jt] = z * pion[jni, jt]

        eion[jni, jt] = 1.5 * opp.KB * temp / mpi
        eele[jni, jt] = z * eion[jni, jt]

        # Compute specific electron entropy as given by the
        # Sackur-Tetrode equation:
        sele[jni,jt] = (opp.KB * opp.NA * z/a) * \
            ( 1.5 * log(opp.KB * temp * opp.ME / (opp.PLANCK * opp.HBAR)) -
              log(z*mpi*nion/a) + 5.0/2.0 )

# Write the EOS data out to a file: "al-idealgas-001.cn6":
opp.writeIonmixFile("al-idealgas-001.cn6", (z,), (1.0,), nions, temps,
                    zbar=zbar, 
                    pion=pion, pele=pele,
                    eion=eion, eele=eele,
                    sele=sele)

# Read the EOS data back in:
opg_imx = opp.OpacIonmix("al-idealgas-001.cn6", mpi, 
                         twot=True, man=True, hassele=True)

for jni, nion in enumerate(nions):
    for jt, temp in enumerate(temps):
        print "%6i %6i %15.6e %15.6e %15.6e" % (jni, jt, nion, temp, opg_imx.sele[jni,jt])
