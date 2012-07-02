import opacplot2 as opp
from opacplot2 import RADCONST as AR, KB
import numpy as np
from matplotlib import pyplot as plt

mpi = 4.002602 / opp.NA
opimx = opp.OpacIonmix("data/al-imx-002.cn4", mpi, twot=True, man=True, verbose=True)

t1 = 0.0
t2 = 20.0
tvals = np.linspace(t1, t2, 1001)

dens = 2.16

# PLOT EELE:
vals =  np.array([opp.interpDT(opimx.eele, opimx.dens, opimx.temps, 
                               dens, t, bctmin=opp.BC_EXTRAP_ZERO, log=False)
                  for t in tvals])
plt.plot(tvals, vals, "r", lw=2)
plt.grid(True)
plt.ylabel("Electron Specific Internal Energy (ergs/g)")
plt.xlabel("Electron Temperature (eV)")

# PLOT CVEELE:
plt.figure()
vals =  np.array([opp.interpDT(opimx.eele, opimx.dens, opimx.temps, 
                               dens, t, bctmin=opp.BC_EXTRAP_ZERO, log=False,
                               lookup=opp.INTERP_DFDT)
                  for t in tvals])
plt.plot(tvals, vals, "r", lw=2)
plt.grid(True)
plt.ylabel("Electron Specific Heat (ergs/g/eV)")
plt.xlabel("Electron Temperature (eV)")
plt.show()
