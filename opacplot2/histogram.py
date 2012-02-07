import numpy as np

def histdata(en, op):
    energies = []
    opacs    = []
    for i in range(len(en)-1):
        energies += [ en[i], en[i+1] ]
        opacs    += [ op[i], op[i] ]
                
    return np.array(energies), np.array(opacs)
