import numpy as np

# Convert Joules to Ergs:
JOULE_TO_ERG = 1.0e+07

# Ergs to Joulse:
ERG_TO_JOULE = 1.0/JOULE_TO_ERG

# Convert GPa to Ergs/cc:
GPA_TO_ERGCC = 1.0e+10

# Convert MJ/kg to Ergs/cc:
MJKG_TO_ERGCC = 1.0e+10

# Convert Kelvin to eV:
KELVIN_TO_EV = 1.0/11604.55

# Avogadros number
NA = 6.0221415e+23

# Boltzmann constant [ergs/eV]
KB = 1.60217653e-12

# Radiation constant [ergs/cm^3/eV^4]
RADCONST = 137.201730837

BC_BOUND       = 0
BC_EXTRAP      = 1
BC_EXTRAP_ZERO = 2

INTERP_FUNC = 0
INTERP_DFDD = 1
INTERP_DFDT = 2

MAX_ELE = 100

ELE_HE = 4
ELE_AL = 13

# Atomic Weights:
ATOMIC_WEIGHTS = -1 * np.ones(MAX_ELE)

ATOMIC_WEIGHTS[ELE_HE] = 4.002602
ATOMIC_WEIGHTS[ELE_AL] = 26.9815386


# Atomic Symbols:
ATOMIC_SYMBOLS = ["H" , "He", "Li", "Be", "B", "C",
                  "N" , "O" ]  
