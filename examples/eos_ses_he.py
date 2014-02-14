#!/usr/python
# -*- coding: utf-8 -*-

# Very simple example to illistrate file convertion to IONMIX format.


import opacplot2 as opp
import numpy as np

# Parse the input file, this part is file dependent.
# For OPAC file format the parser might be opacplot2/opg_tabop.py
# otherwise it would have to be implemented.

eos_data = opp.OpgSesame("./he-ses-5761.txt", opp.OpgSesame.SINGLE, verbose=False)[5761]


# Write EoS/Opacity data, see opg_ionmix.py for writeIonmixFile class
# definion. When writing EoS file in 3T one has to provide eion, eele,
# pion, pele, zbar arguments. For a IONMIX file containing opacity data, 
# following arguments have to be provided: ngroups, opac_bounds, rosseland,
# planck_absorb, planck_emiss
fracs = (1.0,)
filename = 'he-ses-5761.cn4'
numDens = opp.NA * eos_data['ele_dens'] / eos_data["abar"] # convert from g.cm⁻³ to cm⁻³.
opp.writeIonmixFile(filename, (eos_data['zmax'],), fracs,
                        numDens=numDens, temps=eos_data['ele_temps'],
                        eion=eos_data["ioncc_eint"],
                        eele=eos_data["ele_eint"],
                        pion=eos_data["ioncc_pres"],
                        pele=eos_data["ele_pres"],
                        zbar=np.ones(eos_data["ele_pres"].shape)*2 # full ionization, just as an exemple
                        )


# Open generated file to check there are no errors.
ionmix = opp.OpacIonmix(filename, eos_data["abar"]/opp.NA, twot=True, man=True, verbose=False)
