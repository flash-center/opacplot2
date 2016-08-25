import opacplot2 as opp
from opacplot2 import RADCONST as AR, KB, BC_BOUND, BC_EXTRAP_ZERO
import numpy as np
from matplotlib import pyplot as plt


mpi = 4.002602 / opp.NA
heimx = opp.OpacIonmix("data/he-imx-005.cn4", mpi, twot=True, man=True, verbose=True)


mpi = 26.9815386 / opp.NA
alimx = opp.OpacIonmix("data/al-imx-002.cn4", mpi, twot=True, man=True, verbose=True)

def pion(opimx, rho, tion, log=False):
    return opp.interpDT(opimx.pion, opimx.dens, opimx.temps, 
                        rho, tion, log=False, lookup=opp.INTERP_FUNC,
                        bctmin=opp.BC_EXTRAP_ZERO, 
                        bcdmin=opp.BC_EXTRAP_ZERO)

def pele(opimx, rho, tele, log=False):
    return opp.interpDT(opimx.pele, opimx.dens, opimx.temps, 
                        rho, tele, log=False, lookup=opp.INTERP_FUNC,
                        bctmin=opp.BC_EXTRAP_ZERO, 
                        bcdmin=opp.BC_EXTRAP_ZERO)


def eion(opimx, rho, tion, log=False):
    return opp.interpDT(opimx.eion, opimx.dens, opimx.temps, 
                        rho, tion, log=False, lookup=opp.INTERP_FUNC,
                        bctmin=opp.BC_EXTRAP_ZERO)

def eele(opimx, rho, tele, log=False):
    return opp.interpDT(opimx.eele, opimx.dens, opimx.temps, 
                        rho, tele, log=False, lookup=opp.INTERP_FUNC,
                        bctmin=opp.BC_EXTRAP_ZERO)


def deidti(opimx, rho, tion, log=False):
    return opp.interpDT(opimx.eion, opimx.dens, opimx.temps, 
                        rho, tion, log=False, lookup=opp.INTERP_DFDT,
                        bctmin=opp.BC_EXTRAP_ZERO)

def deedte(opimx, rho, tele, log=False):
    return opp.interpDT(opimx.eele, opimx.dens, opimx.temps, 
                        rho, tele, log=False, lookup=opp.INTERP_DFDT,
                        bctmin=opp.BC_EXTRAP_ZERO)


def dpidti(opimx, rho, tion, log=False):
    return opp.interpDT(opimx.pion, opimx.dens, opimx.temps, 
                        rho, tion, log=False, lookup=opp.INTERP_DFDT,
                        bctmin=opp.BC_EXTRAP_ZERO, 
                        bcdmin=opp.BC_EXTRAP_ZERO)

def dpedte(opimx, rho, tele, log=False):
    return opp.interpDT(opimx.pele, opimx.dens, opimx.temps, 
                        rho, tele, log=False, lookup=opp.INTERP_DFDT,
                        bctmin=opp.BC_EXTRAP_ZERO,
                        bcdmin=opp.BC_EXTRAP_ZERO)


def dpidd(opimx, rho, tion, log=False):
    return opp.interpDT(opimx.pion, opimx.dens, opimx.temps, 
                        rho, tion, log=False, lookup=opp.INTERP_DFDD,
                        bctmin=opp.BC_EXTRAP_ZERO, 
                        bcdmin=opp.BC_EXTRAP_ZERO)

def dpedd(opimx, rho, tele, log=False):
    return opp.interpDT(opimx.pele, opimx.dens, opimx.temps, 
                        rho, tele, log=False, lookup=opp.INTERP_DFDD,
                        bctmin=opp.BC_EXTRAP_ZERO, 
                        bcdmin=opp.BC_EXTRAP_ZERO)




def gamc(opimx, rho, tele, tion, trad, 
         bcdmin=BC_BOUND, bcdmax=BC_BOUND, 
         bctmin=BC_BOUND, bctmax=BC_BOUND, 
         log = False, verbose=False):


    # if bcdmin == BC_BOUND and rho  < opimx.dens[0] : rho  = opimx.dens[0]
    # if bctmin == BC_BOUND and tele < opimx.temps[0]: tele = opimx.temps[0]
    # if bctmin == BC_BOUND and tion < opimx.temps[0]: tion = opimx.temps[0]

    # if bcdmax == BC_BOUND and rho  > opimx.dens[-1] : rho  = opimx.dens[-1]
    # if bctmax == BC_BOUND and tele > opimx.temps[-1]: tele = opimx.temps[-1]
    # if bctmax == BC_BOUND and tion > opimx.temps[-1]: tion = opimx.temps[-1]

    # def interpDT(var, r, t, lookup):
    #     return opp.interpDT(var, opimx.dens, opimx.temps, 
    #                         r, t,
    #                         lookup=lookup, log=log,                          
    #                         bcdmin=bcdmin, bcdmax=bcdmax,
    #                         bctmin=bctmin, bctmax=bctmax)

    # Compute needed derivatives:
    deidti_val = deidti(opimx, rho, tion)
    dpidti_val = dpidti(opimx, rho, tion)
    dpidd_val  = dpidd(opimx, rho, tion)

    dpidd_si = tion/(rho**2 * deidti_val) * dpidti_val**2 + dpidd_val    

    pion_val = pion(opimx, rho, tion) 


    deedte_val = deedte(opimx, rho, tele)
    dpedte_val = dpedte(opimx, rho, tele)
    dpedd_val  = dpedd(opimx, rho, tele)

    dpedd_se = tele/(rho**2 * deedte_val) * dpedte_val**2 + dpedd_val    

    pele_val = pele(opimx, rho, tele) 


    #deidti = interpDT(opimx.eion, rho, tion, opp.INTERP_DFDT)
    # dpidti = interpDT(opimx.pion, rho, tion, opp.INTERP_DFDT)
    # dpidd  = interpDT(opimx.pion, rho, tion, opp.INTERP_DFDD)

    # dpidd_si = tion/(rho**2 * deidti) * dpidti**2 + dpidd
    # interpDT(opimx.pion, rho, tion, opp.INTERP_FUNC)


    # deedte = interpDT(opimx.eele, rho, tele, opp.INTERP_DFDT)
    # dpedte = interpDT(opimx.pele, rho, tele, opp.INTERP_DFDT)
    # dpedd  = interpDT(opimx.pele, rho, tele, opp.INTERP_DFDD)
    # dpedd_se = tele/(rho**2 * deedte) * dpedte**2 + dpedd

    # pele = interpDT(opimx.pele, rho, tele, opp.INTERP_FUNC)

    prad = AR*trad**4
    ptot = prad + pion_val + pele_val
    gamc = rho/ptot * (dpedd_se + dpidd_si + (4*prad)/(3*rho))

    if verbose == True:
        print( "rho      = %17.10e" % rho)
        print( "tele     = %17.10e" % tele)
        print( "tion     = %17.10e" % tion)
        print( "deidti   = %17.10e" % deidti_val)
        print( "dpidti   = %17.10e" % dpidti_val)
        print( "dpidd    = %17.10e" % dpidd_val)
        print( "dpidd_si = %17.10e" % dpidd_si)
        print( "pion     = %17.10e" % pion_val)
        
        print( "deedte   = %17.10e" % deedte_val)
        print( "dpedte   = %17.10e" % dpedte_val)
        print( "dpedd    = %17.10e" % dpedd_val)
        print( "dpedd_se = %17.10e" % dpedd_se)
        print( "pele     = %17.10e" % pele_val)
        
        print( "prad     = %17.10e" % prad)
        print( "ptot     = %17.10e" % ptot)
        print( "gamc     = %17.10e" % gamc)
        print( "gamc_ion = %17.10e" % (dpidd_si * rho/pion_val))
        print( "gamc_ele = %17.10e" % (dpedd_se * rho/pele_val))

    return gamc

rho = 1.6377875581e-02
tion =  3.1504757702e-02
tele =  9.7453373671e-01
mfhe =  1.2141027255e-03
mfal =  9.9878591299e-01
eion_tot =  2.7614768640e+09
eele_tot =  4.5507022848e+10
pion_tot =  1.8577990000e+07
pele_tot =  1.0033128000e+08

# rho  = 3.4387408959e-05
# tion = 1.1236895508e+04 / 11604.55
# tele = 1.1360400391e+04 / 11604.55
# mfhe = 0.58300763369
# mfal = 0.41699239612
# eion_tot = 2.3226744832e+11
# eele_tot = 2.1015455334e+11
# pion_tot = 5.175598e+06
# pele_tot = 6.77272125e+05

eion_he = eion(heimx, rho*mfhe, tion)
eion_al = eion(alimx, rho*mfal, tion)
print( "EION:", eion_tot, mfhe*eion_he + mfal*eion_al)

eele_he = eele(heimx, rho*mfhe, tele)
eele_al = eele(alimx, rho*mfal, tele)
print( "EELE:",eele_tot, mfhe*eele_he + mfal*eele_al)

pele_he = pele(heimx, rho*mfhe, tele)
pele_al = pele(alimx, rho*mfal, tele)
print( "PELE:",pele_tot, pele_he + pele_al)

pion_he = pion(heimx, rho*mfhe, tion)
pion_al = pion(alimx, rho*mfal, tion)
print( "PION:",pion_tot, pion_he + pion_al)

print( "\nHELUM GAMC:")
gamc_he = gamc(heimx, rho*mfhe, tele, tion, 0.0, verbose=True)
print( "\nALUMINUM GAMC:")
gamc_al = gamc(alimx, rho*mfal, tele, tion, 0.0, verbose=True)



#         return opp.interpDT(var, opimx.dens, opimx.temps, 
#                             r, t,
#                             lookup=lookup, log=log,                          
#                             bcdmin=bcdmin, bcdmax=bcdmax,
#                             bctmin=bctmin, bctmax=bctmax)



# def gamc(rho, tele, tion, trad, 
#          bcdmin=BC_BOUND, bcdmax=BC_BOUND, 
#          bctmin=BC_BOUND, bctmax=BC_BOUND, 
#          log = False, verbose=False):


#     if bcdmin == BC_BOUND and rho  < opimx.dens[0] : rho  = opimx.dens[0]
#     if bctmin == BC_BOUND and tele < opimx.temps[0]: tele = opimx.temps[0]
#     if bctmin == BC_BOUND and tion < opimx.temps[0]: tion = opimx.temps[0]

#     if bcdmax == BC_BOUND and rho  > opimx.dens[-1] : rho  = opimx.dens[-1]
#     if bctmax == BC_BOUND and tele > opimx.temps[-1]: tele = opimx.temps[-1]
#     if bctmax == BC_BOUND and tion > opimx.temps[-1]: tion = opimx.temps[-1]

#     def interpDT(var, r, t, lookup):
#         return opp.interpDT(var, opimx.dens, opimx.temps, 
#                             r, t,
#                             lookup=lookup, log=log,                          
#                             bcdmin=bcdmin, bcdmax=bcdmax,
#                             bctmin=bctmin, bctmax=bctmax)

#     # Compute needed derivatives:
#     deidti = interpDT(opimx.eion, rho, tion, opp.INTERP_DFDT)
#     dpidti = interpDT(opimx.pion, rho, tion, opp.INTERP_DFDT)
#     dpidd  = interpDT(opimx.pion, rho, tion, opp.INTERP_DFDD)

#     dpidd_si = tion/(rho**2 * deidti) * dpidti**2 + dpidd

#     pion = interpDT(opimx.pion, rho, tion, opp.INTERP_FUNC)

#     deedte = interpDT(opimx.eele, rho, tele, opp.INTERP_DFDT)
#     dpedte = interpDT(opimx.pele, rho, tele, opp.INTERP_DFDT)
#     dpedd  = interpDT(opimx.pele, rho, tele, opp.INTERP_DFDD)
#     dpedd_se = tele/(rho**2 * deedte) * dpedte**2 + dpedd

#     pele = interpDT(opimx.pele, rho, tele, opp.INTERP_FUNC)

#     prad = AR*trad**4
#     ptot = prad + pion + pele
#     gamc = rho/ptot * (dpedd_se + dpidd_si + (4*prad)/(3*rho))

#     if verbose == True:
#         print( "deidti   = %17.10e" % deidti)
#         print( "dpidti   = %17.10e" % dpidti)
#         print( "dpidd    = %17.10e" % dpidd)
#         print( "dpidd_si = %17.10e" % dpidd_si)
#         print( "pion     = %17.10e" % pion)
        
#         print( "deedte   = %17.10e" % deedte)
#         print( "dpedte   = %17.10e" % dpedte)
#         print( "dpedd    = %17.10e" % dpedd)
#         print( "dpedd_se = %17.10e" % dpedd_se)
#         print( "pele     = %17.10e" % pele)
        
#         print( "prad     = %17.10e" % prad)
#         print( "ptot     = %17.10e" % ptot)
#         print( "gamc     = %17.10e" % gamc)

#     return gamc

# print( "\n\n")

# rho  = 1e-05
# tele = 0.025
# tion = 0.025
# trad = 0.0
# print( "\ngamc(rho=%g,tele=%g,tion=%g,trad=%g) = %g" % \
#     (rho, tele, tion, trad, gamc(rho,tele,tion,trad, verbose=True, bctmin=BC_EXTRAP_ZERO))

# print( "\n\n")

# rho  = 2.0e-07
# tele = 4e+05
# tion = 3e+05
# trad = 5.0
# print( "\ngamc(rho=%g,tele=%g,tion=%g,trad=%g) = %g" % \
#     (rho, tele, tion, trad, gamc(rho,tele,tion,trad, verbose=True))

# print( "\n")
# gamc(2e-07,2e+05,2e+05,5.0)
# print( "\n")
# gamc(2e-07,2e+05,3e+05,5.0)

# rho  = 6.646476174630e-08
# tion = 2.0e+05
# pion = opp.interpDT(opimx.pion, opimx.dens, opimx.temps,
#                     rho, tion,
#                     lookup=opp.INTERP_FUNC)
# print( "pion =", pion, rho/mpi * KB * tion)
