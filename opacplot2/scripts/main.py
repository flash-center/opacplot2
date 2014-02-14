#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import os, os.path

import opacplot2 as opp
import re
import numpy as np

def main():

    parser = argparse.ArgumentParser(description= """
    This script is used to automate EoS/opacity tables generation for FLASH.
    """)
    parser.add_argument('-d','--dbdir',
            action="store", type=str,
            default = os.path.abspath(os.curdir),
            help='Path to the database. Default: current directory.')
    parser.add_argument('-t','--dbtype',
            action="store", type=str,
            default='sesame',
            help='Database type. Default: sesame.')
    parser.add_argument('-n','--tablenum',
            action="store", type=int,
            help='Table id', required=True)
    parser.add_argument('-o','--out',
            action="store", type=str,
            help='Ouput filename',
            default=None)
    args = parser.parse_args()
    if args.out is None:
        args.out = '{0}-eos-{1}'.format(args.dbtype, args.tablenum)
    else:
        args.out = os.path.splitext(args.out)[0]

    if args.dbtype == 'sesame':
        print "Parsing sesame input files."
        eos_sesame = opp.OpgSesame(os.path.join(args.dbdir, "xsesame_ascii"),
                                opp.OpgSesame.SINGLE, verbose=False)
        cond_sesame = opp.OpgSesame(os.path.join(args.dbdir, "sescu_ascii"),
                                opp.OpgSesame.DOUBLE, verbose=False)
        eos_data_i = eos_sesame.data[args.tablenum]

        cond_keys = sorted([key for key, val in cond_sesame.data.iteritems()\
                            if val['zmax'] == eos_data_i['zmax']])

        # by default select the last (i.e. newest) table avalable
        cond_data_i = cond_sesame.data[cond_keys[-1]]

        print 'Presets for sesame table', args.tablenum,

        # Merging ele and ion grids
        if args.tablenum in opp.presets.SESAME:
            eos_data = opp.adapt.EosMergeGrids(eos_data_i,
                                    **opp.presets.SESAME[args.tablenum]['merge'])
            fracs = opp.presets.SESAME[args.tablenum]['fracs']
            print 'found'
        else:
            eos_data = opp.adapt.EosMergeGrids(eos_data_i,
                       filter_dens=lambda x: x>0,
                       filter_temps=lambda x: x>0)

            fracs = (1.0,)
            print 'not found'

        output_file = args.out+'.cn4'
        print 'Generating IONMIX file: {0}'.format(output_file)
        numDens = opp.NA * eos_data['ele_dens'] / eos_data["abar"]
        print 'Warning: for now zbar is set to 0.'
        opp.writeIonmixFile(output_file, (eos_data['zmax'],), fracs,
                        numDens=numDens, temps=eos_data['ele_temps'],
                        eion=eos_data["ion_eint"],
                        eele=eos_data["ele_eint"],
                        pion=eos_data["ion_pres"],
                        pele=eos_data["ele_pres"])
        print 'Trying to read if back... ',
        try:
            ionmix = opp.OpacIonmix(output_file, eos_data["abar"]/opp.NA, 
                            twot=True, man=True, verbose=False)
            print '[ok]'
        except:
            print '[failed]'

def opacdump():

    parser = argparse.ArgumentParser(description= """
    This script is used to browse various EoS/Opacity tables formats
    """)
    parser.add_argument('-t','--ftype',
            action="store", type=str,
            choices=['multi', 'propaceos'],
            default='multi',
            help='Database type. Default: sesame.')
    parser.add_argument('-f','--filename',
            action="store", type=str,
            help='Filename', required=True)
    parser.add_argument('--Abar',
            action="store", type=float,
            help='Filename', required=True)
    parser.add_argument('value',
            action="store", type=str,
            help='Parameters to print',
            choices=['grid2prp', 'grid2file'],
            default=None)
    args = parser.parse_args()

    if args.ftype == 'multi':
        from ..opg_multi import MULTI_EXT_FMT, OpgMulti
        FULL_PATH = os.path.abspath(args.filename)
        BASE_PATH, fname = os.path.split(FULL_PATH)
        patrn = re.compile(MULTI_EXT_FMT, re.IGNORECASE)
        fname_base = re.sub(patrn, '', fname) # strip the extension
        tab = OpgMulti.fromfile(BASE_PATH, fname_base)
    from ..opg_propaceos import OpgPropaceosGrid
    from scipy.constants import N_A
    tele = tab['temp']
    nele = tab['rho']*N_A/args.Abar

    if args.value == 'grid2prp':
        print tab.keys()
        print "="*80
        print OpgPropaceosGrid.format_grid1(nele, tele, tab['groups'])
        print "="*80
        print '+'*80
        print '='*80
        print OpgPropaceosGrid.format_grid2(nele, tele)
    elif args.value == 'grid2file':
        outname_base = os.path.join('./', fname_base)
        np.savetxt(outname_base+'.nele-grid.txt', nele,
                header='Electron density grid [1/cc]\nsource: {0}'.format(FULL_PATH))
        np.savetxt(outname_base+'.temp-grid.txt', tele,
                header='Plasma temperature grid [eV]\nsource: {0}'.format(FULL_PATH))
        np.savetxt(outname_base+'.nu-grid.txt', tab['groups'],
                header='Photon energy groups grid [eV]\nsource: {0}'.format(FULL_PATH))


def opac_prp2file():
    from ..opg_propaceos import OpgPropaceosAscii

    parser = argparse.ArgumentParser(description= """
    This script is used to browse various EoS/Opacity tables formats
    """)
    parser.add_argument('-t','--ftype',
            action="store", type=str,
            choices=['multi', 'ionmix'],
            default='multi',
            help='Database type. Default: sesame.')
    parser.add_argument('-f','--filename',
            action="store", type=str,
            help='Filename', required=True)
    args = parser.parse_args()
    op = OpgPropaceosAscii(args.filename)
    basedir, filename = os.path.split(os.path.abspath(args.filename))
    basename, _ = os.path.splitext(filename)
    if args.ftype == 'multi':
        from ..opg_multi import OpgMulti

        op_multi = OpgMulti(**op)
        op_multi.set_id(3717)
        print op_multi.keys()
        op_multi.write(os.path.join(basedir, basename + '-prp'), floor=1e-7)
