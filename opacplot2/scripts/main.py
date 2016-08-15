#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import absolute_import
#from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals

import argparse
import os, os.path

import opacplot2 as opp
import re
import numpy as np
from scipy.constants import N_A
from scipy.constants import physical_constants
          
eV2K_cst = physical_constants['electron volt-kelvin relationship'][0]


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
        print("Parsing sesame input files.")
        eos_sesame = opp.OpgSesame(os.path.join(args.dbdir, "xsesame_ascii"),
                                opp.OpgSesame.SINGLE, verbose=False)
        cond_sesame = opp.OpgSesame(os.path.join(args.dbdir, "sescu_ascii"),
                                opp.OpgSesame.DOUBLE, verbose=False)
        eos_data_i = eos_sesame.data[args.tablenum]

        cond_keys = sorted([key for key, val in cond_sesame.data.iteritems()\
                            if val['zmax'] == eos_data_i['zmax']])

        # by default select the last (i.e. newest) table avalable
        cond_data_i = cond_sesame.data[cond_keys[-1]]

        print('Presets for sesame table', args.tablenum, end='')

        # Merging ele and ion grids
        if args.tablenum in opp.presets.SESAME:
            eos_data = opp.adapt.EosMergeGrids(eos_data_i,
                                    **opp.presets.SESAME[args.tablenum]['merge'])
            fracs = opp.presets.SESAME[args.tablenum]['fracs']
            print('found')
        else:
            eos_data = opp.adapt.EosMergeGrids(eos_data_i,
                       filter_dens=lambda x: x>0,
                       filter_temps=lambda x: x>0)

            fracs = (1.0,)
            print('not found')

        output_file = args.out+'.cn4'
        print('Generating IONMIX file: {0}'.format(output_file))
        numDens = opp.NA * eos_data['ele_dens'] / eos_data["abar"]
        print('Warning: for now zbar is set to 0.')
        opp.writeIonmixFile(output_file, (eos_data['zmax'],), fracs,
                        numDens=numDens, temps=eos_data['ele_temps'],
                        eion=eos_data["ion_eint"],
                        eele=eos_data["ele_eint"],
                        pion=eos_data["ion_pres"],
                        pele=eos_data["ele_pres"])
        print('Trying to read if back... ', end='')
        try:
            ionmix = opp.OpacIonmix(output_file, eos_data["abar"]/opp.NA, 
                            twot=True, man=True, verbose=False)
            print('[ok]')
        except:
            print('[failed]')

def opacdump():

    parser = argparse.ArgumentParser(description= """
    This script is used to browse various EoS/Opacity tables formats
    """)
    parser.add_argument('-t','--ftype',
            action="store", type=str,
            choices=['multi', 'propaceos', 'grid_solid', 'grid_gas', 'grid_log'],
            default='multi',
            help='Database type. Default: sesame.')
    parser.add_argument('-f','--filename',
            action="store", type=str,
            help='Filename')
    parser.add_argument('--Abar',
            action="store", type=float,
            help='Atomic mass of the element')
    parser.add_argument('--rho0',
            action="store", type=float,
            help='Refernce density for a solid [g/cc]')
    parser.add_argument('value',
            action="store", type=str,
            help='Parameters to print',
            choices=['grid2prp', 'grid2file', 'grid2ascii'],
            default=None)
    args = parser.parse_args()

    if args.ftype in ['grid_solid', 'grid_gas', 'grid_log']:
        groups = np.array([])
        if args.ftype == 'grid_solid' and args.rho0 is None:
            raise ValueError('--rho0 argument should be provided with "grid_solid" type!')
        from eospac.base import GridBase
        kind = args.ftype.split('_')[1]
        grid = GridBase(rho_ref=args.rho0, kind=kind)
        rho = grid.rho_grid
        temp = grid.temp_grid/eV2K_cst


    else:
        if args.Abar is None:
            raise ValueError
        if args.filename is None:
            raise ValueError
        if args.ftype == 'multi':
            from ..opg_multi import MULTI_EXT_FMT, OpgMulti
            FULL_PATH = os.path.abspath(args.filename)
            BASE_PATH, fname = os.path.split(FULL_PATH)
            patrn = re.compile(MULTI_EXT_FMT, re.IGNORECASE)
            fname_base = re.sub(patrn, '', fname) # strip the extension
            tab = OpgMulti.open_file(BASE_PATH, fname_base)
        temp = tab['temp']
        rho = tab['rho']
        groups = tab['groups']
    if args.Abar is not None:
        nele = rho*N_A/args.Abar

    if args.value == 'grid2prp':
        #print tab.keys()
        print("="*80)
        from ..opg_propaceos import OpgPropaceosGrid
        temp = np.unique(np.fmax(0.01, temp))
        nele = np.unique(np.fmax(1e10, nele))
        nele = np.unique(np.fmin(5e25, nele))
        print(OpgPropaceosGrid.format_grid1(nele, temp, groups))
        print("="*80)
        print('+'*80)
        print('='*80)
        print(OpgPropaceosGrid.format_grid2(nele, temp))
    elif args.value == 'grid2file':
        outname_base = os.path.join('./', fname_base)
        np.savetxt(outname_base+'.nele-grid.txt', nele,
                header='Electron density grid [1/cc]\nsource: {0}'.format(FULL_PATH))
        np.savetxt(outname_base+'.temp-grid.txt', temp,
                header='Plasma temperature grid [eV]\nsource: {0}'.format(FULL_PATH))
        np.savetxt(outname_base+'.nu-grid.txt', groups,
                header='Photon energy groups grid [eV]\nsource: {0}'.format(FULL_PATH))
    elif args.value == 'grid2ascii':
        def repr_grid(arr, label):
            out = []
            out += ['='*80]
            out += ['     {0}: {1} points'.format(label, len(arr))]
            out += ['='*80]
            out = '\n'.join(out)
            return out + '\n'+ np.array2string(arr, precision=3, separator=',')
        print(repr_grid(rho, "Density grid [g/cc]"))
        print(repr_grid(temp, 'Temperature grid [eV]'))
        if args.Abar is not None:
            print( repr_grid(nele, "Ionic density grid [1/cc]") ) 
        if len(groups):
            print(repr_grid(groups, "Photon energy groups [eV]") )

def opac_checkhdf():
    parser = argparse.ArgumentParser(description= """
    Check consistency of hdf5 opacity file
    """)
    parser.add_argument('input_file',
            action="store", type=str,
            help='Input hdf5 file ')
    args = parser.parse_args()
    from ..opg_hdf5 import OpgHdf5

    op = OpgHdf5.open_file(args.input_file)
    op.run_testsuite('short')

