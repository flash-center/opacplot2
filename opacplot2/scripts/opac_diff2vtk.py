#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse

import os, os.path

import opacplot2 as opp
import re
import numpy as np
from scipy.constants import N_A
from scipy.constants import physical_constants
from ..opg_hdf5 import OpgHdf5

eV2K_cst = physical_constants['electron volt-kelvin relationship'][0]

def opac_diff2vtk():

    parser = argparse.ArgumentParser(description= """
    This script is used to compare multiple hdf5 opacity tables of the same material.
    All radiative grids are projected 
    """)
    parser.add_argument('-o', '--outdir',
            action="store", type=str,
            default='./',
            help='Output directory')
    parser.add_argument('opacity_files',
            action="store", type=list,
            help='List of opacity files. The radiative grid\
of the first file is used for projection\
when applicable.')
    args = parser.parse_args()

    # for 
    basedir, filename = os.path.split(os.path.abspath(args.input_file))
    basename_in, _ = os.path.splitext(filename)



    if args.outfile is not None:
        filename_out = args.outfile
    else:
        filename_out = os.path.join(basedir, basename_in)

    op = OpgHdf5.open_file(os.path.abspath(args.input_file))




    # ====================== Writing to output format ====================
    if args.ftype == 'multi':
        from ..opg_multi import OpgMulti, MULTI_EXT_FMT

        op = OpgMulti(**op)
        op.set_id(3717) # just a random number, shouldn't matter
        op.write(os.path.join(basedir, filename_out), floor=None)
    elif args.ftype == 'ascii':
        outfile = os.path.join(basedir, filename_out)+'.txt'
        def repr_grid(arr, label):
            out = []
            out += ['='*80]
            out += ['     {0}: {1} points'.format(label, len(arr))]
            out += ['='*80]
            out = '\n'.join(out)
            return out + '\n'+ np.array2string(arr, precision=3, separator='')
        print repr_grid(op['dens'][:], "Density grid [g/cc]")
        print repr_grid(op['temp'][:], 'Temperature grid [eV]')
        print repr_grid(op['idens'][:], "Ionic density grid [1/cc]")
        #print repr_grid(op['groups'][:], "Photon energy groups [eV]")
        nu = op['groups'][:]
        nu = 0.5*(nu[1:]+nu[:-1])
        out_op =  np.array([nu]+[op[key+'_mg'][0,0]*op['dens'][0] for key in ['opp', 'opr', 'emp']])
        np.set_printoptions(threshold=1e9)
        print repr_grid(out_op.T, 'Opacity: nu [eV], opp [1/cm], opr[1/cm], emp [??]')



