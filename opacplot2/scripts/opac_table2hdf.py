#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse

import os, os.path

import opacplot2 as opp
import re
import numpy as np
from scipy.constants import N_A
from scipy.constants import physical_constants
          
eV2K_cst = physical_constants['electron volt-kelvin relationship'][0]

def opac_table2hdf():

    parser = argparse.ArgumentParser(description= """
    This script is used to browse various EoS/Opacity tables formats
    """)
    avalable_formats = ['multi', 'ionmix', 'hdf5', 'propaceos']
    parser.add_argument('-i','--inftype',
            action="store", type=str,
            choices=avalable_formats,
            default='multi',
            help='Input filetype. Default: multi')
    parser.add_argument('-o','--outname',
            action="store", type=str,
            default='hdf5',
            help='Output file base name/path')
    parser.add_argument('--Znum',
            action="store", type=str,
            help='Comma separated list of Z num for every component ')
    parser.add_argument('--Xnum',
            action="store", type=str,
            help='Comma separated list of X num fractions for every component')
    parser.add_argument('input_file',
            action="store", type=str,
            help='Input file')
    args = parser.parse_args()

    basedir, filename = os.path.split(os.path.abspath(args.input_file))
    basename_in, _ = os.path.splitext(filename)

    if args.outname is None:
        filename_out, _ = os.path.splitext(os.path.abspath(args.outname))
    else:
        filename_out = os.path.join(basedir, basename_in)
    if args.Znum is not None:
        Znum = [int(el) for el in args.Znum.split(',')]
    else:
        Znum = None
    if args.Xnum is not None:
        Xnum = [float(el) for el in args.Xnum.split(',')]
    else:
        Xnum = None

    # ====================== Parsing input file ==========================
    if args.inftype == 'propaceos':
        from ..opg_propaceos import OpgPropaceosAscii
        op = OpgPropaceosAscii(args.input_file)
        op.write2hdf(filename_out+'.h5')
    elif args.inftype == 'multi':
        if Znum is None:
            raise ValueError('Znum parameter should be provided!')
        if Xnum is None and len(Znum)>1:
            raise ValueError('Xnum parameter should be provided when len(Znum)>1 !')
        from ..opg_multi import OpgMulti, MULTI_EXT_FMT
        patrn = re.compile(MULTI_EXT_FMT, re.IGNORECASE)
        basename_in = re.sub(patrn, '', filename) # strip the extension
        op = OpgMulti.fromfile(basedir, basename_in)
        op.write2hdf(filename_out+'.h5', Znum=Znum, Xnum=Xnum)
