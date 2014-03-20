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

def opac_convert():

    parser = argparse.ArgumentParser(description= """
    This script is used to browse various EoS/Opacity tables formats
    """)
    avalable_formats = ['multi', 'ionmix', 'hdf5', 'propaceos']
    parser.add_argument('-i','--inftype',
            action="store", type=str,
            choices=avalable_formats,
            default='multi',
            help='Input filetype. Default: multi')
    parser.add_argument('-o','--outfype',
            action="store", type=str,
            choices=avalable_formats,
            default='hdf5',
            help='Output filetype. Default: hdf5')
    parser.add_argument('-b', '--basename',
            action="store", type=str,
            help='Output file base name/path')
    parser.add_argument('--Abar',
            action="store", type=float,
            help='Atomic mass')
    parser.add_argument('input_file',
            action="store", type=str,
            help='Input file')
    args = parser.parse_args()

    basedir, filename = os.path.split(os.path.abspath(args.filename))
    basename_in, _ = os.path.splitext(filename)

    # ====================== Parsing input file ==========================
    if args.inftype == 'propaceos':
        from ..opg_propaceos import OpgPropaceosAscii
        op = OpgPropaceosAscii(args.filename)
    elif args.inftype == 'multi':
        from ..opg_multi import OpgMulti
        patrn = re.compile(MULTI_EXT_FMT, re.IGNORECASE)
        basename_in = re.sub(patrn, '', filename) # strip the extension
        op = OpgMulti.fromfile(basedir, basename_in)
    elif agrs.inftype == 'hdf5':
        from ..opg_hdf5 import OpgHdf5






    # ====================================================================
    if args.basename is None:
        filename_out = args.basename
    else:
        filename_out = os.path.join(basedir, basename_in)
    # ====================== Writing to output format ====================
    if args.ftype == 'multi':
        from ..opg_multi import OpgMulti, MULTI_EXT_FMT

        op_multi = OpgMulti(**op)
        op_multi.set_id(3717) # just a random number, shouldn't matter
        print op_multi.keys()
        op_multi.write(os.path.join(basedir, basename + '-prp'), floor=1e-7)
