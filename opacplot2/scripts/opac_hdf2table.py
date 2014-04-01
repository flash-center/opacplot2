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

def opac_hdf2table():

    parser = argparse.ArgumentParser(description= """
    This script is used to browse various EoS/Opacity tables formats
    """)
    avalable_formats = ['multi', 'ionmix']
    parser.add_argument('-t','--ftype',
            action="store", type=str,
            choices=avalable_formats,
            default='multi',
            help='Output filetype. Default: multi')
    parser.add_argument('-o', '--outfile',
            action="store", type=str,
            help='Output file base name/path')
    parser.add_argument('input_file',
            action="store", type=str,
            help='Path to the input hdf5 file')
    args = parser.parse_args()

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
