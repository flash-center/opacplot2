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
import sys
import numpy as np
import matplotlib.pyplot as plt

def opac_diff():
    from ..opg_multi import MULTI_EXT_FMT, OpgMulti
    from ..eos_plotter import plot_diff_mg_opac, plot_2D_map, plot_Zbar
    from hedp.rad import planck
    import itertools

    parser = argparse.ArgumentParser(description= """
    This script is used to browse various EoS/Opacity tables formats
    """)
    parser.add_argument('-f','--force',
            action="store_true",
            help='Overwrite files', required=False)
    parser.add_argument('tables',
            metavar='table',
            action="store", type=str,
            help='Database type. Default: sesame.',
            nargs='+')
    parser.add_argument('-s','--stride',
            action="store",
            default=4,
            type=int,
            help='Stride between plots on the density & temperature grid.', required=False)
    parser.add_argument('-o','--out_dir',
            action="store",
            help='Output directory', required=True)
    args = parser.parse_args()

    if len(args.tables)==1:
        print('One opacity table given only! Nothing to compare with!')
        return
    elif len(args.tables)>3:
        print('More then 2 tables given! It means that error values will be compared between the first 2 tables only!')

    op = []
    alphabet = iter('ABCDEFGH')
    for table_path in args.tables:
        basedir, filename = os.path.split(os.path.abspath(table_path))
        patrn = re.compile(MULTI_EXT_FMT, re.IGNORECASE)
        basename = re.sub(patrn, '', filename) # strip the extension
        op += [OpgMulti.open_file(basedir, basename)]
        op[-1]['label'] = basename
        op[-1]['basedir'] = basedir
        op[-1]['key'] = next(alphabet)


    rho = op[0]['rho']
    temp = op[0]['temp']

    if os.path.isdir(args.out_dir):
        if not args.force:
            print("Directory {0} already exists, use -f flag to overwrite!".format(args.out_dir))
    else:
        os.mkdir(args.out_dir)

        os.mkdir(os.path.join(args.out_dir, '1D_spectra'))

#    for iop in op:
#        print iop['label'],':'
#        for var in ['opp_mg', 'opr_mg', 'eps_mg', 'zbar']:
#            print var, iop[var].min(), iop[var].max()


    with open(os.path.join(args.out_dir, 'README.txt'), 'w') as f:
        w = f.write
        w("This folder contains plots showing difference\n \
           between several multigroup opacity files:\n\n""")

        for el in op:
            w(' {0}: {2:10} => {1}\n'.format(el['key'], el['basedir'], el['label']))

        arr2str = lambda x: np.array2string(x, precision=3, max_line_width=130)

        w('\nDensity grid [g.cm⁻²]: ({0} pts)\n'.format(len(op[0]['rho'])))
        w(arr2str(op[0]['rho']))

        w('\nTemperature grid [eV]: ({0} pts)\n'.format(len(op[0]['temp'])))
        w(arr2str(op[0]['temp']))

        w('\nRadiation grid [eV]: ({0} pts)\n'.format(len(op[0]['groups'])))
        w(arr2str(op[0]['groups']))


    print('Plotting Zbar - done.')
    fig = plt.figure(figsize=(12,4))
    plot_Zbar(fig, op)
    fig.savefig(os.path.join(args.out_dir, 'Zbar.png'), bbox_inches='tight')



    print('Generating spectra ', end='')
    for idx0 in np.arange(0,len(rho), args.stride):
        for idx1 in np.arange(0,len(temp), args.stride):
            fig = plt.figure(figsize=(10,6))
            plot_diff_mg_opac(fig, op, idx=(idx0, idx1))
            fig.savefig(os.path.join(args.out_dir, '1D_spectra','diff_mg_spec_{0}_{1}.png'.format(idx0, idx1)),
                    bbox_inches='tight')
            plt.close(fig)
            sys.stdout.write('.')
            sys.stdout.flush()

    print('')
    print('Generating 2D error maps ')
    for op0, op1 in itertools.combinations(op, 2):
        subdir = '{0}vs{1}'.format(op0['key'], op1['key'])
        subdir_full = os.path.join(args.out_dir, subdir)
        if not os.path.isdir(subdir_full):
            os.mkdir(subdir_full)
        print('  -', subdir, ' ', end='')
        for temp_val in [10, 40, 100, 400, 1000, 2000]:
            fig = plt.figure(figsize=(20,5))
            plot_2D_map(fig, [op0, op1], temp_val)
            fig.savefig(os.path.join(args.out_dir, subdir,'opac_2D_error_{0}_{1}eV.png'.format(
                                          subdir, int(temp_val))), bbox_inches='tight')
            plt.close(fig)
            sys.stdout.write('.')
            sys.stdout.flush()
        print('')











