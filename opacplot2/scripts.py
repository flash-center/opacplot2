#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import os, os.path

import opacplot2 as opp

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




