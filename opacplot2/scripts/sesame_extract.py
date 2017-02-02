#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
sesame_extract.py Extract single table from the SESAME ASCII file and write into a separate file 

Created by JT Laune on Jan 26, 2017 for the Flash Center for Computational Science
  for use in managing the Flash Center's SESAME tables, and for development of opacplot2.

Simple enough to be hopefully functional in both Python 2.x and Python 3.x.

Example usage:
  # Check syntax possibilities
  python sesame_extract.py -h
  
  # Extract water (SESAME element 7150 - see the publicly available PDF) information
  #   from the 160 MB SESAME ASCII file, and save it into "water.ses"
  python sesame_extract.py -o "./water.ses" "$HOME/SESAME/sesame-ec/sesame_ascii" 7150

CHANGELOG:
  2017-01-27 Tested by Scott Feister in Python 2.7.12; successfully extracted H2O from sesame_ascii
             Added lots of comment lines -SKF

TODO:
  * Complementary function would call "Opacplot2" to print a summary of file contents
"""

import argparse

def get_input_data():
    """ Parse command syntax for arguments, and provide help ('-h') option. """
    parser = argparse.ArgumentParser(
                description='This script is used to extract single '
                            'SESAME tables from the SESAME database.')
    parser.add_argument('database', action='store', type=str,
                        help='Database file name.')
    parser.add_argument('tabnum', action='store', type=str,
                        help='SESAME table number.')
    parser.add_argument('-o', '--output', action='store', type=str,
                        help='Output file name.')

    args = parser.parse_args()

    if args.output is None:
        args.output = '{}_{}.ses'.format(args.database[:-4], args.tabnum)
    return args

def find_table(f_out, fhand, tabnum):
    """ Seek file cursor position to the material of interest """
    while True:
        line = fhand.readline()
        if not line: raise Warning('Table not found') # Reached EOF.
        if line[:3] == " 2 ": raise Warning('Table not found') 
        if line.split()[0] == '0':
            if str(line.split()[1]) == tabnum:
                f_out.write(line)
                break # Advanced fhand far enough

def write_entries(f_out, fhand):
    """ Copy lines into output file until another material is reached """
    while True:
        line = fhand.readline()
        if not line: break # Reached EOF.
        if line.split()[0] == '0': break # Reached another material.
        f_out.write(line)

def extract_tables():
    """ Extract single table from the SESAME ASCII file and write into a separate file """ 
    args = get_input_data()
    with open(args.output, 'w') as f_out:
        with open(args.database, 'r') as fhand:
            find_table(f_out, fhand, args.tabnum)
            write_entries(f_out, fhand)


if __name__=='__main__':
    extract_tables()
