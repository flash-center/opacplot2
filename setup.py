#!/usr/bin/env python

import os
import sys
import json
import subprocess 

#from distutils.core import setup
from setuptools import setup


if __name__ == "__main__":
    setup(name="opacplot2",
          version='0.1',
          description='Package for reading, manipulating, and plotting EOS and Opacity data',
          author='Milad',
          author_email='milad@flash.uchicago.edu',
          url='http://flash.uchicago.edu/',
          packages=['opacplot2', 'opacplot2.presets'],
          entry_points = {
                    'console_scripts': ['opacplot2 = opacplot2.scripts:main',
                                        'opac_table2hdf = opacplot2.scripts.opac_table2hdf:opac_table2hdf',
                                        'opac_hdf2table = opacplot2.scripts.opac_hdf2table:opac_hdf2table',
                                        'opacdump = opacplot2.scripts.main:opacdump',
                                        'opacdiff = opacplot2.scripts.opac_diff:opac_diff',
                                        'opac_checkhdf = opacplot2.scripts.main:opac_checkhdf',
                                        ],
                        }
          )
