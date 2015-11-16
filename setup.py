#!/usr/bin/env python

import os

from setuptools import setup, find_packages


if __name__ == "__main__":
    setup(name="opacplot2",
          version='0.9.0',
          description='Package for reading, manipulating, and plotting EOS and Opacity data',
          author='Milad',
          author_email='milad@flash.uchicago.edu',
          url='http://flash.uchicago.edu/',
          packages=find_packages(),
          package_data={'opacplot2': [os.path.join('tests','data', '*')]},
          test_suite="opacplot2.tests.run",

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
