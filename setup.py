#!/usr/bin/env python

import os

from setuptools import setup, find_packages


if __name__ == "__main__":
    setup(name="opacplot2",
          version='1.0.0',
          description='Package for reading, manipulating, and plotting EOS and Opacity data',
          author='Milad',
          author_email='milad@flash.uchicago.edu',
          url='http://flash.uchicago.edu/',
          packages=find_packages(),
          package_data={'opacplot2': [os.path.join('tests','data', '*')]},
          test_suite="opacplot2.tests.run",
          install_requires=[
              "numpy >= 1.6",
              "tables >= 3.0",
              "six >= 1.6",
              "setuptools >= 18.0",
              "periodictable >=1.4.1"
            ],

          entry_points = {
                    'console_scripts': ['opac-convert = opacplot2.scripts.opac_convert:convert_tables',
                                        'opac-error = opacplot2.scripts.opac_error:check_error',
                                        'sesame-extract = opacplot2.scripts.sesame_extract:extract_tables'
                                        ],
                        }
          )
