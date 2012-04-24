#!/usr/bin/env python

import os
import sys
import json
import subprocess 

from distutils.core import setup

if __name__ == "__main__":
    setup(name="opacplot2",
          version='0.1',
          description='Package for reading, manipulating, and plotting EOS and Opacity data',
          author='Milad',
          author_email='milad@flash.uchicago.edu',
          url='http://flash.uchicago.edu/',
          packages=['opacplot2']
          )
