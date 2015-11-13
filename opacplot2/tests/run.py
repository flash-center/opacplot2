#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import os
import nose

_base_dir, _ = os.path.split(__file__)

def run(coverage=False):
    argv=['', '-s', '--where={}'.format(_base_dir), '--verbosity=2']
    if coverage:
        argv += ['--with-coverage', '--cover-package=opacplot2']
    result = nose.run(argv=argv)
    status = int(not result)
    return status


if __name__ == '__main__':
    status = run()
    print('Exit status: {}'.format(status))
    sys.exit(status)




