#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import os
import pytest

_base_dir, _ = os.path.split(__file__)

def run(coverage=False):
    argv=['-x', '{}'.format(_base_dir), '-v']
    if coverage:
        argv += ['--cov-report', 'term-missing', '--cov=opacplot2']
    result = pytest.main(argv)
    status = int(result)
    return status

def run_cli(coverage=False):
    status = run(coverage=coverage)
    print('Exit status: {}'.format(status))
    sys.exit(status)



if __name__ == '__main__':
    run_cli(coverage=True)




