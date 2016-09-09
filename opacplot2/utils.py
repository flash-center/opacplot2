import os
import re
import random

import tables as tb

def munge_h5filename(filename, h5filename):
    """Gets an opacplot hdf5 filename from existing names."""
    if not isinstance(h5filename, basestring):
        if isinstance(filename, basestring): 
            h5filename = filename.rpartition(".")[0] + '.h5'
        else:
            h5filename = "opp.h5"
    return h5filename


def isopacplot(filename, format=None):
    """Determines if a file is in of opacplot format."""
    return (format in [None, "opp"]) and isinstance(filename, basestring) and \
            ((os.path.isfile(filename) and tb.isHDF5File(filename)) or \
              not os.path.exists(filename))


def randomize_ionmix(filename, outfilename):
    """Randomizes the data from an existing ionmix file and rewrites it 
    to the outfile."""
    basepttn = "\d{6}"
    exppttn = "E[+-]\d{2}"
    rand_base = lambda mobj: "{0:06}".format(random.randint(0, 999999))
    rand_exp = lambda mobj: "E{0:+03}".format(random.randint(-99, 99))

    with open(filename) as f:
        lines = f.readlines()

    newlines = []
    for line in lines:
        newline = re.sub(basepttn, rand_base, line)
        newline = re.sub(exppttn, rand_exp, newline)
        newlines.append(newline)

    with open(outfilename, 'w') as f:
        f.write("".join(newlines))

