# The now documented opp_file.py!

import os
import sys

import tables as tb
import numpy as np

from .utils import isopacplot      # determines if file is an opacplot file

# pick filters from pytables
BASIC_FILTERS = tb.Filters(complevel=5, complib='zlib', shuffle=True, fletcher32=False)

# Generate a registrar dictionary to collect file types and respective parsers
FORMAT_PARSERS = {
    'ionmix': ('ionmix', 'parse'),
    }

class OppFile(object):
    """The opacplot file."""

	# determine if function is opacplot file. If file is not opacplot file,
	# generate object 'opp' with _parse (see below), make a properly labelled
	# file with .update, then throw away opp
    def __init__(self, filename, format=None, *args, **kwargs):
        if not isopacplot(filename, format):
            opp = self._parse(filename, format, *args, **kwargs)
            self.__dict__.update(opp.__dict__)
            del opp
            return 
        self._handle = tb.openFile(filename, mode='a', filters=BASIC_FILTERS)
        self.root = self._handle.root
        self.__file__ = filename

	# _parse, passed the name and format of a file and any args and kwargs,
	# determines if the format is in the registrar FORMAT_PARSERS. If the format
	# is in FORMAT_PARSERS, _parse proceeds to call _format_parse (see below)
	# and make object 'opp' assigned to the newly made opacplot file, then returns
	# opp. If the format of the file is not in FORMAT_PARSERS, _parse will attempt
	# to use _format_parse for each possible format in FORMAT_PARCERS. If the
	# format could not be read, _parse gives an error message.
    @classmethod
    def _parse(cls, filename, format, *args, **kwargs):
        if format in FORMAT_PARSERS:
            opp = cls._format_parse(filename, format, *args, **kwargs)
        else:
            for form in FORMAT_PARSERS.keys():
                try:
                    opp = cls._format_parse(filename, form, *args, **kwargs)
                    break
                except:
                    continue
            else:
                msg = "format of {0} could not be determined from available types ({1})"
                msg = msg.format(filename, ", ".join(FORMAT_PARSERS.keys()))
                raise RuntimeError(msg)
        return opp
	

	# _format_parse ensures '.' is in sys.path, whatever that means. It then
	# assigns formatmod and parserfunc to the file's format format and the 
	# format's parse function, respectively. It then assigns object opp to
	# the output of parserfunc and returns opp.
    @classmethod
    def _format_parse(cls, filename, format, *args, **kwargs):
        if '.' not in sys.path:
            sys.path.insert(0, '.')
        formatmod, parserfunc = FORMAT_PARSERS[format]
        formatmod = __import__(formatmod, globals(), locals(), fromlist=[None])
        parserfunc = getattr(formatmod, parserfunc)
        opp = parserfunc(filename, *args, **kwargs)
        return opp

# this does something
    def __del__(self):
        self._handle.close()
