#!/usr/bin/python
# -*- coding: utf-8 -*-
import re
import numpy as np


FLOAT_REGEXP = r'[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'

class OpalComposition(dict):
    def __init__(self, line):
        """
        Parse OPAL table hearder to extract composition
        """
        regexp  = ".*\s+(?P<source>)\s+X=(?<X>{0}) Y=(?<Y>{0}) Z=(?<Z>{0}) dXc=(?<dXc>{0}) dXo=(?<dXo>{0})".format(FLOAT_REGEXP)
        res = re.match(regexp, line)
        if res:
            return res
        else:
            raise ValueError("Couldn't parse line: {0}".format(line))



class OpalTable(dict):
    """
    This is a dictionary derived class that regroups
    """
