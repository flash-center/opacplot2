#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import re
import gzip


def read_multi_eos(path):
    if path.endswith(".gz"):
        fopen = gzip.open
    else:
        fopen = open
    
    eos_tables = []
    FAST_HEADER_FMT = re.compile("^.*\sEOS\s.*$")
    FULL_HEADER_FMT = '\s(?P<matid>\d+)\s+(?P<table>[^0-9]+)(?P<Nr>[0-9E+-.]+)\s(?P<Nt>[0-9E+-.]+)$'
    MINUS_FMT = re.compile("(?<!E)-")

    #def parse

    with fopen(path, 'r') as f:
        for line in f:
            if re.match(FAST_HEADER_FMT, line):
                eos_tables.append({"header_str": line, "data_str": []})
            else:
                eos_tables[-1]['data_str'].append(re.sub(MINUS_FMT, ' -', line)) # replace "-0." by ' -0.'
    # process headers:
    for tab in eos_tables:
        m = re.match(FULL_HEADER_FMT, tab['header_str'])
        if m:
            h = {key: m.group(key).strip() for key in \
                          ['matid', 'table', 'Nr', 'Nt'] }
            tab['matid'], tab['inversed'],  tab['tabid'] = h['matid'][:4], h['matid'][4]=='1',  h['matid'][5:]
            Nr = tab['Nr'] = int(float(h['Nr']))
            Nt = tab['Nt'] = int(float(h['Nt']))
        else:
            raise ValueError('Failed to parse header line: {}'.format(tab['header_str']))
        data_str = ''.join(tab['data_str']).replace('\n','')
        data = tab['data'] = np.fromstring(data_str, sep=' ')
        
        del tab['data_str']
        del tab['header_str']
        tab['rho'] = data[:Nr]
        tab['unk'] = data[Nr:2*Nr] # unknown
        tab['eint'] = data[2*Nr:2*Nr+Nt]
        tab['pres'] = data[2*Nr+Nt:2*Nr+Nt+Nr*Nt].reshape((Nr,Nt), order='F')
        #tab['temp'] = data[2*Nr+Nt+Nr*Nt:].
        #print data[2*Nr+Nt+Nr*Nt:].reshape((Nr,Nt), order='C')

    return eos_tables





if __name__ == "__main__":
    eos_tables = read_multi_eos("/home/rth/luli/codes/multi90/Tables/Aluminium/3717_3eos.multi.gz")
    el =  eos_tables[0]
    Nr = el['Nr']
    Nt = el['Nt']
    Ntot = 2*Nr + Nt + 2*Nr*Nt
    print( el['rho'], el['eint'], el['pres'])

    print( el['data'].shape, Ntot)

    
