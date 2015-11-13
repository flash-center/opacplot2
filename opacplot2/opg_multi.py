#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import absolute_import
#from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals

import os
import os.path
import numpy as np
import re
import gzip
import periodictable as ptab
from scipy.constants import N_A

MULTI_EXT_FMT = r'.(?P<ext>[oe]p[ezprs])(?:\.gz)?'

def get_related_multi_tables(folder, base_name, verbose=False):
    """
    Get all related multi Tables defined by a folder and a base name.

    Parameters:
    -----------
      - folder [str]: folder containing the tables
      - base_name [str]: base name of the table
    """
    patrn = re.compile((r'{0}'+MULTI_EXT_FMT).format(base_name), re.IGNORECASE)
    folder = os.path.abspath(folder)

    res =  {re.match(patrn, fname).group('ext').lower(): os.path.join(folder, fname)\
                for fname in os.listdir(folder)\
                if re.match(patrn, fname)}
    if not res:
        raise ValueError("Couln't find any tables in MULTIv5 format!")
    else:
        if verbose:
            valid_keys = [key for key in ['opp', 'opr', 'opz', 'eps'] if key in res]
            print("Found {0} files of table {1}!".format(valid_keys, base_name))
        if 'ope' in res:
            res['eps'] = res['ope']
            del res['ope']
        return res


FMT = " 13.8E"
PFMT = '%'+FMT


class OpgMulti(dict):
    """
    Can be used either to parse or to write MULIv5 tables.

    Parameters:

     * table - two use cases are possible
       1. (folder, base_name) 
       2. dict
    """

    def __init__(self, *cargs, **vargs):
        self.table_name = {}
        super(OpgMulti, self).__init__(*cargs, **vargs)

    _op_labels = dict(opp='PLANCK M', opr='ROSSELAND M', eps='EPS M ', opz='')

    @classmethod

    def open_file(cls, folder, base_name):
        """
        Parse MULTI format from a file
        Parameters:
        -----------
          - folder
          - base name

        Returns:
        --------
          OpgMulti
        """
        op = cls()
        table =  get_related_multi_tables(folder, base_name)
        print('Parsing opacity tables {0}'.format(re.sub(MULTI_EXT_FMT, '', table.items()[0][1])))
        print(' '.join([' ']*10, end='')
        for tabletype, path in table.items():
            op._parse(path, tabletype.lower())
            print('...',tabletype, end='')
        print('')
        return op

    def set_id(self, _id):
        _id = str(_id)
        self.table_name = {'opz': _id+'2017', 'opp': _id+'3017',
                                'eps': _id+'5017', 'opr': _id+'4017'}

    def _parse(self, path, tabletype):
        """
        Parse MULTIv5 opacity/ionization file
        Parameters:
        -----------
          - path [str]: path to the file
          - tabletype [str]: table type, one of ('opp', 'opr', 'eps', 'opz')

        """
        if os.path.splitext(path)[1] == '.gz':
            f = gzip.open(path)
        else:
            f = open(path)
        idx = []
        i = 0
        # get table name and opacity type (ex: "1232323", "PLANCK M")
        if tabletype == 'opz':
            table_name_pattern = re.compile(r'^\s(\d+)\b\s+')
        else:
            table_name_pattern = re.compile(r'(\d+)\b\s*(\w+\s\w)')
        for line in f:
            #if line.count('PLANCK M') or\
            #        line.count('ROSSE M') or line.count('EPS M'):
            search_result = table_name_pattern.search(line)
            if search_result:
                idx.append(i)
                self.table_name[tabletype] = search_result.groups()[0]
            i += 1
        max_idx = i
        f.seek(0)
        if tabletype == 'opz':
            idx = [0]
        #    self.table_name[tabletype] = ''
        [out, groups, rho, temp] = [[]]*4 # initialize
        for i in np.arange(len(idx)):
            header = f.readline()
            [r_len, T_len] = list(map(int,
                              list(map(float, [header[30:45], header[45:60]]))))
            if tabletype != 'opz':
                group_info = f.readline()

                [l_group, r_group] = list(map(float, [header[30:45], header[45:60]]))
                if not groups:
                    groups = list(map(float, [group_info[:15], group_info[15:30]]))
                else:
                    groups.append(float(group_info[15:30]))
            tmp = []
            if len(idx) > 1:
                max_k = idx[1] - 2
            else:
                max_k = max_idx - 1

            for k in np.arange(max_k):
                cline = f.readline()
                cline = [cline[:15], cline[15:30], cline[30:45], cline[45:60]]
                if cline[0]:
                    # means it's not actually an empty line
                    # this whole thing is very awkward. Has to be
                    # rewritten.
                    if not cline.count("\n") and not cline.count("\r\n"):
                        tmp.extend(list(map(float, cline)))
                    else:
                        if cline.count("\n"):
                            tmp.extend(list(map(float, cline[:cline.index("\n")])))
                        elif cline.count("\r\n"):
                            tmp.extend(list(map(float, cline[:cline.index("\r\n")])))

            out.append(tmp[r_len+T_len:])
            #if not len(rho) and not len(temp):
            rho = np.power(10, tmp[:r_len])
            temp = np.power(10, tmp[r_len:r_len+T_len])

            if 'dens' in self and "temp" in self:
                assert len(rho) == len(self['dens']),\
                "The rho grid is not the same for all op[prez] files"
                assert len(temp) == len(self['temp']),\
                    "The temp grid is not the same for all op[prez] files"
            else:
                self['dens'] = rho
                if tabletype == 'opz':
                    self['temp'] = temp*1e3  # opz is in keV
                else: 
                    self['temp'] = temp       # everything else is in eV

        if tabletype == 'opz':
            self["zbar"] = np.power(10, np.array(out).reshape(
                                    (len(self['temp']), len(self['dens'])))).T
        else:
            if 'groups' in self:
                assert len(groups) == len(self['groups']),\
                    "The group number is not the same for all op[prez] files"
            else:
                self['groups'] = np.array(groups)
            self[tabletype+'_mg'] = np.power(10, np.array(out).reshape(
                                            (len(idx), len(self['temp']), len(self['dens'])))).T
        f.close()

    def write(self, prefix, fmin=None, fmax=None):
        """
        Write multigroup opacities to files specified by a prefix
        """
        for key in ['eps_mg', 'temp', 'dens', 'opp_mg', 'opr_mg', 'emp_mg']:
            if key in self:
                self[key] = self[key][:]
        
        HEADER_FMT2 = "{dim[0]:{f}}{dim[1]:{f}}\n"
        extensions = {'opp_mg': 'opp','opr_mg':'opr','eps_mg': 'eps','zbar':'opz'}
        if 'eps_mg' not in self and 'emp_mg' in self:
            self['eps_mg'] = self['emp_mg']/self['opp_mg']
        if 'zbar' not in self and 'Zf_DT' in self:
            self['zbar'] = self['Zf_DT']

        for opt in filter(lambda k: k in ['opp_mg','opr_mg','eps_mg','zbar'], self):
            ctable =  self[opt]
            ext = extensions[opt]
            f  = gzip.open("{prefix}.{ext}.gz".format(prefix=prefix, ext=ext), 'w')
            if opt == 'zbar':
                HEADER_FMT1 = " {tname:14} 0.60000000E+01"
                f.write((HEADER_FMT1 + HEADER_FMT2).format(tname=self.table_name[ext],
                                                    dim=ctable.shape, f=FMT))
                X = self['dens']
                X = np.append(X, self['temp']*1e-3)
                val = ctable[:,:].T
                if fmin is not None:
                    val = np.fmax(val, fmin)
                X = np.append(X, val)
                X = np.log10(X)
                self._write_vector(f, X)
            else:
                HEADER_FMT1 = " {tname:14}{op_type:14} "
                for n in np.arange(len(self['groups']) - 1):
                    f.write((HEADER_FMT1 + HEADER_FMT2).format(tname=self.table_name[ext],
                            op_type=self._op_labels[ext]+'  ', dim = ctable.shape, f=FMT))
                    f.write("{:{f}}{:{f}}\n".format(self['groups'][n], self['groups'][n+1], f=FMT))

                    X = self['dens']
                    X = np.append(X, self['temp'])
                    val = ctable[:,:,n].T
                    val[np.isnan(val)] = (fmax is not None) and fmax or np.nanmax(val)
                    if fmin is not None:
                        val = np.fmax(val, fmin)
                    if fmax is not None:
                        val = np.fmin(val, fmax)

                    X = np.append(X, val)
                    X = np.log10(X)
                    self._write_vector(f, X)
                f.close()
            f.close()


    def write2hdf(self, filename, Znum=None, Anum=None, Xnum=None):
        """ Convert to hdf5
        Parameters
        """
        import tables
        h5filters = tables.Filters(complib='blosc', complevel=7, shuffle=True,
                least_significant_digit=6)
        f = tables.openFile(filename, 'w', filters=h5filters)

        if Znum is None:
            if "Znum" in self:
                Znum = self['Znum']
            else:
                raise ValueError('Znum Varray should be providied!')
        if type(Znum) is int:
            Znum = [Znum]
        self['Znum'] = np.array(Znum, dtype='int')
        if Anum is None:
            self['Anum'] = np.array([ptab.elements[el].mass for el in self['Znum']])
        else:
            self['Anum'] = np.array(Anum)
        self['Zsymb'] = np.array([ptab.elements[el].symbol for el in self['Znum']], dtype='|S2')
        if Xnum is None:
            if len(Znum) == 1:
                self['Xnum'] = np.array([1.0])
            else:
                raise ValueError('Xnum array should be provided')
        else:
            self['Xnum'] = np.array(Xnum)
        self['Abar'] = np.sum(self['Xnum']*self['Anum'])
        self['Zmax'] = np.sum(self['Xnum']*self['Znum'])
        self['idens'] = self['dens']*N_A/self['Abar']
        
        if "eps_mg" in self:
            self['emp_mg'] = self['opp_mg']*self['eps_mg']
        else:
            print('Warning: looks like this is LTE opacity, no eps file found!')
            print('Setting planck emissivity to be the same as plack opacity...')
            self['emp_mg'] = self['opp_mg'].copy()

        names_dict = {'idens': 'idens',
                      'temp': 'temp',
                      'dens' : 'dens',
                      'zbar': 'Zf_DT',
                      'opp_mg': 'opp_mg',
                      'opr_mg': 'opr_mg',
                      'Znum': 'Znum',
                      'Anum': 'Anum',
                      'Xnum': 'Xnum',
                      'groups': 'groups',}
        if 'eps_mg' in self:
            self['emp_mg'] = self['opp_mg']*self['eps_mg']
            names_dict['emp_mg'] = 'emp_mg'
        for prp_key, h5_key in sorted(names_dict.iteritems()):
            atom = tables.Atom.from_dtype(self[prp_key].dtype)
            ds = f.createCArray(f.root, h5_key, atom, self[prp_key].shape, filters=h5filters)
            ds[:] = self[prp_key]
        # writing attributes
        for attr in ['Abar', 'Zmax']:
            setattr(f.root._v_attrs,attr, self[attr])
        f.close()


    @staticmethod
    def _write_vector(f, X):
        remd = np.mod(len(X), 4)
        if remd:
            np.savetxt(f, X[:-remd].reshape((-1,4)), fmt=PFMT, delimiter='')
            np.savetxt(f, X[-remd:].reshape((1,-1)), fmt=PFMT, delimiter='')
        else:
            np.savetxt(f, X.reshape((-1,4)), fmt=PFMT, delimiter='')
        return f


if __name__ == '__main__':
    BASE_PATH = '/home/rth/multi90/Tables/Aluminium'
    Al = OpgMulti((BASE_PATH, 'AL_mknlte'))

    Al.write('test')
