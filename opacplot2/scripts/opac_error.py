import opacplot2 as opp
import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os.path
import periodictable as ptab
import scipy as sp
plt.rcParams.update({'text.usetex': False})

def get_input_data():

    parser = argparse.ArgumentParser(
                description="This script is used to check error differences"
                            "between two files.")
    
    
    
    parser.add_argument('-v', '--verbose',
                        action='store_const', const=True,
                        help='Verbosity option.')
    
    
    parser.add_argument('input_1',
                        action='store', type=str,
                        help='Input file 1.')
    
    parser.add_argument('input_2',
                        action='store', type=str,
                        help='Input file 2.')
    
    parser.add_argument('-f', '--filetypes',
                        action='store', type=str,
                        help='Input filetypes.')
    
    parser.add_argument('--mpi_1',
                        action='store', type=str,
                        help='Mass per ion for file 1.')
    
    parser.add_argument('--mpi_2',
                        action='store', type=str,
                        help='Mass per ion for file 2.')
    
    parser.add_argument('--Znum_1',
                        action='store', type=str,
                        help='Atomic numbers for file 1.')
    
    parser.add_argument('--Znum_2',
                        action='store', type=str,
                        help='Atomic numbers for file 2.')
    
    parser.add_argument('--Xfracs_1',
                        action='store', type=str,
                        help='Number fractions for file 1.')
    
    parser.add_argument('--Xfracs_2',
                        action='store', type=str,
                        help='Number fractions for file 2.')
    
    parser.add_argument('--filters_1',
                        action='store', type=str,
                        help='dens, temp filter list '
                             'for SESAME for file 1 (g/cm^3, eV).')
    
    parser.add_argument('--filters_2',
                        action='store', type=str,
                        help='dens, temp filter list '
                             'for SESAME for file 2 (g/cm^3, eV).')
    
    parser.add_argument('-p','--plot',
                        action='store_const',
                        const=True, default=False)
    
    parser.add_argument('--writelog',
                        action='store_const',
                        const=True, default=False,
                        help='Write error values to file.')
    
    parser.add_argument('--lin_grid',
                        action='store_const',
                        const=True, default=False,
                        help='Linear values for interpolated grid.')
    
    parser.add_argument('--tabnum_1',
                        action='store', type=str,
                        help='Specify the SESAME table number for file 1.')
    
    parser.add_argument('--tabnum_2',
                        action='store', type=str,
                        help='Specify the SESAME table number for file 2.')
    
    args = parser.parse_args()
    
    # Get the relevant paths and filenames.
    path_in_1 = os.path.abspath(args.input_1)
    path_in_2 = os.path.abspath(args.input_2)
    
    basedir_1, fn_1 = os.path.split(path_in_1)
    basedir_2, fn_2 = os.path.split(path_in_2)
    
    # Split filename twice in case of MULTI files (.opr.gz, etc)
    basename_1 = os.path.splitext(os.path.splitext(fn_1)[0])[0]
    basename_2 = os.path.splitext(os.path.splitext(fn_2)[0])[0]
    
    # Create lists for filetypes.
    if args.filetypes is not None:
        args.filetypes = [typ for typ in args.filetypes.split(',')]
    
    if args.Znum_1 is not None:
        args.Znum_1 = [num for num in args.Znum_1.split(',')]
    if args.Znum_2 is not None:
        args.Znum_2 = [num for num in args.Znum_2.split(',')]
    
    # Convert mpis to float.
    if args.mpi_1 is not None:
        args.mpi_1 = float(args.mpi_1)
    if args.mpi_2 is not None:
        args.mpi_2 = float(args.mpi_2)
    
    # Convert xfracs to float list.
    if args.Xfracs_1 is not None:
        args.Xfracs_1 = [float(x) for x in args.Xfracs_1.split(',')]
    if args.Xfracs_2 is not None:
        args.Xfracs_2 = [float(x) for x in args.Xfracs_2.split(',')]
    
    # Set defaults for SESAME filters.
    if args.filters_1 is not None:
        args.filters_1 = [float(num) for num in args.filters_1.split(',')]
    else:
        args.filters_1 = [0., 0.,]
    if args.filters_2 is not None:
        args.filters_2 = [float(num) for num in args.filters_2.split(',')]
    else:
        args.filters_2 = [0., 0.,]
    
    # Convert tabnum into int.
    if args.tabnum_1 is not None:
        try:
            args.tabnum_1 = int(args.tabnum_1)
        except ValueError:
            raise ValueError('Please provide a valid '
                             'SESAME table number for file 1.')
    if args.tabnum_2 is not None:
        try:
            args.tabnum_2 = int(args.tabnum_2)
        except ValueError:
            raise ValueError('Please provide a valid '
                             'SESAME table number for file 2.')
    
    
    input_data = {'args':args,
                  'basename_1':basename_1,
                  'basename_2':basename_2,
                  'path_in_1':path_in_1,
                  'path_in_2':path_in_2,
                  'basedir_1':basedir_1,
                  'basedir_2':basedir_2,
                  'fn_1':fn_1,
                  'fn_2':fn_2}
    
    
    return input_data

def read_format_ext(args, f_1, f_2):
    # Try to read from the input file extension.
    ext_dict = {'.prp':'propaceos',
                '.eps':'multi',
                '.opp':'multi',
                '.opz':'multi',
                '.opr':'multi',
                '.mexport':'sesame-qeos',
                '.ses':'sesame',
                '.cn4':'ionmix'}
    
    # If the input file is compressed, choose the next extension.
    if os.path.splitext(f_1)[1] == '.gz':
        _, ext_1 = os.path.splitext(os.path.splitext(f_1)[0])
    else:
        _, ext_1 = os.path.splitext(f_1)
    
    if os.path.splitext(f_2)[1] == '.gz':
        _, ext_2 = os.path.splitext(os.path.splitext(f_2)[0])
    else:
        _, ext_2 = os.path.splitext(f_2)
    
    # Choose the correct input type based on extension and set args.input
    # accordingly.
    args.filetypes = []
    if ext_1 in ext_dict.keys():
        args.filetypes = args.filetypes + [ext_dict[ext_1]]
    else:
        raise Warning('Cannot tell filetype from extension {}. Please specify '
                      'input file type with --input.'.format(ext_1))
    if ext_2 in ext_dict.keys():
        args.filetypes = args.filetypes + [ext_dict[ext_2]]
    else:
        raise Warning('Cannot tell filetype from extension {}. Please specify '
                      'input file type with --input.'.format(ext_2))

class Formats_Read(object):
    """
    Reads in a file and returns an object with useful attributes based on the
    corresponding opacplot2 object. This class also includes the naming
    conventions for each format.
    
    
    This procedure is preferred (although it is rather redundant with the
    rest of opacplot2) since it preserves the original structure of the
    opacplot2 object & all calculations & processing is done in this class
    alone. This way, the mechanisms are more transparent for error checking.
    
    
    Returns
    -------
    Formats_Read
        Formats_Read().data is the opacplot2 object corresponding to the input
        file.
        
        Formats_Read().ft is the filetype of the input file.
        
        Formats_Read().common_keys are the "common dictionary style" keys for
        opacplot2.
    """
    
    # Names dictionaries.
    # Convert "common dictionary style" keys -> Format keys
    propaceos_names_dict = {'idens' : 'nion',
                            'temp' : 'temp',
                            'dens' : 'dens',
                            'Zf_DT' : 'zbar',
                            'Ut_DT' : 'eint',
                            'Uec_DT' : 'eele',
                            'Ui_DT' : 'eion',
                            'Pi_DT' : 'pion',
                            'Pec_DT' : 'pele',
                            'opp_mg' : 'opp_mg',
                            'opr_mg' : 'opr_mg',
                            'emp_mg' : 'emp_mg',
                            'opp_int' : 'opp_int',
                            'opr_int' : 'opr_int',
                            'emp_int' : 'emp_int',
                            'Znum' : 'Znum',
                            'Anum' : 'Anum',
                            'Xnum' : 'Xnum',
                            'Anum_prp':  'Anum_prp',
                            'groups' : 'groups',
                            'Zsymb' : 'Zsymb',
                            'BulkMod' : 'BulkMod',
                            'ElemNum' : 'ElemNum',
                            'Abar' : 'Abar',
                            'Zmax' : 'Zmax'}
    
    multi_names_dict = {'idens':'idens',
                        'temp':'temp',
                        'dens':'dens',
                        'Zf_DT':'zbar',
                        'opp_mg':'opp_mg',
                        'opr_mg':'opr_mg',
                        'Znum':'Znum',
                        'Anum':'Anum',
                        'Xnum':'Xnum',
                        'groups':'groups',
                        'Abar':'Abar',
                        'Zmax':'Zmax',
                        'emp_mg':'emp_mg'}
    
    sesame_names_dict = {'idens':'idens',
                         'temp':'ele_temps',
                         'dens':'ele_dens',
                         'Zf_DT':'zbar',
                         'Ut_DT':'total_eint',
                         'Uec_DT':'ele_eint',
                         'Ui_DT':'ioncc_eint',
                         'Pi_DT':'ioncc_pres',
                         'Pec_DT':'ele_pres',
                         'Znum':'Znum',
                         'Xnum':'Xnum',
                         'BulkMod':'bulkmod',
                         'Abar':'abar',
                         'Zmax':'zmax'}

    
    sesame_qeos_names_dict = {'idens' : 'idens',
                              'temp' : 'ele_temps',
                              'dens' : 'ele_dens',
                              'Zf_DT' : 'zbar',
                              'Ut_DT' : 'total_eint',
                              'Uec_DT' : 'ele_eint',
                              'Ui_DT' : 'ion_eint',
                              'Pi_DT' : 'ion_pres',
                              'Pec_DT' : 'ele_pres',
                              'Znum' : 'Znum',
                              'Xnum' : 'Xnum',
                              'BulkMod' : 'bulkmod',
                              'Abar' : 'abar',
                              'Zmax' : 'zmax'}
    
    
    ionmix_names_dict = {'Znum' : 'zvals',
                         'Xnum' : 'fracs',
                         'idens' : 'numDens',
                         'temp' : 'temps',
                         'Zf_DT' : 'zbar',
                         'Pi_DT' : 'pion',
                         'Pec_DT' : 'pele',
                         'Ui_DT' : 'eion',
                         'Uec_DT' : 'eele',
                         'groups' : 'opac_bounds',
                         'opr_mg' : 'rosseland',
                         'opp_mg' : 'planck_absorb',
                         'emp_mg' : 'planck_emiss'}
    
    # Inverted names: Convert format keys -> "common dictionary style" keys
    propaceos_names_dict_inv = {v:k for k, v
                                in propaceos_names_dict.items()}
    multi_names_dict_inv = {v:k for k, v
                            in multi_names_dict.items()}
    sesame_names_dict_inv = {v:k for k, v
                             in sesame_names_dict.items()}
    sesame_qeos_names_dict_inv = {v:k for k, v
                                  in sesame_qeos_names_dict.items()}
    ionmix_names_dict_inv = {v:k for k, v
                             in ionmix_names_dict.items()}
    
    def __init__(self, form, basedir, basename, path_in,
                 mpi=None, znum=None, xnum=None,
                 filters=[0.,0.], verbose=False, tabnum=None):
        # Initialize the dictionary for handling functions.
        self.set_handle_dict()
        
        # Set attributes.
        self.form = form
        self.basedir = basedir
        self.basename = basename
        self.path_in = path_in
        self.mpi = mpi
        self.znum = znum
        self.filters = filters
        self.verbose = verbose
        self.xnum = xnum
        self.tabnum = tabnum
        
        # For SESAME, we need the hedp package to calculate zbar.
        need_hedp_list = ['sesame', 'sesame-qeos']
        if self.form in need_hedp_list:
            try:
                global hedp
                import hedp.eos
            except ImportError:
                raise ImportError('You need the hedp module. You can get it here: '
                                  'https://github.com/luli/hedp.')
        
        # Use handle_dict to create the eos_dict based on the input format.
        try:
            self.data = self.handle_dict[self.form]()
            self.ft = self.form
        except KeyError:
            raise KeyError('{} is not a valid format name!'.format(self.form))
    
    def set_handle_dict(self):
        self.handle_dict = {'propaceos' : self.propaceos_read,
                            'multi' : self.multi_read,
                            'sesame' : self.sesame_read,
                            'sesame-qeos' : self.sesame_qeos_read,
                            'ionmix' : self.ionmix_read}
    
    def propaceos_read(self):
        # If we are unable to find the correct script for opg_propaceos
        # we need to let the user know.
        try:
            import opacplot2.opg_propaceos
            op = opp.opg_propaceos.OpgPropaceosAscii(self.path_in)
            self.common_keys = [self.propaceos_names_dict_inv[key]
                                for key in op.keys()
                                if key in self.propaceos_names_dict_inv.keys()]
            return op
        except ImportError:
            raise ImportError('You do not have the opg_propaceos script.')
    
    def multi_read(self):
        op = opp.OpgMulti.open_file(self.basedir, self.basename)
        
        # Decide if we need Znum, and Xnum, calculate Anum if it is not given.
        if self.znum is None:
            if 'Znum' in op:
                self.znum = op['Znum']
            else:
                raise ValueError('Znum Varray should be provided!')
        if type(self.znum) is int:
            self.znum = [self.znum]
        op['Znum'] = np.array(self.znum, dtype='int')
        op['Anum'] = np.array([ptab.elements[el].mass for el in op['Znum']])
        if self.xnum is None:
            if len(self.znum) == 1:
                op['Xnum'] = np.array([1.0])
            else:
                raise ValueError('Xnum array should be provided')
        else:
            op['Xnum'] = np.array(Xnum)
        
        # Setting more attributes.
        op['Abar'] = np.sum(op['Xnum']*op['Anum'])
        op['Zmax'] = np.sum(op['Xnum']*op['Znum'])
        op['idens'] = op['dens']*opp.NA/op['Abar']
        self.common_keys = [self.multi_names_dict_inv[key]
                            for key in op.keys()
                            if key in self.multi_names_dict_inv.keys()]
        return op
    
    def sesame_read(self):
        # TODO Add options for single vs double
        if self.verbose:
            print('Opening up QEOS SESAME file {}...'.format(self.path_in))
        # Try SINGLE precision and then DOUBLE if that doesn't work.
        try:
            op = opp.OpgSesame(self.path_in, opp.OpgSesame.SINGLE)
        except ValueError:
            op = opp.OpgSesame(self.path_in, opp.OpgSesame.DOUBLE)
        if self.tabnum is not None:
            table_key = self.tabnum
        else:
            if self.verbose:
                print('Selecting the last table available...')
            # Select the last table (newest) table available.
            table_key = sorted(op.data.keys())[-1]
        
        if self.verbose:
            print('Setting the atomic numbers...')
        # Sesame needs Znum.
        if self.znum is None:
            if 'Znum' in op.data[table_key].keys():
                self.znum = op.data[table_key]['Znum']
            else:
                raise ValueError('Znum Varray should be provided!')
        
        op.data[table_key]['Znum'] = np.array(self.znum, dtype='int')
        
        if self.verbose:
            print('Merging the Ion and '
                  'Electron temperature and density grids...')
        
        # We must merge ion_ and ele_ grids for qeos-sesame data.
        # Then we can calculate zbar using hedp module.
        op.data[table_key] = opp.utils.EosMergeGrids(
                            op.data[table_key], intersect=['ele', 'ioncc'],
                            filter_dens=lambda x: (x>self.filters[0]),
                            filter_temps=lambda x: (x>self.filters[1]),
                            qeos=False)
        
        if self.verbose:
            print('Calculating average ionization...')
        dens_arr, temp_arr = np.meshgrid(op.data[table_key]['ele_dens'],
                                         op.data[table_key]['ele_temps'])
        
        zbar = hedp.eos.thomas_fermi_ionization(
                                    dens_arr, temp_arr,
                                    op.data[table_key]['Znum'],
                                    op.data[table_key]['abar']).T
        
        op.data[table_key]['zbar'] = zbar
        
        if self.verbose:
            print('Calculating number densities...')
        # Add in number density key.
        op.data[table_key]['idens'] = ((op.data[table_key]['ele_dens']
                                             * opp.NA)
                                             / op.data[table_key]['abar'])
        
        # Create a list of the "common dictionary format" keys.
        self.common_keys = [self.sesame_names_dict_inv[key]
                            for key in op.data[table_key].keys()
                            if key in self.sesame_names_dict_inv.keys()]
        
        return op
    
    def sesame_qeos_read(self):
        if self.verbose:
            print('Opening up QEOS SESAME file {}...'.format(self.path_in))
        # Try SINGLE precision and then DOUBLE if that doesn't work.
        try:
            op = opp.OpgSesame(self.path_in, opp.OpgSesame.SINGLE)
        except ValueError:
            op = opp.OpgSesame(self.path_in, opp.OpgSesame.DOUBLE)
            
        if self.tabnum is not None:
            table_key = self.tabnum
        else:
            if self.verbose:
                print('Selecting the last table available...')
            # Select the last table (newest) table available.
            table_key = sorted(op.data.keys())[-1]
        
        # Sesame needs Znum.
        if self.znum is None:
            if 'Znum' in op.data[table_key].keys():
                self.znum = op.data[table_key]['Znum']
            else:
                raise ValueError('Znum Varray should be provided!')
        
        op.data[table_key]['Znum'] = np.array(self.znum, dtype='int')
        
        if self.verbose:
            print('Merging the Ion and '
                  'Electron temperature and density grids...')
        # We must merge ion_ and ele_ grids for qeos-sesame data.
        # Then we can calculate zbar using hedp module.
        op.data[table_key] = opp.utils.EosMergeGrids(
                            op.data[table_key], intersect=['ele', 'ion'],
                            filter_dens=lambda x: (x>self.filters[0]),
                            filter_temps=lambda x: (x>self.filters[1]),
                            qeos=True)
        
        if self.verbose:
            print('Calculating average ionization...')
        dens_arr, temp_arr = np.meshgrid(op.data[table_key]['ele_dens'],
                                         op.data[table_key]['ele_temps'])
        
        zbar = hedp.eos.thomas_fermi_ionization(
                                    dens_arr, temp_arr,
                                    op.data[table_key]['Znum'],
                                    op.data[table_key]['abar']).T
        
        op.data[table_key]['zbar'] = zbar
        
        if self.verbose:
            print('Calculating number densities...')
        # Add in number density key.
        op.data[table_key]['idens'] = ((op.data[table_key]['ele_dens']
                                             * opp.NA)
                                             / op.data[table_key]['abar'])
        
        # Create a list of the "common dictionary format" keys.
        self.common_keys = [self.sesame_qeos_names_dict_inv[key]
                            for key in op.data[table_key].keys()
                            if key in self.sesame_qeos_names_dict_inv.keys()]
        
        return op
    
    def ionmix_read(self):
        if self.verbose:
            print('Opening up IONMIX file {}...'.format(self.path_in))
        if self.mpi is None:
            raise Warning('Need mpi for ionmix!')
        else:
            # TODO Add options for man and twot
            op = opp.OpacIonmix(self.path_in, self.mpi, man=True, twot=True)
            self.common_keys = [self.ionmix_names_dict_inv[attr]
                                for attr in dir(op)
                                if attr in self.ionmix_names_dict_inv.keys()]
            return op

class get_eos_array(object):
    """
    Gets an EoS array based on the input opacplot2 object.
    """
    def __init__(self, eos, arr):
        
        # Initialize the dictionary for handling functions.
        self.set_handle_dict()
        
        # Use handle_dict to create the eos_dict based on the input format.
        try:
            self.arr = self.handle_dict[eos.ft](eos, arr)
        # TODO fix this, error handling is not useful since key erros come
        # from everywhere in this class.
        except KeyError:
            raise KeyError('{} is not a valid format name!'.format(eos.ft))
    
    def set_handle_dict(self):
        self.handle_dict = {'propaceos' : self.propaceos,
                            'multi' : self.multi,
                            'sesame' : self.sesame,
                            'sesame-qeos' : self.sesame_qeos,
                            'ionmix' : self.ionmix}
    
    def propaceos(self, eos, arr):
        return eos.data[Formats_Read.propaceos_names_dict[arr]]
    
    def multi(self, eos, arr):
        return eos.data[Formats_Read.multi_names_dict[arr]]
    
    def sesame(self, eos, arr):
        if eos.tabnum is None:
            # Select the last table (newest) table available.
            table_key = sorted(eos.data.data.keys())[-1]
        else:
            table_key = eos.tabnum
        data_dict = eos.data.data[table_key]
        return data_dict[Formats_Read.sesame_names_dict[arr]]
    
    def sesame_qeos(self, eos, arr):
        if eos.tabnum is None:
            # Select the last table (newest) table available.
            table_key = sorted(eos.data.data.keys())[-1]
        else:
            table_key = eos.tabnum
        data_dict = eos.data.data[table_key]
        return data_dict[Formats_Read.sesame_qeos_names_dict[arr]]
    
    def ionmix(self, eos, arr):
        return getattr(eos.data, Formats_Read.ionmix_names_dict[arr])

def compare_eos(eos_1, eos_2, verbose=False,
                plot=False,
                write_log_file=False,
                lin_grid=False):
    
    logfile_name = 'eos_errors.txt'
        
    # Union of all "common dictionary format" keys to do a full error report.
    # Not including 'idens', 'dens', 'temp', 'groups', 'opp_mg', 'opp_int',
    # 'opr_int',  'emp_mg', 'opr_mg', 'emp_int', 'Abar','Zmax', 'Ut_DT', 
    # 'Znum', 'BulkMod', 'Xnum', 'Anum', 'Zsymb',  'ElemNum', 'Anum_prp'.
    # (aka no opacity data currently).
    # TODO add opacity comparison capabilities.
    keys = ['Pec_DT', 'Zf_DT', 'Pi_DT', 'Uec_DT', 'Ui_DT']
    
    shared_keys = [key for key in keys
                   if key in eos_1.common_keys
                   and key in eos_2.common_keys]
    
    if verbose:
        err_report_str = 'Performing error report on:\n'
        for i in range(len(shared_keys)):
            err_report_str += '{}. {}\n'.format((i+1), shared_keys[i])
        print(err_report_str)
    
    # Perform error report using number densities.
    error_report = []
    
    # Freak out if there is no number density.
    if 'idens' not in eos_1.common_keys or 'idens' not in eos_2.common_keys:
        raise Warning('No number density data!')
    
    # Get the temperature and density arrays.
    dens_1 = get_eos_array(eos_1, 'idens').arr
    temp_1 = get_eos_array(eos_1, 'temp').arr
    dens_2 = get_eos_array(eos_2, 'idens').arr
    temp_2 = get_eos_array(eos_2, 'temp').arr
        
    # These will be used for the interpolator function `griddata`.
    d_interp_1, t_interp_1 = np.meshgrid(dens_1, temp_1)
    d_interp_2, t_interp_2 = np.meshgrid(dens_2, temp_2)
    
    # Creating a new grid to interpolate onto.
    d = opp.utils.intersect_1D_sorted_arr(dens_1, dens_2)
    t = opp.utils.intersect_1D_sorted_arr(temp_1, temp_2)
    D_new, T_new = np.meshgrid(d,t)
    
    # These will be used for the interpolator function `griddata`.
    d_interp_1, t_interp_1 = np.meshgrid(dens_1, temp_1)
    d_interp_2, t_interp_2 = np.meshgrid(dens_2, temp_2)
    
    if (d is None) or (t is None):
        raise Warning('Density and temperature arrays must have some overlap!')
    if verbose:
        print('Density range: {:.5E} to {:.5E} #/cm^3.'.format(d[0], d[-1]))
        print('Temperature range: {:.5E} to {:.5E} eV.'.format(t[0], t[-1]))
        print('Generating error report...')
    
    fn_1 = os.path.split(eos_1.path_in)[1]
    fn_2 = os.path.split(eos_2.path_in)[1]
    
    if write_log_file:
        # Append heading for our current grid.
        with open(logfile_name, 'a') as f:
            f.write('Files: {}, {}\n'.format(fn_1, fn_2))
            f.write('[Array, RMS Error, Absolute Error]\n')
    
    # Do analysis on each of the shared keys.
    for key in shared_keys:
        # Get the data.
        data_1 = get_eos_array(eos_1, key).arr
        data_2 = get_eos_array(eos_2, key).arr

        # Use interpolation to account for mismatched grid sizes.
        # `rescale=True` to account for the orders of magnitude difference
        # in the dens/temp grids.
        # `scipy.interpolate.interp2d` was not giving accurate interpolation.
        # I believe this is due to the orders of magnitude difference also,
        # which `griddata` can easily fix. - JT
        # Additionally, `griddata` is much faster than using an interpolator
        # function to fill an empty grid.
        interp_data_1 = sp.interpolate.griddata(
                                (d_interp_1.flatten(), t_interp_1.flatten()),
                                data_1.T.flatten(),
                                (D_new.flatten(), T_new.flatten()), 
                                rescale=True, 
                                method='linear')
        interp_data_2 = sp.interpolate.griddata(
                                (d_interp_2.flatten(), t_interp_2.flatten()),
                                data_2.T.flatten(),(D_new.flatten(), T_new.flatten()),
                                rescale=True,
                                method='linear')
        
        interp_data_1 = interp_data_1.reshape(D_new.shape[0], D_new.shape[1])
        interp_data_2 = interp_data_2.reshape(D_new.shape[0], D_new.shape[1])
        interp_data_1 = interp_data_1.T
        interp_data_2 = interp_data_2.T
                   
        err_1_sqr = np.square((interp_data_1 - interp_data_2)/interp_data_1)
        err_2_sqr = np.square((interp_data_1 - interp_data_2)/interp_data_2)
        
        err_1_rms = np.sqrt(err_1_sqr.mean())
        err_2_rms = np.sqrt(err_2_sqr.mean())
        err_rms = max(err_1_rms, err_2_rms)
        
        err_1_abs = np.sqrt(np.max(err_1_sqr))
        err_2_abs = np.sqrt(np.max(err_2_sqr))
        err_abs = max(err_1_abs, err_2_abs)
        
        fmt='%.0f %%'       
        if plot:
            titles = {'Zf_DT':'Average Ionization',
                      'Pec_DT':'Electron Pressure',
                      'Pi_DT':'Ion Pressure',
                      'Uec_DT':'Electron Energy',
                      'Ui_DT':'Ion Energy'}
            
            fig, axarr = plt.subplots(1,3)
            x, y = np.meshgrid(d, t)
            res_levels = {0:.1, 1:.01, 2:.0001 }
            fig.set_size_inches(21, 6)
                
            for i in range(3):    
                levels = np.linspace(0, res_levels[i], 256)
                cs = axarr[i].contourf(x, y, np.sqrt(err_1_sqr).T, 
                                       levels, extend='max')
                cb = plt.colorbar(cs, ax=axarr[i])
                cb.formatter = matplotlib.ticker.FuncFormatter(lambda x,p: '{:.1e}%'.format(float(x)*100))
                cb.update_ticks()
                if i==2:
                    cb.set_label('% Error')
                if not lin_grid:
                    axarr[i].loglog()
                axarr[i].set_xlim((d[0], d[-1]))
                axarr[i].set_ylim((t[0], t[-1]))
                axarr[i].set_xlabel('rho [#/cm^(3)]')
                axarr[i].set_ylabel('T [eV]')                
            
            fig.tight_layout()
            fig.suptitle('{} % Error for {} vs. {}'.format(titles[key], fn_1, fn_2))
            fig.subplots_adjust(top=0.85)
            fig.savefig('{}.png'.format(key+'_err'))
        
        if write_log_file:
            with open(logfile_name, 'a') as f:
                f.write('{}, {}, {}\n'.format(key, err_rms, err_abs))
        
        print('Error statistics for {}:'.format(key))
        print('RMS % Error: {:.5e}.'.format(err_rms*100))
        print('Max % Absolute Error: {:.5e}.'.format(err_abs*100))

def check_error():
    input_data = get_input_data()
    
    if input_data['args'].filetypes is None:
        read_format_ext(input_data['args'],
                        input_data['fn_1'],
                        input_data['fn_2'])
        
    eos_1 = Formats_Read(input_data['args'].filetypes[0],
                         input_data['basedir_1'],
                         input_data['basename_1'],
                         input_data['path_in_1'],
                         mpi=input_data['args'].mpi_1,
                         znum=input_data['args'].Znum_1,
                         xnum=input_data['args'].Xfracs_1,
                         filters=input_data['args'].filters_1,
                         verbose=input_data['args'].verbose,
                         tabnum=input_data['args'].tabnum_1)
    
    eos_2 = Formats_Read(input_data['args'].filetypes[1],
                         input_data['basedir_2'],
                         input_data['basename_2'],
                         input_data['path_in_2'],
                         mpi=input_data['args'].mpi_2,
                         znum=input_data['args'].Znum_2,
                         xnum=input_data['args'].Xfracs_2,
                         filters=input_data['args'].filters_2,
                         verbose=input_data['args'].verbose,
                         tabnum=input_data['args'].tabnum_2)
    
    compare_eos(eos_1, eos_2, verbose=input_data['args'].verbose,
                plot=input_data['args'].plot,
                write_log_file=input_data['args'].writelog,
                lin_grid=input_data['args'].lin_grid)

if __name__=='__main__':
    check_error()
