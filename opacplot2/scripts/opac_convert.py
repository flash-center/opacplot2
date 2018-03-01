import opacplot2 as opp
import argparse
import os.path

def get_input_data():
    # Available formats.
    avail_output_formats = ['ionmix']
    avail_input_formats = ['propaceos', 'multi', 'sesame', 'sesame-qeos']

    # Creating the argument parser.
    parser = argparse.ArgumentParser(
                description="This script is used to browse various"
                            "EoS/Opacity tables formats.")

    parser.add_argument('-v', '--verbose',
                        action='store_const', const=True,
                        help='Verbosity option.')

    parser.add_argument('--Znum',
                        action='store', type=str,
                        help='Comma separated list of Z'
                             'numbers for every component.')

    parser.add_argument('--Xfracs',
                        action='store', type=str,
                        help='Comma separated list of X'
                             'fractions for every component.')

    parser.add_argument('--outname',
                        action='store', type=str,
                        help='Name for output file without extension.')

    parser.add_argument('-i', '--input',
                        action='store', type=str,
                        choices=avail_input_formats,
                        help='Input filetype.')

    parser.add_argument('-o', '--output',
                        action='store', type=str,
                        choices=avail_output_formats,
                        default='ionmix',
                        help='Output filetype. Default: IONMIX.')

    parser.add_argument('input_file',
                        action='store', type=str,
                        help='Input file.')

    parser.add_argument('--log',
                        action='store', type=str,
                        help='Logarithmic data keys.')

    parser.add_argument('--tabnum',
                        action='store', type=str,
                        help='Specify the SESAME table number.')

    args = parser.parse_args()

    # Get the relevant paths and filenames.
    path_in = os.path.abspath(args.input_file)
    basedir, fn_in = os.path.split(path_in)
    # Split filename twice in case of MULTI files (.opr.gz, etc)
    basename = os.path.splitext(os.path.splitext(fn_in)[0])[0]

    # Adjusting the input.
    # Set the base output name that specified by --outname.
    # Otherwise, set it to the same as the base input name.
    if args.outname is not None:
        args.outname = os.path.splitext(os.path.abspath(args.outname))[0]
    else:
        args.outname = os.path.join(basedir, basename)

    # Create lists out of the strings for Znum, Xfracs, and log if given.
    if args.Znum is not None:
        args.Znum = [int(num) for num in args.Znum.split(',')]
    if args.Xfracs is not None:
        args.Xfracs = [float(num) for num in args.Xfracs.split(',')]
    if args.log is not None:
        args.log = [str(key) for key in args.log.split(',')]

    # Convert tabnum into int.
    if args.tabnum is not None:
        try:
            args.tabnum = int(args.tabnum)
        except ValueError:
            raise ValueError('Please provide a valid SESAME table number.')

    input_data = {'args' : args,
                  'basename' : basename,
                  'path_in' : path_in,
                  'basedir' : basedir,
                  'fn_in' : fn_in}

    return input_data

def read_format_ext(args, fn_in):
    # Try to read from the input file extension.
    ext_dict = {'.prp':'propaceos',
                '.eps':'multi',
                '.opp':'multi',
                '.opz':'multi',
                '.opr':'multi',
                '.mexport':'sesame-qeos',
                '.ses':'sesame'}
    # If the input file is compressed, choose the next extension.
    if os.path.splitext(fn_in)[1] == '.gz':
        _, ext = os.path.splitext(os.path.splitext(fn_in)[0])
    else:
        _, ext = os.path.splitext(fn_in)

    # Choose the correct input type based on extension and set args.input
    # accordingly.
    if ext in ext_dict.keys():
        args.input = ext_dict[ext]
    else:
        raise Warning('Cannot tell filetype from extension. Please specify '
                      'input file type with --input.')

class Formats_toEosDict(object):
    """
    Contains handling functions to convert different types of tables
    into a common EoS dictionary for IONMIX.
    """
    def __init__(self, args, basedir, basename, path_in):
        # Initialize the dictionary for handling functions.
        self.set_handle_dict()

        # Set attributes.
        self.args = args
        self.basedir = basedir
        self.basename = basename
        self.path_in = path_in

        # Use handle_dict to create the eos_dict based on the input format.
        try:
            self.eos_dict = self.handle_dict[args.input]()
        except KeyError:
            raise KeyError('Must use valid format name.')

    def set_handle_dict(self):
        self.handle_dict = {'propaceos' : self.propaceos_toEosDict,
                            'multi' : self.multi_toEosDict,
                            'sesame' : self.sesame_toEosDict,
                            'sesame-qeos' : self.sesame_qeos_toEosDict}

    def propaceos_toEosDict(self):
        # If we are unable to find the correct library for opg_propaceos
        # we need to let the user know.
        try:
            import opacplot2.opg_propaceos
            op = opp.opg_propaceos.OpgPropaceosAscii(self.path_in)
            eos_dict = op.toEosDict(log=self.args.log)
            return eos_dict
        except ImportError:
            raise ImportError('You do not have the opg_propaceos script.')

    def multi_toEosDict(self):
        op = opp.OpgMulti.open_file(self.basedir, self.basename)
        eos_dict = op.toEosDict(Znum=self.args.Znum,
                                Xnum=self.args.Xfracs,
                                log=self.args.log)
        return eos_dict

    def sesame_toEosDict(self):
        try:
            op = opp.OpgSesame(self.path_in, opp.OpgSesame.SINGLE)
        except ValueError:
            op = opp.OpgSesame(self.path_in, opp.OpgSesame.DOUBLE)


        if len(op.data.keys()) > 1:
            raise Warning('More than one material ID found. '
                          'Use sesame-extract to create a file '
                          'with only one material first.')

        if self.args.tabnum is not None:
            eos_dict = op.toEosDict(Znum=self.args.Znum,
                                    Xnum=self.args.Xfracs,
                                    log=self.args.log,
                                    tabnum=self.args.tabnum)
        else:
            eos_dict = op.toEosDict(Znum=self.args.Znum,
                                    Xnum=self.args.Xfracs,
                                    log=self.args.log)
        return eos_dict

    def sesame_qeos_toEosDict(self):
        raise Warning('QEOS-SESAME is not ready yet!')
        try:
            op = opp.OpgSesame(self.path_in, opp.OpgSesame.SINGLE)
        except ValueError:
            op = opp.OpgSesame(self.path_in, opp.OpgSesame.DOUBLE)


        if len(op.data.keys()) > 1:
            raise Warning('More than one material ID found. '
                          'Use sesame-extract to create a file '
                          'with only one material first.')

        if self.args.tabnum is not None:
            eos_dict = op.toEosDict(Znum=self.args.Znum, Xnum=self.args.Xfracs,
                                    qeos=True, log=self.args.log,
                                    tabnum=self.args.tabnum)
        else:
            eos_dict = op.toEosDict(Znum=self.args.Znum, Xnum=self.args.Xfracs,
                                    qeos=True, log=self.args.log)
        return eos_dict

class EosDict_toIonmixFile(object):
    """
    Takes a common EoS dictionary and writes it to the correct output format.
    """
    def __init__(self, args, eos_dict):
        # Initialize the handling function dictionary.
        self.set_handle_dict()

        # Set attributes.
        self.args = args
        self.eos_dict = eos_dict

        # Execute the write function based on output format.
        self.handle_dict[args.output]()

    def set_handle_dict(self):
        self.handle_dict = {'ionmix' : self.eosDict_toIonmix}

    def eosDict_toIonmix(self):
        # These are the naming conventions translated to ionmix arguments.
        imx_conv = {'Znum':'zvals',
                    'Xnum':'fracs',
                    'idens':'numDens',
                    'temp':'temps',
                    'Zf_DT':'zbar',
                    'Pi_DT':'pion',
                    'Pec_DT':'pele',
                    'Ui_DT':'eion',
                    'Uec_DT':'eele',
                    'groups':'opac_bounds',
                    'opr_mg':'rosseland',
                    'opp_mg':'planck_absorb',
                    'emp_mg':'planck_emiss'}

        # Initialize ionmix argument dictionary.
        imx_dict = {}

        # Translating the keys over.
        for key in imx_conv.keys():
            if key in self.eos_dict.keys():
                imx_dict[imx_conv[key]] = self.eos_dict[key]

        # Set ngroups if opacity bounds are present.
        if 'opac_bounds' in imx_dict:
            imx_dict['ngroups'] = len(self.eos_dict['groups']) - 1

        # Check if required FLASH EoS tables are present.
        imx_req_keys = ['zbar', 'eion', 'eele', 'pion', 'pele']
        if not all(key in imx_dict for key in imx_req_keys):
            print("The required EoS data for FLASH is not present...\n"
                  "Aborting the IONMIX file creation...")
            raise Warning("Missing EoS data for IONMIX file to run in FLASH")

        # For verbose flag.
        if self.args.verbose:
            verb_conv = {'zvals':'Atomic numbers',
                         'fracs':'Element fractions',
                         'numDens':'Ion number densities',
                         'temps':'Temperatures',
                         'zbar':'Average ionizations',
                         'pion':'Ion pressure',
                         'pele':'Electron pressure',
                         'eion':'Ion internal energy',
                         'eele':'Electron internal energy',
                         'opac_bounds':'Opacity bounds',
                         'rosseland':'Rosseland mean opacity',
                         'planck_absorb':'Absorption Planck mean opacity',
                         'planck_emiss':'Emission Planck mean opacity',
                         'ngroups':'Number of opacity groups'}

            verb_str = 'Wrote the following data to IONMIX file:\n'
            i = 0
            for key in imx_dict.keys():
                i = i + 1
                if i == len(imx_dict.keys()):
                    verb_str = verb_str + '{}. {}'.format(i, verb_conv[key])
                else:
                    verb_str = verb_str + '{}. {} \n'.format(i, verb_conv[key])
            print(verb_str)

        # Write the ionmix file based on what data is stored in imx_dict.
        opp.writeIonmixFile(self.args.outname + '.cn4', **imx_dict)

def convert_tables():
    # Grab the input data.
    input_data = get_input_data()

    # Read the file extension if the user did not specify an input.
    if input_data['args'].input is None:
        read_format_ext(input_data['args'], input_data['fn_in'])



    # Reading in data and converting it to the common dictionary format.
    eos_dict = Formats_toEosDict(input_data['args'],
                                 input_data['basedir'],
                                 input_data['basename'],
                                 input_data['path_in']).eos_dict

    EosDict_toIonmixFile(input_data['args'], eos_dict)

if __name__=='__main__':
    convert_tables()
