import opacplot2 as opp
import argparse
import os.path
import textwrap

def convert_tables():
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
                        
    args = parser.parse_args()
    
    # Get the relevant paths and filenames.
    path_in = os.path.abspath(args.input_file)
    basedir, fn_in = os.path.split(path_in)
    # Split filename twice in case of MULTI files (.opr.gz, etc)
    basename = os.path.splitext(os.path.splitext(fn_in)[0])[0]
    
    # Adjust input.
    if args.outname is not None:
        fn_out = os.path.splitext(os.path.abspath(args.outname))[0]    
    else:
        fn_out = os.path.join(basedir, basename)
    
    if args.Znum is not None:
        args.Znum = [int(num) for num in args.Znum.split(',')]
    
    if args.Xfracs is not None:
        args.Xfracs = [float(num) for num in args.Xfracs.split(',')]
    
    if args.log is not None:
        args.log = [str(key) for key in args.log.split(',')]
    
    if args.input is None:
        ext_dict = {'.prp':'propaceos',
                    '.eps':'multi',
                    '.opp':'multi',
                    '.opz':'multi',
                    '.opr':'multi',
                    '.mexport':'sesame-qeos',
                    '.ses':'sesame'}
        _, ext = os.path.splitext(fn_in)
        if ext in ext_dict.keys():
            args.input = ext_dict[ext]
        else:
            raise Warning('Cannot tell filetype from extension. Please specify '
                          'input file type with --input.')
    
        
    # Reading in data and converting it to the common dictionary format.
    if args.input == 'propaceos':
        # If we are unable to find the correct library for opg_propaceos
        # we need to let the user know.
        try:
            import opacplot2.opg_propaceos
            op = opp.opg_propaceos.OpgPropaceosAscii(path_in)
            eos_dict = op.toEosDict(log=args.log)
        except ImportError:
            print('Error: You do not have the opg_propaceos script.')
            
    elif args.input == 'multi':
        op = opp.OpgMulti.open_file(basedir, basename)
        eos_dict = op.toEosDict(Znum=args.Znum, 
                                Xnum=args.Xfracs, 
                                log=args.log)
        
    elif args.input == 'sesame':
        op = opp.OpgSesame(path_in, opp.OpgSesame.SINGLE)
        eos_dict = op.toEosDict(Znum=args.Znum, 
                                Xnum=args.Xfracs,
                                log=args.log)
    elif args.input == 'sesame-qeos':
        op = opp.OpgSesame(path_in, opp.OpgSesame.SINGLE)
        eos_dict = op.toEosDict(Znum=args.Znum, Xnum=args.Xfracs, 
                                qeos=True, log=args.log)
    
    if args.output=='ionmix':
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
            if key in eos_dict.keys():
                imx_dict[imx_conv[key]] = eos_dict[key]
        
        # Set ngroups if opacity bounds are present.
        if 'opac_bounds' in imx_dict:
            imx_dict['ngroups'] = len(eos_dict['groups']) - 1                
        
        # For verbose flag.
        if args.verbose:
            verb_conv = {'zvals':'atomic numbers',
                         'fracs':'element fractions',
                         'numDens':'ion number densities',
                         'temps':'temperatures',
                         'zbar':'average ionizations',
                         'pion':'ion pressure',
                         'pele':'electron pressure',
                         'eion':'ion internal energy',
                         'eele':'electron internal energy',
                         'opac_bounds':'opacity bounds',
                         'rosseland':'Rosseland mean opacity',
                         'planck_absorb':'absorption Planck mean opacity',
                         'planck_emiss':'emission Planck mean opacity',
                         'ngroups':'number of opacity groups'}
                
            verb_str = 'Wrote the following data to IONMIX file:'
            i = 0
            for key in imx_dict.keys():
                i = i + 1
                if i == len(imx_dict.keys()):
                    verb_str = verb_str + ' {}.'.format(verb_conv[key])
                else:
                    verb_str = verb_str + ' {},'.format(verb_conv[key])
            print(textwrap.fill(verb_str))
        
        # Write the ionmix file based on what data is stored in imx_dict.
        opp.writeIonmixFile(fn_out+'.cn4', **imx_dict)
                                    
if __name__=='__main__':
    convert_tables()
    