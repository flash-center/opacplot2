'''Convert TOPS opacity table in html format to text file'''

from opacplot2 import tops_html2text
from os.path import splitext

def tops_html2txt():

    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', metavar='INPUT', type=str,
                        help='input filename')
    parser.add_argument('--output', '-o', metavar='OUTPUT', type=str,
                        help='output filename')
    parser.add_argument('-p', action='store_true',
                        help='print txt instead of save')
    args = parser.parse_args()

    text = tops_html2text(args.input)

    if args.p:
        print(text)
    else:
        if args.output is None:
            fn_out = splitext(args.input)[1] + '.tops'
        else:
            fn_out = args.output
        with open(fn_out, 'w') as f:
            f.write(text)


if __name__ == '__main__':
    tops_html2txt()
