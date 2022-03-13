'''Convert TOPS opacity table in html format to text file'''


def tops_html2txt(filename):
    """Convert TOPS opacity table in html format to text

    Parameters
    ----------
    filename : str
        filename of html file for TOPS table

    Returns
    -------
    str
        text of TOPS table
    """

    from bs4 import BeautifulSoup

    with open(filename) as fp:
        soup = BeautifulSoup(fp, 'html.parser')

    code = soup.find_all('code')[0]
    table_text = [' '.join(a.text.split('\xa0')).rstrip() + ' \n'
                  for a in code.contents[0::2]]
    table_text = (' '+''.join(table_text).lstrip()).rstrip()

    return table_text


def main():

    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', metavar='INPUT', type=str,
                        help='input filename')
    parser.add_argument('output', metavar='OUTPUT', type=str,
                        help='input filename')
    parser.add_argument('-p', action='store_true',
                        help='print txt instead of save')
    args = parser.parse_args()

    table_text = tops_html2txt(args.input)

    if args.p:
        print(table_text)
    else:
        with open(args.output, 'w') as f:
            f.write(table_text)


if __name__ == '__main__':
    main()
