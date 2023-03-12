from math import ceil
import numpy as np
from os.path import splitext
from io import StringIO
import periodictable

# Avogadros number
NA = 6.0221415e+23


def tops_html2text(filename):
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
    text = [' '.join(a.text.split('\xa0')).rstrip() + ' \n'
            for a in code.contents[0::2]]
    text = (' '+''.join(text).lstrip()).rstrip()

    return text


class OpgTOPS():
    def __init__(self, filename, ep_max='auto', handle_large='next_group'):
        """
        Parse TOPS Opacities (no unit conversion for this intialization)

        Parameters
        ----------
        fname : str
            Filename of the .tops or .html file
        ep_max : {'log', 'lin', 'auto'} or float, optional
            The upper bound of photon energy (keV) in the last photon group, by
            default 'auto'
            - 'log': use logarithmic extrapolation to calculate `epmax`
            - 'lin': use linear extrapolation to calculate `epmax`
            - 'auto': automatically choose from {'log', 'lin'}
        handle_large : {'no', 'lower_ceiling', 'next_group'}
            Choose how to handle group opacity entries with value 1e10, by
            default 'next_group'
            - 'no': do nothing
            - 'lower_ceiling' : use the largest opacity value throughout the
            table that is below 1e10
            - 'next_group' : use the opacity value of the next photon group at
            the same temperature-density point
        """

        self.filename = filename
        if splitext(filename)[1] == '.html':
            lines = StringIO(tops_html2text(filename)).readlines()
        else:
            with open(filename, 'r') as f:
                lines = f.readlines()

        try:
            assert isinstance(ep_max, float)
        except AssertionError:
            try:
                assert ep_max in {'log', 'lin', 'auto'}
            except AssertionError:
                raise KeyError("epmax should be in {'log', 'lin', 'auto'} or "
                               "float!")
        self.ep_max = ep_max

        try:
            assert handle_large in {'no', 'lower_ceiling', 'next_group'}
        except AssertionError:
            raise KeyError("handle_large shoule be in "
                           "{'no', 'lower_ceiling', 'next_group'}")
        self.handle_large = handle_large

        dats = [int(s) for s in lines[0].split() if s.isdigit()]
        self.NT, self.Nd, self.Nm = dats[0:3]

        # Read elements and number fractions
        for lid, line in enumerate(lines):
            line_split = line.split()
            if line_split[:10] == \
                "No. Fraction Mass Fraction  At. No.  Chem. Sym.  Mat ID.".\
                    split():
                lid_frac = lid
                lid_temp = lid_frac + 1 + self.Nm
                dataio = StringIO(''.join(lines[lid_frac+1:lid_temp]))
                dats = np.loadtxt(dataio,
                                  dtype={'names': ["Xnum", "Massfrac", "Znum",
                                                   "Zsymb",  "MatID"],
                                         'formats': ['f8', 'f8', 'i4',
                                                     'S2', 'i4']
                                         }
                                  )
                dats.reshape(self.Nm)
                # Made some tweaks to ensure Xnum and Znum are iterable
                self.Xnum = dats['Xnum'] if dats['Xnum'].size >1 else [dats['Xnum'].tolist()]
                self.Massfrac = dats['Massfrac']
                self.Znum = dats['Znum'] if dats['Znum'].size >1 else [dats['Znum'].tolist()]
                self.Zsymb = dats['Zsymb']
                self.MatID = dats['MatID']
                self.Zmax = np.average(self.Znum, weights=self.Xnum)
                self.Anum = np.array([periodictable.elements[Znum].mass
                                      for Znum in self.Znum])
                self.Abar = np.average(self.Anum, weights=self.Xnum)

        # Read temperature grid
        line_split = lines[lid_temp].split()
        assert line_split[:5] == 'Temperature grid used the following'.split()
        assert int(line_split[5]) == self.NT
        lid_dens = lid_temp + 1 + int(np.ceil(self.NT/6))
        self.temp = np.fromstring(
            ' '.join(lines[lid_temp+1:lid_dens]), sep=' ')
        self.temp.reshape(self.NT)

        # Read density grid
        line_split = lines[lid_dens].split()
        assert line_split[:5] == 'Density grid used the following'.split()
        assert int(line_split[5]) == self.Nd
        lid_grps = lid_dens + 1 + int(np.ceil(self.Nd/6))
        self.dens = np.fromstring(
            ' '.join(lines[lid_dens+1:lid_grps]), sep=' ')
        self.dens.reshape(self.Nd)
        self.nion = self.dens * NA / self.Abar

        # Read photon grid
        line_split = lines[lid_grps].split()
        if line_split[0] == 'Photon':
            self.Ng = int(line_split[5])
            self.multigroup = True
            lid_opac = lid_grps + 1 + int(np.ceil(self.Ng/6))
            self.grps = np.fromstring(' '.join(lines[lid_grps+1:lid_opac]),
                                      sep=' ')
            self.grps.reshape(self.Ng)
            self.grps = np.concatenate([self.grps, [np.inf]])
            use_fit = 'no'
            if isinstance(ep_max, float):
                self.grps[self.Ng] = ep_max
                use_fit = 'no'
            elif ep_max == 'lin':
                use_fit = 'lin'
            elif ep_max == 'log':
                use_fit = 'log'
            elif ep_max == 'auto':
                corr_log_fit = np.corrcoef(range(self.Ng-1),
                                           np.log(self.grps[:self.Ng-1]))[0, 1]
                corr_lin_fit = np.corrcoef(range(self.Ng-1),
                                           self.grps[:self.Ng-1])[0, 1]
                use_fit = ('log' if corr_log_fit > corr_lin_fit else 'lin')
            if use_fit == 'log':
                fit = np.polyfit(range(self.Ng), np.log(
                    self.grps[:self.Ng]), 1)
                self.grps[self.Ng] = np.exp(fit[1] + fit[0] * (self.Ng))
            elif use_fit == 'lin':
                fit = np.polyfit(range(self.Ng), self.grps[:self.Ng], 1)
                self.grps[self.Ng] = fit[1] + fit[0] * (self.Ng)
        elif line_split[0] == 'Rosseland':
            self.multigroup = False
            lid_opac = lid_grps

        # Read gray opacity and free electron number moments
        self.ross_int = np.zeros((self.NT, self.Nd))
        self.plnk_int = np.zeros((self.NT, self.Nd))
        self.zbar = np.zeros((self.NT, self.Nd))
        self.z2bar = np.zeros((self.NT, self.Nd))

        for t in range(self.NT):
            line_split = lines[lid_opac].split()
            assert line_split[:7] == \
                "Rosseland and Planck opacities and free electrons".split()
            line_split = lines[lid_opac+1].split()
            assert line_split[:11] == \
                "Density Ross opa Planck opa No. Free Av Sq Free T=".split()
            dataio = StringIO("".join(lines[lid_opac+2:lid_opac+2+self.Nd]))
            dats = np.loadtxt(dataio)
            dats.reshape((self.Nd, 5))
            assert (self.dens == dats[:, 0]).all()
            self.ross_int[t] = dats[:, 1]
            self.plnk_int[t] = dats[:, 2]
            self.zbar[t] = dats[:, 3]
            self.z2bar[t] = dats[:, 4]
            lid_opac += 2 + self.Nd

        self.ross_int = self.ross_int.T
        self.plnk_int = self.plnk_int.T
        self.zbar = self.zbar.T
        self.z2bar = self.z2bar.T

        # Read multigroup opacity
        if not self.multigroup:
            return
        self.ross_mg = np.zeros((self.NT, self.Nd, self.Ng))
        self.plnk_mg = np.zeros((self.NT, self.Nd, self.Ng))

        line_split = lines[lid_opac].split()
        assert line_split[:2] == 'Multigroup opacities'.split()
        lid_opac += 1

        for t in range(self.NT):
            for d in range(self.Nd):
                line_split = lines[lid_opac].split()
                assert line_split[:9] == \
                    "Energy Ross mg Planck mg for T, density =".split()
                assert self.temp[t] == float(line_split[9])
                assert self.dens[d] == float(line_split[10])
                dataio = StringIO(
                    "".join(lines[lid_opac+1:lid_opac+1+self.Ng]))
                dats = np.loadtxt(dataio)
                dats.reshape((self.Ng, 3))
                assert (self.grps[:self.Ng] == dats[:, 0]).all()
                self.ross_mg[t, d] = dats[:, 1]
                self.plnk_mg[t, d] = dats[:, 2]
                lid_opac += 1 + self.Ng

        self.ross_mg = self.ross_mg.swapaxes(0, 1)
        self.plnk_mg = self.plnk_mg.swapaxes(0, 1)

        # Handle the entries with unphsically large value 1e10
        if self.handle_large == 'no':
            return
        if self.handle_large == 'lower_ceiling':
            ceiling_value = self.ross_mg[self.ross_mg < 1e10].max()
            self.ross_mg[self.ross_mg == 1e10] = ceiling_value
            ceiling_value = self.plnk_mg[self.plnk_mg < 1e10].max()
            self.plnk_mg[self.plnk_mg == 1e10] = ceiling_value
            return
        if self.handle_large == 'next_group':
            for d in range(self.Nd):
                for t in range(self.NT):
                    for g in range(self.Ng)[::-1]:
                        if self.ross_mg[d, t, g] == 1e10:
                            self.ross_mg[d, t, g] = self.ross_mg[d, t, g+1]
                        if self.plnk_mg[d, t, g] == 1e10:
                            self.plnk_mg[d, t, g] = self.plnk_mg[d, t, g+1]

    def toEosDict(self, fill_eos=False):
        names_dict_req_tops = {
            'idens': 'nion',
            'temp': 'temp',
            'dens': 'dens',
            'Zf_DT': 'zbar',
            'z2bar': 'z2bar',
            'opr_mg': 'ross_mg',
            'opp_mg': 'plnk_mg',
            'emp_mg': 'plnk_mg',
            'opp_int': 'plnk_int',
            'opr_int': 'ross_int',
            'emp_int': 'plnk_int',
            'Xnum': 'Xnum',
            'Massfrac': 'Massfrac',
            'Znum': 'Znum',
            'Zsymb': 'Zsymb',
            'MatID': 'MatID',
            'Zmax': 'Zmax',
            'Anum': 'Anum',
            'Abar': 'Abar',
            'groups': 'grps',
            'Zsymb': 'Zsymb',
            'ElemNum': 'Nm',
        }
        names_dict_mg = {
              'opr_mg': 'ross_mg',
              'opp_mg': 'plnk_mg',
              'emp_mg': 'plnk_mg',
              'groups': 'grps',
        }
        names_list_req_eos = ['Pi_DT', 'Pec_DT', 'Ui_DT', 'Uec_DT']

        eos_dict = {}

        # Data in TOPS table
        for eos_key, tops_key in sorted(names_dict_req_tops.items()):
            eos_dict[eos_key] = getattr(self, tops_key)
        if self.multigroup:
            for eos_key, tops_key in sorted(names_dict_mg.items()):
                eos_dict[eos_key] = getattr(self, tops_key)
        eos_dict['temp'] = eos_dict['temp'] * 1e3
        eos_dict['groups'] = eos_dict['groups'] * 1e3

        # Fill zeros for EOS
        if fill_eos:
            for eos_key in names_list_req_eos:
                eos_dict[eos_key] = np.zeros_like(eos_dict['Zf_DT'])

        return eos_dict
