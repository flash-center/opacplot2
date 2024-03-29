from __future__ import print_function

import re
import random

import math
import numpy as np

from .constants import BC_BOUND, BC_EXTRAP_ZERO
from .constants import INTERP_FUNC, INTERP_DFDD, INTERP_DFDT

import os.path
import opacplot2

import scipy as sp
import scipy.interpolate
import scipy.misc

import copy

# from yt import physical_constants as p
# from yt import units as u
# emis_const = float((60.0 * p.stefan_boltzmann_constant /
#               p.pi**4 * (u.eV/p.kb)**4 * (u.cm**2/u.g)).in_cgs())
emis_const = 633391171028.5317

try:
    import numba
    vectorize = numba.vectorize('float64(float64)', nopython=True)
except:
    try:
        print('numba.vectorize not loaded, will use numpy.vectorize!')
        def vectorize(f): return np.vectorize(f, otypes=[float])
    except:
        print('numpy.vectorize not loaded, will use `list(map(f, x))`!')
        def vectorize(f): return \
            (lambda x: list(map(f, x)) if isinstance(x, list) else f(x))




@vectorize
def planck_int_621(x):
    if x <= 0.0:
        return 0.0
    x2 = x * x
    x3 = x2 * x
    if x < 2.6220995131254904878:
        p = 43867.0/107290978560589824000.0
        p = -3617.0/202741834014720000.0+x2*p
        p = 1.0/1270312243200.0+x2*p
        p = -691.0/19615115520000.0+x2*p
        p = 1.0/622702080.0+x2*p
        p = -1.0/13305600.0+x2*p
        p = 1.0/272160.0+x2*p
        p = -1.0/5040.0+x2*p
        p = 1.0/60.0+x2*p
        p = -1.0/8.0+x*p
        p = 1.0/3.0+x*p
        p = x3*p
        return p
    expmx = math.exp(-x)
    p = 6.4939394022668291491 - expmx*(6.0+6.0*x+3.0*x2+x3)
    expmnx = expmx*expmx
    p = p - expmnx*(0.375+0.75*x+0.75*x2+0.5*x3)
    expmnx = expmnx*expmx
    p = p - expmnx*((2.0/27.0)+(2.0/9.0)*x+(1.0/3.0)*x2+(1.0/3.0)*x3)
    expmnx = expmnx*expmx
    p = p - expmnx*(0.0234375+0.09375*x+0.1875*x2+0.25*x3)
    expmnx = expmnx*expmx
    p = p - expmnx*(0.0096+0.048*x+0.12*x2+0.2*x3)
    expmnx = expmnx*expmx
    p = p - expmnx*((1.0/216.0)+(1.0/36.0)*x+(1.0/12.0)*x2+(1.0/6.0)*x3)
    return p


@vectorize
def planck_int_521(x):
    if x <= 0.0:
        return 0.0
    x2 = x * x
    x3 = x2 * x
    if x < 2.8319431973583671745:
        p = 43867.0/107290978560589824000.0
        p = -3617.0/202741834014720000.0+x2*p
        p = 1.0/1270312243200.0+x2*p
        p = -691.0/19615115520000.0+x2*p
        p = 1.0/622702080.0+x2*p
        p = -1.0/13305600.0+x2*p
        p = 1.0/272160.0+x2*p
        p = -1.0/5040.0+x2*p
        p = 1.0/60.0+x2*p
        p = -1.0/8.0+x*p
        p = 1.0/3.0+x*p
        p = x3*p
        return p
    expmx = math.exp(-x)
    p = 6.4939394022668291491 - expmx*(6.0+6.0*x+3.0*x2+x3)
    expmnx = expmx*expmx
    p = p - expmnx*(0.375+0.75*x+0.75*x2+0.5*x3)
    expmnx = expmnx*expmx
    p = p - expmnx*((2.0/27.0)+(2.0/9.0)*x+(1.0/3.0)*x2+(1.0/3.0)*x3)
    expmnx = expmnx*expmx
    p = p - expmnx*(0.0234375+0.09375*x+0.1875*x2+0.25*x3)
    expmnx = expmnx*expmx
    p = p - expmnx*(0.0096+0.048*x+0.12*x2+0.2*x3)
    return p


@vectorize
def planck_int_421(x):
    if x <= 0.0:
        return 0.0
    x2 = x * x
    x3 = x2 * x
    if x < 3.0912856367280557733:
        p = 43867.0/107290978560589824000.0
        p = -3617.0/202741834014720000.0+x2*p
        p = 1.0/1270312243200.0+x2*p
        p = -691.0/19615115520000.0+x2*p
        p = 1.0/622702080.0+x2*p
        p = -1.0/13305600.0+x2*p
        p = 1.0/272160.0+x2*p
        p = -1.0/5040.0+x2*p
        p = 1.0/60.0+x2*p
        p = -1.0/8.0+x*p
        p = 1.0/3.0+x*p
        p = x3*p
        return p
    expmx = math.exp(-x)
    p = 6.4939394022668291491 - expmx*(6.0+6.0*x+3.0*x2+x3)
    expmnx = expmx*expmx
    p = p - expmnx*(0.375+0.75*x+0.75*x2+0.5*x3)
    expmnx = expmnx*expmx
    p = p - expmnx*((2.0/27.0)+(2.0/9.0)*x+(1.0/3.0)*x2+(1.0/3.0)*x3)
    expmnx = expmnx*expmx
    p = p - expmnx*(0.0234375+0.09375*x+0.1875*x2+0.25*x3)
    return p


@vectorize
def planck_int_321(x):
    if x <= 0.0:
        return 0.0
    x2 = x * x
    x3 = x2 * x
    if x < 3.4233592106975994941:
        p = 43867.0/107290978560589824000.0
        p = -3617.0/202741834014720000.0+x2*p
        p = 1.0/1270312243200.0+x2*p
        p = -691.0/19615115520000.0+x2*p
        p = 1.0/622702080.0+x2*p
        p = -1.0/13305600.0+x2*p
        p = 1.0/272160.0+x2*p
        p = -1.0/5040.0+x2*p
        p = 1.0/60.0+x2*p
        p = -1.0/8.0+x*p
        p = 1.0/3.0+x*p
        p = x3*p
        return p
    expmx = math.exp(-x)
    p = 6.4939394022668291491 - expmx*(6.0+6.0*x+3.0*x2+x3)
    expmnx = expmx*expmx
    p = p - expmnx*(0.375+0.75*x+0.75*x2+0.5*x3)
    expmnx = expmnx*expmx
    p = p - expmnx*((2.0/27.0)+(2.0/9.0)*x+(1.0/3.0)*x2+(1.0/3.0)*x3)
    return p


@vectorize
def planck_int_221(x):
    if x <= 0.0:
        return 0.0
    x2 = x * x
    x3 = x2 * x
    if x < 3.8709108458308740467:
        p = 43867.0/107290978560589824000.0
        p = -3617.0/202741834014720000.0+x2*p
        p = 1.0/1270312243200.0+x2*p
        p = -691.0/19615115520000.0+x2*p
        p = 1.0/622702080.0+x2*p
        p = -1.0/13305600.0+x2*p
        p = 1.0/272160.0+x2*p
        p = -1.0/5040.0+x2*p
        p = 1.0/60.0+x2*p
        p = -1.0/8.0+x*p
        p = 1.0/3.0+x*p
        p = x3*p
        return p
    expmx = math.exp(-x)
    p = 6.4939394022668291491 - expmx*(6.0+6.0*x+3.0*x2+x3)
    p = p - expmx*expmx*(0.375+0.75*x+0.75*x2+0.5*x3)
    return p


@vectorize
def planck_int_217(x):
    if x <= 0.0:
        return 0.0
    x2 = x * x
    x3 = x2 * x
    if x < 3.5772328849933323604:
        p = 1.0/1270312243200.0
        p = -691.0/19615115520000.0+x2*p
        p = 1.0/622702080.0+x2*p
        p = -1.0/13305600.0+x2*p
        p = 1.0/272160.0+x2*p
        p = -1.0/5040.0+x2*p
        p = 1.0/60.0+x2*p
        p = -1.0/8.0+x*p
        p = 1.0/3.0+x*p
        p = x3*p
        return p
    expmx = math.exp(-x)
    p = 6.4939394022668291491 - expmx*(6.0+6.0*x+3.0*x2+x3)
    p = p - expmx*expmx*(0.375+0.75*x+0.75*x2+0.5*x3)
    return p


@vectorize
def planck_int_213(x):
    if x <= 0.0:
        return 0.0
    x2 = x * x
    x3 = x2 * x
    if x < 3.1958913758529258675:
        p = 1.0/622702080.0
        p = -1.0/13305600.0+x2*p
        p = 1.0/272160.0+x2*p
        p = -1.0/5040.0+x2*p
        p = 1.0/60.0+x2*p
        p = -1.0/8.0+x*p
        p = 1.0/3.0+x*p
        p = x3*p
        return p
    expmx = math.exp(-x)
    p = 6.4939394022668291491
    p = p - expmx*(6.0+6.0*x+3.0*x2+x3)
    p = p - expmx*expmx*(0.375+0.75*x+0.75*x2+0.5*x3)
    return p


@vectorize
def planck_int_209(x):
    if x <= 0.0:
        return 0.0
    x2 = x * x
    x3 = x2 * x
    if x < 2.6732813723804727115:
        p = 1.0/272160.0
        p = -1.0/5040.0+x2*p
        p = 1.0/60.0+x2*p
        p = -1.0/8.0+x*p
        p = 1.0/3.0+x*p
        p = x3*p
        return p
    expmx = math.exp(-x)
    p = 6.4939394022668291491
    p = p - expmx*(6.0+6.0*x+3.0*x2+x3)
    p = p - expmx*expmx*(0.375+0.75*x+0.75*x2+0.5*x3)
    return p


@vectorize
def planck_int_205(x):
    if x <= 0.0:
        return 0.0
    x2 = x * x
    x3 = x2 * x
    if x < 1.8785904766714333082:
        p = 1.0/60.0
        p = -1.0/8.0+x*p
        p = 1.0/3.0+x*p
        p = x3*p
        return p
    expmx = math.exp(-x)
    p = 6.4939394022668291491
    p = p - expmx*(6.0+6.0*x+3.0*x2+x3)
    p = p - expmx*expmx*(0.375+0.75*x+0.75*x2+0.5*x3)
    return p


_planck_int_tol_map = [
    [1.4970167008036724492e-02, planck_int_205],
    [1.5564096579167360457e-03, planck_int_209],
    [3.7610575515626060127e-04, planck_int_213],
    [1.3579560624114267097e-04, planck_int_217],
    [6.2367008414830821880e-05, planck_int_221],
    [4.6393333440372567663e-06, planck_int_321],
    [5.4123517973399470964e-07, planck_int_421],
    [8.6098034553907039518e-08, planck_int_521],
    [1.7204616971946109889e-08, planck_int_621],
]
_planck_int = planck_int_213


def planck_int_set_tolerance(tol, verbose=False):
    """Set the orders used in `planck_int` by tolerance

    Parameters
    ----------
    tol : float
        tolerance (maximum allowed error)
    verbose : bool, optional
        print information, by default False
    """
    global _planck_int
    if tol < 0:
        if verbose:
            print('`tol` is negative, `planck_int` not changed')
        return
    for t, f in _planck_int_tol_map:
        if tol >= t:
            if verbose:
                print('`planck_int` set to `{0}`'.format(f.__name__))
            _planck_int = f
            return
    if verbose:
        print('`tol` smaller than {0} is not available, '.format(t)
              + '`planck_int` set to `{0}`'.format(f.__name__))
    _planck_int = f

def planck_int(x):
    """Compute planck integral
       `$\int_{0}^{x}\frac{x^{\prime3}}{\exp(x^{\prime})-1}dx^{\prime}$`
       using `m`-th order expansion as `exp(-x)->0` at large `x` and `n`-th
       order expansion as `x->0` for small `x`. The ordres can be set by
       `planck_int_set_tolerance`.

    Parameters
    ----------
    x : float, list or ndarray
        Input value(s).

    Returns
    -------
    float, list or numpy ndarray
        Output value(s).
    """
    return _planck_int(x)

def randomize_ionmix(filename, outfilename):
    """Randomizes the data from an existing ionmix file and rewrites it
    to the outfile.

    Parameters
    ----------
    filename : str
        Name of file to randomize.
    outfilename : str
        Name of output file.
    """
    basepttn = "\d{6}"
    exppttn = "E[+-]\d{2}"
    rand_base = lambda mobj: "{0:06}".format(random.randint(0, 999999))
    rand_exp = lambda mobj: "E{0:+03}".format(random.randint(-99, 99))

    with open(filename) as f:
        lines = f.readlines()

    newlines = []
    for line in lines:
        newline = re.sub(basepttn, rand_base, line)
        newline = re.sub(exppttn, rand_exp, newline)
        newlines.append(newline)

    with open(outfilename, 'w') as f:
        f.write("".join(newlines))

def interpDT(arr, dens, temps,
             bcdmin=BC_BOUND, bctmin=BC_BOUND,
             lookup=INTERP_FUNC):
    """
    Depending on the choice for lookup, this function returns an interpolation
    function for values in arr, the density derivative of arr, or the
    temperature derivative of arr.

    What ``interpDT()`` returns is dependent upon the ``lookup`` and
    ``bcdmin``/``bctmin`` arguments:

        ``INTERP_FUNC`` will return an interpolation function for any
        density/temperature point within the range.

        ``INTERP_DFDD`` will return a function for the density derivative
        at any density/temperature point within the range.

        ``INTERP_DFDT`` will return a function for the temperature derivative
        at any density/temperature point within the range.

        ``BC_BOUND`` is the default setting. If an input density/temp is smaller
        than the minimum value in the ``dens`` (or ``temps``) array, then it
        will automatically be set to this minimum value.

        ``BC_EXTRAP_ZERO`` will insert zero points into the ``arr`` and either
        ``dens`` or ``temps`` for ``bcdmin``, ``bctmin`` respectively.

    Parameters
    ----------
    arr : numpy.ndarray
        Data array.
    dens : numpy.ndarray
        Density array.
    temps : numpy.ndarray
        Temperature array.
    bcdmin : int
        Boundary conditions for density.
    bctmin : int
        Boundary conditions for temperature.
    lookup : int
        Type of function to return.

    Examples
    --------
    In order to find a function that interpolates between dens/temp points
    for an IONMIX file that we had previously opened as ``imx``,
    we could use::

        >>> import opacplot2 as opp
        >>> f = opp.utils.interpDT(imx.pion, imx.dens, imx.temps,
        ...                        bcdmin=BC_EXTRAP_ZERO,
        ...                        bctmin=BC_EXTRAP_ZERO,
        ...                        lookup=INTERP_FUNC)
        >>> print(f(0,0)) # We added points at zero.
        0
        >>> print(f(123, 456)) # Density of 123 and temperature of 456.
        1234.5678 # Resulting ion pressure at this dens/temp.
    """

    # Adjust for extrapolation to zero.
    if bcdmin == BC_EXTRAP_ZERO and dens[0] != 0:
        # Density arrays should be 1D.
        dens = np.insert(dens, 0, 0)
        arr = np.insert(arr, 0, 0, axis=0)

    if bctmin == BC_EXTRAP_ZERO and temps[0] != 0:
        # Temp arrays should be 1D.
        temps = np.insert(temps, 0, 0)
        arr = np.insert(arr, 0, 0, axis=1)

    # "When on a regular grid with x.size = m and y.size = n,
    # if z.ndim ==2, then z must have shape (n, m)."
    # arr will have shape (dens.size, temps.size), so we must transpose it.
    f = sp.interpolate.interp2d(dens, temps, arr.T, kind='linear')

    if lookup == INTERP_FUNC:
        def interp_func(func):
            def wrapper(d, t):
                # Deal with data being out of range.
                # Adjust for lower bounds.
                if d < dens[0]:
                    d = dens[0]
                if t < temps[0]:
                    t = temps[0]
                # Adjust for upper bounds.
                if d > dens[-1]:
                    d = dens[-1]
                if t > temps[-1]:
                    t = temps[-1]
                return func(d,t)
            return wrapper

        # Return the wrapper of f that takes care of data being out of range.
        return interp_func(f)

    if lookup == INTERP_DFDD:
        def df_dd(func):
            # Input to df_dd(func) will be a (dens, temp) point.
            def outter_wrap(dens, temp):
                # Fix temp point but let dens point range.
                def inner_wrap(x):
                    return func(x, temp)

                return sp.misc.derivative(
                            inner_wrap,
                            dens, # Evaluate derivative at dens point.
                            dx=(dens*1e-12))

            return outter_wrap

        return df_dd(f)

    if lookup == INTERP_DFDT:
        def df_dt(func):
            # Input to df_dt(func) will be a (dens, temp) point.
            def outter_wrap(dens, temp):
                # Fix temp point but let temps point range.
                def inner_wrap(x):
                    return func(dens, x)

                return sp.misc.derivative(
                            inner_wrap,
                            temp, # Evaluate derivative at temp point.
                            dx=(temp*1e-12))

            return outter_wrap

        return df_dt(f)

    raise ValueError("lookup must be INTERP_FUNC, INTERP_DFDD, or INTERP_DFDT")

class fastInterpDT():
    def __init__(self, eosopac, input=None, x='dens', y='temp', g='groups',
                 **kwargs):
        """
        Fast interpolation for values in EoS/Opacity table

        Parameters
        ----------
        eosopac : str or Eos/Opacity table in `opacplot2`
            - First, if `eosopac` is str, `eosopac` is used as the filename of
            the table.
            - Otherwise, if `eosopac` is instance of class`opacplot2.OpacIonmix,
            class`opacplot2.OpgPropaceosAscii, class`opacplot2.OpgSesame,
            or class`opacplot2.OpgTOPS, it is converted to eos_dict
            - Otherwise, if `eosopac` is dict, it is used as dict
        input : str or None, optional
            if eosopac is a filename with unrecognized extension, `input` is
            used to specity input format, by default None
        x : str, optional
            x variable for interpolation, by default 'dens'
        y : str, optional
            y variable for interpolation, by default 'temp'
        g : str, optional
            groups variable for interpolation, by default 'groups'
        **kwargs: passed to class`opacplot2.OpacIonmix or `toEosDict`

        """
        self.x = x
        self.y = y
        self.g = g
        from interpolation.splines import UCGrid, eval_linear
        from interpolation.splines import extrap_options as xto
        self.eval_linear = eval_linear
        self.xto = xto

        if isinstance(eosopac, str):
            ext_dict = {'.cn4':'ionmix',
                        '.prp':'propaceos',
                        '.ses':'sesame',
                        '.html':'tops',
                        '.tops':'tops',
            }

            # If the input file is compressed, choose the next extension.
            if os.path.splitext(eosopac)[1] == '.gz':
                _, ext = os.path.splitext(os.path.splitext(eosopac)[0])
            else:
                _, ext = os.path.splitext(eosopac)

            # Choose the correct input type based on extension and set input
            if ext in ext_dict.keys():
                input = ext_dict[ext]

            if input == 'ionmix':
                op = opacplot2.OpacIonmix(eosopac, **kwargs)
            elif input == 'propaceos':
                try:
                    from opacplot2 import opg_propaceos
                    op = opg_propaceos.OpgPropaceosAscii(eosopac)
                except ImportError:
                    raise ImportError('You do not have opg_propaceos.')
            elif input == 'sesame':
                try:
                    op = opacplot2.OpgSesame(eosopac,
                                             opacplot2.OpgSesame.SINGLE)
                except ValueError:
                    op = opacplot2.OpgSesame(eosopac,
                                             opacplot2.OpgSesame.DOUBLE)
            elif input == 'tops':
                op = opacplot2.OpgTOPS(eosopac)
            else:
                raise ValueError('Unsupported input format')

        elif isinstance(eosopac, (opacplot2.OpacIonmix,
                                  opacplot2.OpgPropaceosAscii,
                                  opacplot2.OpgSesame,
                                  opacplot2.OpgTOPS)):
            op = eosopac

        elif isinstance(eosopac, dict):
            input = 'dict'
        else:
            raise ValueError('Unsupported input format')

        if input == 'dict':
            eos_dict = eosopac
        elif input == 'ionmix':
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
            imx_conv = {value:key for key, value in imx_conv.items()}
            eos_dict = {}
            for key in op.__dict__.keys():
                if key in imx_conv.keys():
                    eos_dict[imx_conv[key]] = op.__dict__[key]
                else:
                    eos_dict[key] = key
        elif input == 'propaceos':
            eos_dict = op.toEosDict(**kwargs)
        elif input == 'sesame':
            eos_dict = op.toEosDict(**kwargs)
        elif input == 'tops':
            eos_dict = op.toEosDict(**kwargs)
        else:
            raise ValueError('Unsupported input format')

        self.ucgrid = UCGrid(eos_dict[x], eos_dict[y])
        self.eos_dict = eos_dict
        self._funcs = {}

    def __getitem__(self, key):
        """
        Get the interpolation function by key

        Parameters
        ----------
        key : str
            Name of the variable for interpolaion. Valid values:
            - key in `self.eos_dict`
            - Absorption coefficient "alphaa_{ig}" (unit:1/cm) where `ig` is the
            group index (from 1 to number of groups)
            - Emission rate "emr_{ig}" (unit:erg/g/s) where `ig` is the group
            index (from 1 to number of groups)

        Returns
        -------
        func : function
            Interpolation function
        """

        if key in self._funcs.keys(): return self._funcs[key]

        if key in self.eos_dict.keys():
            def func(x, y):
                dats = self.eval_linear(self.ucgrid, self.eos_dict[key],
                                        np.array([x, y]).T.copy(),
                                        self.xto.NEAREST)
                return dats
            self._funcs[key] = func
            return self._funcs[key]

        try:
            mainkey, ig = key.split('_')
        except:
            raise KeyError
        if not ig.isdigit():
            raise KeyError
        ig = int(ig)
        if key != "{0}_{1}".format(mainkey,ig):
            key = "{0}_{1}".format(mainkey,ig)
            if key in self._funcs.keys(): return self._funcs[key]

        if mainkey == 'alphaa':
            def func(x, y):
                opac = self.eval_linear(self.ucgrid,
                                        self.eos_dict['opr_mg'][:,:,ig-1],
                                        np.array([x, y]).T.copy(),
                                        self.xto.NEAREST)
                alphaa = opac * np.array(x)
                return alphaa
            self._funcs[key] = func
            return self._funcs[key]

        if mainkey == 'emr':
            def func(x, y):
                t = np.array(y)
                xgl = self.eos_dict[self.g][ig-1] / t
                xgr = self.eos_dict[self.g][ig]   / t
                dp = planck_int(xgr) - planck_int(xgl)
                opac = self.eval_linear(self.ucgrid,
                                        self.eos_dict['emp_mg'][:,:,ig-1],
                                        np.array([x, y]).T.copy(),
                                        self.xto.NEAREST)
                emr = emis_const * opac * t**4 * dp
                return emr
            self._funcs[key] = func
            return self._funcs[key]

        raise KeyError

class EosMergeGrids(dict):
    """This class provides filtering capabilities for the EoS temperature and
    density grids.

    For instance, SESAME tables may have some additional points
    in the ion EoS table, compared to the electron EoS table, and as
    FLASH requires the same density and temperature grid for all species,
    the simplest solution is to remove those extra points.

    Parameters
    ----------
    eos_data : dict
        Dictionary contraining the EoS data.
    intersect : list
        The resulting temperature [eV] and density [g/cm^(-3)]
        grids will be computed as an intersection of grids of all the
        species given in this list. Default: ['ele', 'ioncc']
    filter_dens : function
        A function that takes a grid of densities
        and returns a mask of points we don't wont to keep.
        Defaut: (lamdba x: x>0.) i.e. don't remove anything.
    filter_temps : function
        A function that takes a grid of temperatures
        and returns a mask of points we don't wont to keep.
        Defaut: (lamdba x: x>0.) i.e. don't remove anything.
    thresh : list
        Zero threshold on following keys

    Returns
    -------
    out : dict
        A dictionary with the same keys as eos_data. The species specified by
        ``intersect`` will have equal temperature and density grids.

    Examples
    --------
    >>> eos_sesame = opp.OpgSesame("../sesame/xsesame_ascii",
                                   opp.OpgSesame.SINGLE,verbose=False)
    >>> eos_data  = eos_sesame.data[3720]  # Aluminum
    >>> eos_data_filtered = EosMergeGrids(eos_data,
            intersect=['ele', 'ioncc'],   # Merge ele and ioncc grids
            filter_temps=lamda x: x>1.) # Remove temperatures below 1eV
    """
    def __init__(self, eos_data, filter_dens=lambda x: x>=0.,
                 filter_temps=lambda x: x>=0., intersect=['ele', 'ioncc'],
                 thresh=[], qeos=False):

        user_filter = dict(temps=filter_temps, dens=filter_dens)
        self.origin = eos_data
        self.threshold = thresh
        # Computing the intersection of 'ele' and 'ion' grids.
        i_grids = {}
        for key in ['dens', 'temps']:
            i_grids[key] = eos_data['_'.join((intersect[0],key))]
            for el in intersect[1:]:
                i_grids[key] = np.intersect1d(
                                    i_grids[key],
                                    eos_data['_'.join((el,key))])

        # Defining indexes we want to keep.
        mask = {}
        if qeos:
            species_list = ['ele', 'ion', 'total']
        else:
            species_list = ['ele', 'ion', 'total', 'cc', 'ioncc']
        for species in species_list:
            for var in ['dens', 'temps']:
               key = species + '_' + var
               # mask based on the intersection of 'ele' and 'ion' grids
               mask[key] = np.in1d(eos_data[key],
                                   i_grids[var],
                                   assume_unique=True)
               # mask from user provided parameters
               mask[key] = mask[key]*user_filter[var](eos_data[key])
        self.mask = mask
        # Initalising dictionary.
        for key in eos_data:
            self[key] = None # The actual values returned by __getitem__().
        return

    def _get_mask(self, key):
        if '_dens' in key or '_temps' in key:
            return self.mask[key]
        elif any([key.endswith(prefix) for prefix in ['pres', 'eint', 'free']]):
            species = key.split('_')[0]
            indexes = np.meshgrid(
                        np.nonzero(self.mask[species+'_dens'])[0],\
                        np.nonzero(self.mask[species + '_temps'])[0])
            return indexes[0].T, indexes[1].T

    def __getitem__(self, key):
        if key in self.origin:
            if '_dens' in key or '_temps' in key:
                mask = self.mask[key]
                return self.origin[key][mask]
            elif '_ndens' in key:
                return len(self[key.replace('_n', '_')])
            elif '_ntemp' in key:
                # there is an incoherence between '_dens' -> '_ndens'
                # and '_temps' -> '_ntemps' that should really be fixed
                return len(self[key.replace('_n', '_')+'s'])
            elif any([key.endswith(word) for word in ['pres', 'eint', 'free']]):
                data =  self.origin[key][self._get_mask(key)]
                if key in self.threshold:
                    return np.fmax(data, 0)
                else:
                    return data
            else:
                return self.origin[key]
        else:
            # Now just in case we have added some extra keys in there,
            # reproduce a normal dict's behaviour
            return dict.__getitem__(self, key)

def intersect_1D_sorted_arr(arr_1, arr_2):
    """
    Function to return the venn diagram of two sorted 1D arrays.

    In other words, this function will return the union of all values from
    both arrays that are within the intersection of their ranges. This may be
    used for interpolation schemes.

    Parameters
    ----------
    arr_1 : numpy.ndarray
        1D sorted array number 1.
    arr_2 : numpy.ndarray
        1D sorted array number 2.

    Returns
    -------
    numpy.ndarray
        An array that is the venn diagram of the two input arrays.

    Examples
    --------
    >>> import numpy as np
    >>> import opacplot2 as opp
    >>> a = np.array([1,2,3,4,5])
    >>> b = np.array([3.5,4.5,5.5,6])
    >>> c = opp.utils.intersect_1D_sorted_arr(a,b)
    >>> print(c)
    [3.5 4 4.5 5]
    """
    # Check for some overlap.
    if (arr_1[-1] < arr_2[0]) or (arr_2[-1] < arr_1[0]):
        return None

    # Use interpolation to account for mismatched grid sizes.
    # We will also work with the intersection of the dens/temp grids.
    arr_min = max(arr_1[0], arr_2[0])
    arr_max = min(arr_1[-1], arr_2[-1])

    max_idx_1 = np.argmin(arr_1 <= arr_max)
    # If arr_max is bigger than arr_1[-1] max_idx_1 = 0, so we fix that here.
    if max_idx_1 == 0:
        max_idx_1 = len(arr_1)
    min_idx_1 = np.argmax(arr_1 >= arr_min)

    max_idx_2 = np.argmin(arr_2 <= arr_max)
    if max_idx_2 == 0:
        max_idx_2 = len(arr_2)-1
    min_idx_2 = np.argmax(arr_2 >= arr_min)

    sliced_arr_1 = arr_1[min_idx_1:max_idx_1]
    sliced_arr_2 = arr_2[min_idx_2:max_idx_2]

    merged_arr = np.concatenate((sliced_arr_1, sliced_arr_2))

    return np.unique(merged_arr)

################################################################################
# Functions and classes below this have not been adequately tested nor         #
# documented.                                                                  #
################################################################################

def avgopac(energies_in, opacs_in, trad, ebnds,
            weight="constant", bound="error"):
    """

    Parameters
    ----------
    energies_in :

    opacs_in :

    trad :

    ebnds :

    weight='constant' :

    bound='error' :

    """
    try:
        from scipy.integrate import quad
    except ImportError:
        print('Error: Scipy not installed, cannot caclulate opacity integrals!')
        raise
    except:
        raise
    # Check for errors:

    # Make sure that none of the energy group boundaries is outside of
    # the energy range for this opacity:
    energies = energies_in.copy()
    opacs = opacs_in.copy()

    if bound == "error":
        for en in ebnds:
            if en < energies[0] or en > energies[-1]:
                raise ValueError('Energy outside'
                                 'of range {0} {1}'.format(en, energies[0]))
    elif bound == "continue":
        emin = np.min(ebnds)
        emax = np.max(ebnds)

        energies = np.empty(len(energies)+2)
        energies[0] = emin
        energies[1:-1] = energies_in[:]
        energies[-1] = emax

        opacs = np.empty(len(opacs)+2)
        opacs[0] = opacs_in[0]
        opacs[1:-1] = opacs_in[:]
        opacs[-1] = opacs_in[-1]

    else:
        raise ValueError("Illegal boundary treatment")

    def op(en):

        idx = energies.searchsorted(en)
        # The energy is between index idx and idx-1

        de = energies[idx] - energies[idx-1]
        do = opacs[idx] - opacs[idx-1]

        ans = (en-energies[idx-1]) * do/de + opacs[idx-1]
        return ans

    def weight_constant(en):
        return 1.0

    def weight_planck(en):
        x = en/trad
        if x > 700.0: return 1.0
        return en**3/(1-exp(-x)) * exp(-x)

    # def weight_rosseland(en):



    # Choose the weight function:
    f = lambda en: weight(en)*op(en)
    if weight == "constant":
        weight = weight_constant
    elif weight == "planck":
        weight = weight_planck
    elif weight == "rosseland":
        weight = weight_rossland
        f = lambda en: weight(en)/op(en)

    opavg = np.empty(len(ebnds)-1)

    for i in range(len(opavg)):
        numerator = quad(f, ebnds[i], ebnds[i+1], limit=500, epsrel=0.001)
        denominator = quad(weight, ebnds[i],
                           ebnds[i+1], limit=500,
                           epsrel=1.0e-06)
        opavg[i] = numerator[0]/denominator[0]

    return opavg

def ensure_monotonicity(dens, temp, table_in, axis='dens'):
    table = copy.deepcopy(table_in)
    if axis == "dens":
        X = dens
        Y = temp
    elif axis== "temp":
        X = temp
        Y = dens
        table = table.T
    print("Assuring monotonicity", end='')
    for y_idx in range(1,len(Y)):

        for x_idx in range(1,len(X)):
            df = table[x_idx, y_idx] - table[x_idx-1, y_idx]
            if df<0.0:
                print('.', end='')
                table[x_idx, y_idx] = table[x_idx-1, y_idx] + 1e-9
    print('')
    if axis == "temp":
        table = table.T

    return table

class CheckEosConsistency:
    def __init__(self, eos):
        self.fail = 0
        self.num_tests = 0
        self.eos = eos
        self.check_pos()
        self.check_sound_speed()
        self.check_heat_capacity()
        if not self.fail:
            print('Sucess: passed {0}/{0} tests !'.format(self.num_tests))
        else:
            print('Failure: {0}/{1} tests failed!'.format(self.fail, self.num_tests))

    def check_pos(self):
        res = True
        for spec in ['ele', 'ioncc']:
            for tab_name in ['pres', 'eint']:
                tab = self.eos['_'.join([spec, tab_name])]
                self.num_tests += 1
                if np.any(tab < 0):
                    print("{0}_{1} table has negative values".format(spec, tab_name))
                    bad_idx = np.nonzero(tab<0)
                    print('        dens             temp           {0}_{1}'.format(spec, tab_name))
                    print(np.array([self.eos[spec+'_dens'][bad_idx[0]],
                                   self.eos[spec+'_temps'][bad_idx[1]],
                                   self.eos['_'.join([spec, tab_name])][bad_idx]]).T)
                    self.fail +=1
        return res

    def check_sound_speed(self):
        self._check_deriv('ioncc_pres', "dens")
        self._check_deriv('ele_pres', "dens")
    def check_heat_capacity(self):
        self._check_deriv('ioncc_eint', "temp")
        self._check_deriv('ele_eint', "temp")
    def _check_deriv(self, tab_name, axis):
        tab = self.eos[tab_name]
        ax_idx = {'dens': 0, 'temp': 1}[axis]
        diff = np.diff(tab, axis=ax_idx)
        self.num_tests += 1
        if np.any(diff<0):
            print('d{{{0}}}/d{{{1}}} has negative values'.format(tab_name, axis))
            self.fail +=1

def eint_offset(table):
    if np.any(table<0):
        print("Ensuring eint is positive. Adding offset of {0:.3e}".format(-table.min()))
        table = table + np.abs(table.min()) + 1e-9
    return table

def interp_isochores_1d(eos, table='ele', ref_grid='ioncc'):
    #from scipy.interpolate import interp1d
    eosi = copy.deepcopy(eos)
    temp_mask = np.array([temp_el not in eos[table+'_temps']\
                        for temp_el in eos[ref_grid+'_temps']], dtype='bool')

    eos[table+'_temps'] = eos[ref_grid+'_temps']
    eos[table+'_pres'] = np.zeros((len(eos[table+'_dens']), len(eos[table+'_temps'])))

    eos[table+'_eint'] = np.zeros(eos[table+'_pres'].shape)

    for dens_idx  in range(len(eos[table+'_dens'])):
        for par in ['pres', 'eint']:
            # copy values that are on the same grid
            eos[table+'_'+par][dens_idx, ~temp_mask] =  eosi[table+'_'+par][dens_idx, :]

            # interpolate the couple of values that are not there
            #itp = interp1d(eosi[table+'_temps'], eosi[table+'_'+par][dens_idx])
            #eos[table+'_'+par][dens_idx, temp_mask] =  itp(eos[table+'_temps'][temp_mask])
            eos[table+'_'+par][dens_idx, temp_mask] =  np.interp(eos[table+'_temps'][temp_mask],
                                                        eosi[table+'_temps'],
                                                        eosi[table+'_'+par][dens_idx])
    return eos
