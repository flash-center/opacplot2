import numpy as np
from numbers import Number


class opg_tops():
    def __init__(self, fname, epmax='auto', lowerceiling=False, usenear=False):
        """_summary_

        Parameters
        ----------
        fname : str
            The name of the file to open.
        epmax : str or float, optional
            if epmax == 'log', then use logarithmic extrapolation
            if epmax == 'linear', then use linear extrapolation
        lowerceiling : bool, optional
            _description_, by default False
        usenear : bool, optional
            _description_, by default False
        """
        self.fname = fname
        self.NT = 0
        self.Nd = 0
        self.Ng = 0
        self.dens = []
        self.temp = []
        self.grps = []
        self.ross = []
        self.plnk = []
        self.epmax = epmax
        self.lowerceiling = lowerceiling
        self.usenear = usenear
        self.get_data(fname)

    def get_data(self, fname):

        fid = open(fname, 'r')

        line = fid.readline()
        dats = [int(s) for s in line.split() if s.isdigit()]
        self.NT = dats[0]
        self.Nd = dats[1]

        # Read up to temp grid
        for line in fid:
            if line.split()[0] == 'Temperature':
                break

        # Read temperature grid
        self.temp = np.zeros(self.NT)
        i = 0
        for line in fid:
            lb = line.split()
            if lb[0] == 'Density':
                break
            for t in lb:
                self.temp[i] = float(t)
                i += 1
        assert i == self.NT

        # Read density grid
        self.dens = np.zeros(self.Nd)
        i = 0
        for line in fid:
            lb = line.split()
            if lb[0] == 'Photon':
                dats = [int(s) for s in line.split() if s.isdigit()]
                self.Ng = 1 + dats[0]
                break
            for t in lb:
                self.dens[i] = float(t)
                i += 1
        assert i == self.Nd

        # Read photon grid
        self.grps = np.zeros(self.Ng)
        i = 0
        for line in fid:
            lb = line.split()
            if lb[0] == 'Rosseland':
                break
            for t in lb:
                self.grps[i] = float(t)
                i += 1
        assert i == self.Ng - 1

        # Last group boundary not included in tops table, use epmax
        if isinstance(self.epmax, Number):
            self.grps[Ng-1] = self.epmax
        elif self.epmax == 'auto':
            corr_log = np.corrcoef(range(self.Ng-1),
                                   np.log(self.grps[:self.Ng-1]))[0, 1]
            corr_linear = np.corrcoef(range(self.Ng-1),
                                      self.grps[:self.Ng-1])[0, 1]
            if corr_log > corr_linear:
                self.epmax = 'log'
            else:
                self.epmax = 'linear'

        # Do logrithmic or linear extrapolation if requested
        if self.epmax == 'log':
            fit = np.polyfit(range(self.Ng-1),
                             np.log(self.grps[:self.Ng-1]), 1)
            self.grps[self.Ng-1] = np.exp(fit[1] + fit[0] * (self.Ng-1))
        elif self.epmax == 'linear':
            fit = np.polyfit(range(self.Ng-1), self.grps[:self.Ng-1], 1)
            self.grps[self.Ng-1] = fit[1] + fit[0] * (self.Ng-1)
        else:
            raise ValueError("Invalid value for epmax, use float or one \
string from 'auto', 'log', 'linear'")

        # Skip to multigroup opacities
        self.ross = np.zeros((self.Nd, self.NT, self.Ng-1))
        self.plnk = np.zeros((self.Nd, self.NT, self.Ng-1))
        for line in fid:
            if line.split()[0] == 'Multigroup':
                break
        for t in range(self.NT):
            for d in range(self.Nd):
                line = fid.readline()
                assert line.split()[0] == 'Energy'
                for g in range(self.Ng-1):
                    line = fid.readline()
                    lb = line.split()
                    self.ross[d, t, g] = float(lb[1])
                    self.plnk[d, t, g] = float(lb[2])
        if self.lowerceiling:
            ceilingvalue = self.ross[self.ross < 1e10].max()
            self.ross[self.ross == 1e10] = ceilingvalue
            ceilingvalue = self.plnk[self.plnk < 1e10].max()
            self.plnk[self.plnk == 1e10] = ceilingvalue
        if self.usenear:
            for t in range(self.NT):
                for d in range(self.Nd):
                    for g in range(-2, -self.Ng, -1):
                        if self.ross[d, t, g] == 1e10:
                            self.ross[d, t, g] = self.ross[d, t, g+1]
                        if self.plnk[d, t, g] == 1e10:
                            self.plnk[d, t, g] = self.plnk[d, t, g+1]

        fid.close()
