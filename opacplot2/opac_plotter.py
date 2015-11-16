import numpy as np

class OpacPlotter:
    def __init__(self, fig, opacs, mplkwargs = None):
        from matplotlib import pyplot as plt
        plt.rcParams['font.sans-serif'] = 'FreeSans'
        plt.rcParams['font.size'] = 13
        plt.rcParams['font.weight'] = 'bold'
        plt.rcParams['legend.fontsize'] = 13

        self.fig = fig
        self.opacs = opacs
        self.nopacs = len(opacs)
        self.mplkwargs = mplkwargs if mplkwargs != None else self.nopacs*[{}]
        self.draw()

    def draw(self):
        self.fig.clf()
        ax = self.fig.add_subplot(111)

        for n in range(self.nopacs):
            ax.plot(self.opacs[n][0], self.opacs[n][1], **self.mplkwargs[n])

        ax.loglog()
        ax.set_xlabel("Photon Energy (eV)")
        ax.set_ylabel("Opacity (cm^2/g)")
        ax.grid(True)

