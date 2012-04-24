from matplotlib import pyplot as plt
import numpy as np

def plot_zbar(denss, temps, zbar, zmax, fig):
    ax = fig.add_subplot(111)
    x,y = np.meshgrid(denss, temps)

    contour = ax.contourf(x,y,zbar.T, zmax-1)
    colorbar = fig.colorbar(contour)
    colorbar.set_label("Average Ionization")

    ax.contour(x,y,zbar.T, zmax-1, colors="k")

    ax.loglog()
    ax.set_xlim((denss[0], denss[-1]))
    ax.set_ylim((temps[0], temps[-1]))
    ax.set_xlabel("Mass Density (g/cc)")
    ax.set_ylabel("Electron Temperature (eV)")
    return ax
