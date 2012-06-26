import matplotlib
from matplotlib import pyplot as plt
import numpy as np
plt.rcParams.update({'text.usetex': True})

def plot_zbar(denss, temps, zbar, zmax, fig):
    ax = fig.add_subplot(111)
    x,y = np.meshgrid(denss, temps)

    cs = ax.contourf(x,y, zbar.T, 256)
    cb = plt.colorbar(cs, ticks=matplotlib.ticker.MaxNLocator(nbins=15))
    cb.set_label("Average Ionization")

    cb2 = plt.contour(x, y, zbar.T,
           colors='k', opacity=0.5, linewidths=0.4,
           locator=matplotlib.ticker.MaxNLocator(nbins=15))

    plt.clabel(cb2,fontsize=6, inline=False, inline_spacing=1, fmt='%1.0f',
               rightside_up=True, use_clabeltext=False)

    ax.loglog()
    ax.set_xlim((denss[0], denss[-1]))
    ax.set_ylim((temps[0], temps[-1]))
    plt.xlabel(r'$\rho$ [g.cm$^{-3}$]')
    plt.ylabel(r'$T$ [eV]')
    return ax

def plot_eos_grid(tdata, var):
    """
    Plot EoS ρ and T grids for different species.
    Parameters:
    -----------
     - tdata [dict]: dictionary containg EoS data.
     - var [str] : select which grid to plot 'dens' or 'temps'.
    
    Examples:
    ---------
    see visualization/eos.ipynb
    """
    fig = plt.figure(figsize=(9,2))
    bottom = 0
    ax = plt.subplot(111)
    linthreshx = 1.0
    keys_list = sorted([el for el in tdata.keys() if '_'+var in el])
    for key in keys_list:
       plt.bar(tdata[key][:-1], np.ones(len(tdata[key])-1), width=np.diff(tdata[key]),
             bottom=bottom, color='w', alpha=0.5)
       bottom += 1
       if len(tdata[key])>1:
           linthreshx = min(linthreshx, tdata[key][1])
    ax = plt.gca()
    ax.set_yticks(np.arange(bottom)+0.5)
    ax.set_yticklabels([el.replace('_', '\_') for el in keys_list])
    ax.set_xscale('symlog', linthreshx=0.8*linthreshx)
    if var == 'dens':
        plt.title(r'Density grids: $\rho$ [g.cm$^{-3}$]')
    else:
        plt.title(r'Temperature grids: $T$ [eV]')

def plot_eos_field(eos_data, species, parameter, grad=None):
    if grad is True:
        plot_eos_field(eos_data, species, parameter, grad='rho')
        plot_eos_field(eos_data, species, parameter, grad='T')
        return
    rho, T = eos_data[species+'_dens'], eos_data[species+'_temps']
    rho_g, T_g = np.meshgrid(rho, T)
    Z = eos_data['_'.join([species, parameter])].T
    species_label = dict(ele='Electron', ion='Ion')[species]
    clabel = dict(pres=r'$log(P)$ GPa', eint=r'u J.cm${^{-3}}$')[parameter]
    def idty(x):
        return x
    # Label in the title
    var_title = {'pres': {'rho': '= c_s^2', 'T': '', None: 'p'},
                  'eint':{'rho': '', 'T': '= c_v', None: 'e'}}[parameter]
    # Unit scaling and some custom transformation we want to apply to the
    # plotted variable (for instance plotting c_s instead of c_s²)
    var_transform = { 'pres': {
        'rho': lambda x: 1e-3*(1e-6*0.5*(x +np.abs(x)))**0.5, # to km/s
        'T': idty,
         None: lambda x: x*1e-9}, # to Mbars
                      'eint': {
        'rho': lambda x: (species == 'ele' and -1 or 1)*x, # this derivative seems to negative for ele
        'T': idty,
        None: idty}}[parameter]
    # Label of the colorbar
    var_ctitle = {'pres': {'rho': '$c_s$ [km.s$^{-1}$]',
                           'T': 'log scale',
                           None: '$\log(p)$ [Mbar]'},
                  'eint': {'rho': 'log scale', 'T':'log scale', 
                           None: '$\log(e)$'}}[parameter]
    grad_label = {'rho': r'\rho', 'T': 'T', None: ''}

    if grad is None:
        fig = plt.figure(figsize=(6,4.5))
        ax = plt.subplot(111)
        # Plot the 
        if Z[1:,1:].min() <= 0:
            zmin = 1e-15
        else:
            zmin = Z[1:,1:].min()
        Z = var_transform[grad](Z) # applying a unit scaling
        levs = np.arange(np.floor(np.log10(zmin)-1),
                       np.ceil(np.log10(Z.max())+1), 1)
        cs = plt.contourf(rho, T, np.log10(np.fmax(Z, zmin)),
                          256, alpha=0.9)
        plt.contourf(rho, T, (Z>=0)*1, [0,0.5], colors='w')
        plt.title('{0} EoS - ${1}$'.format(species_label, var_title[grad]))

    elif grad in ['rho', 'T']:
        # computing the partial derivative
        idx =  (grad == 'rho')*1
        if grad == 'rho':
            fig = plt.figure(figsize=(12,4.5))
        ax = plt.subplot(1,2,idx)
        Z = np.gradient(Z)[idx]/np.gradient([T, rho][idx])[idx]
        Z = var_transform[grad](Z) # applying a unit scaling
        cs = plt.contourf(rho,T, np.log10(Z), 256)
        plt.title(r'{0} EoS - $\frac{{\partial {1}}}{{\partial {2}}} {3}$'.format(
                   species_label, var_title[None], grad_label[grad], var_title[grad]))
    cb2 = plt.contour(rho, T, np.log10(Z),
           colors='k', opacity=0.5, linewidths=0.4,
           locator=matplotlib.ticker.MaxNLocator(nbins=15))
    plt.clabel(cb2,fontsize=6, inline=False, inline_spacing=1, fmt='%1.0f',
               rightside_up=True,use_clabeltext=False)
    cb = plt.colorbar(cs, ticks=matplotlib.ticker.MaxNLocator(nbins=15))
    cb.ax.set_ylabel(var_ctitle[grad])
    plt.grid('on')
    plt.xlabel(r'$\rho$ [g.cm$^{-3}$]')
    plt.ylabel(r'$T$ [eV]')

    ax.set_xscale('symlog', linthreshx=0.7*rho[1])
    ax.set_yscale('symlog', linthreshy=0.7*T[1])
