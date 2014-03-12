#!/usr/bin/python
# -*- coding: utf-8 -*-
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
    keys_list = sorted([el for el in tdata.keys() if '_'+var in el and el != 'cc_temps'])
    for key in keys_list:
       if len(tdata[key])>1:
           linthreshx = min(linthreshx, tdata[key][1])
       plt.bar(tdata[key][:-1], np.ones(len(tdata[key])-1), width=np.diff(tdata[key]),
             bottom=bottom, color='w', alpha=0.5)
       bottom += 1
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
    species_label = dict(ele='Electron', ion='Ion', ioncc='Ioncc', total='Total')[species]
    clabel = dict(pres=r'$log(P)$ GPa', eint=r'u J.cm${^{-3}}$')[parameter]
    def idty(x):
        return x
    # Label in the title
    var_title = {'pres': {'rho': '= c_s^2', 'T': '', None: 'p'},
                  'eint':{'rho': '', 'T': '= c_v', None: 'e'} }[parameter]
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
                           None: '$\log(e)$'},
                  }[parameter]
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

def plot_diff_mg_opac(fig, op_list, idx=None):

    plt.rcParams.update({'text.usetex': True})
    ax = [plt.subplot(2,2,i) for i in range(1, 5)] 
    groups = op_list[0]['groups']
    plt.suptitle('Multigroup opacity comparaison: {0:.2e} g.cm$^{{-2}}$, {1:.3f} eV'.format(
                        op_list[0]['rho'][idx[0]], op_list[0]['temp'][idx[1]]))
    for var_idx, var in enumerate(['opp_mg', 'opr_mg', 'eps_mg']):
        for op in op_list:
            ax[var_idx].step(groups[:-1], op[var][idx[0], idx[1]], where='post',
                       label=(var== 'opp_mg' and op['label'].replace('_', '\_')))

    def gnorm(a,b):
        return np.fmin(np.fmax(a/b,b/a) - 1., 9.9)
    def anorm(a,b):
        return np.abs(a-b)
    def mls(var):
        if var == 'opr_mg':
            return 'dashed'
        else:
            return 'solid'
    for var_idx, var in enumerate(['opp_mg', 'opr_mg']):
        err = gnorm(op_list[1][var][idx[0], idx[1]],  op_list[0][var][idx[0], idx[1]])
        ax[3].step(groups[:-1], err, label=var.replace('_', '\_'), where='post',
                       ls=mls(var))
    ax3b = plt.twinx(ax[3])
    err = anorm(op_list[1]['eps_mg'][idx[0], idx[1]],  op_list[0]['eps_mg'][idx[0], idx[1]])
    ax3b.step(groups[:-1], err, 'r-', where='post')
    ax3b.set_ylabel(r'eps\_mg', color='r')
    for tl in ax3b.get_yticklabels():
        tl.set_color('r')

    ax[2].set_ylim(-0.05, 1.05)


    ax[3].legend(loc='best', prop={"size": 10})
    ax[0].legend(loc='best', prop={"size": 10})
    ax[0].set_title('Planck multigroup opacity')
    ax[1].set_title('Rosseland multigroup opacity')
    ax[2].set_title("Deviation from Krichoff's law (=1 if LTE)")
    ax[3].set_title('Relative error')
    ax[3].set_ylim(-0.05)
    for ax_idx, axi in enumerate(ax):
        axi.set_xscale('log')
        axi.set_xlim(groups[0], groups[-1])
        axi.grid()
        if ax_idx <2:
            axi.set_yscale('log')
    ax[2].set_xlabel(r'$\nu$ [eV]')
    ax[3].set_xlabel(r'$\nu$ [eV]')
    ax[0].set_ylabel(r'$\kappa_p$ [cm$^2$.g$^{-1}$]')
    ax[1].set_ylabel(r'$\kappa_R$ [cm$^2$.g$^{-1}$]')
    ax[2].set_ylabel(r'$\kappa_p/\eta_p$')
    return fig

def plot_2D_map(fig, op, temp_val):
    from matplotlib.colors import LogNorm
    from hedp.rad import planck

    ax = [plt.subplot(1,3,idx+1) for idx in range(3)]
    groups = op[0]['groups']
    x = op[0]['rho']
    y = op[0]['temp']
    X, Y = np.meshgrid(x, y, indices='ij')
    Bnu = planck(groups[:-1], temp_val, 'nu_eV')
    Bnu /= np.sum(Bnu)
    Bnu_arr = np.tile(Bnu, (len(x), len(y))).reshape(op[0]['opp_mg'].shape)
    plt.suptitle('Opacity comparaison for a black body radiation at {0:.1f} eV'.format(temp_val))
    def norm(a,b):
      return np.fmax(a/b,b/a) - 1.0
    for var_idx, var in enumerate(['opp_mg', 'opr_mg', 'eps_mg']):
        err = norm(op[0][var], op[1][var])
        err = (err* Bnu_arr).sum(axis=-1)

        cs = ax[var_idx].pcolormesh(X, Y, err, norm=LogNorm(vmin=0.1, vmax=1e4))
        cb = plt.colorbar(cs, ax=ax[var_idx])
        ax[var_idx].set_title('Relative error for {0}'.format(var.replace('_','\_')))

    for axi in ax:
        axi.set_xscale('log')
        axi.set_yscale('log')
        axi.set_xlim(x.min(),x.max())
        axi.set_ylim(y.min(),y.max())
        axi.set_xlabel(r'$\rho$ [g.cm$^{-3}$]')
        axi.set_ylabel(r'Te [eV]')
    return fig

def plot_Zbar(fig, op):
    from matplotlib.colors import LogNorm
    ax = [plt.subplot(121), plt.subplot(122)]
    Zmax = np.ceil(op[0]['zbar'].max())
    contour_lvls = np.arange(0, Zmax, np.ceil(Zmax/10.))
    contour_lvls = np.unique(np.sort(np.concatenate((contour_lvls, np.array([1])))))
    
    linestyles  = iter(('solid', 'dashed', 'dashdot', 'dotted'))
    for op_idx, op_el in enumerate(op): 
        X, Y = np.meshgrid(op_el['rho'], op_el['temp'], indices='ij')
        current_ls = linestyles.next()
        cs = ax[0].contour(X, Y, op_el['zbar'].T, contour_lvls,
                    linestyles=current_ls,
                    vmin=0, vmax=contour_lvls.max())
        if op_idx == 0:
            ax[0].clabel(cs, inline=0, fontsize=10, fmt='%2.0f')
        # plot a slice at 10 mg/cc
        idx_rho = np.argmin(np.abs(X[0] - 1.e-2))
        ax[1].plot(Y[:,idx_rho], op_el['zbar'][:, idx_rho],'k', ls=current_ls,
                label=op_el["label"].replace('_', '\_'))


        
    plt.colorbar(cs, ax=ax[0])
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_xlabel(r'$\rho$ [g.cm$^{-3}$]')
    ax[0].set_ylabel('$T_e$ [eV]')
    ax[0].set_title('Mean ionization')
    ax[1].legend(loc='best')
    ax[1].set_xscale('symlog', linthreshx=1)
    ax[1].set_xlabel(r'Temperature [eV]')
    ax[1].set_ylabel('Zbar')
    ax[1].set_title('slice at $10^{-2}$ g.cm$^{-2}$')
    return fig


#idx = np.argmin(np.abs(temp - 0.025))
#idx_imx = np.argmin(np.abs(imx.temps - 0.025))
#idx_mt = np.argmin(np.abs(mt.Te - 0.025))
#
#ax[1].semilogx(dens, zbar_snop[idx,:], 'k', label='SNOP')
#ax[1].semilogx(dens, zbar_tf[idx,:], 'k--', label='Thomas Fermi')
#ax[1].semilogx(mt.rho, mt.opz[:, idx_mt], 'k:', label='SNOP ref: {0}'.format(
#                        os.path.split(fpath)[1]).replace('_', '\_'))
##ax[1].semilogx([], [], 'k:', label='Ionmix')
#ax[1].legend(loc='best')
#
#
#ax[1].set_xlabel(r'$\rho$ [g.cm$^{-3}$]')
#ax[1].set_ylabel(r'Zbar')



