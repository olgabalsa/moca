import os
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from astropy.time import Time
from matplotlib import rcParams
import scripts.rv_functions as rv_functions # type: ignore
from astropy.timeseries import LombScargle


rcParams['font.size'] = 16.0
rcParams['axes.linewidth'] = 1

rcParams['ytick.major.size'] = 6
rcParams['ytick.major.width'] = 1.5
rcParams['ytick.minor.size'] = 3
rcParams['ytick.minor.width'] = 0.5

rcParams['xtick.major.size'] = 6
rcParams['xtick.major.width'] = 1.5
rcParams['xtick.minor.size'] = 3
rcParams['xtick.minor.width'] = 0.5

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica']


main_colors= ['#DAAF11', '#E7C74F', '#C3CA87', '#9CCA87', '#87CAAA', '#7FB3C4', '#628FB6', '#426BB4']
rgb_colors = [mcolors.hex2color(color) for color in main_colors]
position_colors = np.linspace(0, 1, len(rgb_colors))
cmap_moca = mcolors.LinearSegmentedColormap.from_list('custom_cmap', list(zip(position_colors, rgb_colors)))


def save_ascii(mock_timeseries, dic_inputs):
    tn = Time.now().isot[:16]
    full_dir = f'{dic_inputs["dir"]}/moca_outputs/{dic_inputs["star"]}'
    os.makedirs(full_dir, exist_ok = True)
    
    df_moca_output = pd.DataFrame({'jd': mock_timeseries['jd'], 'rv': mock_timeseries['rv'], 'erv': mock_timeseries['erv']})
    df_moca_output.to_csv(f'{full_dir}/moca_timeseries_{tn}.csv', sep = ',', index = False)


def plot_periodogram(jd, rv, erv, per, Prot, Pmag, ax, title, ipl = None, secax = False):
    
    global cmap_moca 

    mcadence = np.max([0.6, np.min([jd[i+1] - jd[i] for i in range(len(jd)-1)])])

    frequency, power = LombScargle(jd, rv, erv).autopower(minimum_frequency = 1/600.0, maximum_frequency = 1/mcadence)
    periodicity = 1/frequency
    FAP01, FAP005, FAP001 = LombScargle(jd, rv, erv).false_alarm_level([0.1, 0.05, 0.01])

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major', labelsize = 16)
    ax.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')

    if secax:
        secax = ax1.secondary_xaxis('top')
        secax.set_xticks([8.5, 29.7]) #, 183, 367])
        secax.set_xticklabels(['8.5 d', '29.7 d'], fontsize = 15) #, 'Year/2', 'Year']

    ax.plot(periodicity, power/FAP01, c = 'k', lw = 1.5, zorder = 1000)
    ax.axhline(1, color = 'k', lw = 1, ls = ':', zorder = 100) # Horizontal line FAP 10%
    
    ymin, ymax = ax.get_ylim()
    ytext = ymin + 0.9 * (ymax - ymin)
    colors_per = [cmap_moca(xpos) for xpos in np.linspace(0, 1, len(per))]
    if not ipl:
        ipl = np.arange(len(per))
    for i,p in enumerate(per):
        if i not in ipl:
            continue
        ax.axvline(p, color = colors_per[i], lw = 3, alpha = 1, zorder = 1)
        ax.text(p, ytext, f'P{i+1}', rotation = 90, verticalalignment = 'center', horizontalalignment = 'right', fontsize = 14, zorder = 10)
    if Prot:
        ax.axvline(Prot, color = 'crimson', lw = 3, alpha = 1, zorder = 1)
        ax.text(Prot, ytext, r'P$_{\rm rot}$', rotation = 90, verticalalignment = 'center', horizontalalignment = 'right', fontsize = 14, zorder = 10)
    if Pmag:
        ax.axvline(Pmag, color = 'crimson', lw = 3, alpha = 1, zorder = 1)
        ax.text(Pmag, ytext, r'P$_{\rm mag}$', rotation = 90, verticalalignment = 'center', horizontalalignment = 'right', fontsize = 14, zorder = 10)

    ax.set_xscale('log')
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.set_ylabel('Power / FAP 0.1', fontsize = 15)
    ax.tick_params(axis = 'both', labelsize = 16)
    ax.set_title(title, pad = 10, fontsize = 16)


def do_plot(mock_timeseries, dic_inputs):

    global cmap_moca 
    jd = mock_timeseries['jd']
    rv = mock_timeseries['rv']
    erv = mock_timeseries['erv']
    num_pl = dic_inputs['num_pl']
    star = dic_inputs['star']
    K = dic_inputs['K']
    per = dic_inputs['per']
    t0 = dic_inputs['t0']
    ecc = dic_inputs['ecc']
    w = dic_inputs['w']
    c = ecc * np.cos(np.radians(w))
    d = ecc * np.sin(np.radians(w))
    Prot = dic_inputs['Prot']
    Pmag = dic_inputs['Pmag']

    tn = Time.now().isot[:16]
    full_dir = f'{dic_inputs["dir"]}/moca_outputs/{star}'

    fig = plt.figure(figsize = (13.9, 10. + 4 * num_pl))
    gs = gridspec.GridSpec(nrows = num_pl + 1, ncols = 2, hspace = 0.2, wspace = 0.05)

    for i in range(num_pl + 1):
        axi = plt.subplot(gs[i, 0])
    
        axi.xaxis.set_minor_locator(AutoMinorLocator())
        axi.yaxis.set_minor_locator(AutoMinorLocator())
        axi.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major', labelsize = 16)
        axi.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
    
    # RV vs JD
    jd_m = round(np.min(jd), -2) - 100
    ax0 = plt.subplot(gs[0, 0])
    ax0.scatter(jd - jd_m, rv, c = 'k', edgecolor = 'w', lw = 0.5, s = 100)
    ax0.errorbar(jd - jd_m, rv, yerr = erv, c = 'k', lw = 2, alpha = 0.4, linestyle = 'none')
    jd_model = np.linspace(np.min(jd), np.max(jd), 200)
    model_jd = rv_functions.rv_model(jd_model, dic_inputs)
    ax0.plot(jd_model - jd_m, model_jd, color = 'k', alpha = 0.3, lw = 3., ls = '-', zorder = 10)

    ax0.set_title(f'Time [JD $-$ {int(jd_m)}]', pad = 10, fontsize = 16)
    ax0.tick_params('x', top = True, bottom = True, which = 'both')
    ax0.set_ylabel(r'$\Delta$RV [m/s]', labelpad = 10)
    
    ax0p = plt.subplot(gs[0, 1])
    plot_periodogram(jd, rv, erv, per, Prot, Pmag, ax0p, 'GLS full timeseries')

    # RV vs phase
    colors_pl = [cmap_moca(xpos) for xpos in np.linspace(0, 1, num_pl)]
    rvc = rv.copy()
    for pl in range(num_pl):
        per_pl = per[pl]
        t0_pl = t0[pl]
        phase = ((np.array(jd) - t0_pl) % per_pl) / per_pl

        # substract the signal of the other planets
        K_npl = np.delete(K, pl)
        per_npl = np.delete(per, pl)
        t0_npl = np.delete(t0, pl)
        c_npl = np.delete(c, pl)
        d_npl = np.delete(d, pl)
        dic_npl = {'num_pl': num_pl - 1, 'K': K_npl, 'per': per_npl, 't0': t0_npl, 'c': c_npl, 'd': d_npl}
        rvi = rvc - rv_functions.rv_model(jd, dic_npl)
        
        jd_model = np.linspace(t0[pl], t0[pl] + per[pl], 200)
        phase_model = np.linspace(0, 1, 200)
        K_pl = np.array([K[pl]])
        per_pl = np.array([per[pl]])
        t0_pl = np.array([t0[pl]])
        c_pl = np.array([c[pl]])
        d_pl = np.array([d[pl]])
        dic_pl = {'num_pl': 1, 'K': K_pl, 'per': per_pl, 't0': t0_pl, 'c': c_pl, 'd': d_pl}
        model_i = rv_functions.rv_model(jd_model, dic_pl)

        axi = plt.subplot(gs[pl + 1, 0])
        axi.scatter(phase, rvi, c = colors_pl[pl], edgecolor = 'w', lw = 0.5, s = 100, zorder = 100)
        axi.errorbar(phase, rvi, yerr = erv, c = colors_pl[pl], lw = 2, alpha = 0.4, linestyle = 'none', zorder = 10)
        axi.plot(phase_model, model_i, color = colors_pl[pl], alpha = 0.5, lw = 3., ls = '-', zorder = 10)

        axi.set_title(f'P{pl + 1} = {np.round(per_pl[0],2)} d', pad = 12, fontsize = 16)
        axi.set_ylabel(r'$\Delta$RV [m/s]', labelpad = 10)
        axi.set_xlim(0, 1)

        axip = plt.subplot(gs[pl + 1, 1])
        plot_periodogram(jd, rvi, erv, per, Prot, Pmag, axip, f'GLS signal P{pl + 1}', [pl])

        if pl == num_pl -1:
            axi.set_xlabel('Orbital phase')
            axip.set_xlabel(r'$P$ [d]', fontsize = 15)

    plt.savefig(f'{full_dir}/plot_mocaRVs_{tn}.pdf', dpi = 300, bbox_inches = 'tight')