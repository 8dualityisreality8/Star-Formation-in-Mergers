# -*- coding: utf-8 -*-
"""
make plots showing the radial profiles of delta ssfr alongside the delta ssfr
from the nucleus of the paired galaxies. Can split the data by stellar mass
and galaxy separation.
"""
from astropy.io import fits
import numpy as np
import os, sys
import functions as f
import matplotlib.pyplot as plt
import plotstyle as ps
import astropy.stats as astat

""""""""""""
cal = 'ssfr'
bining = ['dm', 'sep']
""""""""""""
# set filepath
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
#T5path = f.T5path()
T5path = '/Volumes/Ext_drive/photoObj8'

# Open drpall
drpallpath = T5path + '/MPL-8/drpall_uniq.fits'
drpall = fits.getdata(drpallpath)
cmd = drpall['CMD']
drmass = np.log10(drpall['nsa_elpetro_mass']/0.7**2)
nsa_z = drpall['nsa_z']
rad = drpall['nsa_elpetro_th50_r']

mngtarg1 = drpall['mngtarg1']
gal = 1 * (((mngtarg1 & 2L**10) != 0) & ((mngtarg1 != 0))) + 2 * (((mngtarg1 & 2L**11) != 0) & ((mngtarg1 != 0)))

yvals = []
stds = []
masses = []
dms = []
seps = []

nuc_yvals = []
nuc_stds = []
nuc_masses = []
nuc_dms = []
nuc_seps = []

gals = []

# Make folder for figures
outpath = filepath + '/14.dexplot_nucplot'
if not os.path.exists(outpath):
    os.makedirs(outpath)
# Open the individual profiles    
sam = fits.getdata(filepath + '/4.samples.fits')

ap = ['nsa', 'sim']
for k in range(len(ap)):
    pairs = sam[ap[k]+'_pair'].astype('bool')
    
    if ap[k]=='nsa' or ap[k]=='sim':    
        stacked = fits.getdata(filepath + '/8.controls_'+ap[k]+'/control_'+ap[k]+'.fits')
        data = stacked[cal+'_dex']
        std = stacked[cal+'_std']
    
    dm = np.abs(sam[ap[k]+'_dm'])
    if ap[k]=='sim':
        dm = np.ones(dm.shape) * -9999
    sep = sam[ap[k]+'_sepkpc']

    mass = sam['all_mass']
    pairs = pairs
    
    gradbin = fits.getdata(filepath + '/3.gradbin/gradbin_nsa.fits')  
    rbins = gradbin['reff'][0]
    
    yvals.append(data[pairs])
    stds.append(std[pairs])
    masses.append(mass[pairs])
    dms.append(dm[pairs])
    seps.append(sep[pairs])
#______________________________________________________________________________
    # nuc files
    nuc = fits.getdata(filepath + '/11.nuc_dex/nuc_dex_'+ap[k]+'.fits')
    
    nucdata = nuc[cal + '_dex'][pairs]
    std = nuc[cal + '_std'][pairs]
    
    if ap[k]=='all':
        savedmass = mass
        savedcal = nucdata
    
    # Group for alls
    nuc_yvals.append(nucdata)
    nuc_stds.append(std)
    nuc_masses.append(mass[pairs])
    nuc_dms.append(dm[pairs])
    nuc_seps.append(sep[pairs])
    
    gals.append(gal[pairs])

data = np.concatenate((yvals[0], yvals[1]))
mass = np.concatenate((nuc_masses[0], nuc_masses[1]))
dm = np.concatenate((nuc_dms[0], nuc_dms[1]))
sep = np.concatenate((nuc_seps[0], nuc_seps[1]))
nuc_data = np.concatenate((nuc_yvals[0], nuc_yvals[1]))
std = np.concatenate((nuc_stds[0], nuc_stds[1]))

gal = np.concatenate((gals[0], gals[1]))
#______________________________________________________________________________

profile, pstd = f.stack(data, error='sem')
med = np.nanmean(nuc_data[~np.isnan(nuc_data)])
#______________________________________________________________________________
"""""""""""""""
Separate Plots
"""""""""""""""
bining = ['mass', 'dm', 'sep']
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
fontsize = 30
for j in range(len(bining)):
    fig, axs = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(16, 8))
    if bining[j] == 'mass':
        label = 'log(M/M$_\odot$) = '
        binsize = 0.5
        mbins = np.arange(9.0, 11.5, binsize) + binsize/2
        z = np.concatenate((masses[0], masses[1]))
    elif bining[j] == 'dm':
        label = '$\Delta$log(M) = '

        binsize = 0.5
        mbins = np.arange(-1.0, 1.0, binsize) + binsize/2
        
        z = np.concatenate((dms[0], dms[1]))        
    elif bining[j] == 'sep':
        label = 'Separation [kpc] = '
        binsize = 12.5
        mbins = np.arange(0.0, 50, binsize) + binsize/2
        z = np.concatenate((seps[0], seps[1]))
        
    dexbin, dexstd = f.stack(data, z=z, zbin=mbins, error='sem')
    dexbinf, dexstdf = f.stack(data, z=z, zbin=mbins, error='sem', trim=0.5)
    
    # Dex plot
    axs[0].plot(rbins, profile, c='k', zorder=9999, linestyle='--', linewidth=2.0)
    for i in range(len(mbins)):
        n = len(data[np.abs(mbins[i]-z) <= binsize/2])
        if len(dexbin[i])>1:
            #axs[0].plot(rbins, dexbin[i], c=colors[i], linestyle='--')
            axs[0].plot(rbins, dexbin[i], c=colors[i], label=label+str(mbins[i]-binsize/2)+'$-$'+str(mbins[i]+binsize/2)+' ('+str(n)+')')
            axs[0].fill_between(rbins, dexbin[i]-dexstd[i], dexbin[i]+dexstd[i], alpha=0.2,
                             antialiased=True)
            
    leg = axs[0].legend(handlelength=0, handletextpad=0, fontsize=fontsize/2.0, loc=3, fancybox=False, edgecolor='k')
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())
    for item in leg.legendHandles:
        item.set_visible(False)
    leg.get_frame().set_linewidth(2.0)
    
    axs[0].plot([-0.1,3.0], [0,0], c='k', linestyle=':') 
    axs[0].axvline(1.5, linestyle='-.', c='k')
    
    axs[0].set_xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)
    axs[0].set_xlim(-0.1,2.7)
    
    axs[0].set_ylabel(r'$\rm \Delta log(sSFR, dex)$', fontsize=fontsize)
    rng = [-11.15, -9.35]

    axs[0].set_ylim(-0.65, 0.65)
    loc1 = 0.5
    loc2 = 0.1
#______________________________________________________________________________
    """ Nuc Data"""
    if bining[j] == 'mass':
        label = 'log(M/M$_\odot$) = '
        binsize = 0.5
        mbins = np.arange(9.0, 11.5, binsize)+binsize/2
        z = mass
    elif bining[j] == 'dm':
        binsize = 0.5
        mbins = np.arange(-1.25, 1.75, binsize) + binsize/2
        #mbins = np.arange(-1.5, 2.0, binsize) + binsize/2
        label = '$\Delta$log(M) = '
        z = dm
    elif bining[j] == 'sep':
        binsize = 10.0
        mbins = np.arange(0.0, 50, binsize) + binsize/2
        label = 'Separation [kpc] = '
        z = sep
#______________________________________________________________________________    
    
    axs[1].errorbar(z[0:dms[0].shape[0]], nuc_data[0:dms[0].shape[0]], yerr=std[0:dms[0].shape[0]], fmt='s', c='deepskyblue', label='Outside-IFU Pairs')
    if bining[j]!= 'dm':
        axs[1].errorbar(z[dms[0].shape[0]::], nuc_data[dms[0].shape[0]::], yerr=std[dms[0].shape[0]::], fmt='s', c='r', label='Inside-IFU Pairs')
        
    yrng = [-1.45, 1.45]
    axs[1].set_ylim(yrng)
    axs[1].set_ylabel(r'$\rm \Delta log(sSFR, dex)$', fontsize=fontsize)
    
    if bining[j] == 'mass':
        axs[1].plot([8.5, 11.75], [med,med], 'k--', linewidth=2.0, zorder=9999)
        axs[1].plot([8.5, 11.75], [0,0], 'k:', linewidth=2.0, zorder=9999)
        ps.ticks(axs[1], xmajor=0.5, ymajor=0.5, xminor=0.1, yminor=0.1)
        axs[1].set_xlim(8.75, 11.75)
        
        axs[1].set_xlabel(r'log(M/M$_\odot$)', fontsize=fontsize)
    elif bining[j] == 'dm':
        axs[1].set_xlim(-1.6, 2.2)
        axs[1].plot([-2.6, 2.6], [med,med], 'k--', linewidth=2.0, zorder=9999)
        axs[1].plot([-2.6, 2.6], [0,0], 'k:', linewidth=2.0, zorder=9999)
        axs[1].plot([0,0], yrng, 'k:', linewidth=2.0, zorder=9999)
        ps.ticks(axs[1], xmajor=1.0, ymajor=0.5, xminor=0.2, yminor=0.1)
        
        axs[1].set_xlabel(r'$\Delta$log(M)', fontsize=fontsize)
    elif bining[j] == 'sep':
        axs[1].set_xlim(-5.0, 55)
        axs[1].plot([-5.0,55], [med,med], 'k--', linewidth=2.0, zorder=9999)
        axs[1].plot([-5.0,55], [0,0], 'k:', linewidth=2.0, zorder=9999)
        ps.ticks(axs[1], xmajor=10.0, ymajor=0.5, xminor=2.0, yminor=0.1)
        
        axs[1].set_xlabel('Separation [kpc]', fontsize=fontsize)
            
    for i in range(len(mbins)):
        ind = (z > mbins[i] - binsize/2.0)&(z < mbins[i] + binsize/2.0)#&(dm!=-9999)
        y = nuc_data[ind]
        yerr = std[ind]
        
        fin = np.isfinite(y)
        y = y[fin]
        yerr = yerr[fin]
        
        ystd = np.nanstd(y) / len(y)**0.5
        
        if len(y) > 0:
            #sig = astat.SigmaClip(sigma=3)
            #y = sig(y)
            #yold = y
            #ystd = astat.biweight_scale(y[np.isfinite(y)], c=6) / len(y)**0.5
            #y = astat.biweight_location(y[np.isfinite(y)], c=6)
            y = np.nanmean(y[np.isfinite(y)])
            x = mbins[i]
            axs[1].errorbar(x, y, xerr=binsize/2.0, yerr=ystd, c='k', ms=10, marker='D', zorder=9999, capsize=5.0, capthick=2.0, elinewidth=2.0)
        
    leg = axs[1].legend(handlelength=0, fontsize=fontsize/2.0, loc=3, fancybox=False, edgecolor='k', borderpad=0.5)
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())
    for item in leg.legendHandles:
        item.set_visible(False)
    leg.get_frame().set_linewidth(2.0)
#______________________________________________________________________________
    
    plt.subplots_adjust(wspace=0.3)
    ps.ticks(axs[0], xmajor=0.5, ymajor=loc1, xminor=0.1, yminor=loc2)
    ps.style(axs, fontsize=fontsize)
    
    plt.savefig(outpath + '/ssfr_'+bining[j]+'.pdf', bbox_inches='tight', overwrite=True)
    plt.close('all')