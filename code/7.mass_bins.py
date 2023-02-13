#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 13:19:24 2019

Plot the gradients for all pairs/controls binned in mass separately for inside
pairs and outside pairs.

Neither of the two samples are large enough by themselves to have statistically
clean profiles.

@author: joshua
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import plotstyle as ps
"""""""""""""""""""""
Which Parameter to Plot
"""""""""""""""""""""
cal = 'ssfr'
#______________________________________________________________________________
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
figpath = '/Users/joshua/Documents/GitHub/Manga_Mergers/'
droppath = '/Users/joshua/Dropbox/manga/figs'

ap = ['nsa', 'sim']
for k in range(len(ap)):
    outpath = filepath + '/7.mass_bins_'+ap[k]
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    
    
    fil = fits.getdata(filepath + '/5.stacking/stacked_trim_'+ap[k]+'.fits')
    
    #binsize=0.5
    #mbins = np.arange(9.0, 11.0, binsize)+binsize/2
    
    mbins = fil['mbins']
    mbins = mbins[mbins!=0]
    binsize = mbins[1] - mbins[0]
    
    rbins = fil['rbins']
    pbin, pstd = fil[cal+'_pbin'].T, fil[cal+'_pstd'].T
    cbin, cstd = fil[cal+'_cbin'].T, fil[cal+'_cstd'].T
    dexbin, dexstd = fil[cal+'_dex'].T, fil[cal+'_dexstd'].T
    
    fil = fits.getdata(filepath + '/5.stacking/stacked_'+ap[k]+'.fits')
    fpbin = fil[cal+'_pbin'].T
    fcbin = fil[cal+'_cbin'].T
    fdexbin = fil[cal+'_dex'].T
#______________________________________________________________________________
    # plotting
    fontsize=30
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=False, figsize=(28, 8))
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    
    # Control plot
    for i in range(len(mbins)):
        if len(cbin[i])>1:
            axs[0].plot(rbins, fcbin[i], linestyle='--', c=colors[i])
    
            axs[0].plot(rbins, cbin[i], c=colors[i], label='log(M/M$_\odot$) = '+str(mbins[i]-binsize/2)+'-'+str(mbins[i]+binsize/2))
            axs[0].fill_between(rbins, cbin[i]-cstd[i], cbin[i]+cstd[i], alpha=0.2,# edgecolor='#1B2ACC', facecolor='#089FFF',
                             antialiased=True)
    
    # Pair plot
    for i in range(len(mbins)):
        if len(np.asarray(pbin[i]).T)>1:
            axs[1].plot(rbins, fpbin[i], c=colors[i], linestyle='--')
            axs[1].plot(rbins, pbin[i], c=colors[i], label='log(M/M$_\odot$) = '+str(mbins[i]-binsize/2)+'-'+str(mbins[i]+binsize/2))
            axs[1].fill_between(rbins, pbin[i]-pstd[i], pbin[i]+pstd[i], alpha=0.2,
                             antialiased=True)
    # Dex plot
    for i in range(len(mbins)):
        if len(dexbin[i])>1:
            axs[2].plot(rbins, fdexbin[i], c=colors[i], linestyle='--')
            axs[2].plot(rbins, dexbin[i], c=colors[i], label='log(M/M$_\odot$) = '+str(mbins[i]-binsize/2)+'-'+str(mbins[i]+binsize/2))
            axs[2].fill_between(rbins, dexbin[i]-dexstd[i], dexbin[i]+dexstd[i], alpha=0.2,
                             antialiased=True)
    
    leg = axs[2].legend(handlelength=0, handletextpad=0, fontsize=fontsize/1.75, loc=4, fancybox=False, edgecolor='k', labelspacing=0.1)
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())
    for item in leg.legendHandles:
        item.set_visible(False)
    leg.get_frame().set_linewidth(2.0)
        
    axs[2].plot([-0.1,3.0], [0,0], c='k', linestyle=':') 
    
    axs[1].set_xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)
    axs[0].set_xlim(-0.1,2.7)
    axs[1].set_xlim(-0.1,2.7)
    axs[2].set_xlim(-0.1,2.7)
    
    if cal=='ewha':
        axs[0].set_ylabel(r'$\rm log(H\alpha$ EW/Ang)', fontsize=fontsize)
        axs[1].set_ylabel(r'$\rm log(H\alpha$ EW/Ang)', fontsize=fontsize)
        axs[2].set_ylabel(r'$\rm \Delta log(H\alpha$ EW/Ang)', fontsize=fontsize)
        if ap[k] == 'nsa':
            rng = [0.75, 1.85]
        elif ap[k] == 'sim':
            rng = [0.75, 1.85]
        axs[0].set_ylim(rng)
        axs[1].set_ylim(rng)
        axs[2].set_ylim(-1*(rng[1]-rng[0])/2, (rng[1]-rng[0])/2)
        loc1 = 0.20
        loc2 = 0.05
    elif cal=='ssfr':
        axs[0].set_ylabel(r'$\rm log(sSFR, yr^{-1})$', fontsize=fontsize)
        axs[1].set_ylabel(r'$\rm log(sSFR, yr^{-1})$', fontsize=fontsize)
        axs[2].set_ylabel(r'$\rm \Delta log(sSFR, yr^{-1})$', fontsize=fontsize)
        
        if ap[k] == 'nsa':
            rng = [-11.15, -9.35]
        elif ap[k] == 'sim':
            rng = [-11.45, -9.65]
        
        axs[0].set_ylim(rng)
        axs[1].set_ylim(rng)
        axs[2].set_ylim(-1*(rng[1]-rng[0])/2, (rng[1]-rng[0])/2)
        loc1 = 0.5
        loc2 = 0.1
    elif cal == 'DOP16':
        axs[0].set_ylabel('12 + log(O/H)', fontsize=fontsize)
        axs[1].set_ylabel('12 + log(O/H)', fontsize=fontsize)
        axs[2].set_ylabel('$\Delta$log(O/H)', fontsize=fontsize)
        rng = [7.95, 9.05]
        axs[0].set_ylim(rng)
        axs[1].set_ylim(rng)
        axs[2].set_ylim(-1*(rng[1]-rng[0])/2, (rng[1]-rng[0])/2)
        axs[0].annotate('D16', xy=(0.05, 0.05), xycoords='axes fraction', fontsize=fontsize, bbox=dict(facecolor='none', edgecolor='k', linewidth=2))
        loc1 = 0.25
        loc2 = 0.05
    elif cal == 'NO':
        axs[0].set_ylabel('log(N/O)', fontsize=fontsize)
        axs[1].set_ylabel('log(N/O)', fontsize=fontsize)
        axs[2].set_ylabel('$\Delta$log(N/O)', fontsize=fontsize)
        rng = [-1.0, 1.0]
        axs[0].set_ylim(rng)
        axs[1].set_ylim(rng)
        axs[2].set_ylim(-1*(rng[1]-rng[0])/2, (rng[1]-rng[0])/2)
        axs[0].annotate('N/O', xy=(0.05, 0.05), xycoords='axes fraction', fontsize=fontsize, bbox=dict(facecolor='none', edgecolor='k', linewidth=2))
        loc1 = 0.2
        loc2 = 0.05
    elif cal == 'mstar':
        axs[0].set_ylabel(r'log($M$, M$_{\odot}$ kpc$^{-2}$)', fontsize=fontsize)
        axs[1].set_ylabel(r'log($M$, M$_{\odot}$ kpc$^{-2}$)', fontsize=fontsize)
        axs[2].set_ylabel(r'$\Delta$log($M$, dex)', fontsize=fontsize)
        if ap[k] == 'nsa':
            rng = [6.5, 9]
        elif ap[k] == 'sim':
            rng = [6.75, 9.75]
        axs[0].set_ylim(rng)
        axs[1].set_ylim(rng)
        axs[2].set_ylim(-1*(rng[1]-rng[0])/2, (rng[1]-rng[0])/2)
        loc1 = 0.5
        loc2 = 0.1
    elif cal == 'sfr':
        axs[0].set_ylabel(r'$\rm{log(SFR, M_{\odot} yr^{-1} kpc^{-2})}$', fontsize=fontsize)
        axs[1].set_ylabel(r'$\rm{log(SFR, M_{\odot} yr^{-1} kpc^{-2})}$', fontsize=fontsize)
        axs[2].set_ylabel(r'$\rm{\Delta log(SFR, dex)}$', fontsize=fontsize)
        if ap[k] == 'nsa':
            rng = [-4.25, -1.25]
        elif ap[k] == 'sim':
            rng = [-3.75, -1.25]
        axs[0].set_ylim(rng)
        axs[1].set_ylim(rng)
        axs[2].set_ylim(-1*(rng[1]-rng[0])/2, (rng[1]-rng[0])/2)
        loc1 = 0.5
        loc2 = 0.1
    elif cal == 'M08':
        axs[0].set_ylabel('12 + log(O/H)', fontsize=fontsize)
        axs[1].set_ylabel('12 + log(O/H)', fontsize=fontsize)
        axs[2].set_ylabel('$\Delta$log(O/H)', fontsize=fontsize)
        #rng = [8.65, 9.15]
        rng = [8.55, 9.25]
        axs[0].set_ylim(rng)
        axs[1].set_ylim(rng)
        axs[2].set_ylim(-1*(rng[1]-rng[0])/2, (rng[1]-rng[0])/2)
        #axs[2].set_ylim(-0.2, 0.2)
        axs[0].annotate('M08', xy=(0.05, 0.05), xycoords='axes fraction', fontsize=fontsize, bbox=dict(facecolor='none', edgecolor='k', linewidth=2))
        loc1 = 0.1
        loc2 = 0.02
    
    axs[0].annotate('(a) Control Sample', xy=(0.05, 0.93), xycoords='axes fraction', fontsize=fontsize)
    axs[1].annotate('(b) Paired Galaxies', xy=(0.05, 0.93), xycoords='axes fraction', fontsize=fontsize)
    axs[2].annotate('(c) Pairs - Controls', xy=(0.05, 0.93), xycoords='axes fraction', fontsize=fontsize)
    
    ps.ticks(axs, xmajor=0.5, ymajor=loc1, xminor=0.1, yminor=loc2)
    ps.style(axs, fontsize=fontsize)
    
    if ap[k] == 'nsa':
        plt.suptitle('Outside Pair Sample', fontsize=fontsize)
    elif ap[k] == 'sim':
        plt.suptitle('Inside Pair Sample', fontsize=fontsize)
    
    plt.subplots_adjust(wspace=0.30)
    plt.savefig(outpath + '/'+cal+'_'+ap[k]+'_dex.pdf', bbox_inches='tight', overwrite=True)
    plt.close('all')