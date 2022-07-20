#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Show the stacking procedure from 5.stacking.py. This creates a figure showing
each mass range with all individual ssfr profiles in that range and the median
profile for the objects. 

"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os, sys
from matplotlib.offsetbox import AnchoredText
import plotstyle as ps
import functions as f
""""""""""""
cal = 'ssfr'
""""""""""""
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
#figpath = '/Users/joshua/Dropbox/manga/figs'

#T5path = f.T5path()
T5path = '/Volumes/Ext_drive/photoObj8'

drpallpath = T5path + '/MPL-8/drpall_uniq.fits'
drpall = fits.getdata(drpallpath)
ba = drpall['nsa_elpetro_ba']
cmd = drpall['cmd']
z = drpall['nsa_z']
m = np.log10(drpall['nsa_elpetro_mass']/0.7**2)

drpallmasterpath = T5path + '/6.master/drpall.fits'
drpall_master = fits.getdata(drpallmasterpath)

simpath = T5path + '/6.master/simard.fits'
sim = fits.getdata(simpath)

gradbin = fits.getdata(filepath + '/3.gradbin/gradbin_nsa.fits')
rbins = gradbin['reff'][0]

#nsa_d16 = gradbin['DOP16']
#nsa_m08 = gradbin['M08']
nsa_ewha = gradbin['ewha']
nsa_ssfr = gradbin['ssfr']
nsa_mstar = gradbin['mstar']
nsa_sfr = gradbin['sfr']

gradbin = fits.getdata(filepath + '/3.gradbin/gradbin_sim.fits')
#sim_d16 = gradbin['DOP16']
#sim_m08 = gradbin['M08']
sim_ewha = gradbin['ewha']
sim_ssfr = gradbin['ssfr']
sim_mstar = gradbin['mstar']
sim_sfr = gradbin['sfr']

CONT = []
PAIR = []
MASS = []
SSFR = []

ap = ['nsa', 'sim']
for k in range(len(ap)):
    # define pair and control samples
    sam = fits.getdata(filepath + '/4.samples.fits')
    controls = sam['all_control'].astype('bool')
    pairs = sam[ap[k]+'_pair'].astype('bool')
    
    oh = nsa_ssfr
    
    if ap[k]=='nsa':
        mass = m
    elif ap[k] == 'sim':
        mass = sam['corrected_mass']
        sim_mendel = sam['sim_mendel'].astype('bool')
        
        oh[sim_mendel] = sim_ssfr[sim_mendel]
    
    #mass = sam['all_mass']
    #gradbin = fits.getdata(filepath + '/3.gradbin/gradbin_'+ap[k]+'.fits')
    
    outpath = filepath + '/6.indiv_grad_'+ap[k]
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        
    #oh = gradbin[cal]
    
    # Open Stacked Gradients
    fil = fits.getdata(filepath + '/5.stacking/stacked_trim_'+ap[k]+'.fits')
    filfull = fits.getdata(filepath + '/5.stacking/stacked_'+ap[k]+'.fits')
    
    binsize=0.5
    mbins = np.arange(9.0, 11.5 + binsize, binsize)+binsize/2

    rbins = fil['rbins']
    pbin, pstd = fil[cal+'_pbin'].T, fil[cal+'_pstd'].T
    
    #pbin, pstd = f.stack(oh[pairs], z=mass[pairs], zbin=mbins, error='sem', trim = 0.7)
    
    cbin, cstd = fil[cal+'_cbin'].T, fil[cal+'_cstd'].T
    dexbin, dexstd = fil[cal+'_dex'].T, fil[cal+'_dexstd'].T
    npair, ncontrol = fil[cal+'_npair'].T, fil[cal+'_ncntrl'].T
    
    fpbin, fpstd = filfull[cal+'_pbin'].T, filfull[cal+'_pstd'].T
    fcbin, fcstd = filfull[cal+'_cbin'].T, filfull[cal+'_cstd'].T
    
    CONT.append(controls)
    PAIR.append(pairs)
    MASS.append(mass)
    SSFR.append(oh)
    
    # Individual Gradient Plotting
    fontsize=30
    fig, axs = plt.subplots(5, 2, sharex=True, sharey=True, figsize=(10, 14))
    
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    for i in range(5):
        mix = abs(mbins[i] - mass) < (mbins[1] - mbins[0])/2
        
        for j in range(len(oh[mix&controls])):
            axs[i,0].plot(rbins, oh[mix&controls][j], c='silver')
            
        axs[i,0].plot(rbins, cbin[i], zorder=9998, c=colors[i])
        axs[i,0].plot(rbins, fcbin[i], zorder=9997, c=colors[i], linestyle='--')
        axs[i,0].fill_between(rbins, cbin[i]-cstd[i], cbin[i]+cstd[i], alpha=0.3,
                     antialiased=True, color=colors[i], zorder=9999)
        
        text = AnchoredText(str(mbins[i]-0.25)+'-'+str(mbins[i]+0.25), loc=4, frameon=False, prop=dict(fontsize=fontsize/1.5))
        axs[i,0].add_artist(text)
        
        text = AnchoredText('N = '+str(len(oh[mix&controls])), loc=3, frameon=False, prop=dict(fontsize=fontsize/1.5))
        axs[i,0].add_artist(text)
        
        for j in range(len(oh[mix&pairs])):
            axs[i,1].plot(rbins, oh[mix&pairs][j], c='silver')
        
        axs[i,1].plot(rbins, pbin[i], zorder=9998, c=colors[i])
        axs[i,1].plot(rbins, fpbin[i], zorder=9997, c=colors[i], linestyle='--')
        axs[i,1].fill_between(rbins, pbin[i]-pstd[i], pbin[i]+pstd[i], alpha=0.3,
                     antialiased=True, color=colors[i], zorder=9999)
    
        text = AnchoredText(str(mbins[i]-0.25)+'-'+str(mbins[i]+0.25), loc=4, frameon=False, prop=dict(fontsize=fontsize/1.5))
        axs[i,1].add_artist(text)
        
        text = AnchoredText('N = '+str(len(oh[mix&pairs])), loc=3, frameon=False, prop=dict(fontsize=fontsize/1.5))
        axs[i,1].add_artist(text)
    
    if cal == 'ewha':
        loc1 = 1.0
        loc2 = 0.25
    elif cal == 'ssfr':
        loc1 = 1.0
        loc2 = 0.2
    elif cal == 'sfr':
        loc1 = 2.0
        loc2 = 0.5
    
    ps.ticks(axs, xmajor=0.5, ymajor=loc1, xminor=0.1, yminor=loc2)
    ps.style(axs[:,0], fontsize=fontsize)
    ps.style(axs[:,1], fontsize=fontsize, labelleft=False)
    
    axs[0,0].set_title('Controls', fontsize=fontsize)
    axs[0,1].set_title('Pairs', fontsize=fontsize)
    
    plt.xlim(-0.2,2.7)
    
    if cal == 'ewha':
        if ap[k] == 'nsa':
            plt.ylim(-0.25, 2.75)
        elif ap[k] == 'sim':
            plt.ylim(-1.5, 2.75)
    elif cal == 'ssfr':
        if ap[k] == 'nsa':
            plt.ylim(-14.2, -8.2)
        elif ap[k] == 'sim':
            plt.ylim(-14.2, -8.2)
    elif cal == 'sfr':
        plt.ylim(-5.5, -0.5)
        
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    if cal == 'ewha':
        plt.ylabel(r'$\rm log(EW(H\alpha))$', fontsize=fontsize, labelpad=10)
    elif cal == 'ssfr':
        plt.ylabel(r'$\rm log(sSFR, yr^{-1})$', fontsize=fontsize, labelpad=20)
    elif cal == 'sfr':
        plt.ylabel(r'$\rm log(SFR, M_{\odot} yr^{-1} kpc^{-2})$', fontsize=fontsize, labelpad=20)
    
    plt.subplots_adjust(wspace=0, hspace=0)
    
    plt.xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)
    
    if ap[k] == 'nsa':
        plt.suptitle('Outside Pair Sample', fontsize=fontsize)
    elif ap[k] == 'sim':
        plt.suptitle('Inside Pair Sample', fontsize=fontsize)
        
    plt.savefig(outpath + '/'+cal+'_'+ap[k]+'_indiv.pdf', bbox_inches='tight', overwrite=True)
    plt.close('all')
#______________________________________________________________________________
    # Plot the number of galaxies in each bin
    fig, axs = plt.subplots(5, 2, sharex=True, sharey=False, figsize=(10, 14))
    
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    for i in range(5):
        mix = abs(mbins[i] - mass) < (mbins[1] - mbins[0])/2
        
        trim = ncontrol * (~np.isnan(cbin)) + -99 * (np.isnan(cbin))
        trim[trim==-99] = np.nan
        axs[i,0].plot(rbins, trim[i], zorder=9998, c=colors[i])
        axs[i,0].plot(rbins, ncontrol[i], zorder=9998, c=colors[i], linestyle='--')
        axs[i,0].fill_between(rbins, trim[i]-np.sqrt(trim[i]), trim[i]+np.sqrt(trim[i]), alpha=0.3,
                         antialiased=True, color=colors[i], zorder=9999)
    
    
        text = AnchoredText(str(mbins[i]-0.25)+'-'+str(mbins[i]+0.25), loc=4, frameon=False, prop=dict(fontsize=fontsize/1.5))
        axs[i,0].add_artist(text)
        
        text = AnchoredText('N = '+str(len(oh[mix&controls])), loc=3, frameon=False, prop=dict(fontsize=fontsize/1.5))
        axs[i,0].add_artist(text)
        
        trim = npair * (~np.isnan(pbin)) + -99 * (np.isnan(pbin))
        trim[trim==-99] = np.nan
        axs[i,1].plot(rbins, trim[i], zorder=9998, c=colors[i])
        axs[i,1].plot(rbins, npair[i], zorder=9998, c=colors[i], linestyle='--')
        axs[i,1].fill_between(rbins, trim[i]-np.sqrt(trim[i]), trim[i]+np.sqrt(trim[i]), alpha=0.3,
                         antialiased=True, color=colors[i], zorder=9999)
        
        text = AnchoredText(str(mbins[i]-0.25)+'-'+str(mbins[i]+0.25), loc=4, frameon=False, prop=dict(fontsize=fontsize/1.5))
        axs[i,1].add_artist(text)
        
        text = AnchoredText('N = '+str(len(oh[mix&pairs])), loc=3, frameon=False, prop=dict(fontsize=fontsize/1.5))
        axs[i,1].add_artist(text)
    
    for i in range(5):
        axs[i,0].set_ylim(0,775)
        axs[i,1].set_ylim(0,55)
    
    if cal == 'ewha' or cal == 'ssfr':
        axs[0,0].annotate(r'EW(H$\alpha$) and sSFR', xy=(0.06, 0.82), xycoords='axes fraction', fontsize=fontsize/1.5, bbox=dict(facecolor='none', edgecolor='k', linewidth=2))
    
    ps.style(axs[:,0], fontsize=fontsize)
    ps.style(axs[:,1], fontsize=fontsize, labelleft=False)
    ps.ticks(axs[:,0], xmajor=0.5, ymajor=200, xminor=0.1, yminor=100)
    ps.ticks(axs[:,1], xmajor=0.5, ymajor=20, xminor=0.1, yminor=10)
    
    axs[0,0].set_title('Controls', fontsize=fontsize)
    axs[0,1].set_title('Pairs', fontsize=fontsize)
    
    plt.xlim(-0.2,2.7)
    
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    
    plt.subplots_adjust(wspace=0.15, hspace=0.0)
    
    plt.ylabel('Number of Galaxies', fontsize = fontsize, labelpad = 20)
    plt.xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)
    
    if ap[k] == 'nsa':
        plt.suptitle('Outside Pair Sample', fontsize=fontsize)
    elif ap[k] == 'sim':
        plt.suptitle('Inside Pair Sample', fontsize=fontsize)
    
    plt.savefig(outpath + '/'+cal+'_'+ap[k]+'_N.pdf', bbox_inches='tight', overwrite=True)
    plt.close('all')
  
# combined samples
MASS = np.asarray(MASS)
SSFR = np.asarray(SSFR)
CONT = np.asarray(CONT)
CONT[1]= False
PAIR = np.asarray(PAIR)
    
outpath = filepath + '/6.indiv_grad_comb'
if not os.path.exists(outpath):
    os.makedirs(outpath)

# Open Stacked Gradients
fil = fits.getdata(filepath + '/5.stacking/stacked_trim_comb.fits')
filfull = fits.getdata(filepath + '/5.stacking/stacked_comb.fits')

binsize=0.5
mbins = np.arange(9.0, 11.5 + binsize, binsize)+binsize/2

rbins = fil['rbins']
pbin, pstd = fil[cal+'_pbin'].T, fil[cal+'_pstd'].T

cbin, cstd = fil[cal+'_cbin'].T, fil[cal+'_cstd'].T

fpbin, fpstd = filfull[cal+'_pbin'].T, filfull[cal+'_pstd'].T
fcbin, fcstd = filfull[cal+'_cbin'].T, filfull[cal+'_cstd'].T

dex = pbin[i] - cbin[i]
dexf = fpbin[i] - fcbin[i]
dexstd = np.sqrt( cstd[i]**2 + pstd[i]**2 )

# Individual Gradient Plotting
fontsize=30
fig, axs = plt.subplots(5, 2, sharex=True, sharey=True, figsize=(10, 14))

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
for i in range(5):
    mix = abs(mbins[i] - MASS) < (mbins[1] - mbins[0])/2
    
    for j in range(len(SSFR[mix&CONT])):
        axs[i,0].plot(rbins, SSFR[mix&CONT][j], c='silver')
        
    axs[i,0].plot(rbins, cbin[i], zorder=9998, c=colors[i])
    axs[i,0].plot(rbins, fcbin[i], zorder=9997, c=colors[i], linestyle='--')
    axs[i,0].fill_between(rbins, cbin[i]-cstd[i], cbin[i]+cstd[i], alpha=0.3,
                 antialiased=True, color=colors[i], zorder=9999)
    
    text = AnchoredText(str(mbins[i]-0.25)+'-'+str(mbins[i]+0.25), loc=4, frameon=False, prop=dict(fontsize=fontsize/1.5))
    axs[i,0].add_artist(text)
    
    text = AnchoredText('N = '+str(len(SSFR[mix&CONT])), loc=3, frameon=False, prop=dict(fontsize=fontsize/1.5))
    axs[i,0].add_artist(text)
    
    for j in range(len(SSFR[mix&PAIR])):
        axs[i,1].plot(rbins, SSFR[mix&PAIR][j], c='silver')
    
    axs[i,1].plot(rbins, pbin[i], zorder=9998, c=colors[i])
    axs[i,1].plot(rbins, fpbin[i], zorder=9997, c=colors[i], linestyle='--')
    axs[i,1].fill_between(rbins, pbin[i]-pstd[i], pbin[i]+pstd[i], alpha=0.3,
                 antialiased=True, color=colors[i], zorder=9999)

    text = AnchoredText(str(mbins[i]-0.25)+'-'+str(mbins[i]+0.25), loc=4, frameon=False, prop=dict(fontsize=fontsize/1.5))
    axs[i,1].add_artist(text)
    
    text = AnchoredText('N = '+str(len(SSFR[mix&PAIR])), loc=3, frameon=False, prop=dict(fontsize=fontsize/1.5))
    axs[i,1].add_artist(text)

if cal == 'ewha':
    loc1 = 1.0
    loc2 = 0.25
elif cal == 'ssfr':
    loc1 = 1.0
    loc2 = 0.2
elif cal == 'sfr':
    loc1 = 2.0
    loc2 = 0.5

ps.ticks(axs, xmajor=0.5, ymajor=loc1, xminor=0.1, yminor=loc2)
ps.style(axs[:,0], fontsize=fontsize)
ps.style(axs[:,1], fontsize=fontsize, labelleft=False)

axs[0,0].set_title('Controls', fontsize=fontsize)
axs[0,1].set_title('Pairs', fontsize=fontsize)

plt.xlim(-0.2,2.7)

if cal == 'ewha':
    plt.ylim(-1.5, 2.75)
elif cal == 'ssfr':
    plt.ylim(-14.2, -8.2)
elif cal == 'sfr':
    plt.ylim(-5.5, -0.5)
    
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.grid(False)
if cal == 'ewha':
    plt.ylabel(r'$\rm log(EW(H\alpha))$', fontsize=fontsize, labelpad=10)
elif cal == 'ssfr':
    plt.ylabel(r'$\rm log(sSFR, yr^{-1})$', fontsize=fontsize, labelpad=20)
elif cal == 'sfr':
    plt.ylabel(r'$\rm log(SFR, M_{\odot} yr^{-1} kpc^{-2})$', fontsize=fontsize, labelpad=20)

plt.subplots_adjust(wspace=0, hspace=0)

plt.xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)

plt.suptitle('Combined Pair Sample', fontsize=fontsize)

plt.savefig(outpath + '/'+cal+'_comb_indiv.pdf', bbox_inches='tight', overwrite=True)
plt.close('all')