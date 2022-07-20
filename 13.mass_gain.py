#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Show the expected mass growth in the paired galaxies due to merger-induced
star formation.
"""
from astropy.io import fits
import numpy as np
import os, sys
import functions as f
import matplotlib.pyplot as plt
import plotstyle as ps

from lmfit import Model
def gaussian(x, amp, cen, wid):
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    cen = 0
    return (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2))
gmodel = Model(gaussian)

""""""""""""
cal = 'ssfr'
bining = ['mass']#, 'dm', 'sep']
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

yvals = []
stds = []
masses = []
dms = []
seps = []

ap = ['nsa', 'sim']
for k in range(len(ap)):
    # Make folder for figures
    outpath = filepath + '/10.dexplots_'+ap[k]
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    # Open the individual profiles    
    sam = fits.getdata(filepath + '/4.samples.fits')
    pairs = sam[ap[k]+'_pair'].astype('bool')
    
    
    if ap[k]=='nsa' or ap[k]=='sim':    
        stacked = fits.getdata(filepath + '/8.controls_'+ap[k]+'/control_'+ap[k]+'.fits')
        data = stacked[cal+'_dex']
        std = stacked[cal+'_std']
    
    dm = sam[ap[k]+'_dm']
    sep = sam[ap[k]+'_sepkpc']
    
    if ap[k]=='nsa':
        mass = drmass
    elif ap[k] == 'sim':
        mass = sam['corrected_mass']
    elif ap[k] == 'all':
        mass = sam['all_mass']
    
    gradbin = fits.getdata(filepath + '/3.gradbin/gradbin_nsa.fits')  
    rbins = gradbin['reff'][0]
    
    yvals.append(data[pairs])
    stds.append(std[pairs])
    masses.append(mass[pairs])
    dms.append(dm[pairs])
    seps.append(sep[pairs])
    

outpath = filepath + '/13.mass_gain'
if not os.path.exists(outpath):
    os.makedirs(outpath)   

#______________________________________________________________________________
for j in range(len(bining)):        
    if bining[j] == 'mass':
        label = 'log(M/M$_\odot$) = '
        binsize = 0.5
        mbins = np.arange(9.0, 11.5, binsize) + binsize/2
        z = np.concatenate((masses[0], masses[1]))
    elif bining[j] == 'dm':
        label = '$\Delta$log(M/M$_\odot$) = '
        binsize = 1.0
        mbins = np.arange(-2.0, 2.0, binsize) + binsize/2
        z = np.concatenate((dms[0], dms[1]))
    elif bining[j] == 'sep':
        label = 'Separation [kpc] = '
        binsize = 12.5
        mbins = np.arange(0.0, 50, binsize) + binsize/2
        z = np.concatenate((seps[0], seps[1]))

    data = np.concatenate((yvals[0], yvals[1]))
        
    dexbin, dexstd = f.stack(data, z=z, zbin=mbins, error='sem')

    dexbin[np.isnan(dexbin)] = 0.0

fontsize=30
linewidth = 1.5

mass_gain = []
mg_err = []
for i in range(len(dexbin)):
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=False, figsize=(18, 6))
    """ Gaussian Fitting """

    result = gmodel.fit(dexbin[i], x=rbins, amp=np.nanmax(dexbin[i]), cen=rbins[(dexbin[i]==np.nanmax(dexbin[i]))][0], wid=0.5)
    
    print(result.fit_report())
    
    axs[0].errorbar(rbins, dexbin[i], dexstd[i], fmt='s', c='grey')
    axs[0].plot(rbins, result.init_fit, 'k:', label='initial fit', zorder=100)
    axs[0].plot(rbins, result.best_fit, 'r-', label='best fit', zorder=99, linewidth=2.0)
    
    leg = axs[0].legend(fontsize=fontsize/1.5, fancybox=False, edgecolor='k')
    leg.get_frame().set_linewidth(linewidth)
    
    """ Controls """
    
    filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
    figpath = '/Users/joshua/Documents/GitHub/Manga_Mergers/'
    droppath = '/Users/joshua/Dropbox/manga/figs'
    
    
    fil = fits.getdata(filepath + '/5.stacking/stacked_comb.fits')
    
    mbins = fil['mbins']
    mbins = mbins[mbins!=0]
    binsize = mbins[1] - mbins[0]
    
    rbins = fil['rbins']
    cbin, cstd = fil[cal+'_cbin'].T, fil[cal+'_cstd'].T
    
    axs[1].errorbar(rbins, cbin[i], yerr=cstd[i], label='Control', linewidth=linewidth, c='k')
    
    fiterr = np.sqrt(result.redchi**2 + cstd[i]**2)
    axs[1].errorbar(rbins, (cbin[i] + result.best_fit), yerr=fiterr, label='Control + Model Dex', linewidth=linewidth, c='r')
    
    axs[1].fill_between(rbins, y1=cbin[i], y2=(cbin[i] + result.best_fit), facecolor='r', alpha=0.1, hatch='\\', edgecolor='r', linewidth=0.0)
    
    leg = axs[1].legend(fontsize=fontsize/1.5, fancybox=False, edgecolor='k')
    leg.get_frame().set_linewidth(linewidth)
    
    fiterr1 = np.abs(10**(cbin[i] + result.best_fit)) * np.log(10) * fiterr * 1*10**9
    axs[2].errorbar(rbins, 10**(cbin[i] + result.best_fit) * 1*10**9, yerr=fiterr1, c='r', linewidth=linewidth, label='Control + Model Dex')
    #axs[2].plot(rbins, 10**(cbin[i] + result.best_fit) * 1*10**9, c='r', linewidth=linewidth, label='Control + Model Dex')
    #axs[2].fill_between(rbins, y1 = 10**(cbin[i] + result.best_fit) *10**9 + fiterr1, y2=10**(cbin[i] + result.best_fit) * 10**9 - fiterr1, alpha=0.1, facecolor='r')
    
    fiterr2 = np.abs(10**cbin[i]) * np.log(10) * cstd[i] * 1*10**9
    axs[2].errorbar(rbins, 10**cbin[i] * 1*10**9, yerr=fiterr2, c='k', linewidth=linewidth, label='Control')
    #axs[2].plot(rbins, 10**cbin[i] * 1*10**9, c='k', linewidth=linewidth, label='Control')
    #axs[2].fill_between(rbins, y1 = 10**cbin[i] * 10**9 + fiterr2, y2=10**cbin[i] * 10**9 - fiterr2, alpha=0.1, facecolor='k')
    
    gain = (10**(cbin[i] + result.best_fit) - 10**cbin[i]) * 0.4*10**9
    fiterr = np.sqrt( fiterr1**2 + fiterr2**2 )
    
    mass_gain.append(gain)
    mg_err.append(fiterr)
    
    axs[2].errorbar(rbins, gain, yerr=fiterr, c='deepskyblue', linewidth=linewidth, label='Difference')
    #axs[2].fill_between(rbins, y1 = gain+fiterr, y2=gain-fiterr, alpha=0.1, facecolor='deepskyblue')
    
    leg = axs[2].legend(fontsize=fontsize/1.5, fancybox=False, edgecolor='k')
    
    leg.get_frame().set_linewidth(linewidth)
    
    axs[0].set_xlim(-0.1, 2.6)
    axs[0].set_xlim(-0.1, 2.6)
    axs[0].set_xlim(-0.1, 2.6)
    
    axs[0].set_ylim(-0.35, 0.55)
    axs[1].set_ylim(-11.11, -9.34)
    axs[2].set_ylim(-0.025, 0.275)

    
    ps.ticks(axs[0], xmajor=0.5, ymajor=0.1, xminor=0.1, yminor=0.02)
    ps.ticks(axs[1], xmajor=0.5, ymajor=0.5, xminor=0.1, yminor=0.1)
    ps.ticks(axs[2], xmajor=0.5, ymajor=0.05, xminor=0.1, yminor=0.01)
    ps.style(axs, fontsize=fontsize)
    
    axs[0].set_xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)
    axs[1].set_xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)
    axs[2].set_xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)
    
    axs[0].set_ylabel(r'$\rm \Delta log(sSFR, yr^{-1})$', fontsize=fontsize)
    axs[1].set_ylabel(r'$\rm log(sSFR, yr^{-1})$', fontsize=fontsize)
    axs[2].set_ylabel(r'$\Delta$M$_{\rm SF}$/M$_{\rm Total}$', fontsize=fontsize)
    
    #axs[0].set_title('sSFR Dex', fontsize=fontsize)
    #axs[1].set_title('Pair - Control', fontsize=fontsize)
    #axs[2].set_title('Mass Ratio', fontsize=fontsize)
    
    axs[0].text(x=0.07, y=0.08, s='log(M/M$_{\odot}$) = '+str(mbins[i]-binsize/2)+'-'+str(mbins[i]+binsize/2), 
       transform=axs[0].transAxes, fontsize=fontsize/1.5, bbox=dict(boxstyle="square", fc='w'))
    
    plt.subplots_adjust(wspace=0.35)
    
    plt.savefig(outpath + '/mass_gain_'+str(i)+'.pdf', bbox_inches='tight', overwrite=True)
    plt.close('all')
    
""" Plot Merger Induced Mass Gain """
fig, axs = plt.subplots(1, 1, sharex=True, sharey=False, figsize=(8, 8))
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

mass_gain = np.asarray(mass_gain)
mg_err = np.asarray(mg_err)

for i in range(len(mass_gain)):
    axs.plot(rbins, mass_gain[i], c=colors[i], label='log(M/M$_{\odot}$) = '+str(mbins[i]-binsize/2)+'-'+str(mbins[i]+binsize/2))
    axs.fill_between(rbins, y1 = mass_gain[i] + mg_err[i], y2= mass_gain[i] - mg_err[i], alpha=0.1, facecolor=colors[i])
    
leg = axs.legend(fontsize=fontsize/1.5, fancybox=False, edgecolor='k')
leg.get_frame().set_linewidth(linewidth)

axs.set_xlim(-0.1, 2.6)
axs.set_ylim(-0.005, 0.085)

ps.ticks(axs, xmajor=0.5, ymajor=0.02, xminor=0.1, yminor=0.005)
ps.style(axs, fontsize=fontsize)

axs.set_xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)
axs.set_ylabel(r'$\Delta$M$_{\rm Merger}$/M$_{\rm Total}$', fontsize=fontsize)
#axs.set_ylabel(r'M$_{\rm New}$/M$_{\rm Original}$', fontsize=fontsize)

#axs.set_title('Stellar Mass Gain From Merger Induced SF', fontsize=fontsize)

plt.savefig(outpath + '/mass_gain.pdf', bbox_inches='tight', overwrite=True)
plt.close('all')