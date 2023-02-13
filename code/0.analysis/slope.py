#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 10:33:50 2019

@author: joshua
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import functions as f
import plotstyle as ps

fontsize=30

filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
figpath = '/Users/joshua/Dropbox/manga/figs'

#drpallpath = '/Volumes/T5/photoObj8/MPL-8/drpall_uniq.fits'
drpallpath = '/Volumes/Ext_drive/Research/metal/outside_pairs/MPL-8/drpall_uniq.fits'

drpall = fits.getdata(drpallpath)
ba = drpall['nsa_elpetro_ba']

gradbin = fits.getdata('/Volumes/Ext_drive/Research/metal/outside_pairs/gradbin.fits')
dm = np.abs(gradbin['dm'])
dv = gradbin['dv']

kpc_da = []
for i in drpall['nsa_z']:
    kpc_da.append(f.cosmo(i))
kpc_da = np.asarray(kpc_da)

sep = gradbin['sep_arcsec'] * kpc_da
sep[sep==0] = -8

#samplepath = '/Volumes/T5/photoObj8/8.samples/ncomp.fits'
samplepath = '/Volumes/Ext_drive/Research/metal/outside_pairs/8.samples/ncomp.fits'

sam = fits.getdata(samplepath)

cal = 'M08'

mass = np.log10(drpall['nsa_elpetro_mass']/0.7**2)
mass = sep

binsize = 10
mbins = np.arange(0, 50+binsize, binsize)

oh = gradbin[cal]
doh = gradbin[cal+'_std']

m = gradbin[cal+'_m'][:,0]
merr = gradbin[cal+'_m'][:,1]

cmd = drpall['cmd']

# define pair and control samples
alls = (cmd=='BC')
pairs = (((sam['nout_manga']>0)|(sam['nout_nsa']>0))&(sam['donotuse']!=1)&(sam['ptype']==1)&(sam['ngood']==1)&alls)
controls = ((sam['ctrl']==1)&alls)

fig, axs = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(16, 8))
#plt.rc('axes', linewidth=2)

axs[0].errorbar(mass[controls], m[controls], fmt='o', c='grey', markersize=5, capsize=0, capthick=1, elinewidth=1, zorder=0)
axs[0].errorbar(mass[pairs], m[pairs], fmt='o', c='grey', markersize=5, capsize=0, capthick=1, elinewidth=1, zorder=1)

axs[0].set_ylim(-0.5, 0.5)
axs[0].set_xlim(-15, 55)

axs[0].set_ylabel(r'Slope [dex/R$_{\rm eff}]$', fontsize=fontsize)

axs[0].set_xlabel('Separation [kpc]', fontsize=fontsize)
axs[1].set_xlabel(r'Slope [dex/R$_{\rm eff}$]', fontsize=fontsize)

step = mbins[1]-mbins[0]
zhig = mbins + step/2.
zlow = mbins - step/2.

axs[0].plot([-15, 55], [0, 0], c='k', linestyle=':')
axs[0].plot([0, 0], [-0.5, 0.5], c='k', linewidth=2)
axs[0].plot([-15, 55], [np.nanmedian(m[controls]), np.nanmedian(m[controls])], c='k', linestyle='--')

# make the bins for the plots
pval = []
pstd = []
for i in range(len(mbins)):       
    mix = (mass >= zlow[i])&(mass < zhig[i])
    # pairs
    med = np.nanmedian(m[pairs&mix])
    std = np.nanstd(m[pairs&mix]) #/ np.sqrt(len(m[controls&mix&out]))
    pval.append(med)
    pstd.append(std)

axs[0].errorbar([-8], [np.nanmedian(m[controls])], yerr=[np.nanstd(m[controls])], fmt='s', markersize=10, c='k', capsize=4, capthick=2, markeredgecolor='k', label='Controls')
axs[0].errorbar(mbins, pval, yerr=np.asarray(pstd), fmt='s', markersize=10, c='r', capsize=4, capthick=2, markeredgecolor='k', label='Pairs')

weights = np.ones_like(m[controls])/len(m[controls])
#axs[1].hist(m[controls], bins=np.arange(-0.5, 0.5, 0.025), weights=weights, alpha=0.5, color='grey', label='Controls')
hist = np.histogram(m[controls], bins=np.arange(-0.5, 0.5, 0.025), weights=weights)

axs[1].plot(hist[1][0:-1] + 0.5*(hist[1][1]-hist[1][0]), hist[0], c='k', linestyle='--', label='Controls')

weights = np.ones_like(m[pairs])/len(m[pairs])
#axs[1].hist(m[pairs], bins=np.arange(-0.5, 0.5, 0.025), weights=weights, alpha=0.5, fill=False, label='Pairs', hatch='//', edgecolor='deepskyblue', linewidth=2)
axs[1].hist(m[pairs], bins=np.arange(-0.5, 0.5, 0.025), alpha=0.3,label='Pairs', linewidth=2, color='r', weights=weights)
val, bins = np.histogram(m[pairs], bins=np.arange(-0.5, 0.5, 0.025), weights=weights)
val[val==0] = np.nan

binc = 0.5*(bins[1:]+bins[:-1])

#axs[1].errorbar(binc, val, yerr=val**0.5, fmt='s', c='r', capsize=4, capthick=2, markeredgecolor='k')

axs[1].axvline(np.nanmedian(m[controls]), linewidth=2, c='k')
axs[1].axvline(np.nanmedian(m[pairs]), linewidth=2, c='crimson')
axs[1].set_ylabel('Normalized Number', fontsize=fontsize)
axs[1].yaxis.set_label_position('right')
axs[1].yaxis.tick_right()
ps.ticks(axs[0], xmajor=10, ymajor=0.2, xminor=2, yminor=0.05)
#ps.ticks(axs[1], xmajor=0.2, ymajor=0.02, xminor=0.05, yminor=0.01)
ps.style(axs, fontsize=fontsize)

xlabels = ['Control', '0', '10', '20', '30', '40', '50']
xticks = [-8, 0, 10, 20, 30, 40, 50]
axs[0].set_xticks(xticks)
axs[0].set_xticklabels(xlabels)

legend = axs[0].legend(fancybox=False, edgecolor='k', fontsize=fontsize/1.5)
frame = legend.get_frame().set_linewidth(2.0)

legend = axs[1].legend(fancybox=False, edgecolor='k', fontsize=fontsize/1.5)
frame = legend.get_frame().set_linewidth(2.0)

plt.subplots_adjust(wspace=0, hspace=0)

plt.savefig(filepath + '/slope2.png', bbox_inches='tight', overwrite=True)
#plt.savefig(figpath + '/slope.pdf', bbox_inches='tight', overwrite=True)
plt.close('all')