#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 15:23:19 2019

Make separation distribution for outside IFU pairs.

@author: joshua
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os, sys
import functions as f
import plotstyle as ps
#______________________________________________________________________________
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
figpath = '/Users/joshua/Dropbox/manga/figs'

#drpallpath = '/Volumes/T5/photoObj8/MPL-8/drpall_uniq.fits'

drpallpath = '/Volumes/Ext_drive/Research/metal/outside_pairs/MPL-8/drpall_uniq.fits'


drpall = fits.getdata(drpallpath)
cmd = drpall['CMD']
z = drpall['nsa_z']

gradbin = fits.getdata('/Volumes/Ext_drive/Research/metal/outside_pairs/gradbin.fits')
sep = gradbin['sep_arcsec']

kpc_da = []
for i in z:
    kpc_da.append(f.cosmo(i))
kpc_da = np.asarray(kpc_da)

rad_kpc = kpc_da * sep

droppath = '/Users/joshua/Dropbox/manga/figs'
#samplepath = '/Volumes/T5/photoObj8/8.samples/ncomp.fits'

samplepath = '/Volumes/Ext_drive/Research/metal/outside_pairs/8.samples/ncomp.fits'

sam = fits.getdata(samplepath)

# define pair and control samples
alls = (cmd=='BC') 
pairs = (((sam['nout_manga']>0)|(sam['nout_nsa']>0))&(sam['donotuse']!=1)&(sam['ptype']==1)&(sam['ngood']==1))
controls = ((sam['ctrl']==1))

binsize = 5
bins = np.arange(0, 50 + binsize, binsize)

fontsize=30
fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, figsize=(8, 16))

axs[0].hist(rad_kpc[pairs], bins = bins, label='All Pairs', edgecolor='k', color='deepskyblue', linewidth=2)
#ax1.hist(rad_kpc[pairs&alls], bins = bins, label='SF Pairs', edgecolor='k', color='r', linewidth=2)
#legend = axs[0].legend(fancybox=False, edgecolor='k', fontsize=fontsize/1.5)
#frame = legend.get_frame().set_linewidth(2.0)
axs[0].annotate('(a) Full Pair Sample', xy=(0.03, 0.90), xycoords='axes fraction', fontsize=fontsize)
axs[1].annotate('(b) SF Pair Sample Fraction', xy=(0.03, 0.90), xycoords='axes fraction', fontsize=fontsize)


ps.ticks(axs[0], xmajor=10, ymajor=20, xminor=5, yminor=5)
ps.ticks(axs[1], xmajor=10, ymajor=0.2, xminor=5, yminor=0.05)
ps.style(axs, fontsize=fontsize)

axs[0].set_xlim(-1, 56)
axs[0].set_ylim(0, 90)

#axs[0].set_xlabel('Separation [kpc]', fontsize=fontsize)
axs[0].set_ylabel('Number of Galaxies', fontsize=fontsize)

apar = np.histogram(rad_kpc[pairs], bins = bins)
sfgs = np.histogram(rad_kpc[pairs&alls], bins=bins)

frac = np.divide(np.float64(sfgs[0]), np.float64(apar[0]))

axs[1].bar(bins[0:-1] + (bins[1]-bins[0])/2.0, frac, width=bins[1]-bins[0], edgecolor='k', color='r', linewidth=2.0, label='Fraction of SFGs')
#legend = axs[1].legend(fontsize=fontsize/1.5, fancybox=False, edgecolor='k', loc=2)
#frame = legend.get_frame().set_linewidth(2.0)
axs[1].set_ylim(0, 1)
axs[1].set_xlabel('Separation [kpc]', fontsize=fontsize)
axs[1].set_ylabel('Fraction', fontsize=fontsize)

plt.subplots_adjust(hspace=0.0)

plt.savefig(filepath + '/sep_dist.png', bbox_inches='tight', overwrite=True)
plt.savefig(figpath + '/sep_dist.pdf', bbox_inches='tight', overwrite=True)
plt.close('all')