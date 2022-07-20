#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 13:59:14 2019

@author: joshua
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import plotstyle as ps

def kewley01(x):
    return ( 0.61 / (x - 0.47) ) + 1.19
def kauffmann03(x):
    return ( 0.61 / (x - 0.05) ) + 1.30

filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
#figpath = '/Users/joshua/Dropbox/manga/figs'

spax = fits.getdata('/Volumes/Ext_drive/Research/metal/outside_pairs/spaxels.fits')

# pair selection
sam = spax['sample']
bpt = spax['bpt']
cmd = spax['cmd']

O3Hb = spax['o3hb']
N2Ha = spax['N2Ha']

# define pair and control samples
fin = (~np.isnan(O3Hb))&(~np.isnan(N2Ha))&(np.isfinite(O3Hb))&(np.isfinite(N2Ha))

sfgs = (((sam==1)|(sam==2))&(cmd==1))
alls = (((sam==1)|(sam==2)))

fig, (ax1) = plt.subplots(1, sharex=False, sharey=False, figsize=(10, 10))
fontsize=30
#x = -1.65, 0.35
#y = -1.25, 1.25

x = -1.65, 0.65
y = -1.75, 1.75
'''
# all spaxels
o3hb = O3Hb[fin&alls]
n2ha = N2Ha[fin&alls] 

hist = np.histogram2d(n2ha, o3hb, bins=[np.arange(x[0], x[1], 0.01), np.arange(y[0], y[1], 0.01)])[0]
hist = hist/hist.max()
hist[hist!=0] = 1
hist[hist==0] = np.nan


im = ax1.imshow(hist.T, origin='lower', cmap='Greys', extent=[x[0], x[1], y[0], y[1]], aspect='auto')
'''
# SFGs
o3hb = O3Hb[fin&sfgs]
n2ha = N2Ha[fin&sfgs] 

hist = np.histogram2d(n2ha, o3hb, bins=[np.arange(x[0], x[1], 0.01), np.arange(y[0], y[1], 0.01)])[0]
hist = hist/hist.max()
hist[hist==0] = np.nan
im = ax1.imshow(hist.T, origin='lower', cmap='Spectral_r', extent=[x[0], x[1], y[0], y[1]], aspect='auto')

cbar_ax = fig.add_axes([0.1625, 0.88, 0.70, 0.05])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')  
ps.cbar_style(cbar, fontsize=fontsize)
#cbar.ax.set_title('Normalized Number Density',  fontsize=fontsize)#, pad=30)
plt.suptitle('Normalized Number Density',  fontsize=fontsize, y=0.9975)
cbar.ax.set_xticklabels([0.2,0.4,0.6,0.8])

xx = np.arange(x[0], 0.45, 0.01)
yy = kewley01(xx)
ax1.plot(xx, yy, c='k', linewidth=2)
xx = np.arange(x[0], 0.0, 0.01)
yy = kauffmann03(xx)
ax1.plot(xx, yy, c='k', linewidth=2, linestyle='--')

ax1.set_xlim(x)
ax1.set_ylim(y)

ax1.set_ylabel(r'log([O III])/H$\beta$', fontsize=fontsize)
ax1.set_xlabel(r'log([N II])/H$\alpha$', fontsize=fontsize)

ps.ticks(ax1, xmajor=0.5, ymajor=0.5, xminor=0.1, yminor=0.1)
ps.style(ax1, fontsize=fontsize)

plt.savefig(filepath + '/bpt_spax.pdf', bbox_inches='tight', overwrite=True)
#plt.savefig(figpath + '/bpt_spax.pdf', bbox_inches='tight', overwrite=True)
plt.close('all')