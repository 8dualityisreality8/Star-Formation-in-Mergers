#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 14:22:41 2019

@author: joshua
"""
import os, sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

import plotstyle as ps

def blue(x):
    b = 3.1682
    m = 0.16
    return b - (x+18) * m
def red(x):
    b = 4.7866
    m = 0.04
    return b - (x+18) * m

T5path = '/Volumes/Ext_drive/photoObj8'

filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
datapath = os.path.abspath(os.path.join(filepath, '..', '..', 'data'))
#droppath = '/Users/joshua/Dropbox/manga/figs'

#figpath = '/Users/joshua/Documents/GitHub/Manga_Mergers/Figures'

drpallpath = '/Volumes/Ext_drive/photoObj8/MPL-8/drpall_uniq.fits'
drpall = fits.getdata(drpallpath)
mngtarg1 = drpall['mngtarg1']
#drmass = np.log10(drpall['nsa_elpetro_mass']/0.7**2)
cmd = drpall['cmd']

samplepath = '/Volumes/Ext_drive/photoObj8/8.samples/ncomp.fits'
ncomp = fits.getdata(samplepath)

samplepath = '/Volumes/Ext_drive/Research/metal/4.samples.fits'
sam = fits.getdata(samplepath)

bptpath = '/Volumes/Ext_drive/Research/data/bpt_r1.3kpc_pars.fits'
bpt = fits.getdata(bptpath)['bpt']

mag = drpall['nsa_elpetro_absmag']

mass = sam['all_mass']

#______________________________________________________________________________
mup = 11.5
gal = (((mngtarg1 & 2L**10 + 2L**11) != 0) & ((mngtarg1 != 0)))
alls = (cmd=='BC')#&(bpt!=0)
alls2 = gal&(mass >= 9.0)&(mass <= mup)#&(mass>11.0)

# Select Inside Pairs
ipairs = ((ncomp['donotuse']!=1)&(ncomp['ptype']==1)&(ncomp['ngood']>=2)&alls2)

#ipairs = sam['sim_pair'].astype('bool')

# Select Outside Pairs
opairs = (((ncomp['nout_manga']>0)|(ncomp['nout_nsa']>0))&(ncomp['donotuse']!=1)&(ncomp['ptype']==1)&(ncomp['ngood']==1)&alls2)

#opairs = sam['nsa_pair'].astype('bool')

pairs = ipairs|opairs

#______________________________________________________________________________

fontsize=25
fig, (ax1) = plt.subplots(1, sharex=False, sharey=False, figsize=(8, 8))

z = ps.gaussian_den(mag[:,5], mag[:,1]-mag[:,4])

im = ax1.scatter(mag[:,5], mag[:,1]-mag[:,4], c=z, s=5, vmax=1.0, vmin=0.0, cmap='gray_r')

cbar_ax = fig.add_axes([0.13, 0.88, 0.76, 0.05])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')  

#ax1.scatter(mag[:,5], mag[:,1]-mag[:,4], c='k', s=5, vmax=0.8, vmin=0.2)

#ax1.scatter(mag[:,5][bpt==0], mag[:,1][bpt==0]-mag[:,4][bpt==0], c='r', s=5, label='Retired')
#ax1.scatter(mag[:,5][(bpt==1)|(bpt==2)], mag[:,1][(bpt==1)|(bpt==2)]-mag[:,4][(bpt==1)|(bpt==2)], c='deepskyblue', s=5, label='SF')
#ax1.scatter(mag[:,5][bpt>2], mag[:,1][bpt>2]-mag[:,4][bpt>2], c='g', s=5, label='AGN')

cbar.ax.set_title('Normalized Number Density', fontsize=fontsize, pad=25)
#plt.suptitle('Normalized Number Density', fontsize=fontsize, y=0.999975)
ps.cbar_style(cbar, fontsize=fontsize)

#ax1.scatter(mag[:,5][ipairs], (mag[:,1]-mag[:,4])[ipairs], c='r', s=30, marker='D', edgecolor='k', linewidth=1, label='Inside-IFU')
#ax1.scatter(mag[:,5][opairs], (mag[:,1]-mag[:,4])[opairs], c='deepskyblue', s=30, marker='s', edgecolor='k', linewidth=1, label='Outside-IFU')

ax1.scatter(mag[:,5][ipairs&((bpt==1)|(bpt==2))], (mag[:,1]-mag[:,4])[ipairs&((bpt==1)|(bpt==2))], c='deepskyblue', s=30, marker='D', edgecolor='k', linewidth=1, label='SF')
ax1.scatter(mag[:,5][ipairs&(bpt==0)], (mag[:,1]-mag[:,4])[ipairs&(bpt==0)], c='r', s=30, marker='D', edgecolor='k', linewidth=1, label='Retired')
ax1.scatter(mag[:,5][ipairs&(bpt>2)], (mag[:,1]-mag[:,4])[ipairs&(bpt>2)], c='lime', s=30, marker='D', edgecolor='k', linewidth=1, label='AGN')

ax1.scatter(mag[:,5][opairs&((bpt==1)|(bpt==2))], (mag[:,1]-mag[:,4])[opairs&((bpt==1)|(bpt==2))], c='deepskyblue', s=30, marker='D', edgecolor='k', linewidth=1)
ax1.scatter(mag[:,5][opairs&(bpt==0)], (mag[:,1]-mag[:,4])[opairs&(bpt==0)], c='r', s=30, marker='D', edgecolor='k', linewidth=1)
ax1.scatter(mag[:,5][opairs&(bpt>2)], (mag[:,1]-mag[:,4])[opairs&(bpt>2)], c='lime', s=30, marker='D', edgecolor='k', linewidth=1)


leg = ax1.legend(edgecolor='k', fancybox=False, fontsize=fontsize/2.)
leg.get_frame().set_linewidth(2.0)

x = np.arange(-24, -16 + 1, 1)
ax1.plot(x, blue(x), c='k', linestyle='--')
ax1.plot(x, red(x), c='k', linestyle='--')

ax1.set_aspect('auto')

ax1.set_xlim(-16, -24)
ax1.set_ylim(1.05,6.5)

ax1.set_xlabel(r'$\mathcal{M}_{\rm i}$', fontsize=fontsize)
ax1.set_ylabel('NUV-r', fontsize=fontsize)

ps.ticks(ax1, xmajor=1, ymajor=1, xminor=0.2, yminor=0.2)
ps.style(ax1, fontsize=fontsize)

plt.savefig(filepath + '/color-mag.pdf', bbox_inches='tight', overwrite=True)

plt.close('all')