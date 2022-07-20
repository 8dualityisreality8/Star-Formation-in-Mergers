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
import csv

def blue(x):
    b = 3.1682
    m = 0.16
    return b - (x+18) * m
def red(x):
    b = 4.7866
    m = 0.04
    return b - (x+18) * m

def extract(path):
    with open(filepath + path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        table = []
        for row in csv_reader:
            table.append(row)
        table.pop(0)
        table = np.asarray(table).astype('float64')
        return table

T5path = '/Volumes/Ext_drive/photoObj8'

filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
datapath = os.path.abspath(os.path.join(filepath, '..', '..', 'data'))
#droppath = '/Users/joshua/Dropbox/manga/figs'

#figpath = '/Users/joshua/Documents/GitHub/Manga_Mergers/Figures'

drpallpath = '/Volumes/Ext_drive/photoObj8/MPL-8/drpall_uniq.fits'
drpall = fits.getdata(drpallpath)
mngtarg1 = drpall['mngtarg1']
drmass = np.log10(drpall['nsa_elpetro_mass']/0.7**2)
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

rs = (cmd=='RS')
gv = (cmd=='GV')
bc = (cmd=='BC')

#im = ax1.scatter(drmass, mag[:,1]-mag[:,5], c=z, s=5, vmax=1.0, vmin=0.0, cmap='gray_r')
ax1.scatter(drmass[rs], mag[rs,1]-mag[rs,5], c=z[rs], cmap='Reds', vmin=0, vmax=0.75)
ax1.scatter(drmass[gv], mag[gv,1]-mag[gv,5], c=z[gv], cmap='Greens', vmin=0, vmax=0.75)
ax1.scatter(drmass[bc], mag[bc,1]-mag[bc,5], c=z[bc], cmap='Blues', vmin=0, vmax=0.75)



RS = extract('/Wake17RS.csv')
RSm = (RS[-1,1] - RS[0,1]) / (RS[-1,0] - RS[0,0])
RSb = (RS[0,1] - (RSm * RS[0,0]))

BC = extract('/Wake17BC.csv')
BCm = (BC[-1,1] - BC[0,1]) / (BC[-1,0] - BC[0,0])
BCb = (BC[0,1] - (BCm * BC[0,0]))

ax1.plot([8.75, 11.75], np.array([8.75, 11.75]) * RSm + RSb, c='crimson', linewidth=2)
ax1.plot([8.75, 11.75], np.array([8.75, 11.75]) * BCm + BCb, c='deepskyblue', linewidth=2)

#cbar_ax = fig.add_axes([0.13, 0.88, 0.76, 0.05])
#cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')  
#cbar.ax.set_title('Normalized Number Density', fontsize=fontsize, pad=25)
#plt.suptitle('Normalized Number Density', fontsize=fontsize, y=0.999975)
#ps.cbar_style(cbar, fontsize=fontsize)


leg = ax1.legend(edgecolor='k', fancybox=False, fontsize=fontsize/2.)
leg.get_frame().set_linewidth(2.0)
'''
x = np.arange(-24, -16 + 1, 1)
ax1.plot(x, blue(x), c='k', linestyle='--')
ax1.plot(x, red(x), c='k', linestyle='--')
'''
ax1.set_aspect('auto')

ax1.set_xlim(8.75, 11.75)
ax1.set_ylim(0.5,7.5)

ax1.set_xlabel('log(M/M$_{\odot}$)', fontsize=fontsize)
ax1.set_ylabel('NUV-i', fontsize=fontsize)

ps.ticks(ax1, xmajor=1, ymajor=1, xminor=0.2, yminor=0.2)
ps.style(ax1, fontsize=fontsize)

plt.savefig(filepath + '/wake17comp.png', bbox_inches='tight', overwrite=True)

plt.close('all')