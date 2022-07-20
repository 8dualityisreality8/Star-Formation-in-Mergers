#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 11:55:46 2020

@author: joshua
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import plotstyle as ps
import functions as f

T5path = f.T5path()

drpallpath = '/Volumes/Ext_drive/Research/data/drpall_gal.fits'
drpall = fits.getdata(drpallpath)

simpath = T5path + '/MPL-8/ancillary/simard11t3.fit'
sim = fits.getdata(simpath)

photosimpath = '/Volumes/Ext_drive/Research/data/simard11/CasJobs/phioffset.fit'
photosim = fits.getdata(photosimpath)
offset = photosim['phioffset_g']

fontsize=30
fig, axs = plt.subplots(1, 3, sharex=False, sharey=False, figsize=(28, 8))
#______________________________________________________________________________
nsa = drpall['NSA_ELPETRO_TH50_R']
simard = sim['RCHL_R']/sim['scale']

axs[0].scatter(nsa, simard, c='grey', s=5)
axs[0].set_xlim(0,25)
axs[0].set_ylim(0, 25)

axs[0].plot([0,25],[0,25], linestyle='--', c='r')

axs[0].set_xlabel('nsa', fontsize=fontsize)
axs[0].set_ylabel('simard', fontsize=fontsize)
axs[0].set_title('Effective Radius [arcsec]', fontsize=fontsize)
#______________________________________________________________________________

nsa = drpall['NSA_ELPETRO_BA']
simard = 1 - sim['e']

axs[1].scatter(nsa, simard, c='grey', s=5)
axs[1].set_xlim(0,1)
axs[1].set_ylim(0,1)

axs[1].plot([0,1], [0,1], c='r', linestyle='--')

axs[1].set_xlabel('nsa', fontsize=fontsize)
axs[1].set_ylabel('simard', fontsize=fontsize)
axs[1].set_title('B/A', fontsize=fontsize)
#______________________________________________________________________________

nsa = drpall['NSA_ELPETRO_PHI']

simard = []
for i in range(len(sim)):
    off = offset[photosim['dr7objid']==sim['objid'][i]]
    if len(off) > 0:
        if len(off) > 1:
            off = off[0]
        simard.append((360 + 90 - (sim['phi'][i] + off)) % 180)
    else:
        simard.append(-9999)
        
simard = np.asarray(simard)

axs[2].scatter(nsa, simard, c='grey', s=5)

axs[2].set_xlim(0,180)
axs[2].set_ylim(0,180)

axs[2].plot([0,180],[0,180], linestyle='--', c='r')

axs[2].set_xlabel('nsa', fontsize=fontsize)
axs[2].set_ylabel('simard', fontsize=fontsize)

axs[2].set_title('Position Angle', fontsize=fontsize)

ps.ticks(axs[0], xmajor=5, ymajor=5, xminor=1, yminor=1)
ps.ticks(axs[1], xmajor=0.2, ymajor=0.2, xminor=0.1, yminor=0.1)
ps.ticks(axs[2], xmajor=40, ymajor=40, xminor=10, yminor=10)
ps.style(axs, fontsize=fontsize)

plt.savefig('mag_comp.png', overwrite=True, bbox_inches='tight')
plt.close('all')

nsamag = drpall['nsa_elpetro_flux'][:,4]
simmag = sim['RGMAG']

plt.scatter(nsamag, simmag)

