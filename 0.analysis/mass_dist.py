#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 14:49:33 2020

@author: joshua
"""

from astropy.io import fits
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import functions as f
import plotstyle as ps

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

sam = fits.getdata('/Volumes/Ext_drive/Research/metal/4.samples.fits')

outside_control = sam['nsa_control'].astype('bool')
outside_pair = sam['nsa_pair'].astype('bool')

inside_control = sam['sim_control'].astype('bool')
inside_pair = sam['sim_pair'].astype('bool')

cormass = sam['corrected_mass']
#______________________________________________________________________________
fontsize=30
fig, axs = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(16, 8))

#axs[0].hist(drmass[outside_control], bins=np.arange(9,11.5,0.5), alpha=0.25, 
#   color='grey', weights = np.ones_like(drmass[outside_control])/float(len(drmass[outside_control])))
axs[0].hist(drmass[outside_pair], bins=np.arange(9, 11.5, 0.5),# alpha=0.25, 
   color='deepskyblue', edgecolor='k', linewidth=2.0)#, weights = np.ones_like(drmass[outside_pair])/float(len(drmass[outside_pair])))
#axs[1].hist(cormass[inside_control], bins=np.arange(9, 11.5, 0.5), alpha=0.25,
#   color='grey', weights = np.ones_like(cormass[inside_control])/float(len(cormass[inside_control])))
axs[1].hist(cormass[inside_pair], bins=np.arange(9, 12.0, 0.5),# alpha=0.25, 
   color='deepskyblue', edgecolor='k', linewidth=2.0)#, weights = np.ones_like(cormass[inside_pair])/float(len(cormass[inside_pair])))

axs[0].set_xlabel('log(M/M$_\odot$)', fontsize=fontsize)
axs[1].set_xlabel('log(M/M$_\odot$)', fontsize=fontsize)

axs[0].set_ylabel('Number of Paired Galaxies', fontsize=fontsize)
axs[1].set_ylabel('Number of Paired Galaxies', fontsize=fontsize)

axs[0].set_title('Outside-IFU Sample', fontsize=fontsize)
axs[1].set_title('Inside-IFU Sample', fontsize=fontsize)

ps.ticks(axs[0], xmajor=0.5, ymajor=10, xminor=0.1, yminor=2)
ps.ticks(axs[1], xmajor=0.5, ymajor=5, xminor=0.1, yminor=1)
ps.style(axs, fontsize=fontsize)

plt.savefig(filepath + '/mass_dist2.png', bbox_inches='tight', overwrite=True)
plt.savefig(filepath + '/mass_dist2.pdf', bbox_inches='tight', overwrite=True)
plt.close('all')