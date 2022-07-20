# -*- coding: utf-8 -*-
"""
Compare the results of this work with previous works.
"""
from astropy.io import fits
import numpy as np
import os, sys
import functions as f
import matplotlib.pyplot as plt
import plotstyle as ps

import csv

def extract(path):
    with open(outpath + path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        table = []
        for row in csv_reader:
            table.append(row)
        table.pop(0)
        table = np.asarray(table).astype('float64')
        return table

""""""""""""
cal = 'ssfr'
bining = ['mass', 'dm', 'sep']
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

SSFR = []
SFR = []
SSFR_stds = []
SFR_stds = []
masses = []
dms = []
seps = []

nuc_yvals = []
nuc_stds = []
nuc_masses = []
nuc_dms = []
nuc_seps = []

# Make folder for figures
outpath = filepath + '/15.prof_comp'
if not os.path.exists(outpath):
    os.makedirs(outpath)
#______________________________________________________________________________
# Other Papers

# Pan19
pan19 = extract('/RadialProfile/Pan19.csv')
# Thorp19
Thorp19 = extract('/RadialProfile/Thorp19.csv')
#______________________________________________________________________________

# Open the individual profiles    
sam = fits.getdata(filepath + '/4.samples.fits')

ap = ['nsa', 'sim']
for k in range(len(ap)):
    pairs = sam[ap[k]+'_pair'].astype('bool')

    if ap[k]=='nsa' or ap[k]=='sim':    
        stacked = fits.getdata(filepath + '/8.controls_'+ap[k]+'/control_'+ap[k]+'.fits')
        ssfr = stacked['ssfr_dex']
        ssfr_std = stacked['ssfr_std']
        sfr = stacked['sfr_dex']
        sfr_std = stacked['sfr_std']
    
    dm = np.abs(sam['nsa_dm'])
    sep = sam[ap[k]+'_sepkpc']
    
    pairs = pairs#&(dm<=1.0)&~np.isnan(dm)
    #

    if ap[k]=='nsa':
        mass = drmass
    elif ap[k] == 'sim':
        mass = sam['corrected_mass']
    elif ap[k] == 'all':
        mass = sam['all_mass']
    
    gradbin = fits.getdata(filepath + '/3.gradbin/gradbin_nsa.fits')  
    rbins = gradbin['reff'][0]
    
    SSFR.append(ssfr[pairs])
    SSFR_stds.append(ssfr_std[pairs])
    SFR.append(sfr[pairs])
    SFR_stds.append(sfr_std[pairs])
    masses.append(mass[pairs])
    dms.append(dm[pairs])
    seps.append(sep[pairs])
#______________________________________________________________________________
    # nuc files
    nuc = fits.getdata(filepath + '/11.nuc_dex/nuc_dex_'+ap[k]+'.fits')
    
    nucdata = nuc[cal + '_dex'][pairs]
    std = nuc[cal + '_std'][pairs]
    
    if ap[k]=='all':
        savedmass = mass
        savedcal = nucdata
    
    # Group for alls
    nuc_yvals.append(nucdata)
    nuc_stds.append(std)
    nuc_masses.append(mass[pairs])
    nuc_dms.append(dm[pairs])
    nuc_seps.append(sep[pairs])

ssfr = np.concatenate((SSFR[0], SSFR[1]))
sfr = np.concatenate((SFR[0], SFR[1]))

mass = np.concatenate((nuc_masses[0], nuc_masses[1]))
dm = np.concatenate((nuc_dms[0], nuc_dms[1]))
sep = np.concatenate((nuc_seps[0], nuc_seps[1]))
nuc_data = np.concatenate((nuc_yvals[0], nuc_yvals[1]))
std = np.concatenate((nuc_stds[0], nuc_stds[1]))
#______________________________________________________________________________
""" Single plot """
fontsize=30
fig, axs = plt.subplots(2, 1, sharex=False, sharey=False, figsize=(6, 15))

#SSFR
profile, pstd = f.stack(ssfr, error='sem')

axs[0].scatter(rbins, profile, c='k', marker='D')
axs[0].plot(rbins, profile, c='k', label='This Work')
axs[0].fill_between(rbins, profile-pstd, profile+pstd, alpha=0.2,
                             color='k', antialiased=True)
axs[0].plot([-0.1, 2.6], [0,0], ':k')

axs[0].scatter(pan19[:,0], pan19[:,1], c='r')
axs[0].plot(pan19[:,0], pan19[:,1], 'r', label='Pan+19')

leg = axs[0].legend(handlelength=0, handletextpad=0, fontsize=fontsize/2.0, loc=1, fancybox=False, edgecolor='k')
for line, text in zip(leg.get_lines(), leg.get_texts()):
    text.set_color(line.get_color())
for item in leg.legendHandles:
    item.set_visible(False)
leg.get_frame().set_linewidth(2.0)

axs[0].set_ylim(-0.15, 0.65)
axs[0].set_xlim(-0.1, 2.6)

axs[0].set_xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)
axs[0].set_ylabel(r'$\rm \Delta log(sSFR, dex)$', fontsize=fontsize)

# SFR
profile, pstd = f.stack(sfr, error='sem')

axs[1].scatter(rbins, profile, c='k', marker='D')
axs[1].plot(rbins, profile, c='k', label='This Work')
axs[1].fill_between(rbins, profile-pstd, profile+pstd, alpha=0.2,
                             color='k', antialiased=True)
axs[1].plot([-0.1, 2.6], [0,0], ':k')

axs[1].scatter(Thorp19[:,0], Thorp19[:,1], c='deepskyblue', marker='^')
axs[1].plot(Thorp19[:,0], Thorp19[:,1], 'deepskyblue', label='Thorp+19')

axs[1].set_ylim(-0.15, 0.65)
axs[1].set_xlim(-0.1, 2.6)

axs[1].set_xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)
axs[1].set_ylabel(r'$\rm \Delta log(SFR, dex)$', fontsize=fontsize)

leg = axs[1].legend(handlelength=0, handletextpad=0, fontsize=fontsize/2.0, loc=1, fancybox=False, edgecolor='k')
for line, text in zip(leg.get_lines(), leg.get_texts()):
    text.set_color(line.get_color())
for item in leg.legendHandles:
    item.set_visible(False)
leg.get_frame().set_linewidth(2.0)

ps.style(axs, fontsize=fontsize)
ps.ticks(axs, xmajor=0.5, ymajor=0.2, xminor=0.1, yminor=0.05)

plt.savefig(outpath + '/prof_comp.pdf', overwrite=True, bbox_inches='tight')
plt.close('all')