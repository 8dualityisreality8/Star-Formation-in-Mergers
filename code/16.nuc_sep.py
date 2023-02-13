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
import astropy.stats as astat

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
cal = 'sfr'
bining = 'sep'
""""""""""""
# set filepath
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))

nuc_sfr = []
nuc__sfr_std = []
nuc_ssfr = []
nuc_ssfr_std = []
nuc_seps = []

# Make folder for figures
outpath = filepath + '/16.nuc_sep'
if not os.path.exists(outpath):
    os.makedirs(outpath)
    
#______________________________________________________________________________
# Other Papers

# Pan19
Ellison08 = extract('/Separation/Ellison08.csv')[0:5]
# Thorp19
Patton13 = extract('/Separation/Patton13.csv')[0:10]
# Thorp19
Scudder12 = extract('/Separation/Scudder12.csv')[0:6]
# Bustamante20
Busta20 = extract('/Separation/Bustamante20.csv')[0:7]
# Li08
Li08 = extract('/Separation/Li2008_fig5_upperleft.csv')[0:8]
Li08[:,1] = Li08[:,1]*1000.
# Patton20
Pat20 = extract('/Separation/Patton20.csv')
Pat20[:,1] = np.log10(Pat20[:,1])
#______________________________________________________________________________
# Open the individual profiles    
sam = fits.getdata(filepath + '/4.samples.fits')

ap = ['nsa', 'sim']
for k in range(len(ap)):
    pairs = sam[ap[k]+'_pair'].astype('bool')
    
    if ap[k]=='nsa' or ap[k]=='sim':    
        stacked = fits.getdata(filepath + '/8.controls_'+ap[k]+'/control_'+ap[k]+'.fits')
        data = stacked[cal+'_dex']
        std = stacked[cal+'_std']
    
    sep = sam[ap[k]+'_sepkpc']
    
    gradbin = fits.getdata(filepath + '/3.gradbin/gradbin_nsa.fits')  
    rbins = gradbin['reff'][0]
#______________________________________________________________________________
    # nuc files
    nuc = fits.getdata(filepath + '/11.nuc_dex/nuc_dex_'+ap[k]+'.fits')
    
    nucdata = nuc['sfr_dex'][pairs]
    std = nuc['sfr_std'][pairs]

    ssfrdata = nuc['ssfr_dex'][pairs]
    ssfrstd = nuc['ssfr_dex'][pairs]
    
    # Group for alls
    nuc_sfr.append(nucdata)
    nuc__sfr_std.append(std)
    nuc_ssfr.append(ssfrdata)
    nuc_ssfr_std.append(ssfrstd)
    nuc_seps.append(sep[pairs])

sep = np.concatenate((nuc_seps[0], nuc_seps[1]))
sfr = np.concatenate((nuc_sfr[0], nuc_sfr[1]))
ssfr = np.concatenate((nuc_ssfr[0], nuc_ssfr[1]))
ssfr_std = np.concatenate((nuc_ssfr_std[0], nuc_ssfr_std[1]))
sfr_std = np.concatenate((nuc__sfr_std[0], nuc__sfr_std[1]))
#______________________________________________________________________________
med = np.nanmedian(sfr[~np.isnan(sfr)])
#______________________________________________________________________________
"""""""""
Combined Samples
"""""""""    
fontsize=30
fig, axs = plt.subplots(2, 1, sharex=False, sharey=False, figsize=(6, 12))
#______________________________________________________________________________
""" Nuc Data"""
binsize = 10.0
mbins = np.arange(0.0, 50, binsize) + binsize/2
label = 'Separation [kpc] = '
z = sep

X = []
Y = []
Ystd = []
Y2 = []
Y2std = []
for i in range(len(mbins)):
    ind = (z > mbins[i] - binsize/2.0)&(z < mbins[i] + binsize/2.0)
    y = sfr[ind]
    yerr = sfr_std[ind]
    fin = np.isfinite(y)
    y = y[fin]
    yerr = yerr[fin]
    ystd = np.nanstd(y) / len(y)**0.5
    
    y2 = ssfr[ind]
    y2err = ssfr_std[ind]
    fin = np.isfinite(y2)
    y2 = y2[fin]
    y2err = y2err[fin]
    y2std = np.nanstd(y2) / len(y2)**0.5
    
    
    if len(y) > 0:
        y = np.nanmean(y[np.isfinite(y)])
        y2 = np.nanmean(y2[np.isfinite(y2)])
        
        #ystd = astat.biweight_scale(y[np.isfinite(y)], c=6) / len(y)**0.5
        #y = astat.biweight_location(y[np.isfinite(y)], c=6)
        
        #y2std = astat.biweight_scale(y2[np.isfinite(y2)], c=6) / len(y2)**0.5
        #y2 = astat.biweight_location(y2[np.isfinite(y2)], c=6)
        
        #print y
        x = mbins[i]

        X.append(x)
        Y.append(y)
        Y2.append(y2)
        Ystd.append(ystd)
        Y2std.append(y2std)

# SFR
axs[0].errorbar(X, Y, xerr=binsize/2.0, yerr=Ystd, c='k', marker='D', 
             zorder=9999, capsize=5.0, capthick=2.0, elinewidth=2.0)
axs[0].plot(X, Y, c='k', label='This Work')

axs[0].scatter(Ellison08[:,0]/.7, np.log10(Ellison08[:,1]), c='r', marker='o')
axs[0].plot(Ellison08[:,0]/.7, np.log10(Ellison08[:,1]), c='r', label='Ellison+08')

axs[0].scatter(Scudder12[:,0]/.7, Scudder12[:,1], c='deepskyblue', marker='^')
axs[0].plot(Scudder12[:,0]/.7, Scudder12[:,1], c='deepskyblue', label='Scudder+12')

axs[0].scatter(Patton13[:,0], np.log10(Patton13[:,1]), c='g', marker='s')
axs[0].plot(Patton13[:,0], np.log10(Patton13[:,1]), c='g', label='Patton+13')

axs[0].scatter(Busta20[:,0], Busta20[:,1], c='purple', marker='v')
axs[0].plot(Busta20[:,0], Busta20[:,1], c='purple', label='Bustamante+20')

axs[0].plot([-5,105], [0,0], ':k')

# sSFR
axs[1].errorbar(X, Y2, xerr=binsize/2.0, yerr=Y2std, c='k', marker='D', 
             zorder=9999, capsize=5.0, capthick=2.0, elinewidth=2.0)
axs[1].plot(X, Y2, c='k', label='This Work')

axs[1].scatter(Li08[:,1], Li08[:,2], c='r', marker='o')
axs[1].plot(Li08[:,1], Li08[:,2], c='r', label='Li+08')

axs[1].scatter(Pat20[:,0], Pat20[:,1], c='deepskyblue', marker='^')
axs[1].plot(Pat20[:,0], Pat20[:,1], c='deepskyblue', label = 'Patton+20')

axs[1].plot([-5,105], [0,0], ':k')

leg = axs[0].legend(handlelength=0, handletextpad=0, fontsize=fontsize/2.0, loc=1, fancybox=False, edgecolor='k')
for line, text in zip(leg.get_lines(), leg.get_texts()):
    text.set_color(line.get_color())
for item in leg.legendHandles:
    item.set_visible(False)
leg.get_frame().set_linewidth(2.0)

leg = axs[1].legend(handlelength=0, handletextpad=0, fontsize=fontsize/2.0, loc=1, fancybox=False, edgecolor='k')
for line, text in zip(leg.get_lines(), leg.get_texts()):
    text.set_color(line.get_color())
for item in leg.legendHandles:
    item.set_visible(False)
leg.get_frame().set_linewidth(2.0)
        
axs[0].set_xlim(-5, 105)
axs[0].set_ylim(-0.05, 0.55)

#axs[0].set_xlabel('Separation [kpc]', fontsize=fontsize)
axs[0].set_ylabel(r'$\rm \Delta log(SFR, dex)$', fontsize=fontsize)

axs[1].set_xlim(-5, 105)
axs[1].set_ylim(-0.05, 0.45)

axs[1].set_xlabel('Separation [kpc]', fontsize=fontsize)
axs[1].set_ylabel(r'$\rm \Delta log(sSFR, dex)$', fontsize=fontsize)

ps.style(axs[0], fontsize=fontsize, labelbottom=False)
ps.style(axs[1], fontsize=fontsize)
ps.ticks(axs, xmajor=25, ymajor=0.1, xminor=5, yminor=0.02)

plt.subplots_adjust(hspace=0.0)

plt.savefig(outpath + '/nuc_sep.pdf', bbox_inches='tight', overwrite=True)
plt.close('all')