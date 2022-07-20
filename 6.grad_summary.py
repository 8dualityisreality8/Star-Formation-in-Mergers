#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

plots the median "stacked" ssfr profiles for pairs and controls and the 
difference between the two samples.

"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import plotstyle as ps
""""""""""""
cal = 'ssfr'
""""""""""""
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))

#T5path = f.T5path()
T5path = '/Volumes/Ext_drive/photoObj8'

drpallpath = T5path + '/MPL-8/drpall_uniq.fits'
drpall = fits.getdata(drpallpath)
ba = drpall['nsa_elpetro_ba']
cmd = drpall['cmd']
z = drpall['nsa_z']
m = np.log10(drpall['nsa_elpetro_mass']/0.7**2)

samplepath = T5path + '/8.samples/ncomp.fits'
sam = fits.getdata(samplepath)

drpallmasterpath = T5path + '/6.master/drpall.fits'
drpall_master = fits.getdata(drpallmasterpath)

simpath = T5path + '/6.master/simard.fits'
sim = fits.getdata(simpath)

outpath = filepath + '/6.indiv_sum'
if not os.path.exists(outpath):
    os.makedirs(outpath)
#______________________________________________________________________________

# define pair and control samples
sam = fits.getdata(filepath + '/4.samples.fits')
controls = sam['all_control'].astype('bool')
pairs = sam['all_pair'].astype('bool')
mass = sam['all_mass']

# Open Stacked Gradients
filf = fits.getdata(filepath + '/5.stacking/stacked_comb.fits')
#fil = fits.getdata(filepath + '/5.stacking/stacked_trim_comb.fits')

rbins = filf['rbins']
#pbin, pstd = fil[cal+'_pbin'].T, fil[cal+'_pstd'].T
#cbin, cstd = fil[cal+'_cbin'].T, fil[cal+'_cstd'].T

pbin, pstd = filf[cal+'_pbin'].T, filf[cal+'_pstd'].T
cbin, cstd = filf[cal+'_cbin'].T, filf[cal+'_cstd'].T
#______________________________________________________________________________
binsize=0.5
mbins = np.arange(9.0, 11.5, binsize)+binsize/2

# Number of pairs/controls in each mass bin
nnp = []
for i in range(len(mbins)):
    nnp.append(len(pairs[pairs&(np.abs(mass-mbins[i])<=binsize/2.)]))
nnc = []
for i in range(len(mbins)):
    nnc.append(len(controls[controls&(np.abs(mass-mbins[i])<=binsize/2.)]))

# Individual Gradient Plotting
fontsize=30
fig, axs = plt.subplots(1, 3, sharex=True, sharey=False, figsize=(20, 6))

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
for i in range(len(mbins)):
    mix = abs(mbins[i] - mass) < (mbins[1] - mbins[0])/2
        
    #axs[0].plot(rbins, cbinf[i], zorder=9998, c=colors[i], linestyle='--')
    
    #axs[0].plot(rbins, cbin[i], zorder=9998, c=colors[i], 
    #   label='log(M/M$_\odot$) = '+str(mbins[i]-binsize/2)+'$-$'+str(mbins[i]+binsize/2))
    
    axs[0].plot(rbins, cbin[i], zorder=9998, c=colors[i], 
       label='log(M/M$_\odot$) = '+str(mbins[i]-binsize/2)+'$-$'+str(mbins[i]+binsize/2)+' ('+str(nnc[i])+')')
    #axs[0].plot(rbins, cbin[i], zorder=9998, c=colors[i], linestyle='--')
    axs[0].fill_between(rbins, cbin[i]-cstd[i], cbin[i]+cstd[i], alpha=0.3,
                 antialiased=True, color=colors[i], zorder=9999)
        
    #axs[1].plot(rbins, pbinf[i], zorder=9998, c=colors[i], linestyle='--')
    axs[1].plot(rbins, pbin[i], zorder=9998, c=colors[i],
       label='log(M/M$_\odot$) = '+str(mbins[i]-binsize/2)+'$-$'+str(mbins[i]+binsize/2)+' ('+str(nnp[i])+')')
    #axs[1].plot(rbins, pbin[i], zorder=9998, c=colors[i], linestyle='--')
    axs[1].fill_between(rbins, pbin[i]-pstd[i], pbin[i]+pstd[i], alpha=0.3,
                 antialiased=True, color=colors[i], zorder=9999)
    
    dex = pbin[i] - cbin[i]
    #dexf = pbinf[i] - cbinf[i]
    dexstd = np.sqrt( cstd[i]**2 + pstd[i]**2 )
    
    axs[2].plot( [-0.2, 2.7] , [0.0, 0.0], c='k', linestyle = ':')
    for j in range(3):
        axs[j].axvline(1.5, linestyle='-.', c='k')
    axs[2].plot(rbins, dex, zorder=9998, c=colors[i])
    #axs[2].plot(rbins, dex, zorder=9998, c=colors[i], linestyle='--')
    axs[2].fill_between(rbins, dex-dexstd, dex+dexstd, alpha=0.3,
                 antialiased=True, color=colors[i], zorder=9999)

dex = np.nanmean(pbin - cbin, axis=0)
axs[2].plot(rbins, dex, c='k', linestyle='--', zorder=9999)

for i in range(2):
    leg = axs[i].legend(handlelength=0, handletextpad=0, fontsize=fontsize/2.25, 
             loc=4, fancybox=False, edgecolor='k')
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())
    for item in leg.legendHandles:
        item.set_visible(False)
    leg.get_frame().set_linewidth(2.0)
    
if cal == 'ewha':
    loc1 = 0.5
    loc2 = 0.1
elif cal == 'ssfr':
    loc1 = 0.5
    loc2 = 0.1
elif cal == 'sfr':
    loc1 = 2.0
    loc2 = 0.5
elif cal == 'mstar':
    loc1 = 2.0
    loc2 = 0.5
else:
    loc1 = 0.5
    loc2 = 0.1

ps.ticks(axs[0:2], xmajor=0.5, ymajor=loc1, xminor=0.1, yminor=loc2)
ps.ticks(axs[2], xmajor=0.5, ymajor=0.5, xminor=0.1, yminor=0.1)
ps.style(axs, fontsize=fontsize)
#ps.style(axs[1], fontsize=fontsize, labelleft=False)

axs[0].set_title('Controls', fontsize=fontsize)
axs[1].set_title('Pairs', fontsize=fontsize)
axs[2].set_title('Pairs - Controls', fontsize=fontsize)

axs[0].set_xlim(-0.2,2.7)

if cal == 'ewha':
    axs[0].set_ylim(0.45, 1.8)
    axs[1].set_ylim(0.45, 1.8)
    axs[2].set_ylim(-0.65, 0.95)
elif cal == 'ssfr':
    axs[0].set_ylim(-11.6, -9.3)     
    axs[1].set_ylim(-11.6, -9.3)
    axs[2].set_ylim(-0.65, 0.95)
elif cal == 'sfr':
    axs[0].set_ylim(-5.5, -0.5)
elif cal == 'mstar':
    axs[0].set_ylim(4.5, 10.5)
    
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
if cal == 'ewha':
    axs[0].set_ylabel(r'$\rm log(EW(H\alpha))$', fontsize=fontsize, labelpad=10)
    axs[1].set_ylabel(r'$\rm log(EW(H\alpha))$', fontsize=fontsize, labelpad=10)
    axs[2].set_ylabel(r'$\rm \Delta log(EW(H\alpha))$', fontsize=fontsize, labelpad=0)
elif cal == 'ssfr':
    axs[0].set_ylabel(r'$\rm log(sSFR, yr^{-1})$', fontsize=fontsize, labelpad=10)
    axs[1].set_ylabel(r'$\rm log(sSFR, yr^{-1})$', fontsize=fontsize, labelpad=10)
    axs[2].set_ylabel(r'$\rm \Delta log(sSFR, dex)$', fontsize=fontsize, labelpad=0)
elif cal == 'sfr':
    axs[0].set_ylabel(r'$\rm log(SFR, M_{\odot} yr^{-1} kpc^{-2})$', fontsize=fontsize, labelpad=20)
elif cal == 'mstar':
    axs[0].set_ylabel(r'$log(M/M_{\odot})$', fontsize=fontsize, labelpad=20)

plt.subplots_adjust(wspace=0.4)
axs[0].set_xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)
axs[1].set_xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)
axs[2].set_xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)


#plt.suptitle('Combined Pair Sample', fontsize=fontsize)
    
plt.savefig(outpath + '/'+cal+'_comb.pdf', bbox_inches='tight', overwrite=True)
plt.close('all')