#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 10:25:31 2020

Select a subsample of control galaxies for each galaxies pair. Select controls
based on stellar mass, redshift, (maybe) radial size, (maybe) local density. 
Select a set of ~## controls for each galaxy pair. 

Stack the individual profiles for the subset of controls and the take the
difference between the pair and stacked control profiles. These dex pair
profiles can then be stacked later by mass, sep, etc.

The code should:
    
    1. Build our full control and pairs samples from before
    2. Iterate over pairs
    3. Have a method for selecting the closest controls based on mass and z
    4. Track how many and which controls were used
    5. Stack the selected controls
    6. Take the difference between the individual pair profiles and the stacked
       set of controls
    7. Save the individual dex pair profiles

LATER: add in something to prevent the recycling of a few controls, maybe set a
flag for a control that is used

For some reason I'm getting the same dex profiles for each param?

@author: joshua
"""
def roundup(x):
    return int(m.ceil(x / 10.0)) * 10

from astropy.io import fits
import math as m
import numpy as np
import os, sys
import functions as f
import matplotlib.pyplot as plt
import plotstyle as ps
import random

# set filepath
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
T5path = f.T5path()
# Open drpall
drpallpath = T5path + '/MPL-8/drpall_uniq.fits'
drpall = fits.getdata(drpallpath)
cmd = drpall['CMD']
drmass = np.log10(drpall['nsa_elpetro_mass']/0.7**2)
z = drpall['nsa_z']
rad = drpall['nsa_elpetro_th50_r']
plateifu = drpall['plateifu']
mngtarg1 = drpall['mngtarg1']
mngtarg2 = drpall['mngtarg2']

# Open ncomp
samplepath = T5path + '/8.samples/ncomp.fits'
sam = fits.getdata(samplepath)

drpallmasterpath = T5path + '/6.master/drpall.fits'
drpall_master = fits.getdata(drpallmasterpath)

simpath = T5path + '/6.master/simard.fits'
sim = fits.getdata(simpath)

ap = ['nsa', 'sim']
sample = ['Oustide', 'Inside']
#for k in range(len(ap)):    

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

gal = 1 * (((mngtarg1 & 2L**10) != 0) & (mngtarg1 != 0)) + 2 * (((mngtarg1 & 2L**11) != 0) & (mngtarg1 != 0))

ncon = [5,10,20,50]
fontsize=20
for k in range(1):
    outpath = filepath + '/8.boostrap'
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    sam = fits.getdata(filepath + '/4.samples.fits')
    controls = sam[ap[k]+'_control'].astype('bool')
    pairs = sam[ap[k]+'_pair'].astype('bool')
    if ap[k]=='nsa':
        mass = drmass
    elif ap[k] == 'sim':
        mass = sam['corrected_mass']
    gradbin = fits.getdata(filepath + '/3.gradbin/gradbin_'+ap[k]+'.fits')
    
    rbins = gradbin['reff'][0]
    # Set the properties to build dex profiles for
    strng = ['ssfr']
    data = [gradbin[i] for i in strng]
    
    std = [gradbin[i+'_std'] for i in strng]
    
    ssfr = np.zeros((len(sam), len(rbins)))

    ssfr_std = np.zeros((len(sam), len(rbins)))
    
    #for i in range(len(sam[pairs])):
    for i in range(1):
        i = np.where(plateifu[pairs] == '8485-3704')[0][0]
        fig, axs = plt.subplots(1, 5, sharex=True, sharey=True, figsize=(64, 8))
        
        # Definse MaNGA Sample
        samp = gal[i]
        
        for y in range(len(ncon)):
            # Number of controls found for a pair
            N = []
            # Number of Interations run
            N_iter = []
            # times a control has been used
            con = np.zeros(len(np.where(controls==True)[0]))
            # Set limits for allowable mass range and redshift range
            lim_mass = 0.1
            #lim_redshift = 0.005
            #lim_radius = 0.05
            # Set minimum number of controls
            min_con = ncon[y]
            max_con = ncon[y]
                    
            cycles = 1000
            test = np.zeros((cycles, len(rbins)))
            # Cycle several times and record median values
            for l in range(cycles):
                n = 0
                n_iter = 0
                lim_m = 0
                #lim_z = 0
                #lim_r = 0.15
                while n < min_con:
                    lim_m += lim_mass
                    #lim_z += lim_redshift
                    #lim_r += lim_radius
                    # param differences
                    del_m = np.abs(mass[pairs][i] - mass[controls])
                    #del_z = np.abs(z[pairs][i]    -    z[controls])
                    #del_r = np.abs(rad[pairs][i]  -  rad[controls])
                    # select controls which fall within in limits
                    subset = (del_m < lim_m)&(gal[controls] == samp)
                    #record number of controls found
                    n = len(np.where(subset==True)[0])
                    n_iter += 1
        
                N.append(n)    
                N_iter.append(n_iter)
                # If control is used add 1
                con+=subset.astype(int)
                
                #radomly select 10 of controls
                selected = np.where(subset==True)[0]
                ran = random.sample(range(len(selected)), max_con)
                
                keep = selected[ran]
                
                newset = np.zeros(subset.shape).astype('bool')
                newset[keep] = True
                if n > 0:
                    # Loop over params to stack
                    pbin, pstd = f.stack(data[0][controls][newset], z=False, zbin=False, error='sem', trim = None)
                    dex = data[0][pairs][i] - pbin
                    dexstd = np.sqrt( std[0][pairs][i]**2 + pstd**2 )
                    ssfr[pairs][i] = dex[0]
                    ssfr_std[pairs][i] = dexstd[0]
                    
                    test[l] = dex
                    axs[y].plot(rbins, dex[0], c='lightgrey')
#______________________________________________________________________________
            # Plot bootstraping
            # median profiles of all cycles
            med = np.nanmedian(test, axis=0)
            err = np.nanstd(test, axis=0)
            
            if ncon[y] == 20:
                save = err
            
            #axs.errorbar(rbins, med, yerr=std, c='r', fmt='s')
            axs[y].plot(rbins, med, c=colors[y], zorder=9999)
            axs[y].fill_between(rbins, med+err, med-err, color=colors[y], zorder=999, alpha=0.25)
        
            pbin, pstd = f.stack(data[0][controls][subset], z=False, zbin=False, error='sem', trim = None)
            dex = data[0][pairs][i] - pbin
            dexstd = np.sqrt( std[0][pairs][i]**2 + pstd**2 )
        
            #axs[y].errorbar(rbins, dex[0], yerr=dexstd[0], linestyle='--', c='k', zorder=99999, capsize=3)
        
            axs[y].set_xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)
            axs[y].set_xlim(-0.1,2.7)
            
            #axs[y].set_ylim(-11.15, -9.35)
            axs[y].set_ylim(-1.0, 1.0)
        
            axs[y].set_title('Number of Controls = '+str(ncon[y]), fontsize=fontsize)
        
            ps.ticks(axs[y], xmajor=0.5, ymajor=0.25, xminor=0.1, yminor=0.05)
            ps.style(axs[y], fontsize=fontsize, labelleft=False)
                
            axs[4].plot(rbins, med, zorder=9999, label='Controls '+str(ncon[y]))
        
        axs[0].text( 0.07, 0.95, 'Control Pool Size = '+str(n), transform=axs[0].transAxes, fontsize=fontsize/1.5)
        axs[0].text( 0.07, 0.92, 'Iterations = '+str(n_iter), transform=axs[0].transAxes, fontsize=fontsize/1.5)
        
        axs[4].legend(fontsize=fontsize, fancybox=False, edgecolor='k')
        axs[4].set_xlabel(r'$\rm R/R_{eff}$', fontsize=fontsize)
        
        ps.ticks(axs[4], xmajor=0.5, ymajor=0.25, xminor=0.1, yminor=0.05)
        ps.style(axs[4], fontsize=fontsize, labelleft=False)
        
        plt.suptitle('Number of Cycles = '+str(cycles), fontsize=fontsize)
        ps.style(axs[0], fontsize=fontsize)
        axs[0].set_ylabel(r'$\rm \Delta log(sSFR, dex)$', fontsize=fontsize)
        plt.subplots_adjust(wspace=0.0)
        plt.savefig(outpath +'/'+ plateifu[pairs][i] +'.pdf', overwrite=True, bbox_inches='tight')
        plt.close('all')
        