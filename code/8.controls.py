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
#T5path = f.T5path()
T5path = '/Volumes/Ext_drive/photoObj8'

SAVE8 = []

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

drpallmasterpath = T5path + '/6.master/drpall.fits'
drpall_master = fits.getdata(drpallmasterpath)

simpath = T5path + '/6.master/simard.fits'
sim = fits.getdata(simpath)

# 1 == primary, 2 == secondary
gal = 1 * (((mngtarg1 & 2L**10) != 0) & (mngtarg1 != 0)) + 2 * (((mngtarg1 & 2L**11) != 0) & (mngtarg1 != 0))# + 3 *(((mngtarg1 & 2L**12) != 0) & (mngtarg1 != 0))

# Error from bootstrapping
berr = 0.0

gradbin = fits.getdata(filepath + '/3.gradbin/gradbin_nsa.fits')
rbins = gradbin['reff'][0]
strng = ['ewha', 'ssfr', 'sfr', 'mstar', 'dop16', 'm08']
nsa_data = np.asarray([gradbin[i] for i in strng])

gradbin = fits.getdata(filepath + '/3.gradbin/gradbin_sim.fits')
rbins = gradbin['reff'][0]

sim_data = np.asarray([gradbin[i] for i in strng])
std = [gradbin[i+'_std'] for i in strng]

ap = ['nsa', 'sim']
sample = ['Oustide', 'Inside']
for k in range(len(ap)):    
    outpath = filepath + '/8.controls_'+ap[k]
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    sam = fits.getdata(filepath + '/4.samples.fits')
    controls = sam[ap[k]+'_control'].astype('bool')
    pairs = sam[ap[k]+'_pair'].astype('bool')
    sim_mendel = sam['sim_mendel'].astype('bool')
    
    data = nsa_data
    if ap[k]=='nsa':
        mass = drmass
    elif ap[k] == 'sim':
        mass = sam['corrected_mass']
        data[:,sim_mendel] = sim_data[:,sim_mendel]

    ewha = np.zeros((len(sam), len(rbins)))
    ssfr = np.zeros((len(sam), len(rbins)))
    sfr = np.zeros((len(sam), len(rbins)))
    mstar = np.zeros((len(sam), len(rbins)))
    dop16 = np.zeros((len(sam), len(rbins)))
    m08 = np.zeros((len(sam), len(rbins)))
    
    ewha_std = np.zeros((len(sam), len(rbins)))
    ssfr_std = np.zeros((len(sam), len(rbins)))
    sfr_std = np.zeros((len(sam), len(rbins)))
    mstar_std = np.zeros((len(sam), len(rbins)))
    dop16_std = np.zeros((len(sam), len(rbins)))
    m08_std = np.zeros((len(sam), len(rbins)))
    
    # Number of controls found for a pair
    N = np.zeros((len(sam)))
    # Number of Interations run
    N_iter = np.zeros((len(sam)))
    # times a control has been used    
    con = np.zeros(len(np.where(controls==True)[0]))
    
    # Set limits for allowable mass range and redshift range
    # limits from Pan+19, mass=0.1, z=0.005, rad=20%
    lim_mass = 0.1
    lim_redshift = 0.025
    lim_radius = 0.05
    # Set minimum number of controls
    min_con = 20
    max_con = 20
    # Loop over all galaxies
    for i in range(len(sam)):
        # Only do operation for pairs
        if pairs[i]==True and rad[i]!=-9999:
            # Iteratively expand limits until minimum number of controls is found
        
            # Definse MaNGA Sample
            samp = gal[i]
            
            # Exclude gals which aren't Primary or Secondary I should do this in
            # sample.py
            if gal[i]!=0 and gal[i] <= 4:
            
                n = 0
                n_iter = 0
                lim_m = 0
                lim_z = 0
                lim_r = 0.15
                while n < min_con:
                    lim_m += lim_mass
                    lim_z += lim_redshift
                    lim_r += lim_radius
                    
                    # param differences
                    del_m = np.abs(mass[i] - mass[controls])
                    del_z = np.abs(z[i]    -    z[controls])
                    del_r = np.abs(rad[i]  -  rad[controls])
                    
                    # select controls which fall within in limits
                    subset = (del_m < lim_m)&(del_z < lim_z)&(del_r < lim_r * rad[i])&(gal[controls] == samp)
                    #subset = (del_m < lim_m)&(gal[controls] == samp)
                    
                    # record number of controls found
                    n = len(np.where(subset==True)[0])
                    n_iter += 1
        
                N[i] = n    
                N_iter[i] = n_iter
                # If control is used add 1
                con+=subset.astype(int)
                
                #radomly select N controls
                '''
                Generalize the error from the randomly selected controls as 
                being 0.1 dex
                
                For random selection, just choose the first 20
                '''

                selected = np.where(subset==True)[0]
                keep = selected[0:max_con]
                
                newset = np.zeros(subset.shape).astype('bool')
                newset[keep] = True
                
                if n > 0:
                    # Loop over params to stack
                    pbin, pstd = f.stack(data[0][controls][newset], z=False, zbin=False, error='sem', trim = None)
                    dex = data[0][i] - pbin
                    dexstd = np.sqrt( std[0][i]**2 + berr**2 )
                    ewha[i] = dex
                    ewha_std[i] = dexstd
                    
                    pbin, pstd = f.stack(data[1][controls][newset], z=False, zbin=False, error='sem', trim = None)
                    dex = data[1][i] - pbin
                    dexstd = np.sqrt( std[1][i]**2 + berr**2 )
                    ssfr[i] = dex
                    ssfr_std[i] = dexstd
                                                       
                    pbin, pstd = f.stack(data[2][controls][newset], z=False, zbin=False, error='sem', trim = None)
                    dex = data[2][i] - pbin
                    dexstd = np.sqrt( std[2][i]**2 + berr**2 )
                    sfr[i] = dex
                    sfr_std[i] = dexstd
                    
                    pbin, pstd = f.stack(data[3][controls][newset], z=False, zbin=False, error='sem', trim = None)
                    dex = data[3][i] - pbin
                    dexstd = np.sqrt( std[3][i]**2 + berr**2 )
                    mstar[i] = dex
                    mstar_std[i] = dexstd
                    
                    pbin, pstd = f.stack(data[4][controls][newset], z=False, zbin=False, error='sem', trim = None)
                    dex = data[4][i] - pbin
                    dexstd = np.sqrt( std[4][i]**2 + berr**2 )
                    dop16[i] = dex
                    dop16_std[i] = dexstd

                    pbin, pstd = f.stack(data[5][controls][newset], z=False, zbin=False, error='sem', trim = None)
                    dex = data[5][i] - pbin
                    dexstd = np.sqrt( std[5][i]**2 + berr**2 )
                    m08[i] = dex
                    m08_std[i] = dexstd

        f.update_progress((i+1.0)/np.float64(len(sam)))       
    
    cols = []
    # add to FITS file    
    cols.append(fits.Column(name='ewha_dex', format=str(sfr.shape[1])+'E', array=ewha))
    cols.append(fits.Column(name='ewha_std', format=str(sfr.shape[1])+'E', array=ewha_std))
    
    cols.append(fits.Column(name='ssfr_dex', format=str(sfr.shape[1])+'E', array=ssfr))
    cols.append(fits.Column(name='ssfr_std', format=str(sfr.shape[1])+'E', array=ssfr_std))
    
    cols.append(fits.Column(name='sfr_dex', format=str(sfr.shape[1])+'E', array=sfr))
    cols.append(fits.Column(name='sfr_std', format=str(sfr.shape[1])+'E', array=sfr_std))
    
    cols.append(fits.Column(name='mstar_dex', format=str(sfr.shape[1])+'E', array=mstar))
    cols.append(fits.Column(name='mstar_std', format=str(sfr.shape[1])+'E', array=mstar_std))

    cols.append(fits.Column(name='dop16_dex', format=str(sfr.shape[1])+'E', array=dop16))
    cols.append(fits.Column(name='dop16_std', format=str(sfr.shape[1])+'E', array=dop16_std))

    cols.append(fits.Column(name='m08_dex', format=str(sfr.shape[1])+'E', array=m08))
    cols.append(fits.Column(name='m08_std', format=str(sfr.shape[1])+'E', array=m08_std))
    
    cols.append(fits.Column(name='N_control', format='1E', array=N))
    cols.append(fits.Column(name='N_iter', format='1E', array=N_iter))
    
    
    cols = fits.ColDefs(cols)
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(outpath+'/control_'+ap[k]+'.fits', overwrite=True)
#______________________________________________________________________________
    # Plotting Statitistic of Controls
    fontsize = 25
    
    N = N[N!=0]
    N_iter = N_iter[N_iter!=0]
    print len(N_iter)
    fig, axs = plt.subplots(1, 3, sharex=False, sharey=False, figsize=(24, 8))
    axs[0].hist(N, np.arange(10, 180, 20), color='deepskyblue', edgecolor='k', linewidth=2.0)
    axs[0].set_xlabel('Number of Controls Found per Pair', fontsize=fontsize/1.25)
    axs[0].set_ylabel('Frequency', fontsize=fontsize/1.25)
    
    axs[1].hist(con, bins = np.arange(-1.0,con.max() + 1, 2), color='deepskyblue', 
       edgecolor='k', linewidth=2.0)
    axs[1].set_yscale('log')
    axs[1].set_xlabel('Number of Times a Given Control is Used', fontsize=fontsize/1.25)
    axs[1].set_ylabel('Frequency', fontsize=fontsize/1.25)    
    
    axs[2].hist(N_iter, np.arange(0.5,max(N_iter) + 0.5 + 1, 1), color='deepskyblue', edgecolor='k', linewidth=2.0)
    axs[2].set_yscale('log')
    axs[2].set_xlabel('Iterations Get '+str(min_con)+' Controls', fontsize=fontsize/1.25)
    axs[2].set_ylabel('Frequency', fontsize=fontsize/1.25)
    
    axs[1].xaxis.set_major_locator(plt.MultipleLocator(2))
    axs[1].xaxis.set_minor_locator(plt.MultipleLocator(1))
    
    ps.ticks(axs[0], xmajor=20, ymajor=5, xminor=5, yminor=1)
    axs[2].xaxis.set_major_locator(plt.MultipleLocator(1))
    ps.style(axs, fontsize=fontsize)
    
    plt.suptitle('Control Sample Statistics for the '+sample[k]+' Pair Sample', fontsize=fontsize)
    
    plt.savefig(outpath + '/control_stats_'+ap[k]+'.pdf', bbox_inches='tight', overwrite=True)
    plt.close('all')