#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Take the difference between the nuclear sSFR between pairs and controls.

"""
from astropy.io import fits
import os, sys
import numpy as np
import functions as f
import astropy.stats as astat

filepath = os.path.abspath(os.path.dirname(sys.argv[0]))

sam = fits.getdata(filepath + '/4.samples.fits')

T5path = '/Volumes/Ext_drive/photoObj8'

drpallpath = T5path + '/MPL-8/drpall_uniq.fits'
drpall = fits.getdata(drpallpath)
mngtarg1 = drpall['mngtarg1']
drmass = np.log10(drpall['nsa_elpetro_mass']/0.7**2)
rad = drpall['nsa_elpetro_th50_r']
z = drpall['nsa_z']

nuc = fits.getdata(filepath + '/10.nuc.fits')
SSFR_nsa = nuc['ssfr_nsa']
ssfr_STD_nsa = nuc['ssfr_std_nsa']
SSFR_sim = nuc['ssfr_sim']
ssfr_STD_sim = nuc['ssfr_std_sim']

SFR_nsa = nuc['sfr_nsa']
sfr_STD_nsa = nuc['sfr_std_nsa']
SFR_sim = nuc['sfr_sim']
sfr_STD_sim = nuc['sfr_std_sim']

EWHA_nsa = nuc['ewha_nsa']
ewha_STD_nsa = nuc['ewha_std_nsa']
EWHA_sim = nuc['ewha_sim']
ewha_STD_sim = nuc['ewha_std_sim']

DOP16_nsa = nuc['dop16_nsa']
DOP16_STD_nsa = nuc['dop16_std_nsa']
DOP16_sim = nuc['dop16_sim']
DOP16_STD_sim = nuc['dop16_std_sim']

M08_nsa = nuc['m08_nsa']
M08_STD_nsa = nuc['m08_std_nsa']
M08_sim = nuc['m08_sim']
M08_STD_sim = nuc['m08_std_sim']

gal = 1 * (((mngtarg1 & 2L**10) != 0) & (mngtarg1 != 0)) + 2 * (((mngtarg1 & 2L**11) != 0) & (mngtarg1 != 0))

ssfr_dex = np.zeros(len(drpall)) * np.nan
ssfr_dexstd = np.zeros(len(drpall)) * np.nan

sfr_dex = np.zeros(len(drpall)) * np.nan
sfr_dexstd = np.zeros(len(drpall)) * np.nan

ewha_dex = np.zeros(len(drpall)) * np.nan
ewha_dexstd = np.zeros(len(drpall)) * np.nan

d16_dex = np.zeros(len(drpall)) * np.nan
d16_dexstd = np.zeros(len(drpall)) * np.nan

m08_dex = np.zeros(len(drpall)) * np.nan
m08_dexstd = np.zeros(len(drpall)) * np.nan

# Set limits for allowable mass range and redshift range
# limits from Pan+19, mass=0.1, z=0.005, rad=20%
lim_mass = 0.1
lim_redshift = 0.025
lim_radius = 0.05

# Set minimum number of controls
min_con = 20
max_con = 20
ap = ['nsa', 'sim']
for k in range(len(ap)): 
    # Loop over all galaxies    
    if ap[k]=='nsa':
        mass = drmass
        controls = sam['nsa_control'].astype('bool')
        pairs = sam['nsa_pair'].astype('bool')
        ssfr = SSFR_nsa
        ssfr_std = ssfr_STD_nsa
        sfr = SFR_nsa
        sfr_std = sfr_STD_nsa
        ewha = EWHA_nsa
        ewha_std = ewha_STD_nsa
        d16 = DOP16_nsa
        d16_std = DOP16_STD_nsa
        m08 = M08_nsa
        m08_std = M08_STD_nsa
    elif ap[k] == 'sim':
        mass = sam['corrected_mass']
        controls = sam['sim_control'].astype('bool')
        pairs = sam['sim_pair'].astype('bool')
        ssfr = SSFR_sim
        ssfr_std = ssfr_STD_sim
        sfr = SFR_sim
        sfr_std = sfr_STD_sim
        ewha = EWHA_sim
        ewha_std = ewha_STD_sim
        d16 = DOP16_sim
        d16_std = DOP16_STD_sim
        m08 = M08_sim
        m08_std = M08_STD_sim
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
                    #subset = (del_m < lim_m)&(gal[controls] == samp)
                    subset = (del_m < lim_m)&(del_z < lim_z)&(del_r < lim_r * rad[i])&(gal[controls] == samp)
                    n = len(np.where(subset==True)[0])    
    
                selected = np.where(subset==True)[0]
                keep = selected[0:max_con]
                
                newset = np.zeros(subset.shape).astype('bool')
                newset[keep] = True
                
                if n > 0:
                    # Loop over params to stack
                    con = ssfr[controls][newset]
                    constd = astat.biweight_scale(con[np.isfinite(con)&~np.isnan(con)], 
                             c=6)/np.sqrt(len(con[np.isfinite(con)&~np.isnan(con)]))
                    ssfr_dex[i] = ssfr[i] - np.nanmedian(con)
                    ssfr_dexstd[i] = np.sqrt(ssfr_std[i]**2  + constd**2)
                    
                    con = sfr[controls][newset]
                    constd = astat.biweight_scale(con[np.isfinite(con)&~np.isnan(con)], 
                             c=6)/np.sqrt(len(con[np.isfinite(con)&~np.isnan(con)])) 
                    sfr_dex[i] = sfr[i] - np.nanmedian(con)
                    sfr_dexstd[i] = np.sqrt(sfr_std[i]**2  + constd**2)
                    
                    con = ewha[controls][newset]
                    constd = astat.biweight_scale(con[np.isfinite(con)&~np.isnan(con)], 
                             c=6)/np.sqrt(len(con[np.isfinite(con)&~np.isnan(con)])) 
                    ewha_dex[i] = ewha[i] - np.nanmedian(con)
                    ewha_dexstd[i] = np.sqrt(ewha_std[i]**2  + constd**2)
                    ''
                    con = d16[controls][newset]
                    constd = astat.biweight_scale(con[np.isfinite(con)&~np.isnan(con)], 
                             c=6)/np.sqrt(len(con[np.isfinite(con)&~np.isnan(con)])) 
                    d16_dex[i] = d16[i] - np.nanmedian(con)
                    d16_dexstd[i] = np.sqrt(d16_std[i]**2  + constd**2)
                    
                    con = m08[controls][newset]
                    constd = astat.biweight_scale(con[np.isfinite(con)&~np.isnan(con)], 
                             c=6)/np.sqrt(len(con[np.isfinite(con)&~np.isnan(con)])) 
                    m08_dex[i] = m08[i] - np.nanmedian(con)
                    m08_dexstd[i] = np.sqrt(m08_std[i]**2  + constd**2)
            #break
        f.update_progress((i+1.0)/np.float64(len(sam)))
            
    cols = []
    
    cols.append(fits.Column(name='ssfr_dex', format='1E', array=ssfr_dex))
    cols.append(fits.Column(name='ssfr_std', format='1E', array=ssfr_dexstd))
    cols.append(fits.Column(name='sfr_dex', format='1E', array=sfr_dex))
    cols.append(fits.Column(name='sfr_std', format='1E', array=sfr_dexstd))
    cols.append(fits.Column(name='ewha_dex', format='1E', array=ewha_dex))
    cols.append(fits.Column(name='ewha_std', format='1E', array=ewha_dexstd))
    
    cols.append(fits.Column(name='DOP16_dex', format='1E', array=d16_dex))
    cols.append(fits.Column(name='DOP16_std', format='1E', array=d16_dexstd))
    cols.append(fits.Column(name='M08_dex', format='1E', array=m08_dex))
    cols.append(fits.Column(name='M08_std', format='1E', array=m08_dexstd))
    
    cols = fits.ColDefs(cols)
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(filepath+'/11.nuc_dex/nuc_dex_'+ap[k]+'.fits', overwrite=True)