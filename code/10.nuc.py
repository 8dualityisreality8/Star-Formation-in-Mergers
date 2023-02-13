#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:23:35 2020

Here I want to extract the nuclear (<0.5Re) values of enhancement take the mean
of the dex within the inner 0.5Re of each paired galaxy and save them into a 
single catalog.

@author: joshua
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

rbins = np.arange(0,2.8, 0.2)

gal = 1 * (((mngtarg1 & 2L**10) != 0) & (mngtarg1 != 0)) + 2 * (((mngtarg1 & 2L**11) != 0) & (mngtarg1 != 0))

outpath = filepath + '/11.nuc_dex'
if not os.path.exists(outpath):
    os.makedirs(outpath)

nuc = 0.5

SSFR_nsa = np.zeros(len(drpall))
STD_ssfr_nsa = np.zeros(len(drpall))
SSFR_sim = np.zeros(len(drpall))
STD_ssfr_sim = np.zeros(len(drpall))

SFR_nsa = np.zeros(len(drpall))
STD_sfr_nsa = np.zeros(len(drpall))
SFR_sim = np.zeros(len(drpall))
STD_sfr_sim = np.zeros(len(drpall))

EWHA_nsa = np.zeros(len(drpall))
STD_EWHA_nsa = np.zeros(len(drpall))
EWHA_sim = np.zeros(len(drpall))
STD_EWHA_sim = np.zeros(len(drpall))

DOP16_nsa = np.zeros(len(drpall))
STD_DOP16_nsa = np.zeros(len(drpall))
DOP16_sim = np.zeros(len(drpall))
STD_DOP16_sim = np.zeros(len(drpall))

M08_nsa = np.zeros(len(drpall))
STD_M08_nsa = np.zeros(len(drpall))
M08_sim = np.zeros(len(drpall))
STD_M08_sim = np.zeros(len(drpall))
for i in range(len(drpall)):
    
    plateifu = drpall['plateifu'][i]
        
    grad = fits.getdata(filepath + '/1.metalgrad/'+plateifu+'.fits')
    
    bpt = grad['bpt'][0]
    ewha = grad['ewha'][0]
    
    snr = grad['snr_ha'][0].astype('bool')& (ewha>=6)#&((bpt==1)|(bpt==2))
    
    lha = grad['LogLha'][0,snr]
    mstar = grad['mstar'][0,snr]
    d16 = grad['DOP16'][0,snr]
    m08 = grad['M08'][0,snr]

    
    sfr = f.SFR_M11(10**lha)
    ssfr = np.log10(f.divide(sfr, 10**mstar))
    sfr = np.log10(sfr)
    
    ewha = np.log10(grad['ewha'][0,snr])
    
    reff_nsa = grad['rad_nsa'][0,snr]
    reff_sim = grad['rad_sim'][0,snr]

    # sSFR
    SSFR_nsa[i] = np.nanmedian(ssfr[reff_nsa<=0.5])
    std = astat.biweight_scale(ssfr[~np.isnan(ssfr)&np.isfinite(ssfr)&(reff_nsa<=0.5)], 
                                    c=6)/np.sqrt(len(ssfr[~np.isnan(ssfr)&np.isfinite(ssfr)&(reff_nsa<=0.5)]))
    STD_ssfr_nsa[i] = std
    
    SSFR_sim[i] = np.nanmedian(ssfr[reff_sim<=0.5])
    std = astat.biweight_scale(ssfr[~np.isnan(ssfr)&np.isfinite(ssfr)&(reff_sim<=0.5)], 
                                c=6)/np.sqrt(len(ssfr[~np.isnan(ssfr)&np.isfinite(ssfr)&(reff_sim<=0.5)]))
    STD_ssfr_sim[i] = std
    
    # SFR
    SFR_nsa[i] = np.nanmedian(sfr[reff_nsa<=0.5])
    std = astat.biweight_scale(sfr[~np.isnan(sfr)&np.isfinite(sfr)&(reff_nsa<=0.5)], 
                                    c=6)/np.sqrt(len(sfr[~np.isnan(sfr)&np.isfinite(sfr)&(reff_nsa<=0.5)]))
    STD_sfr_nsa[i] = std
    
    SFR_sim[i] = np.nanmedian(sfr[reff_sim<=0.5])
    std = astat.biweight_scale(sfr[~np.isnan(sfr)&np.isfinite(sfr)&(reff_sim<=0.5)], 
                                c=6)/np.sqrt(len(sfr[~np.isnan(sfr)&np.isfinite(sfr)&(reff_sim<=0.5)]))
    STD_sfr_sim[i] = std
    
    # EWHA
    EWHA_nsa[i] = np.nanmedian(ewha[reff_nsa<=0.5])
    std = astat.biweight_scale(ewha[~np.isnan(ewha)&np.isfinite(ewha)&(reff_nsa<=0.5)], 
                                    c=6)/np.sqrt(len(ewha[~np.isnan(ewha)&np.isfinite(ewha)&(reff_nsa<=0.5)]))
    STD_EWHA_nsa[i] = std
    
    EWHA_sim[i] = np.nanmedian(ewha[reff_sim<=0.5])
    std = astat.biweight_scale(ewha[~np.isnan(ewha)&np.isfinite(ewha)&(reff_sim<=0.5)], 
                                    c=6)/np.sqrt(len(ewha[~np.isnan(ewha)&np.isfinite(ewha)&(reff_sim<=0.5)]))
    STD_EWHA_sim[i] = std
    
    # DOP16
    DOP16_nsa[i] = np.nanmedian(d16[reff_nsa<=0.5])
    std = astat.biweight_scale(d16[~np.isnan(d16)&np.isfinite(d16)&(reff_nsa<=0.5)], 
                                    c=6)/np.sqrt(len(d16[~np.isnan(d16)&np.isfinite(d16)&(reff_nsa<=0.5)]))
    STD_DOP16_nsa[i] = std
    
    DOP16_sim[i] = np.nanmedian(d16[reff_sim<=0.5])
    std = astat.biweight_scale(d16[~np.isnan(d16)&np.isfinite(d16)&(reff_sim<=0.5)], 
                                    c=6)/np.sqrt(len(d16[~np.isnan(d16)&np.isfinite(d16)&(reff_sim<=0.5)]))
    STD_DOP16_sim[i] = std
    
    # M08
    M08_nsa[i] = np.nanmedian(m08[reff_nsa<=0.5])
    std = astat.biweight_scale(m08[~np.isnan(m08)&np.isfinite(m08)&(reff_nsa<=0.5)], 
                                    c=6)/np.sqrt(len(m08[~np.isnan(m08)&np.isfinite(m08)&(reff_nsa<=0.5)]))
    STD_M08_nsa[i] = std
    
    M08_sim[i] = np.nanmedian(m08[reff_sim<=0.5])
    std = astat.biweight_scale(m08[~np.isnan(m08)&np.isfinite(m08)&(reff_sim<=0.5)], 
                                    c=6)/np.sqrt(len(m08[~np.isnan(m08)&np.isfinite(m08)&(reff_sim<=0.5)]))
    STD_M08_sim[i] = std
    
    f.update_progress((i+1.0)/np.float64(len(drpall)))      
#______________________________________________________________________________
# all
nsa_pair = sam['nsa_pair'].astype('bool')
sim_pair = sam['sim_pair'].astype('bool')

ssfr_data = np.zeros((drpall.shape[0])) * np.nan
ssfr_data[nsa_pair] = SSFR_nsa[nsa_pair]
ssfr_data[sim_pair] = SSFR_sim[sim_pair]

ssfr_std = np.zeros((drpall.shape[0])) * np.nan
ssfr_std[nsa_pair] = STD_ssfr_nsa[nsa_pair]
ssfr_std[sim_pair] = STD_ssfr_sim[sim_pair]

sfr_data = np.zeros((drpall.shape[0])) * np.nan
sfr_data[nsa_pair] = SFR_nsa[nsa_pair]
sfr_data[sim_pair] = SFR_sim[sim_pair]

sfr_std = np.zeros((drpall.shape[0])) * np.nan
sfr_std[nsa_pair] = STD_sfr_nsa[nsa_pair]
sfr_std[sim_pair] = STD_sfr_sim[sim_pair]

ewha_data = np.zeros((drpall.shape[0])) * np.nan
ewha_data[nsa_pair] = EWHA_nsa[nsa_pair]
ewha_data[sim_pair] = EWHA_sim[sim_pair]

ewha_std = np.zeros((drpall.shape[0])) * np.nan
ewha_std[nsa_pair] = STD_EWHA_nsa[nsa_pair]
ewha_std[sim_pair] = STD_EWHA_sim[sim_pair]

d16_data = np.zeros((drpall.shape[0])) * np.nan
d16_data[nsa_pair] = DOP16_nsa[nsa_pair]
d16_data[sim_pair] = DOP16_sim[sim_pair]

d16_std = np.zeros((drpall.shape[0])) * np.nan
d16_std[nsa_pair] = STD_DOP16_nsa[nsa_pair]
d16_std[sim_pair] = STD_DOP16_sim[sim_pair]

m08_data = np.zeros((drpall.shape[0])) * np.nan
m08_data[nsa_pair] = M08_nsa[nsa_pair]
m08_data[sim_pair] = M08_sim[sim_pair]

m08_std = np.zeros((drpall.shape[0])) * np.nan
m08_std[nsa_pair] = STD_M08_nsa[nsa_pair]
m08_std[sim_pair] = STD_M08_sim[sim_pair]

cols = []
cols.append(fits.Column(name='sfr_nsa', format='1E', array=SFR_nsa))
cols.append(fits.Column(name='sfr_std_nsa', format='1E', array=STD_sfr_nsa))
cols.append(fits.Column(name='ssfr_nsa', format='1E', array=SSFR_nsa))
cols.append(fits.Column(name='ssfr_std_nsa', format='1E', array=STD_ssfr_nsa))
cols.append(fits.Column(name='ewha_nsa', format='1E', array=EWHA_nsa))
cols.append(fits.Column(name='ewha_std_nsa', format='1E', array=STD_EWHA_nsa))
cols.append(fits.Column(name='DOP16_nsa', format='1E', array=DOP16_nsa))
cols.append(fits.Column(name='DOP16_std_nsa', format='1E', array=STD_DOP16_nsa))
cols.append(fits.Column(name='M08_nsa', format='1E', array=M08_nsa))
cols.append(fits.Column(name='M08_std_nsa', format='1E', array=STD_M08_nsa))
cols.append(fits.Column(name='sfr_sim', format='1E', array=SFR_sim))
cols.append(fits.Column(name='sfr_std_sim', format='1E', array=STD_sfr_sim))
cols.append(fits.Column(name='ssfr_sim', format='1E', array=SSFR_sim))
cols.append(fits.Column(name='ssfr_std_sim', format='1E', array=STD_ssfr_sim))
cols.append(fits.Column(name='ewha_sim', format='1E', array=EWHA_sim))
cols.append(fits.Column(name='ewha_std_sim', format='1E', array=STD_EWHA_sim))
cols.append(fits.Column(name='DOP16_sim', format='1E', array=DOP16_sim))
cols.append(fits.Column(name='DOP16_std_sim', format='1E', array=STD_DOP16_sim))
cols.append(fits.Column(name='M08_sim', format='1E', array=M08_sim))
cols.append(fits.Column(name='M08_std_sim', format='1E', array=STD_M08_sim))
cols.append(fits.Column(name='sfr_all', format='1E', array=sfr_data))
cols.append(fits.Column(name='sfr_std_all', format='1E', array=sfr_std))
cols.append(fits.Column(name='ssfr_all', format='1E', array=ssfr_data))
cols.append(fits.Column(name='ssfr_std_all', format='1E', array=ssfr_std))
cols.append(fits.Column(name='ewha_all', format='1E', array=ewha_data))
cols.append(fits.Column(name='ewha_std_all', format='1E', array=ewha_std))
cols.append(fits.Column(name='DOP16_all', format='1E', array=d16_data))
cols.append(fits.Column(name='DOP16_std_all', format='1E', array=d16_std))
cols.append(fits.Column(name='M08_all', format='1E', array=m08_data))
cols.append(fits.Column(name='M08_std_all', format='1E', array=m08_std))

cols = fits.ColDefs(cols)
tbhdu = fits.BinTableHDU.from_columns(cols)
tbhdu.writeto(outpath+'/nuc.fits', overwrite=True)