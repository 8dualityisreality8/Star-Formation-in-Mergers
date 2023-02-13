#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:56:29 2019

An overall summary file for the binned gradients of all _uniq galaxies. Also 
calculates SFR and sSFR from luminosity and mass. Use ncomp.fits later to 
select controls and pairs from this file later.

@author: joshua
"""
from astropy.io import fits
import numpy as np
import os, sys
import functions as f
    
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
datapath = os.path.abspath(os.path.join(filepath, '..', 'data'))

T5path = '/Volumes/Ext_drive/photoObj8'

drpallpath = T5path + '/MPL-8/drpall_uniq.fits'
drpall = fits.getdata(drpallpath)

mngtarg1 = drpall['mngtarg1']
gal = 1*(((mngtarg1 & 2L**10)!=0) & ((mngtarg1 != 0)))+2*(((mngtarg1 & 2L**11)!=0) & ((mngtarg1 != 0)))

outpath = filepath + '/3.gradbin'
if not os.path.exists(outpath):
    os.makedirs(outpath)

binsize = 0.20
rbinc = np.arange(0.0, 2.6+binsize, binsize)
rbins = np.append(rbinc - binsize/2, rbinc[-1]+binsize/2)

strings = ['Writing NSA Elliptical Petrosian Profiles', 
           'Writing Simard+11 Profiles']
ap = ['nsa', 'sim']
scale = ['nsa_elpetro_th50_r', 'sim_RCHL_R']
ba = ['nsa_elpetro_ba', 'sim_ba']

for j in range(2):
    print '\n'+strings[j]+'\n'

    good_e = np.zeros((len(drpall), 1))*np.nan

    Reff = np.zeros((len(drpall), len(rbinc)))

    EWHA = np.zeros((len(drpall), len(rbinc)))*np.nan
    EWHA_std = np.zeros((len(drpall), len(rbinc)))*np.nan
    SFR = np.zeros((len(drpall), len(rbinc)))*np.nan
    SFR_std = np.zeros((len(drpall), len(rbinc)))*np.nan
    SSFR = np.zeros((len(drpall), len(rbinc)))*np.nan
    SSFR_std = np.zeros((len(drpall), len(rbinc)))*np.nan
    NSS = np.zeros((len(drpall), len(rbinc)))*np.nan
    
    SSFR_ha = np.zeros((len(drpall), len(rbinc)))*np.nan
    SSFR_std_ha = np.zeros((len(drpall), len(rbinc)))*np.nan
    NSS_ha = np.zeros((len(drpall), len(rbinc)))*np.nan
    
    SSFR_ma = np.zeros((len(drpall), len(rbinc)))*np.nan
    SSFR_std_ma = np.zeros((len(drpall), len(rbinc)))*np.nan
    NSS_ma = np.zeros((len(drpall), len(rbinc)))*np.nan
    
    MSTAR = np.zeros((len(drpall), len(rbinc)))*np.nan
    MSTAR_std = np.zeros((len(drpall), len(rbinc)))*np.nan
    
    DOP16 = np.zeros((len(drpall), len(rbinc)))*np.nan
    DOP16_std = np.zeros((len(drpall), len(rbinc)))*np.nan
    M08 = np.zeros((len(drpall), len(rbinc)))*np.nan
    M08_std = np.zeros((len(drpall), len(rbinc)))*np.nan
    PP04 = np.zeros((len(drpall), len(rbinc)))*np.nan
    PP04_std = np.zeros((len(drpall), len(rbinc)))*np.nan
    
    Pix_kpc = np.zeros((len(drpall), 1))*np.nan

    for i in range(len(drpall)):
        plateifu = drpall['plateifu'][i]
    
        grad = fits.getdata(filepath + '/1.metalgrad/'+plateifu+'.fits')
        
        rad_scale = grad[scale[j]]
        
        if rad_scale==0:
            good_e[i] = 0
        else:
            good_e[i] = 1
            
            Reff[i] = rbinc
                        
            mstar = grad['mstar'][0]
            reff = grad['rad_'+ap[j]][0]
            bpt = grad['bpt'][0]
            ewha = grad['ewha'][0]
            lha = grad['loglha'][0]
            
            dop16 = grad['dop16'][0]
            m08 = grad['m08'][0]
            pp04 = grad['pp04'][0]
            
            mask = grad['snr_ha'][0].astype('bool') & (ewha>=6)
            
            reff = reff[mask]
            ewha = ewha[mask]
            lha = lha[mask]        
            mstar = mstar[mask]
 
            dop16 = dop16[mask]
            m08 = m08[mask]
            pp04 = pp04[mask]
    
            EWHA[i], EWHA_std[i], ne = f.azimuth_profiles(reff, np.log10(ewha), rbins)
            DOP16[i], DOP16_std[i], empty = f.azimuth_profiles(reff, np.log10(dop16), rbins)
            M08[i], M08_std[i], empty = f.azimuth_profiles(reff, np.log10(m08), rbins)
            PP04[i], PP04_std[i], empty = f.azimuth_profiles(reff, np.log10(pp04), rbins)
            
            pixkpc = (f.cosmo(grad[ap[j]+'_z']) * 1.0)**2 / np.sqrt( (grad[ba[j]]**2 - 0.13**2) / (1 - 0.13**2) )
            Pix_kpc[i] = pixkpc
            sfr = f.SFR_M11(10**lha)
            
            ssfr = f.divide(sfr, 10**mstar)
            SSFR[i], SSFR_std[i], NSS[i] = f.azimuth_profiles(reff, np.log10(ssfr), rbins)  
            if pixkpc != 0:
                SFR[i], SFR_std[i], ns = f.azimuth_profiles(reff, np.log10(sfr/pixkpc), rbins)
                MSTAR[i], MSTAR_std[i], nm = f.azimuth_profiles(reff, np.log10(10**mstar/pixkpc), rbins)

        f.update_progress((i+1.0)/np.float64(len(drpall)))

    cols = []
    cols.append(fits.Column(name='plateifu', format='11A', array=drpall['plateifu']))
    cols.append(fits.Column(name='good_e', format='1I', array=good_e))
    cols.append(fits.Column(name='reff', format=str(len(rbinc))+'E', array=Reff))
    cols.append(fits.Column(name='pixkpc', format='1E', array=Pix_kpc))
    cols.append(fits.Column(name='EWHA', format=str(len(rbinc))+'E', array=EWHA))
    cols.append(fits.Column(name='EWHA_std', format=str(len(rbinc))+'E', array=EWHA_std))
    cols.append(fits.Column(name='SFR', format=str(len(rbinc))+'E', array=SFR))
    cols.append(fits.Column(name='SFR_std', format=str(len(rbinc))+'E', array=SFR_std))
    
    cols.append(fits.Column(name='sSFR', format=str(len(rbinc))+'E', array=SSFR))
    cols.append(fits.Column(name='sSFR_std', format=str(len(rbinc))+'E', array=SSFR_std))
    cols.append(fits.Column(name='n_spaxels', format=str(len(rbinc))+'E', array=NSS))

    cols.append(fits.Column(name='Mstar', format=str(len(rbinc))+'E', array=MSTAR))
    cols.append(fits.Column(name='Mstar_std', format=str(len(rbinc))+'E', array=MSTAR_std))

    cols.append(fits.Column(name='DOP16', format=str(len(rbinc))+'E', array=DOP16))
    cols.append(fits.Column(name='DOP16_std', format=str(len(rbinc))+'E', array=DOP16_std))
    cols.append(fits.Column(name='M08', format=str(len(rbinc))+'E', array=M08))
    cols.append(fits.Column(name='M08_std', format=str(len(rbinc))+'E', array=M08_std))
    cols.append(fits.Column(name='PP04', format=str(len(rbinc))+'E', array=PP04))
    cols.append(fits.Column(name='PP04_std', format=str(len(rbinc))+'E', array=PP04_std))
    
    cols = fits.ColDefs(cols)
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(outpath +'/gradbin_'+ap[j]+'.fits', overwrite=True)
    