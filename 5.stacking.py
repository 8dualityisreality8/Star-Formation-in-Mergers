#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Take the ssfr/metallicity gradients from the 3.gradbin files and "stack" them
within 0.5 dex stellar mass bins. Here, stacking is taking the median profile
values of all objects in the sample within separate stellar mass bins.

"""
from astropy.io import fits
import numpy as np
import os, sys
import functions as f

def N(a):
    n = []
    if len(a) > 0:
        for i in range(len(a[0])):
            n.append(len(a[:,i][~np.isnan(a[:,i])&np.isfinite(a[:,i])]))
        return np.asarray(n)
    else:
        return np.zeros(14)
#______________________________________________________________________________
# Paths and Files
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))

#T5path = f.T5path()
T5path = '/Volumes/Ext_drive/photoObj8'

simpath = T5path + '/6.master/simard.fits'
sim = fits.getdata(simpath)

samplepath = T5path + '/8.samples/ncomp.fits'
sam = fits.getdata(samplepath)

drpallpath = T5path + '/MPL-8/drpall_uniq.fits'
drpall = fits.getdata(drpallpath)
cmd = drpall['CMD']
m = np.log10(drpall['nsa_elpetro_mass']/0.7**2)

drpallmasterpath = T5path + '/6.master/drpall.fits'
drpall_master = fits.getdata(drpallmasterpath)

outpath = filepath + '/5.stacking'
if not os.path.exists(outpath):
    os.makedirs(outpath)
#______________________________________________________________________________
gradbin = fits.getdata(filepath + '/3.gradbin/gradbin_nsa.fits')
rbins = gradbin['reff'][0]

nsa_d16 = gradbin['DOP16']
nsa_m08 = gradbin['M08']
nsa_ewha = gradbin['ewha']
nsa_ssfr = gradbin['ssfr']
nsa_mstar = gradbin['mstar']
nsa_sfr = gradbin['sfr']

gradbin = fits.getdata(filepath + '/3.gradbin/gradbin_sim.fits')
sim_d16 = gradbin['DOP16']
sim_m08 = gradbin['M08']
sim_ewha = gradbin['ewha']
sim_ssfr = gradbin['ssfr']
sim_mstar = gradbin['mstar']
sim_sfr = gradbin['sfr']
#______________________________________________________________________________
# Binning parameters
# Loop over both apertures

CONT = []
PAIR = []
MASS = []
SSFR = []
EWHA = []
D16 =  []
M08 =  []
    
ap = ['nsa', 'sim']
for l in range(2):
    sam = fits.getdata(filepath + '/4.samples.fits')
    controls = sam[ap[l]+'_control'].astype('bool')
    pairs = sam[ap[l]+'_pair'].astype('bool')
    
    d16 = nsa_d16
    m08 = nsa_m08
    ewha = nsa_ewha
    ssfr = nsa_ssfr
    mstar = nsa_mstar
    sfr = nsa_sfr
    if ap[l]=='nsa':
        mass = m
        
    elif ap[l] == 'sim':
        mass = sam['corrected_mass']
        """
        Use NSA where there is no simard aperture
        """
        sim_mendel = sam['sim_mendel'].astype('bool')
        
        d16[sim_mendel] = sim_d16[sim_mendel]
        m08[sim_mendel] = sim_m08[sim_mendel]
        ewha[sim_mendel] = sim_ewha[sim_mendel]
        ssfr[sim_mendel] = sim_ssfr[sim_mendel]
        mstar[sim_mendel] = sim_mstar[sim_mendel]
        sfr[sim_mendel] = sim_sfr[sim_mendel]
        
    binsize=0.5
    mbins = np.arange(9.0, 11.5, binsize) + binsize/2.0
#______________________________________________________________________________
    CONT.append(controls)
    PAIR.append(pairs)
    MASS.append(mass)
    SSFR.append(ssfr)
    EWHA.append(ewha)
    D16.append(d16)
    M08.append(m08)
    
    # Stacked Profiles
    data =  [ewha, ssfr, mstar, sfr, d16, m08]
    strng = ['ewha', 'ssfr', 'mstar', 'sfr', 'DOP16', 'M08']
    
    lim = [0, 0.5]
    filnames = ['stacked', 'stacked_trim']
    
    for k in range(2):
        cols = []
        cols.append(fits.Column(name='rbins', format='1E', array=rbins))
        cols.append(fits.Column(name='mbins', format='1E', array=mbins))
        for i in range(len(data)):
            pbin, pstd = f.stack(data[i][pairs], z=mass[pairs], zbin=mbins, error='sem', trim = lim[k])
            cbin, cstd = f.stack(data[i][controls], z=mass[controls], zbin=mbins, error='sem', trim = lim[k])
                        
            #offsets
            dexbin = pbin - cbin
            dexstd = np.sqrt(pstd**2 + cstd**2)
            
            npairs = []
            ncontrols = []
            for j in range(len(mbins)):
                mix = abs(mbins[j] - mass) < (mbins[1] - mbins[0])/2
                npairs.append(N(data[i][pairs&mix]))
                ncontrols.append(N(data[i][controls&mix]))   
            npairs = np.asarray(npairs)
            ncontrols = np.asarray(ncontrols)
                
            # add to FITS file
            cols.append(fits.Column(name=strng[i]+'_pbin', format=str(pbin.shape[0])+'E', array=pbin.T))
            cols.append(fits.Column(name=strng[i]+'_pstd', format=str(pbin.shape[0])+'E', array=pstd.T))
            cols.append(fits.Column(name=strng[i]+'_cbin', format=str(pbin.shape[0])+'E', array=cbin.T))
            cols.append(fits.Column(name=strng[i]+'_cstd', format=str(pbin.shape[0])+'E', array=cstd.T))
            cols.append(fits.Column(name=strng[i]+'_npair', format=str(pbin.shape[0])+'E', array=npairs.T))
            cols.append(fits.Column(name=strng[i]+'_ncntrl', format=str(pbin.shape[0])+'E', array=ncontrols.T))
            cols.append(fits.Column(name=strng[i]+'_dex', format=str(pbin.shape[0])+'E', array=dexbin.T))
            cols.append(fits.Column(name=strng[i]+'_dexstd', format=str(pbin.shape[0])+'E', array=dexstd.T))
            
        cols = fits.ColDefs(cols)
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.writeto(outpath+'/'+filnames[k]+'_'+ap[l]+'.fits', overwrite=True)
#______________________________________________________________________________
"""
Stack the combined sample
"""

MASS = np.asarray(MASS)
SSFR = np.asarray(SSFR)
CONT = np.asarray(CONT)
PAIR = np.asarray(PAIR)
EWHA = np.asarray(EWHA)
D16 =  np.asarray(D16)
M08 =  np.asarray(M08)

lim = [0, 0.75]
filenames = ['stacked_comb', 'stacked_trim_comb']

for k in range(2):

    pbin, pstd = f.stack(SSFR[PAIR], z=MASS[PAIR], zbin=mbins, error='sem', trim = lim[k])
    cbin, cstd = f.stack(SSFR[CONT], z=MASS[CONT], zbin=mbins, error='sem', trim = lim[k])
    
    epbin, epstd = f.stack(EWHA[PAIR], z=MASS[PAIR], zbin=mbins, error='sem', trim = lim[k])
    ecbin, ecstd = f.stack(EWHA[CONT], z=MASS[CONT], zbin=mbins, error='sem', trim = lim[k])
    
    dpbin, dpstd = f.stack(D16[PAIR], z=MASS[PAIR], zbin=mbins, error='sem', trim = lim[k])
    dcbin, dcstd = f.stack(D16[CONT], z=MASS[CONT], zbin=mbins, error='sem', trim = lim[k])
    
    mpbin, mpstd = f.stack(M08[PAIR], z=MASS[PAIR], zbin=mbins, error='sem', trim = lim[k])
    mcbin, mcstd = f.stack(M08[CONT], z=MASS[CONT], zbin=mbins, error='sem', trim = lim[k])
    
    cols = []
    cols.append(fits.Column(name='rbins', format='1E', array=rbins))
    cols.append(fits.Column(name='mbins', format='1E', array=mbins))
    cols.append(fits.Column(name='ssfr_pbin', format=str(pbin.shape[0])+'E', array=pbin.T))
    cols.append(fits.Column(name='ssfr_pstd', format=str(pbin.shape[0])+'E', array=pstd.T))
    cols.append(fits.Column(name='ssfr_cbin', format=str(pbin.shape[0])+'E', array=cbin.T))
    cols.append(fits.Column(name='ssfr_cstd', format=str(pbin.shape[0])+'E', array=cstd.T))
    
    cols.append(fits.Column(name='ewha_pbin', format=str(pbin.shape[0])+'E', array=epbin.T))
    cols.append(fits.Column(name='ewha_pstd', format=str(pbin.shape[0])+'E', array=epstd.T))
    cols.append(fits.Column(name='ewha_cbin', format=str(pbin.shape[0])+'E', array=ecbin.T))
    cols.append(fits.Column(name='ewha_cstd', format=str(pbin.shape[0])+'E', array=ecstd.T))
    
    cols.append(fits.Column(name='DOP16_pbin', format=str(pbin.shape[0])+'E', array=dpbin.T))
    cols.append(fits.Column(name='DOP16_pstd', format=str(pbin.shape[0])+'E', array=dpstd.T))
    cols.append(fits.Column(name='DOP16_cbin', format=str(pbin.shape[0])+'E', array=dcbin.T))
    cols.append(fits.Column(name='DOP16_cstd', format=str(pbin.shape[0])+'E', array=dcstd.T))
    
    cols.append(fits.Column(name='M08_pbin', format=str(pbin.shape[0])+'E', array=mpbin.T))
    cols.append(fits.Column(name='M08_pstd', format=str(pbin.shape[0])+'E', array=mpstd.T))
    cols.append(fits.Column(name='M08_cbin', format=str(pbin.shape[0])+'E', array=mcbin.T))
    cols.append(fits.Column(name='M08_cstd', format=str(pbin.shape[0])+'E', array=mcstd.T))
    
    cols = fits.ColDefs(cols)
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(outpath+'/'+filenames[k]+'.fits', overwrite=True)
