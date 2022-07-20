#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 13:20:54 2019

Compute metallicity gradients in MaNGA pairs where there is an IFU over both
pairs.

@author: joshua
"""
size = 300
from astropy.io import fits
import numpy as np
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
import os, sys
import functions as f
import metallicity as metal
import matplotlib.pyplot as plt
import plotstyle as ps
from PIL import Image
import matplotlib.patches as patches
from matplotlib.patches import Ellipse
from copy import copy
import math as m

T5path = f.T5path()

cosmo = FlatLambdaCDM(H0=70*u.km/u.s/u.Mpc, Om0=0.3)

filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
datapath = os.path.abspath(os.path.join(filepath, '..', '..', 'data'))

drpallpath = T5path + '/6.master/drpall.fits'
drpall = fits.getdata(drpallpath)
ras = drpall['ifura']
decs = drpall['ifudec']
ifus = drpall['ifudesignsize']
plifus = drpall['plateifu']

photosimpath = '/Volumes/Ext_drive/Research/data/simard11/CasJobs/phioffset.fit'
photosim = fits.getdata(photosimpath)
offset = photosim['phioffset_g']

simpath = T5path + '/6.master/simard.fits'
sim = fits.getdata(simpath)

elpetro_mass = np.log10(drpall['nsa_elpetro_mass'])

outpath = filepath + '/1.metalgrad'
if not os.path.exists(outpath):
    os.makedirs(outpath)

plateifu = np.unique(drpall['plateifu'])

plateifu = ['7443-12703']
for i in range(len(plateifu)):
    # Retrieve drpall properties for plateifu
    ind = np.where(drpall['plateifu']==plateifu[i])[0][0]
    ifura = drpall['ifura'][ind]
    ifudec = drpall['ifudec'][ind]
    tra = drpall['objra'][ind]
    tdec = drpall['objdec'][ind]
    nsa_z = drpall['nsa_z'][ind]
    nsa_rad = drpall['nsa_elpetro_th50_r'][ind]
    nsa_ba = drpall['nsa_elpetro_ba'][ind]
    nsa_phi = drpall['nsa_elpetro_phi'][ind]
    ifu = drpall['ifudesignsize'][ind]    
    
    c1 = SkyCoord(tra* u.deg, tdec* u.deg, frame='icrs')
#______________________________________________________________________________     
   # Get ellipse params for Simard+11 with adjustments for position angle offsets
    off = offset[photosim['dr7objid']==sim['objid'][ind]]
    if len(off) > 0:
        if len(off) > 1:
            off = off[0]
        sim_phi = ((360 + 90 - (sim['phi'][ind] + off)) % 180)
        sim_ba = 1 - sim['e'][ind]
        sim_rad = sim['RCHL_R'][ind]/sim['scale'][ind]
    else:
        sim_phi = 0
        sim_ba = 0
        sim_rad = 0
        
    sim_z = sim['z'][ind]
#______________________________________________________________________________
    # SPFIT WCS
    spfit = fits.open(datapath + '/spfit_pars/'+plateifu[i]+'_pars.fits')
    hdr = spfit[2].data['HDR'][0]

    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [float([s for s in hdr[82].split() if s.replace('.','',1).isdigit()][0]), 
                   float([s for s in hdr[83].split() if s.replace('.','',1).isdigit()][0])]
    
    w.wcs.cdelt = [float([s for s in hdr[86].split() if s.replace('.','',1).strip('-').isdigit()][0]), 
                   float([s for s in hdr[87].split() if s.replace('.','',1).strip('-').isdigit()][0])]
    
    w.wcs.crval = [float([s for s in hdr[84].split() if s.replace('.','',1).strip('-').isdigit()][0]), 
                   float([s for s in hdr[85].split() if s.replace('.','',1).strip('-').isdigit()][0])]
    
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ['deg', 'deg']
    w.wcs.set_pv([(2, 1, 45.0)])
#______________________________________________________________________________
    # Extract emission lines and calculate BPT classification and Metallicity
    spfit_data = spfit[1].data
    flux = spfit_data['flux']
    aon = spfit_data['AON']
    binmap = spfit[2].data['binmap'][0]
    bins = binmap.ravel()
    
    mstar = np.log10(np.sum(spfit_data['m_star'][:,0], axis=1))
    
    # Extract Fluxes
    linename = spfit[2].data['linename'][0]
     
    ha = flux[:,0,(linename == 'Ha6565')][:,0]
    hb = flux[:,0,(linename == 'Hb4863')][:,0]
    n2 = flux[:,0,(linename == 'NII6585')][:,0]
    s2_6718 = flux[:,0,(linename == 'SII6718')][:,0]
    s2_6733 = flux[:,0,(linename == 'SII6733')][:,0]
    o3 = flux[:,0,(linename == 'OIII5008')][:,0]
    o2 = flux[:,0,(linename=='OII3730')][:,0]
    
    dha = flux[:,1,(linename == 'Ha6565')][:,0]
    dhb = flux[:,1,(linename == 'Hb4863')][:,0]
    dn2 = flux[:,1,(linename == 'NII6585')][:,0]
    ds2_6718 = flux[:,1,(linename == 'SII6718')][:,0]
    ds2_6733 = flux[:,1,(linename == 'SII6733')][:,0]
    do3 = flux[:,1,(linename == 'OIII5008')][:,0]
    do2 = flux[:,1,(linename=='OII3730')][:,0]

        
    # HA EW 
    ewha = spfit_data['EW'][:,(linename == 'Ha6565')][:,0]
#______________________________________________________________________________
    #Calculate the Ha Luminosity
    ebv = f.EBV(ha, hb)
    #ha_corr = f.ccm_unred([6565.]*len(ha), ha, ebv = ebv)
    r = cosmo.luminosity_distance(nsa_z).value * 10**8 * (3.09*10**16)
    lha = np.log10(ha * 10**-17 * (4 * np.pi * r**2))
    lha = lha.astype('float')
    
    bpt = f.bpt_whan(o3, hb, n2, ha, ewha, unred=False)

    # Get the pixel coordinates of the spaxels and thier deprojected
    # effective radius
    crd = np.array([[tra, tdec]], np.float_)       
    
    # pixel coordinate of object center in SPFIT
    x0, y0 = w.wcs_world2pix(crd, 1)[0]

    #pixkpc = (f.cosmo(grad[ap[j]+'_z']) * 1.0)**2 / np.sqrt( (grad[ba[j]]**2 - 0.13**2) / (1 - 0.13**2) )
    #Pix_kpc[i] = pixkpc
    sfr = f.SFR_M11(10**lha)
    ssfr = np.log10(f.divide(sfr, 10**mstar))
    #SFR[i], SFR_std[i] = f.azimuth_profiles(reff, np.log10(sfr/pixkpc), rbins)
    #ssfr = f.divide(sfr, 10**mstar)

    # Deprojected effective radii using Simard+11 apertures
    ewhamap = np.zeros(binmap.shape)
    for j in range(bins.max() + 1):
        y, x = np.where(binmap == j)
        ewhamap[y,x] = np.log10(sfr[j])
#______________________________________________________________________________
    # Download images from skyserver    
    ind = np.where(drpall['plateifu']==plateifu[i])[0][0]
    ra = ras[ind]
    dec = decs[ind]
    IFU_size = ifus[ind]
    
    IFU_arcsec = np.array(([19, 37, 61, 91, 127], [16.0, 21.0, 26.0, 31.0, 36.0]))
    ifu_scale = IFU_arcsec[1, np.where(IFU_size==IFU_arcsec[0])[0][0]]
    scale = 40.0/size
#______________________________________________________________________________
    # Open the downloaded png for editting
    im = Image.open(datapath +'/png/'+str(plateifu[i])+'.png')
    cent_pix = np.array([im.size[0]/2.0, im.size[1]/2.0])
#______________________________________________________________________________
    # Make a Astropy WCS
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [cent_pix[0], cent_pix[1]]
    w.wcs.cdelt = np.array([-scale/3600.0, -scale/3600.0])
    w.wcs.crval = [ra, dec]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ['deg', 'deg']
    w.wcs.set_pv([(2, 1, 45.0)])
    my_dpi = 96.0
#______________________________________________________________________________ 
    plt.style.use('dark_background')
    fontsize=30
    fig, axs = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(20, 8))
    
    axs[0].imshow(im, aspect='equal')
#______________________________________________________________________________
    # IFU hexagons    
    hexagon = patches.RegularPolygon([size/2.0, size/2.0], 6, 
              radius = size/2.0 * ifu_scale/40.0, ec = 'w', fc = 'none', orientation=m.pi/2, linewidth=2.0)
    axs[0].add_patch(copy(hexagon))

    axs[1].add_patch(copy(hexagon))
    
    # Mask out photo objs
    crd = np.array([[ra, dec]], np.float_)
    cx, cy = w.wcs_world2pix(crd,1)[0]
    hyp = 2.0/scale
    ell = Ellipse(xy=(cx,cy), width=hyp, height=hyp, angle=0.0, ec = 'lime', fc = 'none', linewidth=2.0)
    axs[0].add_patch(copy(ell))
    
#______________________________________________________________________________
    ewhamap[ewhamap==0] = np.nan
    new_len = 300*ifu_scale/40.0
    
    im= axs[1].imshow(ewhamap, origin='lower', cmap='hot', extent=[(300-new_len)/2.,new_len + (300-new_len)/2.,new_len + (300-new_len)/2.,(300-new_len)/2.])
        
    cbar_ax = fig.add_axes([0.585, 0.06, 0.25, 0.05])
    c = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    #c.set_ticks([-0.05, 0.45, 0.95, 1.45, 1.95, 2.45, 2.95])
    #c.ax.set_xticklabels(['Mask', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0'])
    
    ps.cbar_style(c, fontsize=fontsize, labelside='bottom')
    
    
    ps.ticks(axs, xmajor=80.0, ymajor=80.0, xminor=16.0, yminor=16.0)
    ps.style(axs, fontsize, labelbottom=False, labelleft=False)
    
    axs[0].set_title('Mrk 848', fontsize=fontsize)
    axs[1].set_title('log(SFR)', fontsize=fontsize)
    
    plt.subplots_adjust(wspace=0.0)
    plt.savefig('sfrmap.pdf', bbox_inches='tight', overwrite=True)
    plt.close('all')