#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 13:20:54 2019

Compute metallicity, EW(Ha), and luminosity gradients in MaNGA pairs where 
there is an IFU over both pairs.

@author: joshua
"""
from astropy.io import fits
import numpy as np
from astropy import wcs
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
import os, sys
import functions as f
import metallicity as metal

T5path = '/Volumes/Ext_drive/photoObj8/'

cosmo = FlatLambdaCDM(H0=70*u.km/u.s/u.Mpc, Om0=0.3)

filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
datapath = os.path.abspath(os.path.join(filepath, '..', 'data'))

drpallpath = T5path + '/6.master/drpall.fits'
drpall = fits.getdata(drpallpath)

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

#plateifu = ['7443-12703']
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
    
    aha = aon[:,(linename == 'Ha6565')][:,0]
    ahb = aon[:,(linename == 'Hb4863')][:,0]
    an2 = aon[:,(linename == 'NII6585')][:,0]
    as2_6718 = aon[:,(linename == 'SII6718')][:,0]
    as2_6733 = aon[:,(linename == 'SII6733')][:,0]
    ao3 = aon[:,(linename == 'OIII5008')][:,0]
    ao2 = aon[:,(linename=='OII3730')][:,0]
        
    snr = ((aha >= 3)&(ahb >= 3)&(ao3 >= 3)&(an2 >= 3)).astype('int')
    hasnr = (aha >= 3).astype('int')
    
    # HA EW 
    ewha = spfit_data['EW'][:,(linename == 'Ha6565')][:,0]
#______________________________________________________________________________
    #Calculate the Ha Luminosity
    ebv = f.EBV(ha, hb)
    #ha_corr = f.ccm_unred([6565.]*len(ha), ha, ebv = ebv)
    r = cosmo.luminosity_distance(nsa_z).value * 10**8 * (3.09*10**16)
    lha = np.log10(ha * 10**-17 * (4 * np.pi * r**2))
    lha = lha.astype('float')
    
    bpt = f.bpt(o3, hb, n2, ha, unred=True)
#______________________________________________________________________________
    # Metallicity Calibrators
    dop16, ddop16 = metal.DOP16(n2, s2_6718, s2_6733, ha, hb=hb, dn2=dn2, ds2_6718=ds2_6718, ds2_6733=ds2_6733, dha=dha, dhb=dhb)
    m13o3n2, dm13o3n2 = metal.O3N2M13(o3, n2, ha, hb, do3=do3, dn2=dn2, dha=dha, dhb=dhb)
    m13n2, dm13n2 = metal.N2M13(n2, ha, hb=hb, dn2=dn2, dha=dha)
    m08, dm08 = metal.M08(ratio='R23', OII_3727=o2, OIII_5007=o3, NII_6584=n2, 
        NeIII_3869=None, Ha=ha, Hb=hb, unred=True)
    
    n06, dn06 = metal.N06(ratio='R23', OII_3727=o2, OIII_5007=o3, NII_6584=n2, 
                          SII_6720=s2_6718, NeIII_3869=None, Ha=ha, Hb=hb, unred=True)
    pp04, dpp04 = metal.PP04(o3, n2, ha, hb, do3=do3, dn2=dn2, dha=dha, dhb=dhb) 
    
    o3hb = np.log10(np.divide(o3, hb))
    n2ha = np.log10(np.divide(n2, ha))
    
    no = metal.NO(NII_6548=n2*0.34, NII_6584=n2, OIII_5007=o3, OIII_4959=o3*0.35, OII_3727=o2*2.0, Hb=hb)
#______________________________________________________________________________
    # Get the pixel coordinates of the spaxels and thier deprojected
    # effective radius
    crd = np.array([[tra, tdec]], np.float_)       
    
    # pixel coordinate of object center in SPFIT
    x0, y0 = w.wcs_world2pix(crd, 1)[0]

    # Deprojected effective radii using Simard+11 apertures
    # We use q=0.13 from Giovanelli+94
    sim_reff = []
    reffmap = np.zeros(binmap.shape)
    for j in range(bins.max() + 1):
        y, x = np.where(binmap == j)
        sim_reff.append(f.dist_ellipse(x+1, y+1, x0, y0, sim_ba, sim_phi, q=0.13)[0])
        reffmap[y,x] = f.dist_ellipse(x+1, y+1, x0, y0, sim_ba, sim_phi, q=0.13)[0]
    sim_reff = np.asarray(sim_reff)/sim_rad

    # Deprojected effective radii using NSA Elpetro apertures
    nsa_reff = []
    #reffmap = np.zeros(binmap.shape)
    for j in range(bins.max() + 1):
        y, x = np.where(binmap == j)
        nsa_reff.append(f.dist_ellipse(x+1, y+1, x0, y0, nsa_ba, nsa_phi, q=0.13)[0])
        #reffmap[y,x] = f.dist_ellipse(x+1, y+1, x0, y0, ba, phi, q=0.13)[0]
    nsa_reff = np.asarray(nsa_reff)/nsa_rad
#______________________________________________________________________________ 
    
    cols = []
    cols.append(fits.Column(name='tra', format='1D', array=np.array([tra])))
    cols.append(fits.Column(name='tdec', format='1D', array=[tdec]))
    cols.append(fits.Column(name='ifura', format='1D', array=[ifura]))
    cols.append(fits.Column(name='ifudec', format='1D', array=[ifudec]))        
        
    cols.append(fits.Column(name='rad_sim', format='%.fE' % (len(sim_reff)), array=np.expand_dims(sim_reff, 1).T))
    cols.append(fits.Column(name='rad_nsa', format='%.fE' % (len(nsa_reff)), array=np.expand_dims(nsa_reff, 1).T))
    
    cols.append(fits.Column(name='mstar', format='%.fE' % (len(mstar)), array=np.expand_dims(mstar, 1).T))
    cols.append(fits.Column(name='EBV', format='%.fE' % (len(ebv)), array=np.expand_dims(ebv, 1).T))
    cols.append(fits.Column(name='DOP16', format='%.fE' % (len(dop16)), array=np.expand_dims(dop16, 1).T))
    cols.append(fits.Column(name='M13O3N2', format='%.fE' % (len(m13o3n2)), array=np.expand_dims(m13o3n2, 1).T))
    cols.append(fits.Column(name='M13N2', format='%.fE' % (len(m13n2)), array=np.expand_dims(m13n2, 1).T))
    cols.append(fits.Column(name='M08', format='%.fE' % (len(m08)), array=np.expand_dims(m08, 1).T))
    cols.append(fits.Column(name='N06', format='%.fE' % (len(n06)), array=np.expand_dims(n06, 1).T))
    cols.append(fits.Column(name='PP04', format='%.fE' % (len(pp04)), array=np.expand_dims(pp04, 1).T))
    cols.append(fits.Column(name='NO', format='%.fE' % (len(no)), array=np.expand_dims(no, 1).T))
    cols.append(fits.Column(name='dDOP16', format='%.fE' % (len(ddop16)), array=np.expand_dims(ddop16, 1).T))
    cols.append(fits.Column(name='dM13O3N2', format='%.fE' % (len(dm13o3n2)), array=np.expand_dims(dm13o3n2, 1).T))
    cols.append(fits.Column(name='dM13N2', format='%.fE' % (len(dm13n2)), array=np.expand_dims(dm13n2, 1).T))
    cols.append(fits.Column(name='dM08', format='%.fE' % (len(dm08)), array=np.expand_dims(dm08, 1).T))
    cols.append(fits.Column(name='dN06', format='%.fE' % (len(dn06)), array=np.expand_dims(dn06, 1).T))
    cols.append(fits.Column(name='dPP04', format='%.fE' % (len(dpp04)), array=np.expand_dims(dpp04, 1).T))
    cols.append(fits.Column(name='BPT', format='%.fE' % (len(bpt)), array=np.expand_dims(bpt, 1).T))
    cols.append(fits.Column(name='EwHa', format='%.fE' % (len(ewha)), array=np.expand_dims(ewha, 1).T))
    cols.append(fits.Column(name='LogLHa', format='%.fD' % (len(lha)), array=np.expand_dims(lha, 1).T))
    cols.append(fits.Column(name='SNR', format='%.fI' % (len(snr)), array=np.expand_dims(snr, 1).T))
    cols.append(fits.Column(name='SNR_Ha', format='%.fI' % (len(hasnr)), array=np.expand_dims(hasnr, 1).T))
    
    cols.append(fits.Column(name='sim_RCHL_R', format='1D', array=[sim_rad]))
    cols.append(fits.Column(name='sim_ba', format='1D', array=[sim_ba]))
    cols.append(fits.Column(name='sim_phi', format='1D', array=[[sim_phi]]))
    cols.append(fits.Column(name='sim_z', format='1D', array=[sim_z]))
    
    cols.append(fits.Column(name='nsa_elpetro_th50_r', format='1D', array=[nsa_rad]))
    cols.append(fits.Column(name='nsa_elpetro_ba', format='1D', array=[nsa_ba]))
    cols.append(fits.Column(name='nsa_elpetro_phi', format='1D', array=[[nsa_phi]]))
    cols.append(fits.Column(name='nsa_z', format='1D', array=[nsa_z]))
    
    
    cols.append(fits.Column(name='ifudesignsize', format='1K', array=[ifu]))
    cols.append(fits.Column(name='binmap', format='%.fI' % (len(bins)), array=np.expand_dims(bins,1).T))
    
    cols = fits.ColDefs(cols)
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(outpath+'/'+plateifu[i]+'.fits', overwrite=True)
        
    f.update_progress((i+1.0)/np.float64(len(plateifu)))
    