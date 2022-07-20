#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Define inside and outside pair samples and define the control sample of 
isolated galaxies.

1 == in sample
0 == not in sample

@author: joshua
"""
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import os, sys
import functions as f

# set filepath
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))

#T5path = f.T5path()
T5path = '/Volumes/Ext_drive/photoObj8'

# Open drpall
drpallpath = T5path + '/MPL-8/drpall_uniq_all.fits'
drpall_uniq = fits.getdata(drpallpath)
drmass_uniq = np.log10(drpall_uniq['nsa_elpetro_mass']/0.7**2)
drpallpath = T5path + '/MPL-8/drpall_uniq.fits'
drpall = fits.getdata(drpallpath)
cmd = drpall['CMD']
mngtarg1 = drpall['mngtarg1']
mngtarg3 = drpall['mngtarg3']
drmass = np.log10(drpall['nsa_elpetro_mass']/0.7**2)
z = drpall['nsa_z']
rad = drpall['nsa_elpetro_th50_r']
# Open ncomp
samplepath = T5path + '/8.samples/ncomp.fits'
sam = fits.getdata(samplepath)

drpallmasterpath = T5path + '/6.master/drpall.fits'
drpall_master = fits.getdata(drpallmasterpath)

simpath = T5path + '/6.master/simard.fits'
sim = fits.getdata(simpath)

bptpath = '/Volumes/Ext_drive/Research/data/bpt_r1.0r50_pars.fits'
ewha = fits.getdata(bptpath)['ew']
bptpath = '/Volumes/Ext_drive/Research/data/bpt_r1.3kpc_pars.fits'
bpt = fits.getdata(bptpath)['bpt']

mendel = fits.getdata('/Volumes/Ext_drive/Research/data/mendel.fits')
#______________________________________________________________________________
# define outside control sample and outside pair sampele

mup = 11.5

# Select Primary and Secondary galaxies
gal = (((mngtarg1 & 2L**10 + 2L**11) != 0) & ((mngtarg1 != 0)))


alls = (cmd=='BC')&(bpt<3)&(ewha>=6)
alls2 = gal&(drmass >= 9.0)&(drmass <= mup)

outside_pairs = (((sam['nout_manga']>0)|(sam['nout_nsa']>0))&(sam['donotuse']!=1)&(sam['ptype']==1)&(sam['ngood']==1))
line = 'Outside-IFU pairs: ' + str(len(np.where(outside_pairs==True)[0])) + '\n'
outside_pairs = outside_pairs&alls2
line += 'Outside-IFU Subsample+Mass pairs: ' + str(len(np.where(outside_pairs==True)[0])) + '\n'
outside_pairs = outside_pairs&alls
line += 'Outside-IFU SFG+Subsample+Mass pairs: ' + str(len(np.where(outside_pairs==True)[0])) + '\n\n'


controls = ((sam['ctrl']==1))
line += 'Outside-IFU controls: ' + str(len(np.where(controls==True)[0])) + '\n'
controls = controls & alls2
line += 'Outside-IFU Subsample+Mass controls: ' + str(len(np.where(controls==True)[0])) + '\n'
controls = controls & alls
line += 'Outside-IFU SFG+Subsample+Mass controls: ' + str(len(np.where(controls==True)[0])) + '\n\n'


outside_pairs = (((sam['nout_manga']>0)|(sam['nout_nsa']>0))&(sam['donotuse']!=1)&(sam['ptype']==1)&(sam['ngood']==1)&alls&alls2)
controls = ((sam['ctrl']==1)&alls&alls2)
#______________________________________________________________________________
nsapath = T5path + '/7.neighbor/nsa_uniq.fits'
nsa_pair = fits.getdata(nsapath)

nsapath = '/Volumes/Ext_drive/Research/data/nsa_v1_0_1.fits'
nsa = fits.getdata(nsapath)
nsa_mass = np.log10(nsa['elpetro_mass']/0.7**2)

mangapath = T5path + '/7.neighbor/manga_uniq.fits'
manga = fits.getdata(mangapath)

Dm = []
Dv = []
Sep = []
Sepkpc = []

plateifu = nsa_pair['plateifu']
for i in range(len(sam)):
    ind = np.where(drpall['plateifu']==plateifu[i])[0][0]
    tra = drpall['objra'][ind]
    tdec = drpall['objdec'][ind]
    c1 = SkyCoord(tra* u.deg, tdec* u.deg, frame='icrs')
    # outside pair companion positions
    dm = []
    sep = []
    sepkpc = []
    dv = []
    
    # MaNGA companions
    nmanga = sam['nout_manga'][sam['plateifu']==plateifu[i]]
    if nmanga > 0:
        for j in range(nmanga):
            pifus = manga['pifus'][(manga['plateifu']==plateifu[i]), j]
            dv.append(f.relav(drpall['nsa_z'][ind], drpall_uniq['nsa_z'][drpall_uniq['plateifu']==pifus][0]))
            
            dm.append(drmass[ind] - drmass_uniq[drpall_uniq['plateifu']==pifus][0])
            
            c2 = SkyCoord(drpall_uniq['objra'][drpall_uniq['plateifu']==pifus][0] * u.deg, drpall_uniq['objdec'][drpall_uniq['plateifu']==pifus][0] * u.deg, frame='icrs')
            sep.append(c1.separation(c2).arcsecond)
            sepkpc.append(c1.separation(c2).arcsecond * f.cosmo(z[ind]))
            
    # NSA companions
    nnsa = sam['nout_nsa'][sam['plateifu']==plateifu[i]]
    if nnsa > 0:
        for j in range(nnsa):
            idx = nsa_pair['idx'][(nsa_pair['plateifu']==plateifu[i]), j][0]
            dv.append(f.relav(drpall['nsa_z'][ind], nsa['Z'][idx]))
            dm.append(drmass[ind] - nsa_mass[idx])
            
            c2 = SkyCoord(nsa['ra'][idx] * u.deg, nsa['dec'][idx] * u.deg, frame='icrs')
            sep.append(c1.separation(c2).arcsecond)
            sepkpc.append(c1.separation(c2).arcsecond * f.cosmo(z[ind]))

    # If there are no companions, fill with nans
    if len(dm) == 0:
        Dv.append(np.nan)
        Dm.append(np.nan)
        Sep.append(np.nan)
        Sepkpc.append(np.nan)
    else:
        # keep the closest pair
        ix = np.where(np.asarray(sepkpc) == np.asarray(sepkpc).min())[0][0]
        Dm.append(dm[ix])
        Dv.append(dv[ix])
        Sep.append(sep[ix])
        Sepkpc.append(sepkpc[ix])
#______________________________________________________________________________
# Define inside pair and control sample and corrected mass
nsa_mag = 22.5 - 2.5 * np.log10(drpall_master['nsa_sersic_flux'][:,4])

sim_mag = sim['RG2D']

msim = mendel['m14_logm']
plateifu = drpall_master['plateifu']
nosim = (sim_mag==0)

gradbin = fits.getdata(filepath + '/3.gradbin/gradbin_sim.fits')

spfitpars = fits.getdata('/Volumes/Ext_drive/photoObj8/6.master/spftpars.fits')
spclass = fits.getdata('/Volumes/Ext_drive/photoObj8/6.master/spclass.fits')

# inside pair sample.
cor = []
nosimm = []
mag = []
for i in range(len(drpall)):
    ix = np.where(drpall_master['plateifu']==drpall['plateifu'][i])[0][0]
    cor.append(msim[ix])
    nosimm.append(nosim[ix])
    mag.append(sim_mag[ix])
cor = np.asarray(cor)
nosimm = np.asarray(nosimm)
mag = np.asarray(mag)

mass = cor
sim_mendel = ~nosimm&(cor!=0)&(cor <= mup)&(cor >= 9.0)


m = drmass
m[sim_mendel] = mass[sim_mendel]
alls2 = gal&(m>=9.0)&(m<=mup)

inside_pairs = ((sam['donotuse']!=1)&(sam['ptype']==1)&(sam['ngood']>=2))
line += 'Inside-IFU pairs: ' + str(len(np.where(inside_pairs==True)[0])) + '\n'
inside_pairs = inside_pairs&alls2
line += 'Inside-IFU Subsample+Mass pairs: ' + str(len(np.where(inside_pairs==True)[0])) + '\n'
inside_pairs = inside_pairs&alls
line += 'Inside-IFU SFG+Subsample+Mass pairs: ' + str(len(np.where(inside_pairs==True)[0])) + '\n\n'


inside_control = ((sam['ctrl']==1))
line += 'Inside-IFU controls: ' + str(len(np.where(inside_control==True)[0])) + '\n'
inside_control = inside_control&alls2
line += 'Inside-IFU Subsample+Mass controls: ' + str(len(np.where(inside_control==True)[0])) + '\n'
inside_control = inside_control&alls
line += 'Inside-IFU SFG+Subsample+Mass controls: ' + str(len(np.where(inside_control==True)[0])) + '\n\n'


# save output txt
files = open("4.sample.txt","w+")
files.write(line)
files.close() 

inside_pairs = ((sam['donotuse']!=1)&(sam['ptype']==1)&(sam['ngood']>=2)&alls&alls2)
inside_control = ((sam['ctrl']==1)&alls&alls2)

"""
If no Simard ellipse, use NSA ellipse
"""
mass[~sim_mendel] = drmass[~sim_mendel]

in_separc = []
in_sepkpc = []
in_dm = []

Tmstar=[]
Cmstar=[]
for i in range(len(drpall)):
    ix = np.where(plateifu == drpall['plateifu'][i])[0]
    tra = drpall['objra'][i]
    tdec = drpall['objdec'][i]
    c1 = SkyCoord(tra* u.deg, tdec* u.deg, frame='icrs')
    
    if len(ix) > 1:
        
        seps = []
        sepkpcs = []
        for j in range(len(ix) - 1):
            j+=1
            ra = spfitpars['ra'][ix[j]]
            dec = spfitpars['dec'][ix[j]]
            
            c2 = SkyCoord(ra* u.deg, dec* u.deg, frame='icrs')
            seps.append(c1.separation(c2).arcsecond)
            sepkpcs.append(c1.separation(c2).arcsecond * f.cosmo(z[i]))
    
        seps = np.asarray(seps)
        closest = np.where(seps == seps.min())[0][0]
        
        in_separc.append(seps[closest])
        in_sepkpc.append(sepkpcs[closest])
        
        tmstar = np.log10(np.sum(spfitpars['m_star'][ix[0], 0]))
        cmstar = np.log10(np.sum(spfitpars['m_star'][ix[closest + 1], 0]))

        Tmstar.append(tmstar)
        Cmstar.append(cmstar)
        
        in_dm.append(tmstar - cmstar)
    else:
        in_separc.append(0.0)
        in_sepkpc.append(0.0)
        in_dm.append(0.0)
        Tmstar.append(0.0)
        Cmstar.append(0.0)

all_pairs = (outside_pairs | inside_pairs)
all_controls = controls

#arc
all_sep_arc = np.zeros(all_pairs.shape)*np.nan
all_sep_arc[outside_pairs] = np.array(Sep)[outside_pairs]
all_sep_arc[inside_pairs] = np.array(in_separc)[inside_pairs]

#kpc
all_sep = np.zeros(all_pairs.shape)*np.nan
all_sep[outside_pairs] = np.array(Sepkpc)[outside_pairs]
all_sep[inside_pairs] = np.array(in_sepkpc)[inside_pairs]

all_dm = np.zeros(all_pairs.shape)*np.nan
all_dm[outside_pairs] = np.array(Dm)[outside_pairs]
all_dm[inside_pairs] = np.array(in_dm)[inside_pairs]

all_mass = drmass
all_mass[sim_mendel&inside_pairs] = np.array(mass)[sim_mendel&inside_pairs]

#______________________________________________________________________________    
cols = []
# add to FITS file
cols.append(fits.Column(name='nsa_control', format='1I', array=controls.astype(int)))
cols.append(fits.Column(name='nsa_pair', format='1I', array=outside_pairs.astype(int)))
cols.append(fits.Column(name='sim_control', format='1I', array=inside_control.astype(int)))
cols.append(fits.Column(name='sim_pair', format='1I', array=inside_pairs.astype(int)))
cols.append(fits.Column(name='sim_mendel', format='1I', array=sim_mendel.astype(int)))
cols.append(fits.Column(name='corrected_mass', format='1E', array=mass))

cols.append(fits.Column(name='nsa_dm', format='1E', array=np.asarray(Dm)))
cols.append(fits.Column(name='nsa_dv', format='1E', array=np.asarray(Dv)))
cols.append(fits.Column(name='nsa_separc', format='1E', array=np.asarray(Sep)))
cols.append(fits.Column(name='nsa_sepkpc', format='1E', array=np.asarray(Sepkpc)))

cols.append(fits.Column(name='sim_dm', format='1E', array=np.asarray(in_dm)))
#cols.append(fits.Column(name='sim_dv', format='1E', array=np.asarray(dv_sim)))
cols.append(fits.Column(name='sim_separc', format='1E', array=np.asarray(in_separc)))
cols.append(fits.Column(name='sim_sepkpc', format='1E', array=np.asarray(in_sepkpc)))

cols.append(fits.Column(name='all_control', format='1I', array=all_controls.astype(int)))
cols.append(fits.Column(name='all_pair', format='1I', array=all_pairs.astype(int)))
cols.append(fits.Column(name='all_dm', format='1E', array=all_dm))
cols.append(fits.Column(name='all_separc', format='1E', array=all_sep_arc))
cols.append(fits.Column(name='all_sepkpc', format='1E', array=all_sep))
cols.append(fits.Column(name='all_mass', format='1E', array=all_mass))

cols.append(fits.Column(name='sim_RG2D', format='1E', array=mag))

cols.append(fits.Column(name='tmstar', format='1E', array=np.asarray(Tmstar)))
cols.append(fits.Column(name='cmstar', format='1E', array=np.asarray(Cmstar)))

cols = fits.ColDefs(cols)
tbhdu = fits.BinTableHDU.from_columns(cols)
tbhdu.writeto(filepath+'/4.samples.fits', overwrite=True)