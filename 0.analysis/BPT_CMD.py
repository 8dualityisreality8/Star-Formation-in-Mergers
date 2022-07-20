# -*- coding: utf-8 -*-
"""
Created on Thu May 21 15:07:17 2020

@author: Josh
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import functions as f
import plotstyle as ps

def blue(x):
    b = 3.1682
    m = 0.16
    return b - (x+18) * m
def red(x):
    b = 4.7866
    m = 0.04
    return b - (x+18) * m

def kewley01(x):
    return ( 0.61 / (x - 0.47) ) + 1.19
def kauffmann03(x):
    return ( 0.61 / (x - 0.05) ) + 1.30

#def law(x):
#    return (0.359 / (x + 0.032)) + 1.083

#def law2(y):
#    return (-0.479*y**4 - 0.594*y**3 - 0.542*y**2 - 0.056*y - 0.143)

#def seyfert(x):
#    return 1.05*x + 0.45

filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
#figpath = '/Users/joshua/Dropbox/manga/figs'

drpallpath = '/Volumes/Ext_drive/photoObj8/MPL-8/drpall_uniq.fits'
drpall = fits.getdata(drpallpath)
elpetro_mass = np.log10(drpall['nsa_elpetro_mass'])
cmd = drpall['cmd']
mngtarg1 = drpall['mngtarg1']

samplepath = '/Volumes/Ext_drive/photoObj8/8.samples/ncomp.fits'
ncomp = fits.getdata(samplepath)

samplepath = '/Volumes/Ext_drive/Research/metal/4.samples.fits'
sam = fits.getdata(samplepath)
mass = sam['all_mass']

mangapath = '/Volumes/Ext_drive/photoObj8/7.neighbor/manga_uniq.fits'
manga = fits.getdata(mangapath)

mag = drpall['nsa_elpetro_absmag']

# Select Pairs
#pairs = (((sam['nout_manga']>0)|(sam['nout_nsa']>0))&(sam['donotuse']!=1)&(sam['ptype']==1)&(sam['ngood']==1))
#controls = ((sam['ctrl']==1))

bpt = fits.getdata('/Volumes/Ext_drive/Research/data/bpt_r1.3kpc_pars.fits')
o3hb = np.log10(bpt['OIII5008']/bpt['hb'])
n2ha = np.log10(bpt['NII6585']/bpt['ha'])
ewha = bpt['ew']
#BPT = f.bpt(bpt['OIII5008'], bpt['hb'], bpt['NII6585'], bpt['ha'])
BPT = bpt['bpt']

mask = (~np.isnan(n2ha))&(np.isfinite(n2ha))&(~np.isnan(o3hb))&(np.isfinite(o3hb))

o3hb = o3hb[mask]
n2ha = n2ha[mask]
ewha = ewha[mask]
BPT = BPT[mask]
mag = mag[mask]

#pars = (pairs)&(cmd=='BC')&(ewha>=6)

#______________________________________________________________________________
mup = 11.5
gal = (((mngtarg1 & 2L**10 + 2L**11) != 0) & ((mngtarg1 != 0)))
alls = (cmd=='BC')
alls2 = gal&(mass >= 9.0)&(mass <= mup)#&(ewha>=6)

# Select Inside Pairs
ipairs = ((ncomp['donotuse']!=1)&(ncomp['ptype']==1)&(ncomp['ngood']>=2)&alls2)

#ipairs = sam['sim_pair'].astype('bool')

# Select Outside Pairs
opairs = (((ncomp['nout_manga']>0)|(ncomp['nout_nsa']>0))&(ncomp['donotuse']!=1)&(ncomp['ptype']==1)&(ncomp['ngood']==1)&alls2)

#opairs = sam['nsa_pair'].astype('bool')

pairs = ipairs|opairs

pairs = pairs[mask]

sf = BPT==1
comp = BPT==2
agn = BPT>2
ret = BPT<=0

#______________________________________________________________________________


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=False, sharey=False, figsize=(30, 10))
fontsize=30.0

x = -1.25, 0.75
y = -1.45, 1.45

z = ps.gaussian_den(n2ha, o3hb)

im = ax2.scatter(n2ha, o3hb, c=z, s=5, vmax=1.0, vmin=0.0, cmap='gray_r')

ax2.scatter(n2ha[pairs&ret], o3hb[pairs&ret], c='purple', label='Retired', marker='D', s=30, edgecolor='k', linewidth=1)
ax2.scatter(n2ha[pairs&sf], o3hb[pairs&sf], c='deepskyblue', label='Starforming', marker='D', s=30, edgecolor='k', linewidth=1)
ax2.scatter(n2ha[pairs&comp], o3hb[pairs&comp], c='lime', label='Composite', marker='D', s=30, edgecolor='k', linewidth=1)
ax2.scatter(n2ha[pairs&agn], o3hb[pairs&agn], c='crimson', label='AGN', marker='D', s=30, edgecolor='k', linewidth=1)

leg = ax2.legend(fontsize=fontsize/2.0, loc=3, fancybox=False, edgecolor='k')
leg.get_frame().set_linewidth(2.0)

xx = np.arange(x[0], 0.45, 0.01)
yy = kewley01(xx)
ax2.plot(xx, yy, c='k', linewidth=2)
xx = np.arange(x[0], 0.0, 0.01)
yy = kauffmann03(xx)
ax2.plot(xx, yy, c='k', linewidth=2, linestyle='--')

ax2.set_xlim(x)
ax2.set_ylim(y)

ax2.set_ylabel(r'log([O III]/H$\beta$)', fontsize=fontsize)
ax2.set_xlabel(r'log([N II]/H$\alpha$)', fontsize=fontsize)

ps.ticks(ax2, xmajor=0.5, ymajor=0.5, xminor=0.1, yminor=0.1)
ps.style(ax2, fontsize=fontsize)

ax2.annotate('SF', (-1, -0.25), fontsize=fontsize/1.5)
ax2.annotate('Comp.', (-0.15, -0.9), fontsize=fontsize/1.5)
ax2.annotate('AGN', (-0.5, 0.75), fontsize=fontsize/1.5)

#______________________________________________________________________________
"""WHAN Diagram"""
z = ps.gaussian_den(n2ha, np.log10(ewha))

im = ax3.scatter(n2ha, np.log10(ewha), c=z, s=5, vmax=1.0, vmin=0.0, cmap='gray_r')

ax3.scatter(n2ha[pairs&ret], np.log10(ewha[pairs&ret]), c='purple', label='Retired', marker='D', s=30, edgecolor='k', linewidth=1)
ax3.scatter(n2ha[pairs&sf], np.log10(ewha[pairs&sf]), c='deepskyblue', label='Starforming', marker='D', s=30, edgecolor='k', linewidth=1)
ax3.scatter(n2ha[pairs&comp], np.log10(ewha[pairs&comp]), c='lime', label='Composite', marker='D', s=30, edgecolor='k', linewidth=1)
ax3.scatter(n2ha[pairs&agn], np.log10(ewha[pairs&agn]), c='crimson', label='AGN', marker='D', s=30, edgecolor='k', linewidth=1)

ax3.axhline(np.log10(3), c='k')
ax3.plot([-0.4, 0.75], [np.log10(6), np.log10(6)], c='k', linestyle=':')
ax3.plot([-0.4, -0.4], [np.log10(3), 2.75], c='k', linestyle='--')

ax3.set_xlim(-1.5, 0.75)
ax3.set_ylim(-0.75, 2.75)

ax3.set_ylabel(r'log(EW(H$\alpha$)', fontsize=fontsize)
ax3.set_xlabel(r'log([N II]/H$\alpha$)', fontsize=fontsize)

ps.ticks(ax3, xmajor=0.5, ymajor=0.5, xminor=0.1, yminor=0.1)
ps.style(ax3, fontsize=fontsize)

ax3.annotate('SF', (-1, 0.75), fontsize=fontsize/1.5)
ax3.annotate('sAGN', (0.2, 1.0), fontsize=fontsize/1.5)
ax3.annotate('wAGN', (0.4, 0.6), fontsize=fontsize/1.5)
ax3.annotate('RG', (-1, -0.25), fontsize=fontsize/1.5)
#______________________________________________________________________________
mag = drpall['nsa_elpetro_absmag']

ewha = bpt['ew']
BPT = bpt['bpt']

#mask = (~np.isnan(n2ha))&(np.isfinite(n2ha))&(~np.isnan(o3hb))&(np.isfinite(o3hb))

ewha = ewha#[mask]
BPT = BPT#[mask]
mag = mag#[mask]

pairs = ipairs|opairs

sf = BPT==1
comp = BPT==2
agn = BPT>2
ret = BPT<=0

z = ps.gaussian_den(mag[:,5], mag[:,1]-mag[:,4])

im = ax1.scatter(mag[:,5], mag[:,1]-mag[:,4], c=z, s=5, vmax=1.0, vmin=0.0, cmap='gray_r')

cbar_ax = fig.add_axes([0.405, 0.88, 0.215, 0.05])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')  
cbar.ax.set_title('Normalized Number Density', fontsize=fontsize, pad=30)
#ax1.set_title('Normalized Number Density', fontsize=fontsize, pad=25)
ps.cbar_style(cbar, fontsize=fontsize)

ax1.scatter(mag[pairs&ret,5], (mag[pairs&ret,1]-mag[pairs&ret,4]), c='purple', s=30, marker='D', edgecolor='k', linewidth=1, label='Retired')
ax1.scatter(mag[pairs&sf,5], (mag[pairs&sf,1]-mag[pairs&sf,4]), c='deepskyblue', s=30, marker='D', edgecolor='k', linewidth=1, label='SF')
ax1.scatter(mag[pairs&comp,5], (mag[pairs&comp,1]-mag[pairs&comp,4]), c='lime', s=30, marker='D', edgecolor='k', linewidth=1, label='Composite')
ax1.scatter(mag[pairs&agn,5], (mag[pairs&agn,1]-mag[pairs&agn,4]), c='crimson', s=30, marker='D', edgecolor='k', linewidth=1, label='AGN')

x = np.arange(-24, -16 + 1, 1)
ax1.plot(x, blue(x), c='k', linestyle='--')
ax1.plot(x, red(x), c='k', linestyle='--')

ax1.set_aspect('auto')

ax1.set_xlim(-16, -24)
ax1.set_ylim(1.05,6.5)

ax1.set_xlabel(r'$\rm M_i$', fontsize=fontsize)
ax1.set_ylabel('NUV-r', fontsize=fontsize)

ps.ticks(ax1, xmajor=1, ymajor=1, xminor=0.2, yminor=0.2)
ps.style(ax1, fontsize=fontsize)

ax1.annotate('Red Sequence', (-16.5, 6), fontsize=fontsize/1.5)
ax1.annotate('Green\nValley', (-16.25, 3.5), fontsize=fontsize/1.5)
ax1.annotate('Blue Cloud', (-21, 1.5), fontsize=fontsize/1.5)

plt.savefig(filepath + '/bpt-cmd.pdf', bbox_inches='tight', overwrite=True)

plt.close('all')