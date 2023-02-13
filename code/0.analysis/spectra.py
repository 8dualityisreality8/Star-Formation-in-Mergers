#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 12:22:30 2020

@author: jlsteffen
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import plotstyle as ps

fil = fits.getdata('/Volumes/Ext_drive/Research/data/spfit1arc/7443-12703.fits')

lam = 10**fil['log10lam'][0]

gal = fil['galaxy'][0]
best = fil['best'][0]
neat = fil['neat'][0]

ha = fil['lambda'][0][np.where(fil['name'][0] == 'Ha')[0][0]]

fontsize=20
fig, axs = plt.subplots(1, 1, figsize=(32, 8))

axs.plot(lam, gal, 'k')
axs.plot(lam, best, 'r--')

#step = 50
#axs.axvspan(ha-step, ha+step, color='b', alpha=0.25)

axs.axvline(6450, c='b')
axs.axvline(6800, c='b')

#axs.plot(lam, neat, 'r--')
axs.set_yscale('log')
axs.set_xlabel('Wavelength [$\AA$]', fontsize=fontsize)
axs.set_ylabel('Flux [10$^{-17}$ erg/s/cm$^2$/$\AA$]', fontsize=fontsize)

axs.set_xlim(lam.min(), lam.max())

ps.style(axs, fontsize=fontsize)

plt.savefig('spectra.pdf', bbox_inches='tight', overwrite=True)
plt.close('all')

#______________________________________________________________________________
fontsize=20
fig, axs = plt.subplots(1, 1, figsize=(8, 8))

axs.plot(lam, gal, 'k')
axs.plot(lam, best, 'r--')
#axs.plot(lam, neat, 'r--')
axs.set_yscale('log')
axs.set_xlabel('Wavelength [$\AA$]', fontsize=fontsize)
axs.set_ylabel('Flux [10$^{-17}$ erg/s/cm$^2$/$\AA$]', fontsize=fontsize)

step = 7
axs.axvspan(ha-step, ha+step, color='b', alpha=0.25)
#axs.axvline(ha, c='b')

axs.set_xlim(6450, 6800)

ps.style(axs, fontsize=fontsize)

plt.savefig('halpha.pdf', bbox_inches='tight', overwrite=True)
plt.close('all')