#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 15:08:46 2019

@author: joshua
"""

import numpy as np
from functions import ccm_unred, EBV, divide, log10

"""Line Ratios"""

def R23(OII_3727, OIII_5007, Hb, OIII_4959=None):
    '''
    if OIII_4959 is None:
        OIII_4959 = 0.35 * OIII_5007
    '''
    O = OII_3727 + OIII_5007
    return divide(O, Hb)
def N2Ha(NII_6584, Ha):
    return divide(NII_6584, Ha)
def O3N2(OIII_5007, NII_6584):
    return divide(OIII_5007, NII_6584)
def N2O2(NII_6584, OII_3727):
    return divide(NII_6584, OII_3727)
def N2S2(NII_6584, SII_6720):
    return divide(NII_6584, SII_6720)
def O3O2(OIII_5007, OII_3727):
    return divide(OIII_5007, OII_3727)
def Ne3O2(NeIII_3869, OII_3727):
    return divide(NeIII_3869, OII_3727)
def HaN2S2(Ha, NII_6584, SII_6720, NII_6548=None):
    if NII_6548 is None:
        NII_6548 = 0.34 * NII_6584
    return divide((Ha+NII_6548+NII_6584), SII_6720)
def O3Hb(OIII_5007, Hb):
    return divide(OIII_5007, Hb)
def O2Hb(OII_3727, Hb):
    return divide(OII_3727, Hb)

"""Metallicity Calibrators"""
def T04(o2, o3, hb, ha=None, do2=None, do3=None, dhb=None):
    # R23 with the Tremonti et al. 2004 calibration
    if ha is not None:
        ebv = EBV(ha, hb)
        o2 = ccm_unred([3730.]*len(o2), o2, ebv = ebv)
        o3 = ccm_unred([5008.]*len(o3), o3, ebv = ebv)
        hb = ccm_unred([4863.]*len(hb), hb, ebv = ebv)
    r = divide((2*o2 + 1.35*o3), hb)
    R23 = log10(r)
    oh = 9.185 - 0.313*R23 - 0.264*R23**2 - 0.321*R23**3
    if do2 is not None and do3 is not None and dhb is not None:
        if ha is not None:
            do2 = ccm_unred([3730.]*len(do2), do2*o2, ebv = ebv)
            do3 = ccm_unred([5008.]*len(do3), do3*o3, ebv = ebv)
            dhb = ccm_unred([4863.]*len(dhb), dhb*hb, ebv = ebv)
        dR23 = np.sqrt(np.abs(divide(do2,(o2*np.log(10))))**2 + np.abs(divide(do3,(o3*np.log(10))))**2 +\
                       np.abs(divide(dhb,(hb*np.log(10))))**2)
        
        #re1 = divide(np.sqrt(2**2 * do2**2 + 1.35**2 * do3**2), (2*o2 + 1.35*o3))
        #re2 = divide(dhb, hb)
        #re = np.sqrt(re1**2 + re2**2) * abs(r)
        #R23e = divide(re, np.log(10)*r)
        doh = 0.313*dR23 + 2*0.264*R23*dR23 + 3*0.321*R23**2 * dR23
    else:
        doh = np.zeros_like(oh)
    return oh, doh
def O3N2M13(o3, n2, ha, hb, do3=None, dn2=None, dha=None, dhb=None):
    # O3N2 with the Marino et al. 2013 calibration
    #o3 = o3 * 1.35 # correct for doublet
    
    ebv = EBV(ha, hb)
    o3 = ccm_unred([5008.]*len(o3), o3, ebv = ebv)
    n2 = ccm_unred([6585.]*len(n2), n2, ebv = ebv)
    ha = ccm_unred([6565.]*len(ha), ha, ebv = ebv)
    hb = ccm_unred([4863.]*len(hb), hb, ebv = ebv)
    
    o3n2 = log10(divide(o3, hb) * divide(ha, n2))
    
    oh = 8.533 - 0.214 * o3n2
    if do3 is not None and dn2 is not None and dha is not None and dhb is not None:
        do3 = ccm_unred([5008.]*len(o3), do3*o3, ebv = ebv)
        dn2 = ccm_unred([6585.]*len(ha), dn2*n2, ebv = ebv)
        dha = ccm_unred([6565.]*len(ha), dha*ha, ebv = ebv)
        dhb = ccm_unred([4863.]*len(hb), dhb*hb, ebv = ebv)

        do3n2 = np.sqrt(np.abs(divide(do3,(o3*np.log(10))))**2 + np.abs(divide(dhb,(hb*np.log(10))))**2 +\
                np.abs(divide(dha,(ha*np.log(10))))**2 + np.abs(divide(dn2,(n2*np.log(10))))**2)
        doh = np.abs(0.214*do3n2)
    else:
        doh = np.zeros_like(oh)
    return oh, doh
def N2M13(n2, ha, hb=None, dn2=None, dha=None):
    # N2 with the Marino et al. 2013 calibraiton
    if hb is not None:
        ebv = EBV(ha, hb)
        n2 = ccm_unred([6585.]*len(n2), n2, ebv = ebv)
        ha = ccm_unred([6565.]*len(ha), ha, ebv = ebv)
    N2 = log10(divide(n2,ha))
    oh = 8.743 + 0.462*N2
    if dn2 is not None and dha is not None:
        dn2 = ccm_unred([6585.]*len(ha), dn2*n2, ebv = ebv)
        dha = ccm_unred([6565.]*len(ha), dha*ha, ebv = ebv)
        
        dN2 = np.sqrt(np.abs(divide(dn2,(n2*np.log(10))))**2 + np.abs(divide(dha,(ha*np.log(10))))**2)
        doh = np.abs(0.462 * dN2)
    else:
        doh = np.asarray([0]*len(oh))
    return oh, doh
def PP04(o3, n2, ha, hb, do3=None, dn2=None, dha=None, dhb=None):
    '''
    O3N2 with the Pettini & Pegal 2004 calibration
    
    O3N2 is only valid in the regime O3N2 < 2.0
    '''
    ebv = EBV(ha, hb)
    o3 = ccm_unred([5008.]*len(o3), o3, ebv = ebv)
    n2 = ccm_unred([6585.]*len(n2), n2, ebv = ebv)
    ha = ccm_unred([6565.]*len(ha), ha, ebv = ebv)
    hb = ccm_unred([4863.]*len(hb), hb, ebv = ebv)
    
    o3n2 = log10(divide(o3, hb) * divide(ha, n2))
    
    oh = 8.73 - 0.32*o3n2
    #oh = oh * (o3n2 < 2.0) + 0 * (o3n2 >= 2.0)
    oh[o3n2>=2.0] = np.nan
    
    if do3 is not None and dn2 is not None and dha is not None and dhb is not None:
        do3 = ccm_unred([5008.]*len(o3), do3*o3, ebv = ebv)
        dn2 = ccm_unred([6585.]*len(ha), dn2*n2, ebv = ebv)
        dha = ccm_unred([6565.]*len(ha), dha*ha, ebv = ebv)
        dhb = ccm_unred([4863.]*len(hb), dhb*hb, ebv = ebv)
        
        do3n2 = np.sqrt(np.abs(divide(do3,(o3*np.log(10))))**2 + np.abs(divide(dhb,(hb*np.log(10))))**2 +\
                np.abs(divide(dha,(ha*np.log(10))))**2 + np.abs(divide(dn2,(n2*np.log(10))))**2)
        doh = np.abs(0.32*do3n2)
        doh = doh * (o3n2 < 2.0) + 0 * (o3n2 >= 2.0)
    else:
        doh = np.zeros_like(oh)
    return oh, doh

def ionization(oh, y):
    # Tables from Kewley and Dopita 2002
    k = np.array([[7.39167, 7.46218, 7.57817, 7.73013], 
                  [0.667891, 0.685835, 0.739315, 0.843125], 
                  [0.0680367, 0.0866086, 0.0843640, 0.118166]])
    
    kq = np.array([[-27.0004, -31.2133, -36.0239, -40.9994, -44.7026, -46.1589, -45.6075],
                   [6.03910, 7.15810, 8.44804, 9.78396, 10.8052, 11.2557, 11.2074], 
                   [-0.327006, -0.399343, -0.483762, -0.571551, -0.640113, -0.672731, -0.674460]])
    
    zl = np.array([0.05, 0.1, 0.2, 0.5])  # Columns from Table 2
    ql = np.array([5*10**6, 1*10**7, 2*10**7, 4*10**7, 8*10**7, 1.5*10**8, 3*10**8]) # Columns from Table 3
    
    z = np.array([oh]) - 12 - (-3.07) # Solar abundance of Oxygen
    k0 = []
    k1 = []
    k2 = []
    for i in range(len(z)):
        if ~np.isnan(z[i]):
            ix = np.where(abs(z[i] - zl) == abs(z[i] - zl).min())[0][0]
            k0.append(k[0,ix])
            k1.append(k[1,ix])
            k2.append(k[2,ix])
        else:
            k0.append(np.nan)
            k1.append(np.nan)
            k2.append(np.nan)
    k0 = np.asarray(k0)
    k1 = np.asarray(k1)
    k2 = np.asarray(k2)
    
    q = 10**(k0 + k1*y + k2*y**2)
    
    k0 = []
    k1 = []
    k2 = []
    for i in range(len(q)):
        if ~np.isnan(q[i]):
            ix = np.where(abs(q[i] - ql) == abs(q[i] - ql).min())[0][0]
            k0.append(kq[0,ix])
            k1.append(kq[1,ix])
            k2.append(kq[2,ix])
        else:
            k0.append(np.nan)
            k1.append(np.nan)
            k2.append(np.nan)
    k0 = np.asarray(k0)
    k1 = np.asarray(k1)
    k2 = np.asarray(k2)
    return k0, k1, k2

def KD02(o2, o3, n2, s2_6718, s2_6733, hb, ha=None, do2=None, do3=None, dn2=None, ds2_6718=None, ds2_6733=None):
    '''
    Kewley & Dopita 2002 recomended method
    
    1. if 12+logO/H >= 8.6
        Use the N2O2 method
    2. if 8.6 > 12+logO/H >= 8.5
        Use the average of the Z94 and upper branch of the M91 calibrators
    3. Else if 12+logO/H < 8.5
        Use the average of the C01 and their R23 calibrators where the R23
        calibrator is calculated with the ionization parameter q. To get q,
        use C01 as an initial guess for 12+logO/H and use OII/OIII with the 
        coefficients from Table 2. The use q to get the coefficients from 
        Table 3 to solve for the R23 calibrator. Now use the 12+logO/H as the
        initial guess for q and repeat the procces until the solution converges.
        Typically this just needs to loop once. Use the average of the R23 
        12+logO/H and C01 as the metallicity.    
    
    WARNING: This pipeline utilizes the OII and OIII doublets in some cases and 
    just a single emission line in others.
    '''
    
    # Deredden the emission lines if have ha and hb
    if ha is not None and hb is not None:
        ebv = EBV(ha, hb)
        o2 = ccm_unred([3730.]*len(o2), o2, ebv=ebv)
        n2 = ccm_unred([6585.]*len(n2), n2, ebv=ebv)
        o3 = ccm_unred([5008.]*len(o3), o3, ebv = ebv)
        hb = ccm_unred([4863.]*len(hb), hb, ebv = ebv)
        s2 = ccm_unred([6733.]*len(s2_6733), s2_6733, ebv = ebv) + ccm_unred([6718.]*len(s2_6718), s2_6718, ebv = ebv)
    else:
        s2 = s2_6718 + s2_6733

    oh = np.zeros(o2.size)
    for i in range(len(o2)):
        # Upper Branch, O/H > 8.6
        n2o2 = log10(divide(n2[i], 2*o2[i]))  # uses [OII]lamlam 3726, 3729
        ohup = log10(1.54020 + 1.26602 * n2o2 + 0.167977*n2o2**2) + 8.93
        
        if ohup >= 8.6:
            oh[i] = ohup
        else:
            # Middle Branch, 8.6 > O/H > 8.5
            R23 = log10(divide((o2[i] + 1.35*o3[i]), hb[i])) # uses the doublet for OIII but not OII
            y = log10(divide(1.35*o3[i], o2[i]))  # uses [OIII]lamlam 4959, 5007 and [OII]lam 3727
            Z94 = 9.265 - 0.33 * R23 - 0.202 * R23**2 - 0.207 * R23**3 - 0.333 * R23**4 # Zaritsky et al. 1994
            M91 = 12. - 4.944 + 0.767 * R23 + 0.602 * R23**2 - y * (0.29 + 0.332 * R23 - 0.311 * R23**2) #McGaugh 1991
            ohmid = np.nanmean((Z94, M91), axis=0)

            if ohmid >= 8.5:
                oh[i] = ohmid
                print ohmid
            else:
                # Lower Branch, O/H < 8.5
                C01 = log10(5.09*10**-4 * (divide(o2[i], o3[i])/1.5)**0.17 * (divide(n2[i], s2[i])/0.85)**1.17) + 12
                
                R = log10(divide(o3[i], o2[i]))
                k0, k1, k2 = ionization(C01, R)
                ohloop = divide((-k1 + np.sqrt(k1**2 - 4*k2*(k0 - R23))) , 2*k2)
                ohlast_iter = C01
                while abs(ohlast_iter-ohloop)/ohloop > 0.01:
                    ohlast_iter = ohloop
                    k0, k1, k2 = ionization(ohloop, R)
                    ohloop = divide((-k1 + np.sqrt(k1**2 - 4*k2*(k0 - R23))) , 2*k2)    
                ohlow = np.nanmean((ohloop, C01), axis=0)
                oh[i] = ohlow

    oh = np.asarray(oh)
    #oh = ohmid
    #print ohlow[~np.isnan(ohlow)]
    '''
    NEED to update the error propagation
    '''
    if do2 is not None and dn2 is not None:
        dn2 = ccm_unred([6585.]*len(n2), dn2*n2, ebv = ebv)
        do2 = ccm_unred([3730.]*len(o2), do2*o2, ebv = ebv)
        do3 = ccm_unred([5008.]*len(o3), do3*o3, ebv = ebv)
        ds2 = np.sqrt(ccm_unred([6733.]*len(ha), ds2_6733*s2_6733, ebv = ebv)**2 + ccm_unred([6718.]*len(ha), ds2_6718*s2_6718, ebv = ebv)**2)
        
        dn2o2 = np.sqrt(divide(do2,(o2*np.log(10)))**2 + divide(dn2,(n2*np.log(10)))**2)
        dy = np.sqrt(1.26602**2 * dn2o2**2 + ((2.0*0.167977*n2o2)**2*dn2o2**2))
        doh = np.abs(divide(dy, y*np.log(10)))
    else:
        doh = np.zeros_like(oh)
    return oh, doh

def KD02KE08(o2, o3, n2, s2_6718, s2_6733, hb, ha=None, do2=None, do3=None, dn2=None, ds2_6718=None, ds2_6733=None):
    '''
    KD02 as updated in Kewley & Ellison 2008.
    '''
    
    # Deredden the emission lines if have ha and hb
    if ha is not None and hb is not None:
        ebv = EBV(ha, hb)
        o2 = ccm_unred([3730.]*len(o2), o2, ebv=ebv)
        n2 = ccm_unred([6585.]*len(n2), n2, ebv=ebv)
        o3 = ccm_unred([5008.]*len(o3), o3, ebv = ebv)
        hb = ccm_unred([4863.]*len(hb), hb, ebv = ebv)
        s2 = ccm_unred([6733.]*len(s2_6733), s2_6733, ebv = ebv) + ccm_unred([6718.]*len(s2_6718), s2_6718, ebv = ebv)
    else:
        s2 = s2_6718 + s2_6733

    oh = np.zeros(o2.size)
    
    total_errors = 0
    for i in range(len(o2)):
        # Upper Branch, O/H > 8.6
        n2o2 = log10(divide(n2[i], 2*o2[i]))  # uses [OII]lamlam 3726, 3729
        
        if n2o2 >= -1.2:
            '''
            solve the roots of log(N2/O2) = 1106.8660 - 532.15451 Z + 96.373260 Z^2 - 7.8106123 Z^3 + 0.23928247 Z^4
            '''
            #coeffs = [1106.8660 - n2o2, -532.15451, 96.373260, -7.8106123, 0.23928247]
            coeffs = np.array([0.23928247, -7.8106123, 96.373260, -532.15451, 1106.8660 - n2o2])
            root = np.roots(coeffs)
            oh[i] = root[~np.iscomplex(root)].real.max()
            if oh[i] <=8.4:
                print 'error'
                break
            
        else:
            R23 = log10(divide((o2[i] + 1.35*o3[i]), hb[i])) # uses the doublet for OIII but not OII
            y = log10(divide(1.35*o3[i], o2[i]))  # uses [OIII]lamlam 4959, 5007 and [OII]lam 3727
            M91 = 12. - 4.944 + 0.767 * R23 + 0.602 * R23**2 - y * (0.29 + 0.332 * R23 - 0.311 * R23**2) #McGaugh 1991
            
            #KK04
            y = log10(o3[i]/o2[i])
            
            logoh = 8.2 # for the lower branch, for the upper branch use 8.7
            
            logq = (32.81 - 1.153 * y**2 + logoh*(-3.396-0.025*y+0.1444*y**2))*(4.603-0.3119*y-0.163*y**2+logoh*(-0.48+0.0271*y+0.02037*y**2))**-1
            
            KK04 = 9.40 + 4.65*R23 - 3.17*R23**2 - logq*(0.272 + 0.547*R23-0.513*R23**2)
            
            #KK04 = 9.72 -0.777*R23 - 0.951*R23**2 -0.072*R23**3 -0.811*R23**4 -logq*(0.0737-0.0713*R23-0.141*R23**2 +0.0373*R23**3 -0.058*R23**4)
            
            # now use the use KK04 to solve for logq and iterate the logq and
            # KK04 unitl the solution converges
            
            error = 0
            itera = []
            while abs((KK04-logoh)/logoh) > 0.00001:
                logoh = KK04
                logq = (32.81 - 1.153 * y**2 + logoh*(-3.396-0.025*y+0.1444*y**2))*(4.603-0.3119*y-0.163*y**2+logoh*(-0.48+0.0271*y+0.02037*y**2))**-1
                KK04 = 9.40 + 4.65*R23 - 3.17*R23**2 - logq*(0.272 + 0.547*R23-0.513*R23**2)
                #KK04 = 9.72 -0.777*R23 - 0.951*R23**2 -0.072*R23**3 -0.811*R23**4 -logq*(0.0737-0.0713*R23-0.141*R23**2 +0.0373*R23**3 -0.058*R23**4)
                itera.append(KK04)
                error += 1
                if error > 100:
                    total_errors += 1
                    
                    KK04 = np.nan
                    #print i
                    break
                #print abs(abs(KK04-logoh)/KK04) 
            if len(itera) > 100:
                break
            oh[i] = np.mean((M91, KK04), axis=0)
            #update_progress((i+1.0)/np.float64(len(o2)))
    #print str(total_errors) + ' iteration failures'
    oh = np.asarray(oh)
    #oh = ohmid
    #print ohlow[~np.isnan(ohlow)]
    '''
    NEED to update the error propagation
    '''
    '''
    if do2 is not None and dn2 is not None:
        dn2 = ccm_unred([6585.]*len(n2), dn2*n2, ebv = ebv)
        do2 = ccm_unred([3730.]*len(o2), do2*o2, ebv = ebv)
        do3 = ccm_unred([5008.]*len(o3), do3*o3, ebv = ebv)
        ds2 = np.sqrt(ccm_unred([6733.]*len(ha), ds2_6733*s2_6733, ebv = ebv)**2 + ccm_unred([6718.]*len(ha), ds2_6718*s2_6718, ebv = ebv)**2)
        
        dn2o2 = np.sqrt(divide(do2,(o2*np.log(10)))**2 + divide(dn2,(n2*np.log(10)))**2)
        dy = np.sqrt(1.26602**2 * dn2o2**2 + ((2.0*0.167977*n2o2)**2*dn2o2**2))
        doh = np.abs(divide(dy, y*np.log(10)))
    else:
        doh = np.zeros_like(oh)
    '''
    doh = np.zeros_like(oh)
    return oh, doh#, itera

def test(o2, o3, n2, s2_6718, s2_6733, hb, ha=None, do2=None, do3=None, dn2=None, ds2_6718=None, ds2_6733=None):
    '''
    KD02 as updated in Kewley & Ellison 2008.
    '''
    
    # Deredden the emission lines if have ha and hb
    if ha is not None and hb is not None:
        ebv = EBV(ha, hb)
        o2 = ccm_unred([3730.]*len(o2), o2, ebv=ebv)
        n2 = ccm_unred([6585.]*len(n2), n2, ebv=ebv)
        o3 = ccm_unred([5008.]*len(o3), o3, ebv = ebv)
        hb = ccm_unred([4863.]*len(hb), hb, ebv = ebv)
        s2 = ccm_unred([6733.]*len(s2_6733), s2_6733, ebv = ebv) + ccm_unred([6718.]*len(s2_6718), s2_6718, ebv = ebv)
    else:
        s2 = s2_6718 + s2_6733

    oh = np.zeros(o2.size)
    n = np.zeros(o2.size)
    for i in range(len(o2)):
        # Upper Branch, O/H > 8.6
        n2o2 = log10(divide(n2[i], 2*o2[i]))  # uses [OII]lamlam 3726, 3729
        n[i] = n2o2
        '''
        solve the roots of log(N2/O2) = 1106.8660 - 532.15451 Z + 96.373260 Z^2 - 7.8106123 Z^3 + 0.23928247 Z^4
        '''
        coeffs = [1106.8660 - n2o2, -532.15451, 96.373260, -7.8106123, 0.23928247]
        coeffs = np.array([0.23928247, -7.8106123, 96.373260, -532.15451, 1106.8660 - n2o2])
        root = np.roots(coeffs)
            
        if any(~np.iscomplex(root)):
            
            oh[i] = root[~np.iscomplex(root)].real.min()
            
            
        else:
            oh[i] = np.nan


    oh = np.asarray(oh)
    #oh = ohmid
    #print ohlow[~np.isnan(ohlow)]
    '''
    NEED to update the error propagation
    '''
    if do2 is not None and dn2 is not None:
        dn2 = ccm_unred([6585.]*len(n2), dn2*n2, ebv = ebv)
        do2 = ccm_unred([3730.]*len(o2), do2*o2, ebv = ebv)
        do3 = ccm_unred([5008.]*len(o3), do3*o3, ebv = ebv)
        ds2 = np.sqrt(ccm_unred([6733.]*len(ha), ds2_6733*s2_6733, ebv = ebv)**2 + ccm_unred([6718.]*len(ha), ds2_6718*s2_6718, ebv = ebv)**2)
        
        dn2o2 = np.sqrt(divide(do2,(o2*np.log(10)))**2 + divide(dn2,(n2*np.log(10)))**2)
        dy = np.sqrt(1.26602**2 * dn2o2**2 + ((2.0*0.167977*n2o2)**2*dn2o2**2))
        doh = np.abs(divide(dy, y*np.log(10)))
    else:
        doh = np.zeros_like(oh)
    return oh, n

def DOP16(n2, s2_6718, s2_6733, ha, hb=None, dn2=None, ds2_6718=None, ds2_6733=None, dha=None, dhb=None):    
    # N2S2 from Dopita et al. 2016
    if hb is not None:
        ebv = EBV(ha, hb)
        n2c = ccm_unred([6585.]*len(ha), n2, ebv = ebv)
        s2c = ccm_unred([6733.]*len(ha), s2_6733, ebv = ebv) + ccm_unred([6718.]*len(ha), s2_6718, ebv = ebv)
        hac = ccm_unred([6565.]*len(ha), ha, ebv = ebv)
    
        ha = ha * (ebv<0) + hac * (ebv>=0)
        s2 = (s2_6718+s2_6733) * (ebv<0) + s2c * (ebv>=0)
        n2 = n2 * (ebv<0) + n2c * (ebv>=0) 
    else:
        s2 = s2_6733 + s2_6718
    #y = 1.264*log10(n2) - log10(s2) - 0.264 * log10(ha)
    y = log10(divide(n2, s2)) + 0.264 * log10(divide(n2, ha))
    oh = 8.77 + y + 0.45*(y + 0.3)**5
    if dn2 is not None and ds2_6718 is not None and dha is not None:
        if hb is not None:
            dn2 = ccm_unred([6585.]*len(ha), dn2*n2, ebv = ebv)
            ds2 = np.sqrt(ccm_unred([6733.]*len(ha), ds2_6733*s2_6733, ebv = ebv)**2 + ccm_unred([6718.]*len(ha), ds2_6718*s2_6718, ebv = ebv)**2)
            dha = ccm_unred([6565.]*len(ha), dha*ha, ebv = ebv)
        else:
            ds2 = np.sqrt((ds2_6733*s2_6733)**2 + (ds2_6718*s2_6718)**2)
        dy = np.sqrt(np.abs(divide(1.264*dn2,(n2*np.log(10))))**2 + np.abs(divide(ds2,(s2*np.log(10))))**2 + np.abs(divide(0.264*dha,(ha*np.log(10))))**2)
        doh = np.abs(dy + (5*0.45)*(y+0.3)**4 * dy)
    else:
        doh = np.zeros_like(oh)
    return oh, doh
'''
# Issue with ONS where there is a large gap in probed metallicities
def ONS(obj):
    # From Pilyugin et al. 2010
    o2 = pars['OII3730'][:,obj]
    o3 = pars['OIII5008'][:,obj]
    n2 = pars['NII6585'][:,obj]
    s2 = pars['SII6733'][:,obj] + pars['SII6718'][:,obj]
    hb = pars['Hb4863'][:,obj]
    R2 = np.divide(2*o2, hb, out=np.zeros(o2.shape), where=hb!=0)
    N2 = np.divide(1.34*n2, hb, out=np.zeros(n2.shape), where=hb!=0)
    S2 = np.divide(s2, hb, out=np.zeros(s2.shape), where=hb!=0)
    R3 = np.divide(1.35*o3, hb, out=np.zeros(o3.shape), where=hb!=0)
    P = np.divide(R3, (R3+R2), out=np.zeros(R3.shape), where=(R3+R2)!=0)
    N2R2 = np.divide(N2, R2, out=np.zeros(N2.shape), where=R2!=0)
    S2R2 = np.divide(S2, R2, out=np.zeros(S2.shape), where=R2!=0)
    N2S2 = np.divide(N2, S2, out=np.zeros(N2.shape), where=S2!=0)
    
    oh = []
    for i in range(len(o2)):
        if np.log10(N2[i]) >= -0.1:
            oh.append(8.277 + 0.657*P[i] - 0.399*np.log10(R3[i]) - 0.061*np.log10(N2R2[i]) + 0.005*np.log10(S2R2[i]))
        elif np.log10(N2[i]) < -0.1 and np.log10(N2S2[i]) >= -0.25:
            oh.append(8.816 - 0.733*P[i]-0.454*np.log10(R3[i]) + 0.71*np.log10(N2R2[i]) - 0.337*np.log10(S2R2[i]))
        elif np.log10(N2[i]) < -0.1 and np.log10(N2S2[i]) < -0.25:
            oh.append(8.774 - 1.855*P[i] + 1.517*np.log10(R3[i])+0.304*np.log10(N2R2[i]) + 0.328*np.log10(S2R2[i]))
        else:
            oh.append(np.nan)
    oh = np.asarray(oh)
    return oh
'''

def poly_N06(logR, a):
    '''
    Find the roots of polynomials for Nagao et al. 2006
    '''
    
    roots = np.roots([a[3], a[2], a[1], a[0] - logR])
    # range of Nagao is 12+log(O/H) = 7 - 9.5
    if any(~np.iscomplex(roots)&(roots>=7.05)&(roots <= 9.25)):
        return roots[~np.iscomplex(roots)&(roots>=7.05)&(roots <= 9.25)].real
    else:
        return np.array([np.nan])
    
def rootap_N06(ratio, cal, oh, coeffs):
    R = log10(cal)
    for i in range(len(R)):
        Z = poly_N06(R[i], coeffs[ratio])
        oh.append(Z.min())

def N06(ratio, OII_3727=None, OIII_5007=None, NII_6584=None, 
        SII_6720=None, NeIII_3869=None, Ha=None, Hb=None, unred=False):
    '''
    Line ratios available: R23, N2Ha, O3N2, N2O2, N2S2, O3O2, Ne3O2, Ha+N2/S2, 
    O3Hb, or O2Hb
    '''
    
    lines = [OII_3727, OIII_5007, NII_6584, SII_6720, NeIII_3869, Ha, Hb]
    wave = [3730, 5008, 6585, 6720, 3870, 6565, 4863]
    
    # polynomial coefficients from Tables 6 and 10 from Nagao et al. 2006   
    coeffs = {"R23":[1.2299e0, -4.1926e0, 1.0246e0, -6.3169e-2], 
              "N2Ha":[9.6641e1, -3.9941e1, 5.2227e0, -2.2040e-1],
              "O3N2":[-2.3218e2, 8.4423e1, -9.9330e0, 3.7941e-1],
              "N2O2":[-1.2894e2, 4.8818e1, -6.2759e0, 2.7101e-1],
              "N2S2":[-8.0632e1, 3.1323e1, -4.1010e0, 1.7963e-1],
              "O3O2":[-1.4089e0, 1.3745e0, -1.4359e-1, 0],
              "Ne3O2":[-8.2202e1, 3.2014e1, -4.0706e0, 1.6784e-1],
              "HaN2S2":[-1.4789e2, 5.9974e1, -7.8963e0, 3.4069e-1],
              "O3Hb":[-5.7906e1, 1.5437e1, -1.0872e0, 9.1853e-3],
              "O2Hb":[1.5253e2, -6.1922e1, 8.2370e0, -3.5951e-1]}

    # perform de-reddening for lines 
    if unred == True and Ha is not None and Hb is not None:
        ebv = EBV(Ha, Hb)
        for i in range(len(lines)):
            if lines[i] is not None:
                lines[i] = ccm_unred([wave[i]]*len(lines[i]), lines[i], ebv=ebv)
    
    oh = []
    if ratio == 'R23':
        # separate upper branch from lower branch with O3/O2
        R = log10(R23(OII_3727, OIII_5007, Hb))
        for i in range(len(R)):
            if OIII_5007[i]/OII_3727[i] <= 2:
                Z = poly_N06(R[i], coeffs[ratio]).max()
            else:
                Z = poly_N06(R[i], coeffs[ratio]).min()
            oh.append(Z)
            
    elif ratio =='N2Ha':
        rootap_N06(ratio, N2Ha(NII_6584, Ha), oh, coeffs)
    elif ratio=='O3N2':
        rootap_N06(ratio, O3N2(OIII_5007, NII_6584), oh, coeffs)
    elif ratio=='N2O2':
        rootap_N06(ratio, N2O2(NII_6584, OII_3727), oh, coeffs)
    elif ratio=='N2S2':
        rootap_N06(ratio, N2S2(NII_6584, SII_6720), oh, coeffs)
    elif ratio=='O3O2':
        rootap_N06(ratio, O3O2(OIII_5007, OII_3727), oh, coeffs)
    elif ratio=='Ne3O2':
        rootap_N06(ratio, Ne3O2(NeIII_3869, OII_3727), oh, coeffs)
    elif ratio=='HaN2S2':
        rootap_N06(ratio, HaN2S2(Ha, NII_6584, SII_6720), oh, coeffs)
    elif ratio=='O3Hb':
        rootap_N06(ratio, O3Hb(OIII_5007, Hb), oh, coeffs)
    elif ratio=='O2Hb':
        rootap_N06(ratio, O2Hb(OII_3727, Hb), oh, coeffs)
            
    return np.asarray(oh), np.asarray(oh)

def poly_M08(logR, a):
    '''
    Find the roots of polynomials 
    '''
    
    roots = np.roots([a[4], a[3], a[2], a[1], a[0] - logR])
    # range of M08 is 12+log(O/H) = 7 - 9.5
    if any(~np.iscomplex(roots)):
        return np.real(roots[~np.iscomplex(roots)])
    else:
        return np.array([np.nan])

def M08(ratio, OII_3727=None, OIII_5007=None, NII_6584=None, 
        NeIII_3869=None, Ha=None, Hb=None, unred=False):
    '''
    Line ratios available: R23, N2Ha, O3N2, O3O2, Ne3O2, 
    O3Hb, or O2Hb
    '''
    lines = [OII_3727, OIII_5007, NII_6584, NeIII_3869, Ha, Hb]
    wave = [3730, 5008, 6585, 3870, 6565, 4863]
    
    # polynomial coefficients from Tables 6 and 10 from Maiolino et al. 2008   
    coeffs = {"R23":[0.7462, -0.7149, -0.9401, -0.6154, -0.2524], 
              "N2Ha":[-0.7732, 1.2357, -0.2811, -0.7201, -0.3330],
              "O3Hb":[0.1549, -1.5031, -0.9790, -0.0297, 0],
              "O3O2":[-0.2839, -1.3881, -0.3172, 0, 0],
              "O2Hb":[0.5603, 0.0450, -1.8017, -1.8434, -0.6549],
              "O3N2":[0.4520, -2.6096, -0.7170, 0.1347, 0],
              "Ne3O2":[-1.2608, -1.0861, -0.1470, 0, 0]
              }
    '''
    # perform de-reddening for lines 
    if unred == True and Ha is not None and Hb is not None:
        ebv = EBV(Ha, Hb)
        for i in range(len(lines)):
            if lines[i] is not None:
                lines[i] = ccm_unred([wave[i]]*len(lines[i]), lines[i], ebv=ebv)
    '''
    
    ebv = EBV(Ha, Hb)
    
    OII_3727c = ccm_unred([3730.]*len(OII_3727), OII_3727, ebv=ebv)
    NII_6584c = ccm_unred([6585.]*len(NII_6584), NII_6584, ebv=ebv)
    OIII_5007c = ccm_unred([5008.]*len(OIII_5007), OIII_5007, ebv = ebv)
    Hbc = ccm_unred([4863.]*len(Hb), Hb, ebv = ebv)
    Hac = ccm_unred([6565.]*len(Ha), Ha, ebv = ebv)
    
    OII_3727 = OII_3727 * (ebv<0) + OII_3727c * (ebv>=0)
    OIII_5007 = OIII_5007 * (ebv<0) + OIII_5007c * (ebv>=0)
    NII_6584 = NII_6584 * (ebv<0) + NII_6584c * (ebv>=0)
    Ha = Ha * (ebv<0) + Hac * (ebv>=0)
    Hb = Hb * (ebv<0) + Hbc * (ebv>=0)
    
    # Account for Doublets
    OII_3727 = OII_3727 * 2.0
    OIII_5007 = OIII_5007 * 1.35
    
    oh = []
    if ratio == 'R23':
        # separate upper branch from lower branch with O3/O2
        R = log10(R23(OII_3727, OIII_5007, Hb))
        for i in range(len(R)):
            # Use the NII/Ha ratio to break the double valued R23
            # log(NII/Ha) = -1.39 at 12 + log(O/H) = 8.08
            #if OIII_5007[i]/OII_3727[i] <= 2:
            
            if log10(divide(NII_6584[i], Ha[i])) > -1.39:
                Z = poly_M08(R[i], coeffs[ratio]).max() + 8.69
            else:
                Z = poly_M08(R[i], coeffs[ratio]).min() + 8.69
            
            if Z > 9.5 or Z < 7:
                oh.append(np.nan)
            else:
                oh.append(Z)
    '''   
    elif ratio =='N2Ha':
        rootap(ratio, N2Ha(NII_6584, Ha), oh, coeffs)
    elif ratio=='O3N2':
        rootap(ratio, O3N2(OIII_5007, NII_6584), oh, coeffs)
    elif ratio=='O3O2':
        rootap(ratio, O3O2(OIII_5007, OII_3727), oh, coeffs)
    elif ratio=='Ne3O2':
        rootap(ratio, Ne3O2(NeIII_3869, OII_3727), oh, coeffs)
    elif ratio=='O3Hb':
        rootap(ratio, O3Hb(OIII_5007, Hb), oh, coeffs)
    elif ratio=='O2Hb':
        rootap(ratio, O2Hb(OII_3727, Hb), oh, coeffs)
    '''     
    return np.asarray(oh), np.asarray(oh)

def NO(NII_6548, NII_6584, OIII_5007, OIII_4959, OII_3727, Hb):
    '''
    The log(N+/O+) ratio.
    '''
    r23 = R23(OII_3727, OIII_5007, Hb)
    
    '''[NII] Temperature'''
    tII = 6065. + 1600. * log10(r23) * 1878. * (log10(r23))**2 + 2803. * (log10(r23))**3
    
    n2o2 = log10(divide((NII_6548 + NII_6584), OII_3727))
    
    return n2o2 + 0.307 - 0.02 * log10(tII) - 0.726 / tII