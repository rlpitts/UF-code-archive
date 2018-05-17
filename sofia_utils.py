#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 19:21:50 2018

@author: rlpitts
"""

import numpy as np
import astropy.io.fits as pf
from astropy.wcs import WCS
import matplotlib as mpl
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit as curvf

#imcctran to gal(J2000) coords
#convolve blue channel images to red channel resolution (FWHM=10")
#regrid blue channel images to red channel templates (2"/px)
#integrate over lines

#oi63=pf.open('workfiles/F0313_FI_IFS_04006131_BLU_WXY_OI63um.fits')
#oi145=pf.open('workfiles/F0313_FI_IFS_04006132_RED_WXY_OI145um.fits')
#oiii=pf.open('workfiles/F0313_FI_IFS_04006132_BLU_WXY_OIII88um.fits')
#cii=pf.open('workfiles/F0313_FI_IFS_04006131_RED_WXY_CII158um.fits')

#ext0 has no data, only the complete header - copy to flux?
#ext1-10 are the following:
#1 'FLUX': nx x ny x nw cube of flux values,
#2 'ERROR': associated error values on the flux (also nx x ny x nw),
#3 'UNCORRECTED_FLUX': nx x ny x nw cube of flux values that have not
#   been corrected for atmospheric transmission,
#4 'UNCORRECTED_ERROR': associated error values on the uncorrected flux
#5 'WAVELENGTH': wavelength values associated with each plane of the cube (nw),
#6 'X': The x-coordinates of the data, in arcsecond offsets from the base position (nx),
#7 'Y': The y-coordinates of the data, in arcsecond offsets from the base position (ny),
#8 'TRANSMISSION': The atmospheric transmission model (nw),
#9 'RESPONSE': The instrumental response curve (nw),
#10 'EXPOSURE_MAP': The exposure map (nx x ny x nw); basically shows dither positions &
#   what parts have overlapping coverage


def fitslice(fname):
    if type(fname) is pf.hdu.hdulist.HDUList:
        hdu = fname
    elif type(fname) is str:
        hdu=pf.open(fname)
        mainhd = hdu[0].header
    
    for i in xrange(len(hdu)):
        if i==0:
            pass
        else:
            hdr=hdu[i].header
            hdr.set('CHANNEL',mainhd['CHANNEL'],'Detector channel',after='EXTNAME')
            if 'BLUE' in mainhd['CHANNEL'] and 'CROTA2' in hdr.tostring():
                hdr.set('BMAJ',5/3600.,'beam major axis (deg)',after='CROTA2')
                hdr.set('BMIN',5/3600.,'beam minor axis (deg)',after='BMAJ')
            elif 'RED' in mainhd['CHANNEL'] and 'CROTA2' in hdr.tostring():
                hdr.set('BMAJ',10/3600.,'beam major axis (deg)',after='CROTA2')
                hdr.set('BMIN',10/3600.,'beam minor axis (deg)',after='BMAJ') 
                
            if 'OI63um' in fname:
                obase='fifi_oi_63um_'
            elif 'OI145um' in fname:
                obase='fifi_oi_145um_'
            elif 'OIII88um' in fname:
                obase='fifi_oiii_88um_'
            elif 'CII158um' in fname:
                obase='fifi_cii_158um_'
                
            oname = 'workfiles/'+obase+hdr['EXTNAME']+'.fits'
            newhdu = pf.PrimaryHDU(data=hdu[i].data,header=hdr)
            newhdu.writeto(oname,clobber=True)
            print 'wrote {} extension to {}'.format(hdr['EXTNAME'],oname)
    hdu.close()
    return None

def gauss(x,a,m,s):
    return a * np.exp( (-(x-m)**2) / (2 * s**2) )

def lorentz(x,a,m,hwhm):
    return a / ( 1 + ( (x - m)/hwhm )**2 )

def absorbed_gauss(x,a1,m1,s1,a2,m2,s2):
    return gauss(x,a1,m1,s1)+gauss(x,a2,m2,s2)

def pltspec(cube,wvl,x,y,ecube=None,clr='k',mrkr=None,mec=None,mfc=None,ms=None,
            ls=None,lw=0.75,ds='steps-mid',gfit=False,dgfit=False,lorfit=False):
    spec=cube[:,x,y]
    if ecube is not None:
        uspec=ecube[:,x,y]
        pl.errorbar(wvl,spec,yerr=uspec,color=clr,elinewidth=lw,capsize=ms,
                    marker=mrkr,mec=mec,mfc=mfc,ms=ms,ls=ls,lw=lw,drawstyle=ds)
    else:
        pl.plot(wvl,spec,color=clr,marker=mrkr,mec=mec,mfc=mfc,ms=ms,
                ls=ls,lw=lw,drawstyle=ds)
    pl.xlabel('Wavelength ($\mu$m)')
    pl.ylabel('Flux Density (Jy)')
    if gfit is True:
        p,pcov=curvf(gauss,wvl[~np.isnan(spec)],spec[~np.isnan(spec)],
                     sigma=uspec[~np.isnan(spec)] if uspec is not None else None,
                     p0=[np.nanmax(spec),wvl[np.argmax(spec)],0.05])
        up=np.sqrt(np.diag(pcov))
        print p, pcov, up
        wgp=np.linspace(np.nanmin(wvl),np.nanmax(wvl),100)
        specgp=gauss(wgp,*p)
        wfwhm=p[2]*np.sqrt(8*np.log(2))
#        lblp = ['$\pm$'.join(['{:.2}'.format(p[0]),'{:.1}'.format(up[0])]),
#                '$\pm$'.join(['{:.3}'.format(p[1]),'{:.1}'.format(up[1])]),
#                '$\pm$'.join(['{:.3}'.format(p[2]),'{:.1}'.format(up[2])])]
        pl.plot(wgp,specgp,'r-',
                label='FWHM={}'.format(wfwhm)+'$\mu$m')
    if dgfit is True:
        q,qcov=curvf(absorbed_gauss,wvl[~np.isnan(spec)],spec[~np.isnan(spec)],
                     sigma=uspec[~np.isnan(spec)] if uspec is not None else None,
                     p0=[5.,np.nanmean(wvl),0.05,2.,np.nanmean(wvl),0.01])
        uq=np.sqrt(np.diag(qcov))
        print q, qcov, uq
        wgq=np.linspace(np.nanmin(wvl),np.nanmax(wvl),100)
        specgq=gauss(wgq,*q)
        lblq = ['$\pm$'.join(['{:.2}'.format(q[0]),'{:.1}'.format(uq[0])]),
                '$\pm$'.join(['{:.3}'.format(q[1]),'{:.1}'.format(uq[1])]),
                '$\pm$'.join(['{:.3}'.format(q[2]),'{:.1}'.format(uq[2])])]
        pl.plot(wgq,specgq,'g-',label='{}{}({},{})'.format(lblq[0],'$N$',lblq[1],lblq[2]))
    if lorfit is True:
        l,lcov=curvf(lorentz,wvl[~np.isnan(spec)],spec[~np.isnan(spec)],
                     sigma=uspec[~np.isnan(spec)] if uspec is not None else None,
                     p0=[np.nanmax(spec),wvl[np.argmax(spec)],0.5])#,
                     #bounds=[[0,wvl[0],0],[np.nanmax(spec)+1,wvl[1],wvl[1]-wvl[0]]])
        ul=np.sqrt(np.diag(lcov))
        print l, lcov, ul
        wl=np.linspace(np.nanmin(wvl),np.nanmax(wvl),100)
        specl=lorentz(wl,l[0],l[1],abs(l[2]))
        fwhm=2*abs(l[2])
#        lblp = ['$\pm$'.join(['{:.2}'.format(p[0]),'{:.1}'.format(up[0])]),
#                '$\pm$'.join(['{:.3}'.format(p[1]),'{:.1}'.format(up[1])]),
#                '$\pm$'.join(['{:.3}'.format(p[2]),'{:.1}'.format(up[2])])]
        pl.plot(wl,specl,'r-',
                label='FWHM={:.3}'.format(fwhm)+'$\mu$m')

    pl.legend(loc=0)
    pl.show()
    return None

