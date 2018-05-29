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
import scipy.integrate as spi
import warnings
from scipy.optimize import OptimizeWarning

warnings.simplefilter("error", OptimizeWarning)

#imcctran to gal(J2000) coords
#convolve blue channel images to red channel resolution (FWHM=10")
#regrid blue channel images to red channel templates (2"/px)
#integrate over lines

oi63=pf.open('workfiles/F0313_FI_IFS_04006131_BLU_WXY_OI63um.fits')
oi145=pf.open('workfiles/F0313_FI_IFS_04006132_RED_WXY_OI145um.fits')
oiii=pf.open('workfiles/F0313_FI_IFS_04006132_BLU_WXY_OIII88um.fits')
cii=pf.open('workfiles/F0313_FI_IFS_04006131_RED_WXY_CII158um.fits')

lw={'oi63':63.185,'oi145':145.53,'oiii':88.36,'cii':157.74}

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
            newhdu.writeto(oname,overwrite=True)
            print 'wrote {} extension to {}'.format(hdr['EXTNAME'],oname)
    hdu.close()
    return None

def gauss(x,a,m,s,y,c):
    return ( a/np.sqrt(2*np.pi*s) ) * np.exp( (-(x-m)**2) / (2 * s**2) ) + y*x+c

def lorentz(x,a,m,hwhm,y,c):
    return ( a/(np.pi*hwhm) ) / ( 1 + ( (x - m)/hwhm )**2 ) + y*x+c

def lreg(x,m,b):
    return m*x+b

def pltspec(cube,wvl,x,y,ecube=None,clr='k',mrkr=None,mec=None,mfc=None,ms=None,
            ls=None,lw=0.75,ds='steps-mid',gfit=False,lorfit=False):
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
                     #p0=[np.nanmax(spec),wvl[np.argmax(spec)],0.05,0.6],
                     bounds = [[0.01,min(wvl),np.mean(wvl[1:]-wvl[:-1]),-10,0],
                               [50.,max(wvl),0.5,10,10.0]])
        up=np.sqrt(np.diag(pcov))
        print p, pcov, up
        wgp=np.linspace(np.nanmin(wvl),np.nanmax(wvl),100)
        specgp=gauss(wgp,*p)
        wfwhm=p[2]*np.sqrt(8*np.log(2))
        pl.plot(wgp,specgp,'g-',label='FWHM={:.3}'.format(wfwhm)+'$\mu$m')
    if lorfit is True:
        l,lcov=curvf(lorentz,wvl[~np.isnan(spec)],spec[~np.isnan(spec)],
                     sigma=uspec[~np.isnan(spec)] if uspec is not None else None,
                     #p0=[5.,wvl[np.argmax(spec)],0.05,0.06],
                     bounds = [[0.01,min(wvl),np.mean(wvl[1:]-wvl[:-1]),-10,0.],
                               [50.,max(wvl),0.5,10,10.0]])#,
                     #bounds=[[0,wvl[0],0],[np.nanmax(spec)+1,wvl[1],wvl[1]-wvl[0]]])
        ul=np.sqrt(np.diag(lcov))
        print l, lcov, ul
        wl=np.linspace(np.nanmin(wvl),np.nanmax(wvl),100)
        specl=lorentz(wl,l[0],l[1],abs(l[2]),l[3])
        fwhm=2*abs(l[2])
        pl.plot(wl,specl,'r-',label='FWHM={:.3}'.format(fwhm)+'$\mu$m')

    pl.legend(loc=0)
    pl.show()
    return None

def splitflux(fden, wavelen=None,ferr=None,writeto=None,hdr=None,
            fitfn='gauss'):
    if type(fden) is str:
        hdu=pf.open(fden)
        if len(hdu)>1:
            cube=hdu[1].data
            ecube=hdu[2].data
            wvl=hdu[5].data
            hdr=hdu[0].header
        else:
            cube=hdu[0].data
            hdr=hdu[0].header
    else:
        cube=fden
    if ferr is not None:
        ecube=pf.open(ferr)[0].data if type(ferr) is str else ferr
    if wavelen is not None:
        wvl=pf.open(wavelen)[0].data if type(wavelen) is str else wavelen
    else:
        raise ValueError("If wavelength array isn't defined by now, you goofed")
    for key,val in lw.iteritems():
        if wvl[0] <= val <= wvl[-1]:
            lc = val
            ltxt = key
            print ltxt,lc
            if any(c.isdigit() for c in key):
                hdr.set('ION','[OI]',after='CHANNEL')
            else:
                hdr.set('ION','['+key.upper()+']',after='CHANNEL')
            hdr.set('LINECNTR',val,after='ION')
            break
    if fitfn not in ['gauss','lorentz']:
        raise IOError("Fitting function name must be 'lorentz' or 'gauss'")
    cubedims=np.shape(cube)
    frames=np.array([(n,w) for n,w in enumerate(wvl) if (lc-0.2<w<lc+0.2)])
    #3 sigma is about 27% larger than 2 fwhm
    inds,wls=np.transpose(frames)
    fdencube=np.zeros((cubedims[1],cubedims[2]))
    fduncube=np.zeros((cubedims[1],cubedims[2]))
    spindex=np.zeros((cubedims[1],cubedims[2]))
    spundex=np.zeros((cubedims[1],cubedims[2]))
    lincube=np.zeros((len(frames),cubedims[1],cubedims[2]))
    elincube=np.zeros((len(frames),cubedims[1],cubedims[2]))
    concube=np.zeros((cubedims))
    econcube=np.zeros((cubedims))
    for x in xrange(cubedims[1]):
        for y in xrange(cubedims[2]):
            spec = cube[:,x,y]
            uspec = ecube[:,x,y]
            try:
                'fitting...'
                p,pcov=curvf(globals()[fitfn],wvl[~np.isnan(spec)],
                             spec[~np.isnan(spec)],
                             sigma=uspec[~np.isnan(spec)],
                             bounds = [[0.01,min(wvl),np.mean(wvl[1:]-wvl[:-1]),-10,0.],
                                       [50.,max(wvl),0.4,10,10.0]])
                fwhm = 2*abs(p[2]) if fitfn=='lorentz' else p[2]*np.sqrt(8*np.log(2))
                print 'p = ',p
                goodlc = lc-0.05<p[1]<lc+0.05
                goodcov = np.inf not in pcov
                print goodlc,goodcov
            except (ValueError,RuntimeError,NameError) as e:
                print e
                pass
                
            print ('fwhm' in locals() and locals()['fwhm']<0.16)
            print ('goodlc' in globals() and globals()['goodlc'] is True)
            print ('goodcov' in locals() and locals()['goodcov'] is True)
            
            if (('fwhm' in locals() and locals()['fwhm']<0.16) and 
                ('goodlc' in globals() and globals()['goodlc'] is True) and
                ('goodcov' in locals() and locals()['goodcov'] is True)):
                print 'all good'
                stdp=np.sqrt(np.diag(pcov))
                cvw=p[-2]*frames[:,1]+p[-1]
                lincube[:,x,y]=spec[inds.astype(int)]-cvw
                elincube[:,x,y]=np.sqrt(uspec[inds.astype(int)]**2+
                                        stdp[-2]**2+stdp[-1]**2)
                
                lvw=gauss(wvl,p[0],p[1],p[2],0,0)
                concube[:,x,y]=spec-lvw
                econcube[:,x,y]=np.sqrt(uspec**2+stdp[0]**2+
                                        stdp[1]**2+stdp[2]**2)
                spindex[x,y]=p[-2]
                spundex[x,y]=stdp[-2]
                fdencube[x,y]=p[0]
                fduncube[x,y]=stdp[0]
            else:
                try:
                    s=spec[~inds.astype(int)]
                    u=uspec[~inds.astype(int)]
                    q,qcov=curvf(lreg,wls[~np.isnan(s)],s[~np.isnan(s)],
                                 sigma=u[~np.isnan(s)],bounds = [[-10,0.],[10,10.0]])
                    r,rcov=curvf(globals()[fitfn], wvl[inds.astype(int)],
                                 spec[inds.astype(int)],
                                 sigma=uspec[inds.astype(int)],
                                 bounds = [[0.01,min(wvl),np.mean(wvl[1:]-wvl[:-1]),-10,0.],
                                            [50.,max(wvl),0.4,10,10.0]])
                    rfwhm = 2*abs(p[2]) if fitfn=='lorentz' else p[2]*np.sqrt(8*np.log(2))
                    goodlc2 = lc-0.05<r[1]<lc+0.05
                    goodcov2 = np.inf not in pcov
                    print rfwhm, goodlc2, goodcov2
                except (ValueError,RuntimeError,NameError) as e2:
                    print e2
                
                if (('rfwhm' in locals() and locals()['rfwhm']<0.16) and 
                    ('goodlc2' in globals() and globals()['goodlc2']==True) and
                    ('goodcov2' in locals() and locals()['goodcov2'] is True)):
                    print 'backup: fwhm = ',rfwhm
                    stdr=np.sqrt(np.diag(rcov))
                    stdq=np.sqrt(np.diag(qcov))
                    
                    lvw=gauss(wvl,r[0],r[1],r[2],0,0)
                    concube[:,x,y]=spec-lvw
                    econcube[:,x,y]=np.sqrt(uspec**2+stdr[0]**2+
                                    stdr[1]**2+stdr[2]**2)
                    cvw=q[0]*frames[:,1]+q[1]
                    lincube[:,x,y]=spec[inds.astype(int)]-cvw
                    elincube[:,x,y]=np.sqrt(uspec[inds.astype(int)]**2+
                                            stdq[-2]**2+stdq[-1]**2)
                                            
                    spindex[x,y]=q[0]
                    spundex[x,y]=stdq[0]
                    fdencube[x,y]=r[0]
                    fduncube[x,y]=stdr[0]
                else:
                    print 'line fit failed, attempting continuum fit'
                    fdencube[x,y]=np.NaN
                    fduncube[x,y]=np.NaN
                    lincube[:,x,y]=np.NaN
                    elincube[:,x,y]=np.NaN 
                    try:
                        q,qcov=curvf(lreg,wvl[~np.isnan(spec)],spec[~np.isnan(spec)],
                                     sigma=uspec[~np.isnan(spec)],bounds = [[-10,0.],[10,10.0]])
                        stdq = np.sqrt(np.diag(qcov))
                        concube[:,x,y]=spec
                        econcube[:,x,y]=uspec
                        spindex[x,y]=q[0]
                        spundex[x,y]=stdq[0] 
                    except (ValueError,RuntimeError,NameError,OptimizeWarning):
                        print 'all fits failed'
                        concube[:,x,y]=spec
                        econcube[:,x,y]=uspec
                        spindex[x,y]=np.NaN
                        spundex[x,y]=np.NaN
    lhdr=hdr.copy()
    lhdr.set('EMTYPE','Line emission',before='ION')
    lhdr.append('HISTORY','continuum removed with sofia_utils.splitflux()')
    lhdr['CRVAL3']=wls[0]
    lhdu1 = pf.PrimaryHDU(data=lincube,header=lhdr)
    lhdu2 = pf.ImageHDU(data=elincube,name='error')
    lhdu3 = pf.ImageHDU(data=fdencube,name='flux_density_fit')
    lhdu4 = pf.ImageHDU(data=fduncube,name='flux_density_err')
    lhdu5 = pf.ImageHDU(data=frames,name='frames')
    lhdu0 = pf.HDUList([lhdu1,lhdu2,lhdu3,lhdu4,lhdu5])
    lhdu0.writeto(writeto[0] if type(writeto) is list else 'workfiles/fifils-'+ltxt+'-lineonly.fits',
                  overwrite=True)
    chdr=hdr.copy()
    chdr.set('EMTYPE','Continuum emission',before='ION')
    chdr.append('HISTORY',hdr['ION']+' line emission removed with sofia_utils.splitflux()')
    chdu1 = pf.PrimaryHDU(data=concube,header=lhdr)
    chdu2 = pf.ImageHDU(data=econcube,name='error')
    chdu3 = pf.ImageHDU(data=spindex,name='spectral_index')
    chdu4 = pf.ImageHDU(data=spundex,name='spectral_index_error')    
    chdu0 = pf.HDUList([chdu1,chdu2,chdu3,chdu4])
    chdu0.writeto(writeto[1] if type(writeto) is list else 'workfiles/fifils-'+ltxt+'-continuum.fits',
                  overwrite=True)
    d = {'lflux':lincube,'cflux':concube,'lerror':elincube,'cerror':econcube}
    return d
    
    
def fluxint(fden, wavelen=None,ferr=None,method='trapz',writeto=None,hdr=None,
            intgauss=False,intlorentz=False,out_area='px',full_out=False):
    if type(fden) is str:
        hdu=pf.open(fden)
        if len(hdu)>1:
            cube=hdu[1].data
            ecube=hdu[2].data
            wvl=hdu[5].data
            hdr=hdu[0].header
        else:
            cube=hdu[0].data
            hdr=hdu[0].header
    else:
        cube=fden
    if ferr is not None:
        ecube=pf.open(ferr)[0].data if type(ferr) is str else ferr
    if wavelen is not None:
        wvl=pf.open(wavelen)[0].data if type(wvl) is str else wvl
    else:
        raise ValueError("If wavelength array isn't defined by now, you goofed")
    cubedims = np.shape(cube)
    freq=(2.99792458*10**8)/wvl
    icube = getattr(spi,method)(np.ma.masked_invalid(np.flip(cube,axis=0)*10**-26),
                    x=freq,axis=0)
    if ferr is not None:
        esqr = np.sqrt(np.nansum(np.flip(ecube,axis=0)**2,axis=0))*10**-26
    if intgauss is True or intlorentz is True:
        fn = 'lorentz' if intlorentz is True else 'gauss'
        fitsqr=np.zeros((cubedims[1],cubedims[2]))
        resid=np.zeros((cubedims[1],cubedims[2]))
        parcube=np.zeros((5,cubedims[1],cubedims[2]))
        uparcube=np.zeros((5,cubedims[1],cubedims[2]))
        for x in xrange(cubedims[1]):
            for y in xrange(cubedims[2]):
                spec = cube[:,x,y]
                if ferr is not None:
                    uspec = ecube[:,x,y]
                    u=uspec[~np.isnan(spec)]
                try: 
                    p,pcov=curvf(globals()[fn],
                             wvl[~np.isnan(spec)],spec[~np.isnan(spec)],sigma=u,
                             bounds = [[0.01,min(wvl),np.mean(wvl[1:]-wvl[:-1]),0],
                                        [50.,max(wvl),0.3,5.0]])
                    parcube[:,x,y]=p
                    uparcube[:,x,y]=np.sqrt(np.diag(pcov))
                    fitsqr[x,y],err=spi.quadrature(globals()[fn],
                          min(wvl),max(wvl),args=p)
                    resid[x,y]=fitsqr[x,y]-icube[x,y]
                except Exception:
                    fitsqr[x,y]=np.NaN
                    resid[x,y]=np.NaN
                    parcube[:,x,y]=np.NaN
                    uparcube[:,x,y]=np.NaN
        fit_dict={'p':parcube,'p_u':uparcube,'intfit':fitsqr,'resid':resid}
    if writeto is not None:
        if hdr is None:
            raise ValueError("I need a header to write to file")
        hdr['BUNIT']='W m^-2 sr^-1'
        newhdu = pf.PrimaryHDU(data=icube,header=hdr)
        if ferr is not None:
            ehdu = pf.ImageHDU(data=esqr)
            nhdu = pf.HDUList([newhdu, ehdu])
            nhdu.writeto(writeto,overwrite=True)
        else:
            newhdu.writeto(writeto,overwrite=True)
    if (intgauss is True or intlorentz is True) and (full_out is True):
        return icube,esqr,fit_dict
    else:
        return icube,esqr
