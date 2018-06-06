#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 19:21:50 2018

@author: rlpitts
"""

import numpy as np
import astropy.io.fits as pf
from astropy.wcs import WCS
#import matplotlib as mpl
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

oi63=pf.open('workfiles/F0313_FI_IFS_04006131_BLU_WXY_OI63um_L4.fits')
oi145=pf.open('workfiles/F0313_FI_IFS_04006132_RED_WXY_OI145um_L4.fits')
oiii=pf.open('workfiles/F0313_FI_IFS_04006132_BLU_WXY_OIII88um_L4.fits')
cii=pf.open('workfiles/F0313_FI_IFS_04006131_RED_WXY_CII158um_L4.fits')

lw={'oi63':63.185,'oi145':145.53,'oiii':88.36,'cii':157.74}
flxs = {'oi63':oi63[3].data,'oi145':oi145[1].data,'oiii':oiii[1].data,
        'cii':cii[1].data}
errs = {'oi63':oi63[4].data,'oi145':oi145[2].data,'oiii':oiii[2].data,
        'cii':cii[2].data}
wlas = {'oi63':oi63[5].data,'oi145':oi145[5].data,'oiii':oiii[5].data,
        'cii':cii[5].data}

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
        if hdu[i].data is not None:
            if len(hdu)>=10:
                hdr=hdu[i].header
            elif len(hdu)==3:
                hdr=mainhd
                try:
                    del hdr['TRACMODE']
                    del hdr['TRACERR']
                    del hdr['NODDING']
                    del hdr['CHOPPING']
                except KeyError:
                    pass
#                if '146' in fname:
#                    del hdr[141:145]
#                else:
#                    del hdr[90:94]
                if 'XTENSION' not in hdr:
                    hdr.set('XTENSION','IMAGE', 'IMAGE extension',after='BITPIX')                          
                if i==1:
                    hdr.set('EXTNAME','ERROR',after='XTENSION')
                elif i==2:
                    hdr.set('EXTNAME','FLAG',after='XTENSION')
                
            hdr.set('CHANNEL',mainhd['CHANNEL'],'Detector channel',after='EXTNAME')
            if hdr['CHANNEL'] in ['BLUE','BLU','blue','blu']:
                hdr.set('BMAJ',6/3600.,'beam major axis (deg)',after='CROTA2')
                hdr.set('BMIN',6/3600.,'beam minor axis (deg)',after='BMAJ')
            elif hdr['CHANNEL'] in ['RED','red']:
                hdr.set('BMAJ',12/3600.,'beam major axis (deg)',after='CROTA2')
                hdr.set('BMIN',12/3600.,'beam minor axis (deg)',after='BMAJ') 
                
            if 'OI63' in fname:
                obase='fifils_OI63_'
            elif 'OI145' in fname:
                obase='fifils_OI145_'
            elif 'OI146' in fname:
                obase='fifils_OI146_'              
            elif 'OIII88' in fname:
                obase='fifils_OIII88_'
            elif 'CII158' in fname:
                obase='fifils_CII158_'
#            try:    
            if 'Cont' in fname:
                obase=obase+'cont_'
            elif 'Flux' in fname:
                obase=obase+'line_'
            oname = 'workfiles/'+obase+hdr['EXTNAME']+'.fits'
#            except KeyError:
#                if i==1:
#                    oname = 'workfiles/'+obase+'ERR.fits'
#                elif i==2:
#                    oname = 'workfiles/'+obase+'flag.fits'
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
        specl=lorentz(wl,l[0],l[1],abs(l[2]),l[3],l[4])
        fwhm=2*abs(l[2])
        pl.plot(wl,specl,'r-',label='FWHM={:.3}'.format(fwhm)+'$\mu$m')

    pl.legend(loc=0)
    pl.show()
    return None

def splitflux(fden, wavelen=None,ferr=None,writeto=None,hdr=None,
            fitfn='gauss',clobber=True):
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
    print lc
    if fitfn not in ['gauss','lorentz']:
        raise IOError("Fitting function name must be 'lorentz' or 'gauss'")
    cubedims=np.shape(cube)
    frames=np.array([(n,w) for n,w in enumerate(wvl) if (lc-0.2<w<lc+0.2)])
    #3 sigma is about 27% larger than 2 fwhm
    inds,fwls=np.transpose(frames)
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
                             bounds = [[0.01,lc-0.05,np.mean(wvl[1:]-wvl[:-1]),-10,0.],
                                       [50.,lc+0.05,0.4,10,10.0]])
                fwhm = 2*abs(p[2]) if fitfn=='lorentz' else p[2]*np.sqrt(8*np.log(2))
                #goodlc = lc-0.05<p[1]<lc+0.05 #why is this True here...
                goodcov = np.inf not in pcov
                print '1st try:  ', fwhm,goodcov#, goodlc
            except (ValueError,RuntimeError,NameError) as e:
                print e
                pass
            print ('fwhm' in locals() and locals()['fwhm']<0.16) 
            #print ('goodlc' in globals() and globals()['goodlc'] is True) #...but False here?
            print ('goodcov' in locals() and locals()['goodcov'] is True)
            
            if (('fwhm' in locals() and locals()['fwhm']<0.16) and 
                #('goodlc' in globals() and globals()['goodlc'] is True) and
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
                    w=wvl[~inds.astype(int)]
                    q,qcov=curvf(lreg,w[~np.isnan(s)],s[~np.isnan(s)],
                                 sigma=u[~np.isnan(s)],bounds = [[-10,0.],[10,10.0]])
                    r,rcov=curvf(globals()[fitfn], fwls,
                                 spec[inds.astype(int)],
                                 sigma=uspec[inds.astype(int)],
                                 bounds = [[0.01,lc-0.05,np.mean(wvl[1:]-wvl[:-1]),-10,0.],
                                            [50.,lc+0.05,0.4,10,10.0]])
                    rfwhm = 2*abs(p[2]) if fitfn=='lorentz' else p[2]*np.sqrt(8*np.log(2))
                    #goodlc2 = lc-0.05<r[1]<lc+0.05
                    goodcov2 = np.inf not in pcov
                    print '2nd try:  ', rfwhm, goodcov2#, goodlc2
                except (ValueError,RuntimeError,NameError) as e2:
                    print e2
                print ('rfwhm' in locals() and locals()['rfwhm']<0.16)
                #print ('goodlc2' in globals() and globals()['goodlc2'] is True)
                print ('goodcov2' in locals() and locals()['goodcov2'] is True)
                
                if (('rfwhm' in locals() and locals()['rfwhm']<0.16) and 
                    #('goodlc2' in globals() and globals()['goodlc2'] is True) and
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
                    print 'line fit failed, attempting pure continuum fit'
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
    lhdr['CRVAL3']=fwls[0]
    lhdu1 = pf.PrimaryHDU(data=lincube,header=lhdr)
    lhdu2 = pf.ImageHDU(data=elincube,name='error')
    lhdu3 = pf.ImageHDU(data=fdencube,name='flux_density_fit')
    lhdu4 = pf.ImageHDU(data=fduncube,name='flux_density_err')
    lhdu5 = pf.ImageHDU(data=frames,name='frames')
    lhdu0 = pf.HDUList([lhdu1,lhdu2,lhdu3,lhdu4,lhdu5])
    lhdu0.writeto(writeto[0] if type(writeto) is list else 'workfiles/fifils-'+ltxt+'-lineonly.fits',
                  overwrite=clobber)
    chdr=hdr.copy()
    chdr.set('EMTYPE','Continuum emission',before='ION')
    chdr.append('HISTORY',hdr['ION']+' line emission removed with sofia_utils.splitflux()')
    chdu1 = pf.PrimaryHDU(data=concube,header=lhdr)
    chdu2 = pf.ImageHDU(data=econcube,name='error')
    chdu3 = pf.ImageHDU(data=spindex,name='spectral_index')
    chdu4 = pf.ImageHDU(data=spundex,name='spectral_index_error')    
    chdu0 = pf.HDUList([chdu1,chdu2,chdu3,chdu4])
    chdu0.writeto(writeto[1] if type(writeto) is list else 'workfiles/fifils-'+ltxt+'-continuum.fits',
                  overwrite=clobber)
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
        wvl=pf.open(wavelen)[0].data if type(wavelen) is str else wavelen
    else:
        raise ValueError("If wavelength array isn't defined by now, you goofed")
    cubedims = np.shape(cube)
    freq=np.flip((2.99792458*10**8)/wvl,axis=0)
    icube = getattr(spi,method)(np.ma.masked_invalid(np.flip(cube,axis=0)),
                    x=freq,axis=0)
    if ferr is not None:
        esqr = np.sqrt(np.nansum(np.flip(ecube,axis=0)**2,axis=0))
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

pdrwebs={'ts':pf.open('workfiles/tsweb.fits')[0],
         'cii':pf.open('workfiles/cpweb.fits')[0],
         'fir':pf.open('workfiles/firweb.fits')[0],
         'oioi':pf.open('workfiles/oioiweb.fits')[0], #[O I] 145 Micron/[O I] 63 Micron
         'o145cii':pf.open('workfiles/o145ciiweb.fits')[0],
         'o63cii':pf.open('workfiles/oicpweb.fits')[0]}
#def matchlines(rats,urats,smaps,smhds):
    
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def match_lines3d(arr1,arr2,arr3, errmaps=None, fliprats=[0,0,0],
                  scanmaps=[pdrwebs['oioi'],pdrwebs['o145cii'],pdrwebs['o145cii']]):
    #add coversion to/from cgs
    if type(arr1) is str:
        hdu=pf.open(arr1)
        if len(hdu)>1:
            im1=hdu[0].data
            uim1=hdu[1].data
            hdr1=hdu[0].header
        else:
            im1=hdu[0].data
            hdr1=hdu[0].header
            uim1=None
    else:
        im1=arr1
        uim1=None
    if type(arr2) is str:
        hdu2=pf.open(arr2)
        if len(hdu2)>1:
            im2=hdu2[0].data
            uim2=hdu2[1].data
            hdr2=hdu2[0].header
        else:
            im2=hdu2[0].data
            hdr2=hdu2[0].header
            uim2=None
    else:
        im2=arr2
        uim2=None
    if type(arr3) is str:
        hdu3=pf.open(arr3)
        if len(hdu3)>1:
            im3=hdu3[0].data
            uim3=hdu3[1].data
            hdr3=hdu3[0].header
        else:
            im3=hdu2[0].data
            hdr3=hdu2[0].header
            uim3=None
    else:
        im3=arr3
        uim3=None
        
    if errmaps is not None:
        uim1=pf.open(errmaps[0])[0].data if uim1 is None else uim1
        uim2=pf.open(errmaps[1])[0].data if uim2 is None else uim2
        uim3=pf.open(errmaps[1])[0].data if uim3 is None else uim3
    #this will raise an index error if they're not all the same size
    if 'hdr1' in globals() and 'hdr2' in globals() and 'hdr3' in globals():
        if not hdr1['CDELT2']==hdr2['CDELT2']==hdr3['CDELT2']:
            raise IndexError('disimilar coordinate systems.'+
                             ' Check header[CDELTi] cards & regrid where needed')
    #MUST REGRID EVERYTHING TO CII MAPS - it won't work otherwise
    ratio1 = im2/im1 if fliprats[0]==0 else im1/im2 #o145/cii
    ratio2 = im3/im1 if fliprats[1]==0 else im1/im3 #o63/cii
    ratio3 = im2/im3 if fliprats[2]==0 else im3/im2 #o145/o63
    urat1=ratio1*np.sqrt((uim1/im1)**2+(uim2/im2)**2)
    urat2=ratio1*np.sqrt((uim1/im1)**2+(uim3/im3)**2)
    urat3=ratio1*np.sqrt((uim2/im2)**2+(uim3/im3)**2)
    #uratiof = np.sqrt((uim1/im1)**2+(uim2/im2)**2+(uim3/im3)**2)
    
    web1,web2,web3 = [sm.data for sm in scanmaps]
    #whd1,whd2,whd3 = [sm.header for sm in scanmaps]
    w1,w2,w3 = [WCS(sm) for sm in scanmaps]
    
    ratshape = np.shape(ratio1)
    narr = np.zeros((ratshape[0],ratshape[1],2))
    garr = np.zeros((ratshape[0],ratshape[1],2))
    for x in xrange(ratshape[0]):
        for y in xrange(ratshape[1]):
            inds1=np.where(np.logical_and(web1>ratio1[x,y]-urat1[x,y],
                                          web1<ratio1[x,y]+urat1[x,y]))
            inds2=np.where(np.logical_and(web2>ratio2[x,y]-urat2[x,y],
                                          web2<ratio2[x,y]+urat2[x,y]))
            inds3=np.where(np.logical_and(web3>ratio3[x,y]-urat3[x,y],
                                          web3<ratio3[x,y]+urat3[x,y]))
            #remember, switch x & y and add 1 to get the proper coordinates
            #c = np.concatenate((web1[inds1],web2[inds2],web3[inds3]))
            n1,g01 = w1.all_pix2world(inds1[0],inds1[1],1)
            n2,g02 = w2.all_pix2world(inds2[0],inds2[1],1)
            n3,g03 = w3.all_pix2world(inds3[0],inds3[1],1)
            logn = set.intersection(set(n1),set(n2),set(n3))
            logg0 = set.intersection(set(g01),set(g02),set(g03))
            narr[x,y,0],narr[x,y,1] = np.nanmin(logn),np.nanmax(logn)
            garr[x,y,0],narr[x,y,1] = np.nanmin(logg0),np.nanmax(logg0)
    return narr,garr 