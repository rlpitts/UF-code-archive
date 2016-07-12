# -*- coding: utf-8 -*-
"""
Created on Wed Jul 01 22:49:06 2015

Just a dump heap for functions I needed to make my files usable

@author: Owner
"""
import glob, re, math, os, sys
import numpy as np
#import WISE2MJySr as wise_fcon
from astropy import wcs
from scipy.optimize import fsolve
from astropy.io import fits as pf
import coordconv as ccon
from matplotlib.path import Path

ctab = np.genfromtxt('NantenTab.dat',dtype=None,skip_header=2, usecols=(0,1,2,6))
#need to replace with comparison to coordinates of all BYF clumps
dirlist = ['R01','R02to03','R04','R05','R06','R07','R08','R09','R10','R11',
           'R12to15','R16to18','R17','R19','R20','R21','R22','R23','R24','R25',
           'R26','R27','byf01','byf11','byf38','byf57','byf123']
cdic={d:[] for d in dirlist}
for row in ctab:
    ndig = sum(c.isdigit() for c in str(row[3]))
    if row[0] in [1,11,38,57,123]:
        cdic['byf01' if row[0]==1 else 'byf{}'.format(row[0])].append([row[1],row[2]])
    else:
        if ndig==1:
            if any(dig in str(row[3]) for dig in ('2','3')):
                cdic['R02to03'].append([row[1],row[2]])
            else:
                cdic['R0{}'.format(row[3])].append([row[1],row[2]])
        elif ndig==2:
            if row[3] in ('12','13','14','15'):
                cdic['R12to15'].append([row[1],row[2]])
            elif row[3] in ('16','18'):
                cdic['R16to18'].append([row[1],row[2]])
            else:
                cdic['R{}'.format(row[3])].append([row[1],row[2]])


srclist = ['ATLASGAL','MIPS','GLIMPSE','WISE','MSX','Herschel']
dest = '/astro/data/rhocnc22/pitts/Regions/'
mlist = open('Regions/molist','r').read().splitlines()
mdic = {}
for m in mlist:
    hdu = pf.open(os.path.join('Regions',m))
    hdr = hdu[0].header
    hdr['NAXIS']=2
    del hdr['CRVAL3']
    del hdr['CRPIX3']
    del hdr['CDELT3']
    del hdr['CTYPE3']
    #^does not permanently delete axis 3 info, just removes them from working copy of hdr
    foot = wcs.WCS.calc_footprint(wcs.WCS(hdr))
    hdu.close()
    r = m.split('/')[-2]
    mdic[r] = Path(np.vstack((np.copy(foot),np.copy(foot)[0])),closed=True)
    #only close path AFTER appending a copy of the starting point to the end
    #.cleaned() is a poor substitute.
    
def symclobber(f1,f2):
    try:
        return os.symlink(f1,f2)
    except OSError:
        os.remove(f2)
        return os.symlink(f1,f2)

#use only if necessary - prefer CDELT, CROTA        
def rotmat(c11,c22,cr,hdr):
    if type(c11) is not float:
        float(c11)
    if type(c22) is not float:
        float(c22)
    if type(cr) is not float:
        float(cr)
    cd1_1=c11*math.cos(math.radians(cr))
    cd1_2=-c22*math.sin(math.radians(cr))
    cd2_1=c11*math.sin(math.radians(cr))
    cd2_2=c22*math.cos(math.radians(cr))
    hdr.set('CD1_1',cd1_1,'converted from CDELTi,CROTA2',after='CDELT2')
    hdr.set('CD1_2',cd1_2,'converted from CDELTi,CROTA2',after='CD1_1')
    hdr.set('CD2_1',cd2_1,'converted from CDELTi,CROTA2',after='CD1_2')
    hdr.set('CD2_2',cd2_2,'converted from CDELTi,CROTA2',after='CD2_1')
    return None

def fcc(cr1,cr2,ct1,ct2,foot):
#DO NOT CLOSE PATHS W/O COPYING STARTING POINT TO END OF VERTICES
    if 'GLON' in ct1 or 'GLAT' in ct2:
        box = Path(np.copy(foot)).cleaned()
        if len(box)==5:
            box=Path(np.vstack((box.vertices[:-1],box.vertices[0])),closed=True)
        l,b = cr1,cr2
    elif 'RA' in ct1 or 'DEC' in ct2:
        l,b = np.transpose(ccon.eq2gal(cr1,cr2))[0]
        a,d = np.transpose(np.copy(foot))
        verts=np.transpose(np.reshape(ccon.eq2gal(a,d),(2,4)))
        box = Path(np.vstack((verts,verts[0])),closed=True)
    return l,b,box

def setbeam(hd,fname,d): #must be run after sorting into correct directory
    fwhms = {'ATLASGAL':19.2, 'PLW':35.2, 'PMW':23.9, 'PSW':17.6, 'MG29':6.25, 'msx':18.3,
             'w4':12.0, 'w3':6.5, 'w2':6.4, 'w1':6.1, 'I4':1.98, 'I3':1.88,
             'I2':1.72, 'I1':1.66}
    if 'hpacs' in fname:
    	pacsbm = np.genfromtxt('Regions/resolns.tbl',dtype=None,names=True,missing_values='-',filling_values=0)
    	pacspa = np.genfromtxt('Regions/PACS_PAsByFile.tbl',dtype=None,names=True,missing_values='nan')
    	#go look at http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord.position_angle
        om = hd['INSTMODE'] if 'Para' in hd['INSTMODE'] else hd['SCANSP']
        if 'MAPR' in fname:
            band='R'
        elif 'MAPB' in fname:
            band='B'
        elif 'MAPG' in fname:
            band='G'
        try:
            x = [i for i in xrange(len(pacsbm)) if (pacsbm['band'][i]==band and pacsbm['mode'][i] in om)][0]
            y = 'pa_'+str(pacsbm['pa'][x]).replace('.','o')
	    #names including '.' are legal but may cause issues later; genfromtxt uses NameValidator()
        except IndexError:
            print fname, ': ', om, ' mode not recognized. Aborting.'
            sys.exit()
        try:
            z = [j for j in xrange(len(pacspa)) if pacspa['MopraMap'][j] in os.listdir('Regions/'+d)][0]
            pa = pacspa[y][z]
        except IndexError: #just assume fast scan/PacsSpireParallel mode - probably won't even do these regions
            if d=='R04':
                pa=26.7
            elif d=='R17':
                pa=44.0
            elif d=='R19' or d=='R20':
                pa=45.49
            elif d=='R22':
                pa=52.4
            elif d=='R24':
                pa=53.93
            elif d=='R25':
                pa=53.67
            elif d=='R27':
                pa=55.8
        hd.set('BMAJ',pacsbm['bmaj'][x]/3600,'FWHM in deg')
        hd.set('BMIN',pacsbm['bmin'][x]/3600,'FWHM in deg',after='BMAJ')
        hd.set('BPA',round(pa,1),'in (l,b) for {0}; was {1:0.1f} in fk5'.format(d,pacsbm['pa'][x]),after='BMIN')
	#numbers in {} [before ':'] are indices in args of format
	#(required if there are alphanumeric codes after ':' for at least 1 arg in .format)
        hd.set('PAUNIT','deg    ',after='BPA')
    else:
        for k,v in fwhms.items():
            if k in fname:
                if 'BMAJ' not in hd.keys():
                    hd.set('BMAJ',v/3600,'FWHM in deg')
                    hd.set('BMIN',v/3600,'FWHM in deg', after='BMAJ')
                    hd.set('BPA',0.0,'beam PA (always 0 if symmetric)', after='BMIN')
                    hd.set('PAUNIT','deg    ',after='BPA')
                elif hd['BMAJ']>0.012:
                    #replaces values not already in deg
                    hd.set('BMAJ',v/3600,'FWHM in deg')
                    hd.set('BMIN',v/3600,'FWHM in deg', after='BMAJ')
                if 'BPA' not in hd.keys(): #catch anything missed by above if-elif statement
                    hd.set('BPA',0.0,'beam PA (always 0 if symmetric)', after='BMIN')
                    hd.set('PAUNIT','deg    ',after='BPA')
    return 'bmaj: {0:0.8f}, bmin: {1:0.8f}, bpa:  {2:0.1f}'.format(hd['BMAJ'],hd['BMIN'],hd['BPA']) #headers are modified in-place regardless of what I do.

def flux_norm(hd,dat,omega):
    '''omega = angular area unit
    use sr to output in MJy/sr, and arcsec to output in Jy/arcsec^2'''
    #see http://irsa.ipac.caltech.edu/data/SPITZER/docs/spitzermission/missionoverview/spitzertelescopehandbook/19/
    # and maybe http://newton.cx/~peter/2011/12/reference-the-ultimate-resolved-source-cheatsheet/
    bunit = hd['BUNIT']
    print bunit
    if 'Jy/arcsec^2' in bunit:
        print 'flux already normalized'
        return hd,dat
    elif 'MJy/sr' in bunit:
        sf=2.350443*(10**(-5))
    elif any(s in bunit for s in ('Jy/pixel', 'Jy/pix', 'Jy/px')):
        if 'CDELT' in hd.tostring():
            x,y=hd['CDELT1'],hd['CDELT2'] #x,y in deg
        else: #x,y should be in deg by the end
            if hd['CD2_1']!=0:
                theta=math.atan(float(hd['CD2_1'])/float(hd['CD1_1']))
                x=math.degrees((float(hd['CD2_1'])-float(hd['CD1_1']))/(math.cos(theta)-math.sin(theta)))
                y=math.degrees((float(hd['CD1_2'])+float(hd['CD2_2']))/(math.cos(theta)-math.sin(theta)))
            else:
                x=float(hd['CD1_1'])
                y=float(hd['CD2_2'])
        pa = abs(float(x)*float(y))
        sf = 1./(pa*(3600**2))
    elif 'Jy/beam' in bunit: #MUST BE RUN AFTER SETBEAM()
        #after setbeam(), bmaj, bmin should be in degrees
        #Fix according to http://herschel.esac.esa.int/hcss-doc-13.0/load/spire_drg/html/ch06s09.html ?
        ba = (np.pi*float(hd['BMAJ'])*float(hd['BMIN']))/(4*math.log(2))
        sf = 1./(ba*(3600**2))
    elif 'DN' in bunit: #this is WISE
	dn_to_jy = {1:1.9350E-06,2:2.7048E-06,3:1.8326e-06,4:5.2269E-05} #DN are per pixel
	beam_area = {1:  6.08   * 5.60 * np.pi / (4*np.log(2)),
	    	     2:  6.84   * 6.12 * np.pi / (4*np.log(2)),
	    	     3:  7.36   * 6.08 * np.pi / (4*np.log(2)),
	    	     4:  11.99  * 11.65* np.pi / (4*np.log(2))} #sq arcsec
	pix_area = hd['PXSCAL1']*hd['PXSCAL2'] #1.375"/pixel at the center
	bmaj_min_pa = {1:(6.08,5.60,3),2:(6.84,6.12,15),3:(7.36,6.08,6),4:(11.99,11.65,0)}
	#bmaj/bmin in arcsec, bpa in deg
    	#PSFs in .csh code assume azimuthally averaged FWHMs, but that should be OK - the 2D PSFs are square anyway
	band = hd['BAND']
	hd.set('OMEGABM',beam_area[band],'(arcsec^2)')
	hd.set('DNTOJY',dn_to_jy[band],'from http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/sec2_3f.html')
	sf = dn_to_jy[band] / pix_area
    elif 'W/m^2-sr' in bunit: #this is MSX
	band = hd['BAND']
	Wm2_to_jy = {'A':7.133E+12,'C':2.863E+13,'D':3.216E+13,'E':2.476E+13}
	passband = {'A':'6.8 - 10.8','C':'11.1 - 13.2','D':'13.5 - 15.9','E':'18.2 - 25.1'}
	wl={'A':8.28E-06,'C':1.213E-05,'D':1.465E-05,'E':2.134E-05}
	hd.set('WAVELENG',wl[band],'central wavelength', after='BAND')
	hd.set('BANDFWHM',passband[band],'upper and lower FWHM limits', after='WAVELENG')
	hd.set('WM2TOJY',Wm2_to_jy[band],'W/m^2sr to Jy/sr conversion factor (see comments)',after='BUNIT')
	hd.set('COMMENT','radiance to flux conversion relies on isophotal assumption; adjust based on S_nu')
	sf = Wm2_to_jy[band]/4.254517E10 # 4.254517E10 as^2/sr
    else:
        raise ValueError('Unit not recognized')
        sys.exit()
    if 'sr' in omega:
        newdat=dat*sf/(2.350443*(10**(-5)))
        hd.set('BUNIT','MJy/sr','converted from {}'.format(bunit))
    elif 'arcsec' in omega:
        newdat=dat*sf
        hd.set('BUNIT','Jy/arcsec^2','converted from {}'.format(bunit))
    return hd,newdat

def make_hdu(dat,hdr,sname,fname,d):
    dat0=np.copy(dat)
    hdr0=hdr.copy()
    setbeam(hdr0,sname,d) #hdr0 modified in-place
    #d = k = key in mdic = region folder name; sname = source name
    hdr1,dat1 = flux_norm(hdr0,dat0,'arcsec')
    if (np.array(dat1)==np.array(dat0)).all() and 'wise' not in fname:
	raise IOError('unit conversion failed to execute')
	sys.exit()
    else:
        new_hdu = pf.PrimaryHDU(data = dat1,header = hdr1)
        pf.HDUList(new_hdu).writeto(fname,clobber=True)
        return 'New header & data written to {}'.format(fname)

def make_workfile(fi):
    hdu = pf.open(fi)
    hdr = hdu['image'].header.copy() if 'Herschel' in fi else hdu[0].header
    apath = os.path.abspath(fi)
    hdr['RAWFNAME']=apath #tracer
    if 'CDELT' in hdr.tostring() and 'CD1' not in hdr.tostring():
        rotmat(hdr['CDELT1'],hdr['CDELT2'],hdr['CROTA2'] if 'CROTA2' in hdr else 0.000,hdr)
        #add rotation matrix just in case anything requires it
        #prefer CDELTs, but the full matrix did come up somewhere...
    foot = wcs.WCS.calc_footprint(wcs.WCS(hdr.tostring() if 'Herschel' in fi else fi))
    c1,c2 = float(hdr['CRVAL1']),float(hdr['CRVAL2'])
    ct1,ct2 = hdr['CTYPE1'],hdr['CTYPE2']
    #get dimensions for search box, convert crvals and box edges to galactic coords.
    l,b,box = fcc(c1,c2,ct1,ct2,foot)
    imdat = hdu['image'].data if 'Herschel' in fi else hdu[0].data
    if 'Herschel' in fi:
        phdr = hdu[0].header #not used outside this loop
        hdr.set('INSTRUME',phdr['INSTRUME'],'Instrument attached to this product',after='META_0')
        hdr.set('CUSMODE',phdr['CUSMODE'],'Common Uplink System observation mode',after='INSTRUME')
        hdr.set('INSTMODE',phdr['INSTMODE'],'Instrument Mode',after='INSTRUME')
        if 'hpacs' in fi and 'Para' not in phdr['INSTMODE']:
            hdr.set('SCANSP',phdr['SCANSP'],phdr.comments['SCANSP'],after='INSTMODE')
            hdr.set('SCANVELO',phdr['SCANVELO'],phdr.comments['SCANVELO'],after='INSTMODE')
        try:
            edat = hdu['error'].data
        except KeyError:
            try:
                edat = hdu['stDev'].data
            except KeyError:
                edat = None
    else:
        edat=None
    hdu.close() #copies of hdr & data now in memory, don't need file open anymore
    seq = fi.rsplit('/')
    froot='{}/workfiles/'.format(seq[0])
    used=0 #toggles work-file replacement off after 1st time
    klist=[]
    for k in cdic.iterkeys():
        if any(box.contains_points(cdic[k])):
            lroot='{}/{}/'.format(dest,k)
            #Now working only on files in 'workfiles' subdirectory of arg
            if 'Herschel' in apath:
                cmap = [('level','lvl'),('1342',''),('HPPPMAP','php'),('HPPJSMAP','jsm'),('HPPUNIMAP','uni'),
                        ('extd','x'),('psrc','p'),('PLW','R'),('PMW','G'),('PSW','B')]
                nm = '_'.join(['hpacs' if 'pacs' in fi else 'hspire',seq[1][3:],'id'+seq[-4],seq[-2]+'.fits'])
                for old,new in cmap:
                    nm = nm.replace(old,new) if old in nm else nm
                erfnm = os.path.join(froot,nm[:-5]+'err.fits')
                erlnm = os.path.join(lroot,nm[:-5]+'err.fits')
                if not os.path.isfile(erfnm) or used==0:
                    make_hdu(edat,hdr,fi,erfnm,k)
                symclobber(os.path.abspath(erfnm), erlnm)
            elif 'WISE' in apath:
                nm='wise'+seq[-1]
            elif 'GLM' in apath or 'VELACAR' in apath:
                nm=seq[-1][:-5]+'_0.6'+seq[-1][-5:] if '0.6' in apath else seq[-1][:-5]+'_1.2'+seq[-1][-5:]
            else:
                nm=seq[-1]
            fullnm = os.path.join(froot,nm) #technically not a full name, just full relative path from cwd
            lnm = os.path.join(lroot,nm)
            if os.path.isfile(fullnm)==True and used>0:
                print 'work file already made'
            else:
                make_hdu(imdat,hdr,fi,fullnm,k)
                used+=1
            symclobber(os.path.abspath(fullnm), lnm)
            klist.append(k)
    if len(klist)==0:
        unused=fi
    else:
        print '{} found in region(s) {}'.format(fi,klist)
        unused=None
    return unused

print 'Choose source folder {} or type "esc" to quit _'.format(srclist)
def mover(arg):
    unused=[]
    if arg == 'esc':
        return 'Quitting'
    elif arg not in srclist:
        return 'Please retry with a valid option.'
    else:
        if arg == 'ATLASGAL':
            for nm in glob.glob('ATLASGAL/ATLASGAL*.fits'):
                u = make_workfile(nm)
                if u is not None:
                    unused.append(u)

        elif arg == 'Herschel':
            for nm in glob.glob('Herschel/lvl2-*/*/level*/*/*.fits.gz'):
                #readymade=os.listdir('Herschel/workfiles/')
                #seq=nm.rsplit('/')
                #cmap = [('level','lvl'),('HPPPMAP',''),('HPPJSMAP',''),('HPPUNIMAP',''),
                #        ('extd','x'),('prsc','p'),('PLW','R'),('PMW','G'),('PSW','B')]
                #newnm = '_'.join(['hpacs' if 'pacs' in nm else 'hspire',seq[1][3:],'id'+seq[-4][-5:],seq[-2]+'.fits'])
                #for old,new in cmap:
                #    newnm = newnm.replace(old,new) if old in newnm else newnm
                if 'diag' in nm or 'browse' in nm:
                    continue
                #elif newnm in readymade:
                #    print newnm, ' exists; skipping'
                #    continue
                else:
                    u = make_workfile(nm)
                if u is not None:
                    unused.append(u)

        elif arg == 'WISE':
            for nm in glob.glob('WISE/*.fits'):
                u = make_workfile(nm)
                if u is not None:
                    unused.append(u)

        elif arg == 'MSX':
            for nm in glob.glob('MSX/*.fits'):
                u = make_workfile(nm)
                if u is not None:
                    unused.append(u)

        elif arg == 'MIPS':
            for nm in glob.glob('MIPS/*.fits'):
                if 'cube' in nm:
                    continue
                u = make_workfile(nm)
                if u is not None:
                    unused.append(u)

        elif arg == 'GLIMPSE':
            for root, dirs, files in os.walk('GLIMPSE'):
                for nm in files:
                    if nm.endswith('.fits'):
                        u = make_workfile(os.path.join(root, nm))
                        if u is not None:
                            unused.append(u)
    return 'Soft-linking complete. The following {} files were unused: {}'.format(len(unused),unused)

def addbeam():
    pacslist=glob.glob('Herschel/workfiles/hpacs*.fits')
    pacsbm = np.genfromtxt('Regions/resolns.tbl',dtype=None,names=True,missing_values='-',filling_values=0)
    patab = np.genfromtxt('Regions/PAsByRegion.tbl',dtype=None,names=True,missing_values='nan')
    for f in pacslist:#just add BPAs under HISTORY cards with "ADDBPA.PY:" prefix
        trunk=f.rsplit('/')[-1]
        print 'working on ',trunk
        hdu=pf.open(f,mode='update')
        hdr=hdu[0].header
        rlist=[r.rsplit('/')[-2] for r in glob.glob('Regions/*/{}'.format(trunk))]
        bpalist=[]
        for d in rlist:
            om = hdr['INSTMODE'] if 'Para' in hdr['INSTMODE'] else hdr['SCANSP']
            band=trunk.split('.')[0][-4] if 'err' in trunk else trunk.split('.')[0][-1]
            x = [i for i in xrange(len(pacsbm)) if (pacsbm['band'][i]==band and pacsbm['mode'][i] in om)][0]
            y = 'eq_'+str(pacsbm['pa'][x]).replace('.','o')
            z = [j for j in xrange(len(patab)) if patab['Region'][j]==d][0]
            pa = patab[y][z]
            bpalist+=[float(pa)]
            print 'ADDBEAM.PY: BPA = {0:05.1f} in (l,b) for {1}'.format(pa,d)
            hdr.set('HISTORY','ADDBEAM.PY: BPA = {0:05.1f} in (l,b) for {1}'.format(pa,d))
            #:05.1f = single-decimal-place float with left-hand 0 padding to 5 total characters
        avgpa=np.median(bpalist)
        hdr.set('BPA',round(pa,1),'median BPA in (l,b); was {0:0.1f} in fk5'.format(pacsbm['pa'][x]))
        if f==pacslist[0]:
            print hdr
            key=raw_input('everything look OK? If not, type "x" to kill_')
            if key=='x':
                sys.exit()
        hdu.flush()
        hdu.close()
    return None #and done

def msx_errmap(fi):
    band=fi.split('msxs3')[1][0]
    hdu1=pf.open(fi)
    nfi = fi.split('-')[0]+'err-'+fi.split('-')[1]
    edat = hdu1[0].data*0.17 if band=='E' else hdu1[0].data*0.16
    hdu2=pf.writeto(nfi,edat,hdu1[0].header, clobber=True)
    hdu1.close()
    return "Naive MSX error map written to ", nfi
#Forget this, it's making everything worse.
#def wise_offset(fi):
#    band = fi.split('wise')[1][:2]
#    r = fi.split('/')[1]
#    fs = glob.glob('Regions/{}/wise*{}*int*.fits'.format(r,band))
#    if len(fs)>1:
#        meds = np.zeros(len(fs))
#        for i,f in enumerate(fs):
#            hdu = pf.open(os.path.realpath(f))
#            # np.sort is much faster than sorted()
#            decile = np.sort(hdu[0].data, axis=None)[0:len(hdu[0].data)/10]
#            # np.sort w/o axis kwarg should flatten array before sorting
#            # guess this version won't unless I specify axis=None
#            meds[i] = np.median(decile)
#            hdu.close()
#        med = np.mean(meds) #np.average allows weights, but don't need that.
#    elif len(fs)==1:
#        hdu = pf.open(os.path.realpath(fs[0]))
#        decile = np.sort(hdu[0].data, axis=None)[0:len(hdu[0].data)/10]
#        med = np.median(decile)
#        hdu.close()
#    else:
#        raise ValueError('files not found')
#    #something's making difference too negative
#    hdu1=pf.open(fi)
#    hdu1[0].header.set('OFFSET',med,'estimated ZL correction factor')
#    nfi = fi[:fi.find('.mr')]+'-ZL'+fi[fi.find('.mr'):]
#    hdu2=pf.writeto(nfi,hdu1[0].data - med,hdu1[0].header, clobber=True)
#    hdu1.close()
#    return 'ZL-corrected data written to {}'.format(nfi)


        
            
    

