import numpy as np
import astropy.io.fits as pf
import astropy.convolution as acon
import pylab as plt
#from matplotlib.colors import LogNorm

def SAM(ff,fwhm=None,maskrad=None,preview=True,ovmin=None,ovmax=None,nvmin=None,nvmax=None):
    hdu = pf.open(ff)
    hdr = hdu[0].header
    odat = hdu[0].data

    cd = np.abs(hdr['CDELT1']) if 'CDELT' in hdr.tostring() else np.abs(hdr['CD1_1'])
    if fwhm is None:
	fwhm = hdr['BMAJ'] #if hdr['BMAJ']==hdr['BMIN'] else [hdr['BMAJ'],hdr['BMIN']]
    else:
	fwhm = cd*fwhm/3600.
    hdu.close()
    sigma = fwhm/(cd*np.sqrt(8.0*np.log(2.0)))
    print cd, fwhm, sigma
    kern = acon.Gaussian2DKernel(sigma) #if type(sigma) is float else ... figure it out later if need be

    omed = np.median(odat)
    om = np.mean(odat)
    ostd = np.std(odat)
    ndat = acon.convolve_fft(odat,kern,interpolate_nan=True)
    nstd = np.std(ndat)
    noiselvl = np.std(np.sort(np.ravel(odat))[:-len(np.ravel(odat))/10])
    print noiselvl
    #for images that SAM is good for, we can only be reasonably sure the top 10% or so is signal
    nm = np.mean(ndat)
    if maskrad is None:
	maskrad = 3.0
    samd=np.ma.masked_where(ndat<(nm + maskrad*nstd),odat)

    if preview==True:
	fig, (ax1, ax2) = plt.subplots(ncols=2)
	ovmin=0.0 if not ovmin else ovmin
	ovmax=odat.max() if not ovmax else ovmax
	nvmin=ndat.min() if not nvmin else nvmin
	nvmax=ndat.max() if not nvmax else nvmax
	plt.get_cmap().set_bad(color = 'w', alpha = 1.)
	ax1.pcolormesh(odat,vmin=ovmin,vmax=ovmax)
	ax1.set_title('original image')
	ax1.set_xlabel('x [pixels]')
	ax1.set_ylabel('y [pixels]')

	ax2.pcolormesh(samd,vmin=nvmin,vmax=nvmax)
	ax2.set_title('SAMed image')
	ax2.set_xlabel('x [pixels]')
	ax2.set_ylabel('y [pixels]')

	plt.show()
	return samd
    elif preview==False:
	hdr['history']='SAM: FWHM={:.4} deg, threshold={:.4} sigma above mean noise level'.format(fwhm,maskrad)
	return hdr,samd
	
