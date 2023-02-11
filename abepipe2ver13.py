# -*- coding: utf-8 -*-
"""
Created on Sun Feb 15 23:47:17 2015

@author: Owner
"""

from astropy.io import fits as pf
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as pl
import os, sys, time#, gc
from collections import defaultdict
#from scipy.optimize import curve_fit
#import slidewin as slide #must have in cwd
#import trashdump

#gc.set_debug(gc.DEBUG_LEAK)
#************************************************
def wtd_stdev(dat,wts):
    '''for use in the next function'''
    av=np.average(dat,weights=wts)
    var=np.average((dat-av)**2, weights=wts)
    return np.sqrt(var)
def noise_scrambler(a,ea,fa):
    '''
    function to rearrange noise in a spectrum with the following parameters:
    a = (dtype=array, 1D) flux array
    ea = (dtype=array, 1D) error array
    fa = (dtype=array, 1D) flag array
        
    Returns:
    ra = array of fluxes with noise added
    '''
    ra = np.zeros(np.shape(a))
    wt=np.ma.masked_invalid(1/ea).filled(0)
    netmask=reduce(np.logical_or,[np.ma.getmask(np.ma.masked_outside(fa,-1,0)),
                                  np.ma.getmask(np.ma.masked_invalid(ea)),
                                  np.ma.getmask(np.ma.masked_outside(ea,0,2*wtd_stdev(a,wt))),
                                  np.ma.getmask(np.ma.masked_less_equal(a,0.))])
    errs=np.ma.array(ea,mask=netmask)
    rms=errs.filled(np.ma.median(errs[~errs.mask]))
    ra = [np.random.normal(q,rms[p]) for p,q in enumerate(a)]
    return ra

def detrend(indepv,depv):
    '''
    Removes effects of residual sky absorption before cross-correlating
    (See section 2.3.2)
    '''
#    try:
#        cy=depv.compressed()
#        cx=indepv.compressed()
#    except AttributeError:
#        pass
    if np.ma.is_masked(depv)==True:
        cy=depv.compressed()
    else:
        cy=depv
    if np.ma.is_masked(indepv)==True:
        cx=indepv.compressed()
    else:
        cx=indepv
    p = np.polyfit(cx,cy,2)
    flat=p[0]*cx**2+p[1]*cx+p[2]
    return cy-flat

def pointsum(d):
    '''Co-adds cross-correlations'''
    tups=[]
    for v, A in d.iteritems():
        tups+=[[int(v),sum(A)]]
    return sorted(tups)
    
def gauss_fit(X,U,S,A,offset):
    return A*np.exp(-0.5*((X-U)/S)**2) - offset
    
#def logit(X,M,S,A,B):
#    z = (X-M)/S
#    return A*np.exp(-z)/(S*(1+np.exp(-z))**2) - B
#************************************************
#Prologue: Make dictionary of stars, epochs, and fluxes:
start_time = time.time()
flist = []
#vlist = []
glist = []
epsqr = dict()
clrlist=['r','g','b']
klist = ['hdr','lmda','flux','fluxerr','flag','qual','dvbr11','bc']
hduno = [0,4,1,2,3]
#keep unpack=False for the following 2 matrices
dvpseps=np.genfromtxt('Br11PeakSeps.dat', dtype=None, delimiter='\t')
#probably won't ever need the 3rd column, but it's there just in case
br11lines=np.genfromtxt('BrackettSeries.dat', dtype=None, usecols=(0,1))
#at least 10 stars have epochs with ruined B-bands
# --> skip until I can isolate the individual bad epochs
skip=['2M20162816+3703229',
      '2M18211606-1301256',
      '2M03464087+3217247',
      '2M06014161+2224036',
      '2M06154017+0603582',
      '2M18043735+0155085',
      '2M18072725+2506165',
      '2M19120326+0237212',
      '2M20234436+3728351',
      '2M21103095+4741321',
      '2M03292627+4656162'] #<--older data wasn't ruined (as badly)
#test stars only:
rdvpseps=[['2M03221892+7846028',5],
          ['2M07103563+0626568',5],
          ['2M14030949+2935508',50],
          ['2M06355721+0519000',210],
          ['2M07375513+2116094',5],
          ['2M03534224+8002056',5],
          ['2M06362588+0604595',200],
          ['2M07441014+2205066',5]]

c = 2.99792*10**5 #in km/s
dvmin=c/22500

for file in os.listdir('ABeNewMEFs'):
    if file.startswith("2M") and file.endswith('fit') and file[:18] not in skip:
        meff = pf.open('ABeNewMEFs/'+file, mode='readonly', ignore_missing_end=True, memmap=False)
#for file in os.listdir('RVbinaries'):
#    if file.startswith("r2M14030949") and file.endswith('fit') and file[:18] not in skip:
#        meff = pf.open('RVbinaries/'+file, mode='readonly', ignore_missing_end=True, memmap=False)
        hdr = meff[0].header
        oid = str(hdr["OBJID"])
        if "HJD" in hdr:
            if oid in epsqr:
                ##epsqr already has key oid
                epsqr[oid].append(float(hdr["HJD"]))
            else:
                ##epsqr does not have key oid
                epsqr[oid] = [float(hdr["HJD"])]
        else:
            ##I assume the try/except block for the KeyError was to catch HJD as not a key in hdr
            meff.close()
            continue
        mjd = int(hdr["MJD5"]) #int for keys, for consistency in memory
        try:#qual for picking best spectrum as template
            qual=float(hdr["SNR"])
        except KeyError:#sometimes SNR doesn't exist; use 1/chisq as proxy
            qual=1/float(hdr["CHISQ"])
        try:
            dvp=dvpseps[np.where([dvpseps[nm][0]==oid for nm in range(0,len(dvpseps))])][0][1]
        except IndexError:
            dvp=rdvpseps[np.where([rdvpseps[nm][0]==oid for nm in range(0,len(rdvpseps))])[0]][1]
        bc=float(hdr["BC"])        
        for ind,no in enumerate(hduno):
            if ind==0:
                flist += [[oid, int(mjd), klist[ind], meff[no].header]]
            else:
                flist += [[oid, int(mjd), klist[ind], meff[no].data]]
        flist += [[oid, int(mjd), klist[-3], qual],
                  [oid, int(mjd), klist[-2], dvp],
                  [oid, int(mjd), klist[-1], bc]]
        meff.close()
        #gc.collect()
dcube = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for o,m,g,itsfullofstars in flist:
    dcube[o][m][g].append(itsfullofstars) #:3
dcube.default_factory=None
print '--- Prologue took {} sec ---'.format(time.time() - start_time)
start_time = time.time() #reset
#sys.exit()    
#***************************************
#Main loop:
for star in dcube.iterkeys():
    tot_time = time.time() #to calculate time spent on each star
    epochs=np.array(sorted(dcube[star].iterkeys()))
    if len(epochs) < 3:
        continue
    print star, len(epochs)
    tser=dict.fromkeys(epochs)
    psep=dcube[star][epochs[0]]['dvbr11'][0]
    if psep==-1:
        psep=270 #might be better not to attempt stars with undefined pseps
        #~avg. psep for ABE stars (not including sigma Ori E candidates)
    print psep
    #1. Build up a mask of brackett 11-20 & NaNs; @ the end, need to
    # combine them so only points unmasked for all epochs are
    # cross-correlated
    for ep in epochs:
        lmda=dcube[star][ep]['lmda'][0]
        flux=dcube[star][ep]['flux'][0]
        ferr=dcube[star][ep]['fluxerr'][0]
        flags=dcube[star][ep]['flag'][0]
        if len(lmda[0])<4096: #up-sample the low-res images
            newx=np.arange(0,4096)
            oldx=np.arange(0,4096,2)
            lmda1=[np.interp(newx,oldx,lmda[j]) for j in [0,1,2]]
            flux1=[np.interp(newx,oldx,flux[j]) for j in [0,1,2]]
            ferr1=[np.interp(newx,oldx,ferr[j]) for j in [0,1,2]]
            flag1=[np.floor(np.interp(newx,oldx,flags[j])) for j in [0,1,2]]
            lmda=lmda1
            flux=flux1
            ferr=ferr1
            flags=flag1
        hw=psep*1.5 #cover the Brackett line 504gs
        masklims = []
        masklims = [masklims+[ln[1]/(1+hw/c),ln[1]/(1-hw/c)] for ln in br11lines]
        #maxrv=1130 for all known HMXBs -> ~120A for the entire spectrum
        #in 2-tuples of blue limit, red limit
        brmaskset=[np.ma.masked_inside(lmda, red, blu) for [red,blu] in masklims]
        brmask=reduce(np.logical_or, [np.ma.getmask(m) for m in brmaskset])
        endmask=[[True]*32+[False]*4032+[True]*32]*3
        nanmask=np.ma.getmask(np.ma.masked_invalid(flux))
        fullmask=reduce(np.logical_or, [brmask,endmask,nanmask])
        invalmask=reduce(np.logical_or, [endmask,nanmask])
        #Based on the typical width of the Br11 masked segments, I choose 504 px as
        # my window size since it's the closest number that divides evenly into 4032
        tser[ep]={'x':np.ma.asanyarray(lmda),'mask':fullmask,'mask_nobr':invalmask,'barc':dcube[star][ep]['bc'][0],
                  'fx':[np.ma.asanyarray(flux),np.ma.asanyarray(ferr),np.ma.asanyarray(flags)]}       
    #2. Prep for cross-correlation
        #2.1. Merge masks
    unimask=reduce(np.logical_or, [tser[ep]['mask'] for ep in epochs])
    ubmask=reduce(np.logical_or, [tser[ep]['mask_nobr'] for ep in epochs])
    ums={'unimask':unimask, 'ubmask':ubmask}
        #2.2. Choose epoch with highest S/N; pair everything else off
    dummy=[(ep, dcube[star][ep]['qual'][0]) for ep in epochs]
    bestep=epochs[np.where(epochs==dummy[np.argmax([tup[1] for tup in dummy])][0])][0]
    pairlist = zip([bestep]*len(epochs[:-1]),np.delete(epochs, np.where(epochs==bestep)[0][0]))
    print '--- Masking and pairing took {} sec ---'.format(time.time() - start_time)
    start_time = time.time() #reset
    try:
        bhep = [hjd for hjd in epsqr[star] if str(bestep) in str(hjd)][0]
    except IndexError:
        try:            
            bhep = [hjd for hjd in epsqr[star] if str(bestep+1) in str(hjd)][0]
        except IndexError:
            bhep = [hjd for hjd in epsqr[star] if str(bestep-1) in str(hjd)][0]
        #2.3. Create output files
#    ofi=open('RVbinaries/{}.out'.format(star),'w')
#    ofi.write('#Star ID: {}\n'.format(star))
#    ofi.write('#If Br11 listed as -1, 270 km/s average was used instead\n')
#    ofi.write('#RVerr = standard deviation of RV shifts in noisy Monte-Carlo-generated spectra\n')
#    ofi.write('#Ref Epoch (RHJD): {}, BC: {}\n'.format(bhep,tser[bestep]['barc']))
#    ofi.write('#Epoch    BC    Br11psep    RV    RVerr    BrMasked\n')
#    ofi.write('#(RHJD)    (km/s)    (km/s)    (km/s)    (km/s)    (Y/N)\n')
#    ofi2=open('RVbinaries/{}_noisy.out'.format(star),'w')
#    ofi2.write('#Monte-Carlo noise-added cross correlations\n')
    #3. Cross-correlate only points that appear in all epochs
    nits = 0 #control number of noise_scrambler sims
    for ep in np.delete(epochs, np.where(epochs==bestep)[0][0]):
        bc0=tser[bestep]['barc']
        bc1=tser[ep]['barc']
#        bc0,bc1=0,0
        xcors=None
        xcors_B=None
        Nxcors=[]
        Nxcors_B=[]         
        w=[slice(x,504+x) for x in xrange(0,4032-(504/2),504/2)]
#        baseshift=tser[ep]['x']-tser[bestep]['x']
#        if np.all(baseshift) > 0:
#            tser[ep]['x']-=baseshift
#        elif np.all(baseshift) < 0:
#            tser[ep]['x']+=baseshift
        #3.1. Apply barycenter corrections
        tser0x=[tser[bestep]['x'][j]*(1+bc0/c) for j in xrange(3)]
        tser1x=[tser[ep]['x'][j]*(1+bc1/c) for j in xrange(3)]
        tser0f=tser[bestep]['fx']
        tser1f=tser[ep]['fx']
        maxp = [[[c*(np.max(tser1x[j][pc])-np.min(tser0x[j][pc]))/np.min(tser0x[j][pc]) for pc in w] for j in [0,1,2]]][0]
        maxn = [[[c*(np.min(tser1x[j][pc])-np.max(tser0x[j][pc]))/np.max(tser0x[j][pc]) for pc in w] for j in [0,1,2]]][0]
        for key in ums.iterkeys():
            for i in range(0, 1+nits):
                xcor_i=[]
                for j in [0,1,2]: #
                    for p, pc in enumerate(w):
                        if np.sum(ums[key][j][pc]) > 252:
                            #print pair, j, pc
                            continue
        #3.2. Do fake Monte-Carlos
                        if i>0:
                            t0f=np.ma.array(noise_scrambler(tser0f[0][j][pc],tser0f[1][j][pc],tser0f[2][j][pc]),
                                            mask=ums[key][j][pc])
                            t1f=np.ma.array(noise_scrambler(tser1f[0][j][pc],tser1f[1][j][pc],tser1f[2][j][pc]),
                                            mask=ums[key][j][pc])
                        else:
                            t0f=np.ma.array(tser0f[0][j][pc],mask=ums[key][j][pc])
                            t1f=np.ma.array(tser1f[0][j][pc],mask=ums[key][j][pc])
        #3.3. De-trend, get rid of residual broad sky absorption
                        fx0=detrend(np.ma.array(tser0x[j][pc],mask=ums[key][j][pc]),t0f)
                        fx1=detrend(np.ma.array(tser1x[j][pc],mask=ums[key][j][pc]),t1f)
                        xcor=np.correlate(fx0,fx1,mode='full')
        #3.4. Resample to *1 km/s*-wide bins
                        if len(xcor)%2==0:
                            vshifts=np.concatenate((np.linspace(maxn[j][p],0,num=len(xcor)/2,endpoint=False),
                                                    np.linspace(0,maxp[j][p],num=(len(xcor)/2))))
                        else:
                            vshifts=np.concatenate((np.linspace(maxn[j][p],0,num=len(xcor)/2,endpoint=False),
                                                    np.linspace(0,maxp[j][p],num=(1+len(xcor)/2))))
                        vbins=np.arange(np.floor(maxn[j][p]), np.ceil(maxp[j][p])+1)
                        xcor_i+=zip(vbins, np.interp(vbins,vshifts,xcor))
                if i==0:
                    if key=='unimask':
                        xcors=xcor_i
                    elif key=='ubmask':
                        xcors_B=xcor_i
                elif i>0:
                    if key=='unimask':
                        Nxcors+=[xcor_i]
                    elif key=='ubmask':
                        Nxcors_B+=[xcor_i]
        print '--- Cross-correlation took {} sec ---'.format(time.time() - start_time)
        start_time = time.time() #reset
    #4. Co-Add cross-correlations (DO NOT UNINDENT)
        dsum=dict()
        dsumB=dict()
        for (vs, amps) in xcors:
            if vs in dsum:
                dsum[vs].append(amps)
            else:
                dsum[vs]=[amps]
        for (vsb, ambs) in xcors_B:
            if vsb in dsumB:
                dsumB[vsb].append(ambs)
            else:
                dsumB[vsb]=[ambs]
        ccc=[[h,g] for h,g in pointsum(dsum)] 
        ccc_B=[[h,g] for h,g in pointsum(dsumB)]
#        Nccc=len(Nxcors)*[[]]
#        Nccc_B=len(Nxcors_B)*[[]]
#        for i,arr in enumerate(Nxcors):
#            dNsum=dict()
#            for (vs, amps) in arr:
#                if vs in dNsum:
#                    dNsum[vs].append(amps)
#                else:
#                    dNsum[vs]=[amps]
#            Nccc[i]=[[h,g] for h,g in pointsum(dNsum)] 
#        for i,arr in enumerate(Nxcors_B):
#            dNsumB=dict()
#            for (vsb, ambs) in arr:
#                if vsb in dNsumB:
#                    dNsumB[vsb].append(ambs)
#                else:
#                    dNsumB[vsb]=[ambs]
#            Nccc_B[i]=[[h,g] for h,g in pointsum(dNsumB)]
        print '--- Co-addition took {} sec ---'.format(time.time() - start_time)
        start_time = time.time() #reset
    #5. Extract & output peak velocities with estimated errors
        vels,power=np.transpose(ccc)
        rel_rv=vels[np.argmax(power)]
        #5.1. Gaussian-fit cross-correlation curves with Brackett Lines included
        # block out central blip b/c it's probably due to telluric contamination
        velsb,powerb=np.transpose(ccc_B)
#        zpt=np.argmax(powerb)#This will probably be contamination
#        fitdom=np.concatenate((velsb[zpt-200:zpt-10],velsb[zpt+10:zpt+200]))
#        fitran=np.concatenate((powerb[zpt-200:zpt-10],powerb[zpt+10:zpt+200]))
#        try:
#            bps,bcov = curve_fit(gauss_fit,fitdom,fitran,p0=[velsb[zpt],50.,np.max(powerb)-np.min(powerb),0.4])
#        except RuntimeError:
#            try:
#                bps,bcov = curve_fit(gauss_fit,fitdom,fitran, p0=[velsb[zpt],100.,np.max(powerb)-np.min(powerb),0.4])
#            except RuntimeError:
#                bps=4*[-9999.99]
#        newdom=velsb[zpt-200:zpt+200]
#        psfit=gauss_fit(newdom,bps[0],bps[1],bps[2],bps[3])
        rel_rv_B=velsb[np.argmax(powerb)]#[bps[0],bps[1]]
#        ofi2.write('#epoch: {}, ref_ep: {}\n'.format(ep,bestep))
#        ofi2.write('#original RV shift: {} with Brackett mask, {} without\n'.format(rel_rv,rel_rv_B[0]))
#        ofi2.write('#run    maxRVshift    Brmask(Y/N)\n')
#        Nrel_rv=len(Nccc)*[[]]
        try:
            hep = [hjd for hjd in epsqr[star] if str(ep) in str(hjd)][0]
        except IndexError:
            try:            
                hep = [hjd for hjd in epsqr[star] if str(ep+1) in str(hjd)][0]
            except IndexError:
                hep = [hjd for hjd in epsqr[star] if str(ep-1) in str(hjd)][0]
        runs=1
#        for k,v in enumerate(Nccc):
#            vels,power=np.transpose(v)
#            Nrel_rv[k]=vels[np.argmax(power)]
##            ofi2.write('{}    {}    Y\n'.format(runs,Nrel_rv[k]))
#            runs+=1
#        rms_rel_rv=round(np.std(Nrel_rv[:]),3)
##        ofi2.write('#    stdev: {}\n'.format(rms_rel_rv))
#        runs=1
#        Nrel_rv_B=len(Nccc_B)*[[]]
#        for k,v in enumerate(Nccc_B):
#            vels,power=np.transpose(v)
#            Nrel_rv_B[k]=vels[np.argmax(power)]
##            ofi2.write('{}    {}    N\n'.format(runs,Nrel_rv_B[k]))
#            runs+=1
#        rms_rel_rv_B=round(np.std(Nrel_rv_B[:]),3)
#        ofi2.write('#    stdev: {}\n\n'.format(rms_rel_rv_B))
        #Write these results to file:
#        print 'RV, no brackett lines: ', rel_rv, ' +/- ', rms_rel_rv
#        print 'RV, with brackett lines: ', rel_rv_B, ' +/- ', rms_rel_rv_B
#        ofi.write('{:.7}    {}    {}    {}    {}    {}\n'.format(hep,bc1,
#                  psep,round(rel_rv,3),
#                  rms_rel_rv, 'Y'))        
#        ofi.write('{:.7}    {}    {}    {}    {}    {}\n'.format(hep,bc1,
#                  psep,round(rel_rv_B[0],3),
#                  rms_rel_rv_B, 'N'))
                
        pl.figure(int(star[-5:]))
        figure = pl.gcf()
        figure.set_size_inches(8, 6)
        pl.rcParams.update({'font.size':8})
        pl.subplot(3,1,1)
        x,y=np.transpose(ccc)
        pl.plot(x,y/(max(y)-min(y)),label='ep {}: v={:.2}'.format(ep,rel_rv))#$\pm${:.2},rms_rel_rv
        pl.axvline(linewidth=2,linestyle=':',color='r')
        pl.xlim(-400,400)
        pl.title('{}: Co-Added Cross-Correlation, Brackett Series masked'.format(star))
        pl.ylabel('Cross-correlation')
        pl.legend(fontsize=6).get_frame().set_alpha(0.5)
        ##################
        pl.subplot(3,1,2)
        xb,yb=np.transpose(ccc_B)
        pl.plot(xb,yb/(max(yb)-min(yb)),label='ep {}: v={:.2}'.format(ep,rel_rv_B))#$\pm${:.2},rms_rel_rv_B
#        pl.plot(newdom,psfit/(max(psfit)-min(psfit)),linestyle='--',
#                label='ep {}: $\mu$={:.2}, $\sigma$={:.2}'.format(ep,bps[0],bps[1]))
        pl.axvline(linewidth=2,linestyle=':',color='r')
        pl.xlabel('RV shift relative to epoch {} (km/s)'.format(bestep))
        pl.ylabel('Cross-correlation')
        pl.title('...Brackett Series included')
        pl.xlim(-400,400)
        pl.legend(fontsize=6).get_frame().set_alpha(0.5)
#        pl.subplot(3,1,3)
#        if np.abs(rms_rel_rv)<0.1:
#            rms_rel_rv=0.1
#        else:
#            rms_rel_rv=rms_rel_rv*(nits/(nits-1))
#        if np.abs(rms_rel_rv_B)<0.1:
#            rms_rel_rv_B=0.1
#        else:
#            rms_rel_rv=rms_rel_rv*(nits/(nits-1))
#        pl.errorbar(ep, rel_rv, yerr=rms_rel_rv, marker='o', color='m')
#        pl.errorbar(ep, rel_rv_B, yerr=rms_rel_rv_B, marker='o', color='b')
#        pl.axhline(color='k', linestyle='--')
#        pl.xlabel('Date (MJD5)')
#        pl.ylabel('RV shift (km/s)')
#        pl.tight_layout()
        #pl.show()
        print '--- Plotting and output took {} sec ---'.format(time.time() - start_time)
    pl.subplot(3,1,3)
    visit,rv,erv=np.loadtxt('Outputs/'+star+'.out', usecols=(0,3,4), unpack=True)  
    pl.errorbar(visit[::2], rv[::2], yerr=erv[::2], capsize=4, fmt='.', color='m', label='RV, Br-masked')
    pl.errorbar(visit[1::2], rv[1::2], yerr=erv[1::2], capsize=3, fmt='.', color='b', label='RV, Br included')
    pl.axhline(color='k', linestyle='--')
    pl.xlabel('Date (MJD5)')
    pl.ylabel('RV shift (km/s)')
    pl.legend(loc=8)
    pl.tight_layout()
    pl.savefig('ABeOutputs/{}_xcor-Redo.png'.format(star), dpi=150)
    pl.close()
#    ofi2.close()
#    ofi.close() #There is a sys.exit() below to comment/uncomment depending on
    print '-*-*-Total time for star {}: {} s-*-*-'.format(star,time.time() - tot_time)
    # whether or not you want to stop after 1 star
    #trashdump.dump_garbage(save=True)
    #sys.exit() #comment to stop breaking after 1 star
        
