#module taken from www.johnvinyard.com/blog/?p=268 "Efficient Overlapping Windows with Numpy." Zounds: Chronicling the Development of Zounds, a Python Library for Machine Listening. Aug. 30, 2012

import numpy as np
from numpy.lib.stride_tricks import as_strided as ast
#from itertools import product 
#import exceptions

def sliding_window(a,ws,ss = None,flatten = True):
    '''
    Return a sliding window over a in any number of dimensions
     
    Parameters:
        a  - an n-dimensional numpy array
        ws - an int (a is 1D) or tuple (a is 2D or greater) representing the size of each dimension of the window
        ss - an int (a is 1D) or tuple (a is 2D or greater) representing the amount to slide the window in each dimension. If not specified, it defaults to ws.
        flatten - if True, all slices are flattened, otherwise, there is an extra dimension for each dimension of the input.
     
    Returns
        an array containing each n-dimensional window from parameter a
    '''
    def norm_shape(shape):
        '''
        Normalize numpy array shapes so they're always expressed as a tuple,
        even for one-dimensional shapes.

        Parameters
            shape - an int, or a tuple of ints

        Returns
            a shape tuple
        '''
        try:
            i = int(shape)
            return (i,)
        except TypeError:
            # shape was not a number
            pass

        try:
            t = tuple(shape)
            return t
        except TypeError:
            # shape was not iterable
            pass

        raise TypeError('shape must be an int, or a tuple of ints')
 
    if None is ss:
        # ss was not provided. the windows will not overlap in any direction.
        ss = ws
    ws = norm_shape(ws)
    ss = norm_shape(ss)
     
    # convert ws, ss, and a.shape to numpy arrays so that we can do math in every
    # dimension at once.
    ws = np.array(ws)
    ss = np.array(ss)
    shape = np.array(a.shape)
     
     
    # ensure that ws, ss, and a.shape all have the same number of dimensions
    ls = [len(shape),len(ws),len(ss)]
    if 1 != len(set(ls)):
        raise ValueError(\
        'a.shape, ws and ss must all have the same length. They were %s' % str(ls))
     
    # ensure that ws is smaller than a in every dimension
    if np.any(ws > shape):
        raise ValueError(\
        'ws cannot be larger than a in any dimension.\
 a.shape was %s and ws was %s' % (str(a.shape),str(ws)))
     
    # how many slices will there be in each dimension?
    newshape = norm_shape(((shape - ws) // ss) + 1)
    # the shape of the strided array will be the number of slices in each dimension
    # plus the shape of the window (tuple addition)
    newshape += norm_shape(ws)
    # the strides tuple will be the array's strides multiplied by step size, plus
    # the array's strides (tuple addition)
    newstrides = norm_shape(np.array(a.strides) * ss) + a.strides
    strided = ast(a,shape = newshape,strides = newstrides)
    if not flatten:
        return strided
     
    # Collapse strided so that it has one more dimension than the window.  I.e.,
    # the new array is a flat list of slices.
    meat = len(ws) if ws.shape else 0
    firstdim = (np.product(newshape[:-meat]),) if ws.shape else ()
    dim = firstdim + (newshape[-meat:])
    # remove any dimensions with size 1
    dim = filter(lambda i : i != 1,dim)
    return strided.reshape(dim)

#-----------------------------------------------------------------------------

def gen_iterwin(bgn, end, maxlen):
    '''
    A generator function. Not currently generalized to N-D (will work on it)
    
    Returns a list of iterators over a fixed-length sub-array whose position
    within & size relative to the parent array is not known a priori.
    Truncates automatically to fit within array boundaries.
    (Basically I made this to avoid the IndexError exception that formerly
    required a 3-tiered try statement or a 3-case if statement to resolve)
    
    Paramaters:
    bgn - first index of the sub-array; scalar input & output 
    end - 1 + last index of the sub-array; scalar input & output
    maxlen - length of the parent array (recommend use of len())
    
    Be careful about entering floats as bgn & end: gen_iterwin won't return an
    error so unless you explicitly force the numbers therein to ints, they 
    may not work as indices of a different array.
    '''
    if end >= maxlen:
        yield [ndx for ndx in np.arange(bgn,maxlen)]
    elif bgn <= 0:
        yield [ndx for ndx in np.arange(0,end)]
    else:
        yield [ndx for ndx in np.arange(bgn,end)]
        
#-----------------------------------------------------------------------------

def win_wash(winarr2d, overlap, flags2d=np.array([None]),
             noflag=None, nanflag=None, fixedflag=None):
    '''
    Returns the 1D recombination of a set of overlapping sliding windows.
    If flags are supplied, data in the overlapping regions may be averaged,
    substituted, or thrown out according to the flagging hierarchy, & win_wash
    will return an array containing the 1D recombination and a separate array
    of flags corresponding to the data retained.
    
    Notes:
    Designed to undo the work of sliding_window (above)
    May in the future be generalized to 3D
    
    Parameters:
    winarr2d : array_like
            Input array of shape (# of windows, # of entries per window)
    overlap: int
            No. of points shared by 2 windows in winarr3d; must be > 0
    flags2d : array_like, optional
            Input array of same shape as winarr2d; defaults to None; if None,
            all subsequent **kwargs are ignored
    noflag : int, optonal
            Value in flags2d corresponding to unaltered data in winarr2d;
            defaults to 0
    nanflag : int, optional
            Value in flags2d corresponding to an instance of Inf or NaN in
            winarr2d; defaults to -99999. Does not currently distinguish Inf
            from NaN, but that may be changed. If inf or NaN occurs in winarr2d
            that was not flagged in flags2d, a nanflag will be assigned in the
            output
    fixedflag : int, optional
            Value in flags corresponding to an interpolated or artificial
            value; defaults to noflag if noflag is supplied, 0 otherwise
            
    Returns:
    shutters(, shutterflag) : 1D recombination of winarr2d, retained or
            averaged according to flags2d if supplied; (optional) 1D
            recombination of flags2d
                    
    
    '''
    span=np.arange(0,len(winarr2d)-1)
    window=0
    
    if flags2d.any() is None:
        shutters=winarr2d[0][:-overlap]
        while window in span:
            panes1=np.split(winarr2d[window],[-overlap])
            panes2=np.split(winarr2d[window+1],[overlap])
            for index in range(0,overlap):
                wv1=panes1[1][index]
                wv2=panes2[0][index]
                if wv1==wv2:
                    panes2[0][index]=wv1
                else:
                    panes2[0][index]=np.average((wv1,wv2))
            shutters=np.concatenate((shutters,panes2[0]))
            if window==span[-1]:
                shutters=np.concatenate((shutters,panes2[1]))
                break
            else:
                window+=1
        return shutters
        
    else:
        shutters=winarr2d[0][:-overlap]
        shutterflag=flags2d[0][:-overlap]
        if noflag is None:
            noflag=0
        if np.NaN in winarr2d and nanflag is None:
            nanflag=-99999
        if fixedflag is None:
            fixedflag=noflag
        flagrank=[noflag, fixedflag, nanflag]
        
        while window in span:
            panes1=np.split(winarr2d[window],[-overlap])
            panes2=np.split(winarr2d[window+1],[overlap])
            squeegee1=np.split(flags2d[window],[-overlap])
            squeegee2=np.split(flags2d[window+1],[overlap])
            
            for index in range(0,overlap):
                wv1=panes1[1][index]
                wv2=panes2[0][index]
                flg1=squeegee1[1][index]
                flg2=squeegee2[0][index]
                
                if (flg1==flg2!=nanflag and np.isnan(wv1)==np.isnan(wv2)==False
                    and np.inf not in [wv1,wv2]):
                    panes2[0][index]=np.average((wv1,wv2))
                    squeegee2[0][index]=flg2
                elif np.isnan(wv1)==np.isnan(wv2)==True or wv1==wv2==np.inf or flg1==flg2==nanflag:
                    panes2[0][index]=np.NaN
                    squeegee2[0][index]=nanflag
                elif np.isnan(wv1)==True and np.isnan(wv2)==False:
                    panes2[0][index]=wv2
                    squeegee2[0][index]=flg2                    
                elif np.isnan(wv2)==True and np.isnan(wv1)==False:
                    panes2[0][index]=wv1
                    squeegee2[0][index]=flg1                  
                else:
                    rank1=flagrank.index(flg1)
                    rank2=flagrank.index(flg2)
                    #less is more
                    if rank1<rank2:
                        panes2[0][index]=wv1
                        squeegee2[0][index]=flg1
                    elif rank1>rank2:
                        panes2[0][index]=wv2
                        squeegee2[0][index]=flg2
                    else:
                        panes2[0][index]=np.average((wv1,wv2))
                        squeegee2[0][index]=flg2                
            
            shutters=np.concatenate((shutters,panes2[0]))
            shutterflag=np.concatenate((shutterflag,squeegee2[0]))
            if window==span[-1]:
                shutters=np.concatenate((shutters,panes2[1]))
                shutterflag=np.concatenate((shutterflag,squeegee2[1]))
                break
            else:
                window+=1
        return shutters, shutterflag

#-----------------------------------------------------------------------------

def slide_diff(a1, a2, scale_by=None, fill_val=None, normed=True, masked=False):
    '''
    Finds the pointwise shift between 2 1D arrays, where both values exist, for
    all possible relative shifts between them.
    
    Designed for use alongside np.correlate when cross-correlating 2 signals
    whose digital or pixel spacing maps to some other phase space (e.g. 
    frequency, wavelength, or velocity) nonlinearly or with irregular sampling
    intervals.
    If you have to do anything like this in higher dimensions, you might have
    better luck with template matching using scipy.ndimage
    
    Parameters:
    -----------
    a1 : 1d array or masked array
        Reference array; must be the same length or longer than a2
    a2 : 1d array or masked array
        Array to "slide across" a1; must be no longer than a1
    scale_by : int or float, optional
        Additional scalar to multiply the differences by; default is 1
    fill_val : int or float, optional
        Fill value to assign to points where a2 does not overlap a1; if None,
        defaults to NaN
    normed : bool, optional
        If True (default), divides the difference by the value of a1 at each
        point; if False, difference is left alone
    masked : bool, optional
        If True, masks any entry in diffs equal to fill_val; default is False
        *Note: if you want other masking options, you can either set 'masked'
         to False & make the output a masked array later, or you can output
         diffs with a mask & then add to, alter, or replace the mask
         
    Returns:
    --------
    diffs : 2d array with potentially masked values
        Rows are all possible shifts of a2 relative to a1 where a2 & a1
        have at least 1 overlapping point
        Columns are the pointwise differences between a1 & a2 for the given
        shift, where the 2 overlap, using a1 as the reference frame
        
    Raises:
    -------
    IndexError
        Raised if a2 is longer than a1

    Example:
    >>> a = np.arange(2,18,2)
    >>> a
    array([ 2,  4,  6,  8, 10, 12, 14, 16])
    >>> b = np.arange(0,15,3)
    >>> b
    array([ 0,  3,  6,  9, 12])
    >>> slide_diff(a,b, normed=False)
array([[-5.        ,         nan,         nan,         nan,         nan,
                nan,         nan,         nan],
       [-3.5       , -2.        ,         nan,         nan,         nan,
                nan,         nan,         nan],
       [-2.        , -1.25      , -1.        ,         nan,         nan,
                nan,         nan,         nan],
       [-0.5       , -0.5       , -0.5       , -0.5       ,         nan,
                nan,         nan,         nan],
       [ 1.        ,  0.25      ,  0.        , -0.125     , -0.2       ,
                nan,         nan,         nan],
       [        nan,  1.        ,  0.5       ,  0.25      ,  0.1       ,
         0.        ,         nan,         nan],
       [        nan,         nan,  1.        ,  0.625     ,  0.4       ,
         0.25      ,  0.14285714,         nan],
       [        nan,         nan,         nan,  1.        ,  0.7       ,
         0.5       ,  0.35714286,  0.25      ],
       [        nan,         nan,         nan,         nan,  1.        ,
         0.75      ,  0.57142857,  0.4375    ],
       [        nan,         nan,         nan,         nan,         nan,
         1.        ,  0.78571429,  0.625     ],
       [        nan,         nan,         nan,         nan,         nan,
                nan,  1.        ,  0.8125    ],
       [        nan,         nan,         nan,         nan,         nan,
                nan,         nan,  1.        ]])


    >>> slide_diff(a,b, fill_val=0., normed=False, masked=True)
    masked_array(data =
     [[-10.0 -- -- -- -- -- -- --]
     [-7.0 -8.0 -- -- -- -- -- --]
     [-4.0 -5.0 -6.0 -- -- -- -- --]
     [-1.0 -2.0 -3.0 -4.0 -- -- -- --]
     [2.0 1.0 -- -1.0 -2.0 -- -- --]
     [-- 4.0 3.0 2.0 1.0 -- -- --]
     [-- -- 6.0 5.0 4.0 3.0 2.0 --]
     [-- -- -- 8.0 7.0 6.0 5.0 4.0]
     [-- -- -- -- 10.0 9.0 8.0 7.0]
     [-- -- -- -- -- 12.0 11.0 10.0]
     [-- -- -- -- -- -- 14.0 13.0]
     [-- -- -- -- -- -- -- 16.0]],
                 mask =
     [[False  True  True  True  True  True  True  True]
     [False False  True  True  True  True  True  True]
     [False False False  True  True  True  True  True]
     [False False False False  True  True  True  True]
     [False False  True False False  True  True  True]
     [ True False False False False  True  True  True]
     [ True  True False False False False False  True]
     [ True  True  True False False False False False]
     [ True  True  True  True False False False False]
     [ True  True  True  True  True False False False]
     [ True  True  True  True  True  True False False]
     [ True  True  True  True  True  True  True False]],
           fill_value = 0.0)

    '''
    if len(a2)>len(a1):
        raise IndexError("a2 cannot be longer than a1")
    
    a2=a2[::-1]
    diffs=np.empty((len(a1),(len(a1)+len(a2)-1)))
    if fill_val!=None:
        diffs.fill(fill_val)
    else:
        diffs.fill(np.NaN)
    if scale_by!=None:
        s=scale_by
    else:
        s=1.
    for p in range(0,len(a1)):
        for q in range(0,len(a2)):
            if a1[p] is np.ma.masked or a2[q] is np.ma.masked:
                diffs[p][q+p]=np.NaN
            if not normed: #note: when I first wrote this, a2[q] and a1[p] were switched
            #for abepipe2, it worked better after reversing them
            #changes are not reflected in the examples, but they're basically the same - just rotated 180 degrees
                diffs[p][q+p]=s*(a2[q]-a1[p])
            else:
                diffs[p][q+p]=s*(a2[q]-a1[p])/a1[p]
    if masked and fill_val==None:
        return np.ma.masked_invalid(np.transpose(diffs))
    elif masked and fill_val!=None:
        return np.ma.masked_values(np.transpose(diffs), fill_val)
    else:
        return np.transpose(diffs)
            



                