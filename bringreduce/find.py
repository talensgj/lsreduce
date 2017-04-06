# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 16:25:04 2016

@author: talens
"""

import numpy as np

from . import misc


def peaks(image, sigma=25.):
    """ Find local maxima in an image. """

    from skimage.feature import peak_local_max
    from scipy.ndimage.filters import median_filter, gaussian_filter 
    
    ny, nx = image.shape    
    
    # Prepare the image.
    simage = image - median_filter(image, size=11)
    simage = gaussian_filter(simage, sigma=2.5) 
         
    # Detect local maxima.
    m0, m1 = misc.sigma_clip(simage, sigma=3.)    
    t_abs = m0 + sigma*m1
    tmp = peak_local_max(simage, indices=True, threshold_abs=t_abs, threshold_rel=0)
    
    if tmp.ndim != 2:
        return np.array([]), np.array([])
    else:
        xi, yi = tmp[:,1], tmp[:,0]    
    
    xi, yi = xi.astype('float'), yi.astype('float')    
    
    return xi, yi
    
def gaussfit(image, x, y, nhalf=15):
    """ Fit 2D Guassians to sub-images. """    
    
    from astropy.modeling import models, fitting    
    
    ny, nx = image.shape
    nstars = len(x)
    nbox = 2*nhalf + 1    

    # Seperate the coordinates in pixel center and deviation.
    xi, yi = np.around(x), np.around(y)
    dx, dy = x - xi, y - yi

    # Set up the 2D Gaussian model.
    tmp = np.arange(nbox) - nhalf
    xx, yy = np.meshgrid(tmp, tmp)

    t_init = models.Gaussian2D(x_mean=0., y_mean=0., x_stddev=2.5, y_stddev=2.5, theta=0.) + models.Const2D()
    fit_t = fitting.LevMarLSQFitter()
    
    flag = np.zeros(nstars, dtype='uint8')
    for i in range(nstars):
        
        # Extract the sub-image.
        lx = int(xi[i] - nhalf)
        ux = int(xi[i] + nhalf + 1)
        ly = int(yi[i] - nhalf)
        uy = int(yi[i] + nhalf + 1)
    
        if (lx < 0) | (ux > nx) | (ly < 0) | (uy > ny):
            flag[i] += 1
            continue
    
        subim = np.copy(image[ly:uy,lx:ux])
        
        # Set initial values and perform the fit.
        
        t_init.amplitude_0 = np.maximum(500, subim[nhalf,nhalf] - np.median(subim))
        t_init.amplitude_1 = np.median(subim)
        
        t = fit_t(t_init, xx, yy, subim)
        
        if (t.amplitude_0.value <= 0):
            flag[i] += 2
            
        if (t.amplitude_1.value <= 0):
            flag[i] += 4
        
        # Extract the best fit coordinates.
        dx[i] = t.x_mean_0.value
        dy[i] = t.y_mean_0.value
        
        if (np.abs(dx[i]) > nhalf) | (np.abs(dy[i]) > nhalf):
            flag[i] += 8
            
    return xi + dx, yi + dy, flag
    
def gaussfit_mp(image, x, y, nhalf=15, nproc=4):

    import multiprocessing as mp

    idx = np.arange(len(x))

    p = mp.Pool(nproc)
    results = [(idx[i::nproc], p.apply_async(gaussfit, (image, x[i::nproc], y[i::nproc], nhalf))) for i in range(nproc)]
    p.close()
    p.join()

    xr = np.zeros(x.shape, x.dtype)
    yr = np.zeros(y.shape, y.dtype)
    flag = np.zeros(x.shape, dtype='uint8')

    for idx, result in results:

        xr[idx], yr[idx], flag[idx] = result.get()

    return xr, yr, flag 

def find(image, sigma=25., margin=50, nhalf=15, method='gaussfit'):
    """ Detect stars in an image. """    
    
    ny, nx = image.shape    
    
    # Obtain positions of local maxima.
    xi, yi = peaks(image, sigma)
    
    # Refine the positions.
    if (method == 'gaussfit'):
        x, y, flag = gaussfit_mp(image, xi, yi, nhalf) 
    else:
        print 'Unknown method: {}'.format(method)
        exit()
        
    # Remove the margin.
    mask = misc.clip_rectangle(x, y, nx, ny, margin)  
    x, y, flag = x[mask], y[mask], flag[mask] 
    
    return x, y, flag
    