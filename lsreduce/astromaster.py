# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 14:12:13 2016

@author: talens
"""

import os

import numpy as np

from astropy import wcs
from astropy.io import fits

import matplotlib.pyplot as plt

from . import io, misc
from . import find, astrometry

import logging 

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('lsreduce')

def initial_wcs(pars, scale, lst):
    
    # Extract parameters.
    ha0, dec0 = pars[0], pars[1]
    x0, y0 = pars[2], pars[3]
    th0 = pars[4]
    
    # Convert HA to RA. 
    ra0 = np.mod(15.*lst - ha0, 360.)
    
    # Compute PC matrix.
    sn = np.sin(th0*np.pi/180.)
    cs = np.cos(th0*np.pi/180.)        
    pc0 = np.array([[cs, -sn], [sn, cs]])

    # Convert scale to degrees.
    scale = np.rad2deg(scale)    
    
    # Generate WCS object.
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [x0, y0]
    w.wcs.cdelt = [scale, scale]
    w.wcs.crval = [ra0, dec0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.pc = pc0
    
    return w    
    
def match(xccd, yccd, xwcs, ywcs, chunksize=None):  
    
    if chunksize is None:
        chunksize = len(xccd)/10 + 1
    
    #  Create arrays.
    idx1 = np.arange(len(xccd), dtype='int')
    idx2 = np.zeros(len(xccd), dtype='int')
    dist_sq = np.zeros(len(xccd))
    
    xccd, yccd = xccd[:,np.newaxis], yccd[:,np.newaxis]    
    xwcs, ywcs = xwcs[np.newaxis,:], ywcs[np.newaxis,:]
    
    # Match stars.
    for i in range(0, len(xccd), chunksize):
        
        dtmp = (xccd[i:i+chunksize] - xwcs)**2 + (yccd[i:i+chunksize] - ywcs)**2
        
        args = np.argmin(dtmp, axis=1)
        
        idx2[i:i+chunksize] = args
        dist_sq[i:i+chunksize] = np.diag(dtmp[:,args])    
    
    return idx1, idx2, np.sqrt(dist_sq)
    
def lnlike_wcs(pars, xccd, yccd, ra, dec, scale, lst):
    
    # Predict positions using proposed solution.
    w = initial_wcs(pars, scale, lst)     
    xwcs, ywcs = w.wcs_world2pix(ra, dec, 0)         
          
    # Match detected positions to prediction.
    idx1, idx2, dist = match(xccd, yccd, xwcs, ywcs)       
      
    # Reject a fixed number of the worst matches.
    mask = dist < np.percentile(dist, 90)    
    dist = dist[mask]       
    
    return -np.sum(dist**2)

def blind_solution(image, header, ra, dec, pars0, sigma=20., margin=-100, scale=9e-3/24., order=5, maxiter=20, debug=False):

    import emcee

    ny, nx = image.shape   
    
    # Thin out the catalogue.
    print 'Selecting catalogue stars based on initial parameters.'    

    ra, dec = astrometry.j2000_to_equinox(ra, dec, header['jd'])
    
    w = initial_wcs(pars0, scale, header['lst'])
    xwcs, ywcs = w.wcs_world2pix(ra, dec, 0)
    
    mask = np.isfinite(xwcs) & misc.clip_rectangle(xwcs, ywcs, nx, ny, margin)
    ra, dec = ra[mask], dec[mask]              
    
    # Detect stars in the image.
    print 'Performing blind search for stars.'
    
    xccd, yccd, flags = find.find(image, sigma=sigma)
    mask = (flags == 0)
    xccd, yccd = xccd[mask], yccd[mask]    

    print 'Detected {} stars out of {} expected stars.'.format(len(xccd), len(ra))
    
    if debug:
        
        xwcs, ywcs = w.wcs_world2pix(ra, dec, 0)        
        
        vmin = np.nanpercentile(image, 1)
        vmax = np.nanpercentile(image, 99)    
        
        plt.imshow(image, origin='lower', interpolation='None', vmin=vmin, vmax=vmax, cmap=plt.cm.Greys)
        
        plt.plot(xwcs, ywcs, 'ro', label='WCS')
        plt.plot(xccd, yccd, 'go', label='CCD')
        plt.xlim(-0.5, nx-0.5)
        plt.ylim(-0.5, ny-0.5)        
        
        plt.legend()
        
        plt.tight_layout()
        plt.show()      
        plt.close()
    
    print 'Guessing WCS parameters:'    
    
    # Set up the MCMC to guess the parameters.
    ndim, nwalkers = 5, 1000 # Large number of walkers.
    dpars = np.array([2., 2., 50., 50., 5.]) # Large space to sample.
    pars = [pars0 + dpars*np.random.randn(ndim) for i in range(nwalkers)]
    
    # Run the MCMC.
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike_wcs, args=(xccd, yccd, ra, dec, scale, header['lst']), threads=4)
    sampler.run_mcmc(pars, 1) # Single step.

    # Get the best-fit parameters.
    samples = sampler.flatchain
    lnprob = sampler.flatlnprobability
    
    arg = np.argmax(lnprob)
    pars = samples[arg] 
    
    print '  HA = {:.2f}, Dec = {:.2f}, x = {:.2f}, y = {:.2f}, theta = {:.2f}'.format(*pars)    
    
    print 'Refining WCS parameters:'    
    
    # Set up the MCMC to refine the parameters.
    ndim, nwalkers = 5, 50 # Small number of walkers.
    dpars = np.array([5./60, 5./60, 5., 5., 5./60]) # 5 arcmin (i.e. pixel) space to sample. 
    pars = [pars + dpars*np.random.randn(ndim) for i in range(nwalkers)]
    
    # Run the MCMC.
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike_wcs, args=(xccd, yccd, ra, dec, scale, header['lst']), threads=4)
    sampler.run_mcmc(pars, 250) # Many steps.   
    
    # Get the best fit parameters.
    samples = sampler.flatchain
    lnprob = sampler.flatlnprobability 
    
    arg = np.argmax(lnprob)
    pars = samples[arg]  
    
    print '  HA = {:.2f}, Dec = {:.2f}, x = {:.2f}, y = {:.2f}, theta = {:.2f}'.format(*pars)       
    
    # Create a WCS object from the best fit.
    w = initial_wcs(pars, scale, header['lst'])   
         
    if debug:

        plt.subplot(321)            
        
        plt.plot(sampler.lnprobability.T, c='k')
        
        plt.xlabel('Step')
        plt.ylabel(r'$\ln L$')
       
        plt.subplot(322)            
        
        plt.plot(sampler.chain[:,:,4].T, c='k')
        
        plt.xlabel('Step')
        plt.ylabel(r'$\theta_0$')
        
        plt.subplot(323)            
        
        plt.plot(sampler.chain[:,:,0].T, c='k')
        
        plt.xlabel('Step')
        plt.ylabel(r'$\rm{HA}_0$')
        
        plt.subplot(324)            
        
        plt.plot(sampler.chain[:,:,1].T, c='k')
        
        plt.xlabel('Step')
        plt.ylabel(r'$\rm{Dec}_0$')
        
        plt.subplot(325)            
        
        plt.plot(sampler.chain[:,:,2].T, c='k')
        
        plt.xlabel('Step')
        plt.ylabel(r'$x_0$')
        
        plt.subplot(326)            
        
        plt.plot(sampler.chain[:,:,3].T, c='k')
        
        plt.xlabel('Step')
        plt.ylabel(r'$y_0$')
        
        plt.tight_layout()
        plt.show()
        plt.close()
        
        xwcs, ywcs = w.wcs_world2pix(ra, dec, 0)          
        
        vmin = np.nanpercentile(image, 1)
        vmax = np.nanpercentile(image, 99)    
        
        plt.imshow(image, origin='lower', interpolation='None', vmin=vmin, vmax=vmax, cmap=plt.cm.Greys)
        
        plt.plot(xwcs, ywcs, 'ro', label='WCS')
        plt.plot(xccd, yccd, 'go', label='CCD')
        plt.xlim(-0.5, nx-0.5)
        plt.ylim(-0.5, ny-0.5)        
        
        plt.legend()
        
        plt.tight_layout()
        plt.show()  
        plt.close()
        
    print 'Creating wcspars and polpars dictionaries.'        
        
    wcspars = dict()
    wcspars['lst'] = header['lst']
    wcspars['crpix'] = w.wcs.crpix
    wcspars['crval'] = w.wcs.crval
    wcspars['cdelt'] = w.wcs.cdelt
    wcspars['pc'] = w.wcs.pc  
    
    polpars = dict()
    polpars['order'] = order 
    polpars['nx'] = nx
    polpars['ny'] = ny
    polpars['x_pix2wcs'] = np.zeros((order+1, order+1))
    polpars['y_pix2wcs'] = np.zeros((order+1, order+1))
    polpars['x_wcs2pix'] = np.zeros((order+1, order+1))
    polpars['y_wcs2pix'] = np.zeros((order+1, order+1))    
        
    print 'Iteratively refining the parameters.'        
        
    for niter in range(maxiter):
        
        print 'Iteration {} of {}:'.format(niter+1, maxiter)       
        
        # Evaluate the solution.
        xwcs, ywcs = astrometry.world2wcs(wcspars, ra, dec)
        xpix, ypix = astrometry.wcs2pix(polpars, xwcs, ywcs)
        
        # Match positions.
        idx1, idx2, dist = match(xccd, yccd, xpix, ypix)

        # Reject objectively bad matches.
        mask = dist < 20.
        idx1, idx2, dist = idx1[mask], idx2[mask], dist[mask]
        
        # Reject based on quality of remaining matches.
        mask = dist < 7.5*np.median(dist) 
        idx1, idx2, dist = idx1[mask], idx2[mask], dist[mask]
        
        print '  N = {:04d}, MAD = {:.3f}'.format(len(idx1), np.median(dist))
        
        # Solve the WCS parameters.
        wcspars, dx, dy, dr = astrometry.wcs_solve(wcspars, header['lst'], xccd[idx1], yccd[idx1], ra[idx2], dec[idx2])
        print '  WCSonly: dx = {:.3f}, dy = {:.3f}, dr = {:.3f}'.format(dx, dy, dr)
        
        # Solve the distortions.
        polpars, dx, dy, dr = astrometry.pol_solve(wcspars, polpars, xccd[idx1], yccd[idx1], ra[idx2], dec[idx2])
        print '  WCS+POL: dx = {:.3f}, dy = {:.3f}, dr = {:.3f}'.format(dx, dy, dr)
        
    if debug:

        xwcs, ywcs = astrometry.world2wcs(wcspars, ra, dec)
        xpix, ypix = astrometry.wcs2pix(polpars, xwcs, ywcs)
        
        plt.subplot(121, aspect='equal')
        
        plt.plot(xccd[idx1] - xwcs[idx2], yccd[idx1] - ywcs[idx2], 'k.')
        plt.xlim(-2, 2)
        plt.ylim(-2, 2)        
        
        plt.subplot(122, aspect='equal')
        
        plt.plot(xccd[idx1] - xpix[idx2], yccd[idx1] - ypix[idx2], 'k.')
        plt.xlim(-2, 2)
        plt.ylim(-2, 2)
         
        plt.tight_layout()
        plt.show()
        plt.close()
        
        vmin = np.nanpercentile(image, 1)
        vmax = np.nanpercentile(image, 99)  
        
        plt.imshow(image, origin='lower', interpolation='None', vmin=vmin, vmax=vmax, cmap=plt.cm.Greys)
        
        plt.plot(xwcs, ywcs, 'ro', label='WCS')
        plt.plot(xccd, yccd, 'go', label='CCD')
        plt.plot(xpix, ypix, 'bo', label='PIX')
        plt.xlim(-0.5, nx-0.5)
        plt.ylim(-0.5, ny-0.5)
        
        plt.legend()
        
        plt.tight_layout()
        plt.show() 
        plt.close()
    
    return wcspars, polpars

def astromaster(imagefile, pars0, cat, astromask=None, darkfile=None, outpath='.'):    
       
    # Read the image.
    image, header = fits.getdata(imagefile, header=True)
    
    # Check for the LST.
    try:
        header['lst']
    except KeyError:
        raise
        
    # Check for the lstsequence.
    try:
        header['lstseq']
    except KeyError:
        header['lstseq'] = -1    
    
    # Remove the overscan region.
    try:
        lx = header['X0']
        ux = lx + header['XSIZE']
        ly = header['Y0']
        uy = ly + header['YSIZE'] 
        
        image = image[ly:uy,lx:ux]    
    except:
        pass

    # Convert to float. 
    image = image.astype('float32')    
    
    # Subtract the dark.
    if darkfile is not None:
        
        dark = fits.getdata(darkfile)
        image = image - dark   
    
    # Read the catalogue.
    starcat = io.read_catalogue(cat)
    
    # Select suitably bright and isolated stars.
    mask = (starcat['vmag'] < 7.5) & (starcat['vmag'] > 4.) & (starcat['dist8'] > 10.)    
    ra, dec = starcat['ra'][mask], starcat['dec'][mask]    

    # Create a master solution.
    wcspars, polpars = blind_solution(image, header, ra, dec, pars0, debug=True)
    
    # Refine the mastersolution.
    astro = astrometry.Astrometry(wcspars, polpars, astromask)
    aflag, astrosol = astro.solve(image, header, ra, dec)
    
    # Plot the result.
    xpix, ypix, mask = astro.world2pix(header['lst'], ra, dec, jd=header['jd'])    
    
    ny, nx = image.shape
    vmin = np.nanpercentile(image, 1)
    vmax = np.nanpercentile(image, 99)  
    
    plt.imshow(image, origin='lower', interpolation='None', vmin=vmin, vmax=vmax, cmap=plt.cm.Greys)
    plt.plot(xpix, ypix, 'bo')
    plt.xlim(-0.5, nx-0.5)
    plt.ylim(-0.5, ny-0.5)
    plt.show()       
    
    # Write the result to file.
    wcspars = astro.get_wcspars()
    polpars = astro.get_polpars()
    
    head, tail = os.path.split(imagefile)    
    filename = tail.rsplit('.')[0] + '.hdf5'
    filename = os.path.join(outpath, filename)    
    
    io.write_astromaster(filename, wcspars, polpars, astromask) 
    
    return
