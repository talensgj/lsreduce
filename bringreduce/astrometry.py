# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 17:59:53 2016

@author: talens
"""

import logging

import numpy as np

from astropy import wcs

from . import misc
from . import find


###############################################################################
### Functions for determining the position of the sun and moon.
###############################################################################

def sun_position(siteinfo, t=None):
    """ Compute sun position. """

    import datetime
    from astropy.time import Time
    from astropy.coordinates import EarthLocation, get_body, AltAz

    # If no time is given use system time.
    if t is None:
        t = Time(datetime.datetime.utcnow(), scale='utc')
    else:
        t = Time(t, scale='utc', format='jd')
    
    # Location of observatory and Alt-Az coordinate frame.
    loc = EarthLocation.from_geodetic(siteinfo['lon'], siteinfo['lat'], siteinfo['height'])  
    frame = AltAz(obstime=t, location=loc)       
    
    # Coordinates of the sun.
    sun = get_body('sun', t)
    ra, dec = sun.ra.value, sun.dec.value
    sun = sun.transform_to(frame)
    alt = sun.alt.value

    return ra, dec, alt
    
def moon_position(siteinfo, t=None):
    """ Compute moon position. """    
    
    import datetime
    from astropy.time import Time
    from astropy.coordinates import EarthLocation, get_body, AltAz

    # If no time is given use system time.
    if t is None:
        t = Time(datetime.datetime.utcnow(), scale='utc')
    else:
        t = Time(t, scale='utc', format='jd')
    
    # Location of observatory and Alt-Az coordinate frame.
    loc = EarthLocation.from_geodetic(siteinfo['lon'], siteinfo['lat'], siteinfo['height'])  
    frame = AltAz(obstime=t, location=loc)       
    
    # Coordinates of the sun.
    moon = get_body('moon', t)
    ra, dec = moon.ra.value, moon.dec.value
    moon = moon.transform_to(frame)
    alt = moon.alt.value
    
    return ra, dec, alt
    
def last_sunset(siteinfo, t=None):    
    """ Compute time of previous sunset. """    
    
    import datetime
    import ephem    
    
    site = ephem.Observer()
    site.lat = siteinfo['lat']*np.pi/180.
    site.lon = siteinfo['lon']*np.pi/180.
    site.elevation = siteinfo['height']
    site.horizon = '-10'

    if t is None:    
        site.date = datetime.datetime.utcnow()
    else:
        site.date = t
    
    sunset = site.previous_setting(ephem.Sun(), use_center=True).datetime()
    
    return sunset   
    
def next_sunrise(siteinfo, t=None):
    """ Compute time of next sunrise. """ 
    
    import datetime
    import ephem    
    
    site = ephem.Observer()
    site.lat = siteinfo['lat']*np.pi/180.
    site.lon = siteinfo['lon']*np.pi/180.
    site.elevation = siteinfo['height']
    site.horizon = '-10'    
        
    if t is None:    
        site.date = datetime.datetime.utcnow()
    else:
        site.date = t 
    
    sunrise = site.next_rising(ephem.Sun(), use_center=True).datetime()
    
    return sunrise    
    
###############################################################################
### 2-D Polynomial helper functions.
###############################################################################    
    
def poly2d_mat(x, y, order=6):
    """ Create a matrix for fitting 2-D polynomials of the given order. """     
    
    mat1 = np.vander(x, N=order+1, increasing=True)
    mat2 = np.vander(y, N=order+1, increasing=True)
        
    idx1, idx2 = np.indices((order+1, order+1))
    mask =((idx1 + idx2) <= order)
    idx1, idx2 = idx1[mask], idx2[mask]        
        
    mat = mat1[:,idx1]*mat2[:,idx2]
    
    return mat, idx1, idx2

def poly2d_eval(x, y, a, b, order=6):
    """ Evaluate 2-D polynomials. """
    
    mat, idx1, idx2 = poly2d_mat(x, y, order)
    
    xres = np.dot(mat, a[idx1, idx2])
    yres = np.dot(mat, b[idx1, idx2])
    
    return xres, yres   
    
def poly2d_solve(x1, y1, x2, y2, order=6):
    """ Solve for the coefficients of 2-D polynomials. """    
    
    mat, idx1, idx2 = poly2d_mat(x1, y1, order)  
    
    xtmp = np.linalg.lstsq(mat, x2)[0]
    ytmp = np.linalg.lstsq(mat, y2)[0]    
    
    a = np.zeros((order+1, order+1))
    b = np.zeros((order+1, order+1))  

    a[idx1, idx2] = xtmp
    b[idx1, idx2] = ytmp    
    
    return a, b
    
###############################################################################
### Coordinate transormations.
###############################################################################

def ra2ha(ra, lst):
    """ Convert Right Ascension to Hour Angle. """
    
    ha = np.mod(ra - lst*15., 360.)

    return ha

def ha2ra(ha, lst):
    """ Convert Hour Angle to Right Ascension. """    
    
    ra = np.mod(ha + lst*15., 360.)

    return ra

def create_wcs(wcspars, lst=None):
    """ Create and astropy WCS instance from a dictionary of parameters. """    
    
    if lst is not None:
        ra0, dec0 = wcspars['crval']
        ha0 = ra2ha(ra0, wcspars['lst'])
        ra0 = ha2ra(ha0, lst)
        crval = np.array([ra0, dec0])
        
    else:
        crval = wcspars['crval']
            
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = wcspars['crpix']
    w.wcs.cdelt = wcspars['cdelt']
    w.wcs.crval = crval
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.pc = wcspars['pc']
    
    return w
    
def world2wcs(wcspars, ra, dec, lst=None):
    """ Convert world coordinates to WCS-only coordinates. """    
        
    w = create_wcs(wcspars, lst)
    xwcs, ywcs = w.wcs_world2pix(ra, dec, 0)
    
    return xwcs, ywcs
    
def wcs2world(wcspars, xwcs, ywcs, lst=None):
    """ Convert WCS-only coordinates to world coordinates. """    
    
    w = create_wcs(wcspars, lst)
    ra, dec = w.wcs_pix2world(xwcs, ywcs, 0)
    
    return ra, dec    

def wcs2pix(polpars, xwcs, ywcs, crpix):
    """ Convert WCS-only coordinates to pixel coordinates. """     
    
    U, V = (xwcs - crpix[0]), (ywcs - crpix[1])

    dU, dV = poly2d_eval(U/polpars['nx'], V/polpars['ny'], polpars['x_wcs2pix'], polpars['y_wcs2pix'], polpars['order'])    
    
    u, v = U + dU, V + dV
    
    xpix, ypix = u + crpix[0], v + crpix[1]
    
    return xpix, ypix
    
def pix2wcs(polpars, xpix, ypix, crpix):
    """ Convert pixel coordinates to WCS-only coordinates. """    
    
    u, v = xpix - crpix[0], ypix - crpix[1]
    
    du, dv = poly2d_eval(u/polpars['nx'], v/polpars['ny'], polpars['x_pix2wcs'], polpars['y_pix2wcs'], polpars['order']) 
    
    U, V = u + du, v + dv
    
    xwcs, ywcs = U + crpix[0], V + crpix[1] 
    
    return xwcs, ywcs

###############################################################################
### Functions to solve for transformation parameters.
###############################################################################

def wcs_param(crval, cdelt, xpix, ypix, ra, dec):
    """ Compute the optimal crpix and pc values given crval and cdelt. """    
    
    crval[0] = np.mod(crval[0], 360.)    
    
    wcspars = dict()
    wcspars['crval'] = crval
    wcspars['crpix'] = np.zeros(2)    
    wcspars['cdelt'] = cdelt
    wcspars['pc'] = np.eye(2)
    
    w = create_wcs(wcspars)
    xwcs, ywcs = w.wcs_world2pix(ra, dec, 0)    
    
    a, b = poly2d_solve(xwcs, ywcs, xpix, ypix, 1)    
  
    wcspars['crpix'] = np.array([a[0,0], b[0,0]])  
  
    pci = np.zeros((2,2))    
    pci[0,0] = a[1,0]
    pci[0,1] = a[0,1]
    pci[1,0] = b[1,0]
    pci[1,1] = b[0,1]    
    
    wcspars['pc'] = np.linalg.inv(pci)    
    
    return wcspars

def wcs_chisq(crval, cdelt, xpix, ypix, ra, dec):
    """ Compute the chi-square statistic given crval and cdelt. """    
    
    if (crval[1] < -90) | (crval[1] > 90):
        return 1e24 
    
    try:
        wcspars = wcs_param(crval, cdelt, xpix, ypix, ra, dec)    
    except:
        return 1e24

    try:
        xwcs, ywcs = world2wcs(wcspars, ra, dec)
    except:
        return 1e24  
    
    chisq = np.sum((xpix - xwcs)**2 + (ypix - ywcs)**2)
    
    return chisq

def wcs_solve(wcspars, lst, xpix, ypix, ra, dec):
    """ Find the best-fit WCS parameters. """    
    
    from scipy.optimize import minimize    
    
    # Compute initial values for crval.
    cdelt = wcspars['cdelt']
    ra0, dec0 = wcspars['crval']
    ha0 = ra2ha(ra0, wcspars['lst'])
    ra0 = ha2ra(ha0, lst)
    crval = np.array([ra0, dec0])
        
    # Find the best fit values for crval.
    res = minimize(wcs_chisq, crval, args=(cdelt, xpix, ypix, ra, dec), bounds=[(None, None), (-90, 90)])
    
    # Get the new WCS parameters and add the LST.
    wcspars = wcs_param(res.x, cdelt, xpix, ypix, ra, dec)    
    wcspars['lst'] = lst    
    
    # Compute quality of the solution.
    xwcs, ywcs = world2wcs(wcspars, ra, dec, lst=lst)    
    
    dx = np.std(xpix - xwcs)
    dy = np.std(ypix - ywcs)
    dr = np.sqrt(np.mean((xpix - xwcs)**2 + (ypix - ywcs)**2))      
    
    return wcspars, dx, dy, dr
    
def pol_solve(wcspars, polpars, xpix, ypix, ra, dec):   
    """ Find the best-fit polynomial distortions. """    
    
    # Compute wcs-only coordinates.
    xwcs, ywcs = world2wcs(wcspars, ra, dec)
    
    # Compute intermediate coordinates.
    crpix = wcspars['crpix']
    U, V = xwcs - crpix[0], ywcs - crpix[1]
    u, v = xpix - crpix[0], ypix - crpix[1]
    
    # Solve for the pix2wcs and wcs2pix transformations.
    x_pix2wcs, y_pix2wcs = poly2d_solve(u/polpars['nx'], v/polpars['ny'], U - u, V - v, polpars['order'])   
    x_wcs2pix, y_wcs2pix = poly2d_solve(U/polpars['nx'], V/polpars['ny'], u - U, v - V, polpars['order'])
    
    polpars['x_pix2wcs'] = x_pix2wcs
    polpars['y_pix2wcs'] = y_pix2wcs
    polpars['x_wcs2pix'] = x_wcs2pix
    polpars['y_wcs2pix'] = y_wcs2pix    
    
    # Compute quality of the solution.
    xpix1, ypix1 = wcs2pix(polpars, xwcs, ywcs, crpix)    
    
    dx = np.std(xpix - xpix1)
    dy = np.std(ypix - ypix1)
    dr = np.sqrt(np.mean((xpix - xpix1)**2 + (ypix - ypix1)**2))    
    
    return polpars, dx, dy, dr
 
###############################################################################
### Class for using the astrometry starting from a master solution.
###############################################################################
   
class Astrometry(object):
    
    def __init__(self, wcspars, polpars):
        
        self.wcspars = wcspars
        self.polpars = polpars      
        
        self._backup_solution()        
        
        return
        
    def get_wcspars(self):
        
        return self.wcspars
        
    def get_polpars(self):
        
        return self.polpars
        
    def _backup_solution(self):
        
        self.wcspars_bak = self.wcspars.copy()
        self.polpars_bak = self.polpars.copy()
        
        return
        
    def _reset_solution(self):
        
        self.wcspars = self.wcspars_bak.copy()
        self.polpars = self.polpars_bak.copy()
        
        return
        
    def world2pix(self, lst, ra, dec, nx=4008, ny=2672, margin=50):
        
        # Compute wcs-only pixel coordinates.
        xwcs, ywcs = world2wcs(self.wcspars, ra, dec, lst)        
        
        # Remove coordinates not within the margins.
        mask = np.isfinite(xwcs) & misc.clip_rectangle(xwcs, ywcs, nx, ny, margin)
        xwcs, ywcs = xwcs[mask], ywcs[mask]        
        
        # Convert to actual pixel coordinates.
        xpix, ypix = wcs2pix(self.polpars, xwcs, ywcs, self.wcspars['crpix'])
        
        return xpix, ypix, mask
    
    def pix2world(self, lst, xpix, ypix):
        
        # Convert to wcs-only pixel coordinates.
        xwcs, ywcs = pix2wcs(self.polpars, xpix, ypix, self.wcspars['crpix'])
        
        # Convert to world coordinates.
        ra, dec = wcs2world(self.wcspars, xwcs, ywcs, lst)
        
        return ra, dec
        
    def solve(self, image, header, ra, dec, drate=.6, maxiter=10, margin=668, ntries=1):
        """ Try to create a new solution from the give image. 
        
            Flags:
            0 : All went well.
            1 : The image contained bright elements, but a solution was made.
            2 : The image is too bright no new solution was made.
            4 : Detection rate fell below the threshold, no solution was made.
        """
        log = logging.getLogger('bringreduce').getChild('astrometry')
        aflag = 0        
        
        # Check for bright elements.
        if (np.nanpercentile(image, 95) > 9500) | (np.nanstd(image) > 2000):
            
            log.info('Image may be bright, attempting a solution...')                       
            aflag += 1        
        
        # Check for extreme amounts of bright elements.
        if (np.nanmedian(image) > 20000):
            
            log.info('Image is bright, terminating...')                        
            aflag += 2
            
            return aflag, dict()
        
        # Evaluate the solution.
        xpix, ypix, mask = self.world2pix(header['lst'], ra, dec)
        ra, dec = ra[mask], dec[mask]

        fstars = len(ra)
        
        mask = (np.abs(xpix - 2004) < margin) & (np.abs(ypix - 1336) < margin)
        fstars_c = len(ra[mask])        
        
        for niter in range(maxiter):        
            
            log.info('Solving the astrometry iteration {} of {}.'.format(niter+1, maxiter))        
            
            # Search for stars.
            xpix, ypix, mask = self.world2pix(header['lst'], ra, dec)
            ra_, dec_ = ra[mask], dec[mask]            
        
            xccd, yccd, flag = find.gaussfit_mp(image, xpix, ypix)  
           
            mask = (flag == 0) 
            ustars = sum(mask)

            mask_c = (flag == 0) & (np.abs(xpix - 2004) < margin) & (np.abs(ypix - 1336) < margin)
            ustars_c = np.sum(mask_c)

            log.info('Found {} of {} expected stars.'.format(ustars, fstars))          
            
            # Check the detection rate.
            if (ustars >= drate*fstars):
                
                log.info('Current detection rate is {:.2f}'.format(float(ustars)/fstars))
                
            elif (ustars_c >= drate*fstars_c) & (ntries > 0): 
                
                log.info('Central detection rate is {:.2f}. Attempting a solution...'.format(float(ustars_c)/fstars_c))
                ntries -= 1
                
            else:
            
                log.info('Too few stars detected, terminating...')                
                aflag += 4
                
                self._reset_solution()
                
                return aflag, dict()
                    
            # Solve the WCS parameters.
            log.info('Solving for the WCS transformation.')
            wcspars, dx, dy, dr = wcs_solve(self.wcspars, header['lst'], xccd[mask], yccd[mask], ra_[mask], dec_[mask])
            log.info('Best fit WCS transformation gives dx={:.3f}, dy={:.3f}, dr={:.3f}'.format(dx, dy, dr))

            # Solve the polynomial coefficients.
            log.info('Solving for the polynomial distortions.')
            polpars, dx, dy, dr = pol_solve(wcspars, self.polpars, xccd[mask], yccd[mask], ra_[mask], dec_[mask])
            log.info('Best fit WCS+pol transformation gives dx={:.3f}, dy={:.3f}, dr={:.3f}'.format(dx, dy, dr))     
            
            self.wcspars = wcspars
            self.polpars = polpars            
            
        # Succesfully created a new solution.
        self._backup_solution()

        # Dictionary describing the solution.
        astrosol = dict()
        
        # Parameters.
        astrosol['lstseq'] = header['lstseq']
        astrosol['fstars'] = fstars
        astrosol['ustars'] = ustars
        astrosol['crval'] = wcspars['crval']
        astrosol['crpix'] = wcspars['crpix']
        astrosol['pc'] = wcspars['pc']
        astrosol['x_wcs2pix'] = polpars['x_wcs2pix']
        astrosol['y_wcs2pix'] = polpars['y_wcs2pix']
        astrosol['x_pix2wcs'] = polpars['x_pix2wcs']
        astrosol['y_pix2wcs'] = polpars['y_pix2wcs']
        
        # Scatter between predicted and detected coordinates.
        astrosol['dx'] = dx
        astrosol['dy'] = dy
        astrosol['dr'] = dr

        # World coordinates of fixed points on the CCD.
        xpix = np.array([1002., 1002., 2004., 3006., 3006.])
        ypix = np.array([668., 2004., 1336., 668., 2004.])           
        ra, dec = self.pix2world(header['lst'], xpix, ypix)            
        astrosol['ra0'] = ra[0]
        astrosol['ra1'] = ra[1]
        astrosol['ra2'] = ra[2]
        astrosol['ra3'] = ra[3]
        astrosol['ra4'] = ra[4]
        astrosol['dec0'] = dec[0]
        astrosol['dec1'] = dec[1]
        astrosol['dec2'] = dec[2]
        astrosol['dec3'] = dec[3]
        astrosol['dec4'] = dec[4]
        
        return aflag, astrosol
        