# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 15:54:55 2017

@author: talens
"""

import os
import logging

import h5py
import numpy as np

from lsreduce import astrometry, io, misc
        
class GetSystematics(object):

    def __init__(self, filename):
        
        f = io.SysFile(filename)
    
        ascc, nobs, mag, sigma = f.read_magnitudes()    
        self.magnitudes = dict()
        self.magnitudes['ascc'] = ascc
        self.magnitudes['nobs'] = nobs
        self.magnitudes['mag'] = mag
        self.magnitudes['sigma'] = sigma
    
        pg, nobs, trans = f.read_transmission()
        self.pgcam = pg
        self.transmission = dict()
        self.transmission['nobs'] = nobs
        self.transmission['trans'] = trans
        
        pg, nobs, amplitudes = f.read_intrapix()  
        self.pgipx = pg
        self.intrapix = dict()
        self.intrapix['nobs'] = nobs
        self.intrapix['amplitudes'] = amplitudes
        
        hg, nobs, clouds, sigma, lstmin, lstmax = f.read_clouds()        
        self.hg = hg
        self.clouds = dict()
        self.clouds['nobs'] = nobs
        self.clouds['clouds'] = clouds
        self.clouds['sigma'] = sigma
        self.lstmin = lstmin        
        
        return
        
    def get_magnitudes(self, ascc):

        mask = self.magnitudes['ascc'] == ascc
        mag = self.magnitudes['mag'][mask]
        nobs = self.magnitudes['nobs'][mask]
        sigma = self.magnitudes['sigma'][mask]

        return mag, nobs, sigma        
        
    def get_transmission(self, ra, dec, lst):

        ha = astrometry.ra2ha(ra, lst)
        dec = np.repeat(dec, len(lst))

        idx1, idx2 = self.pgcam.radec2idx(ha, dec)

        trans = self.transmission['trans'][idx1,idx2]
        nobs = self.transmission['nobs'][idx1,idx2]

        return trans, nobs
        
    def get_intrapix(self, ra, dec, lst, x, y):
        
        ha = astrometry.ra2ha(ra, lst)
        dec = np.repeat(dec, len(lst))

        idx1, idx2 = self.pgipx.radec2idx(ha, dec)
        
        ipx_x = self.intrapix['amplitudes'][idx1,idx2,0]*np.sin(2*np.pi*x) + self.intrapix['amplitudes'][idx1,idx2,1]*np.cos(2*np.pi*x)
        ipx_y = self.intrapix['amplitudes'][idx1,idx2,2]*np.sin(2*np.pi*y) + self.intrapix['amplitudes'][idx1,idx2,3]*np.cos(2*np.pi*y)        
        nobs = self.intrapix['nobs'][idx1,idx2]        
        
        return ipx_x + ipx_y, nobs
        
    def get_clouds(self, ra, dec, lstseq):
      
        idx1 = self.hg.radec2idx(ra, dec)
        idx2 = lstseq - self.lstmin
        
        clouds = self.clouds['clouds'][idx1,idx2]
        nobs = self.clouds['nobs'][idx1,idx2]
        sigma = self.clouds['sigma'][idx1,idx2]
      
        return clouds, nobs, sigma
        
    def get_systematics(self, ascc, ra, dec, lstseq, lst, x, y):

        flag = np.zeros(len(lstseq), dtype='uint8')
        
        mag, nobs, sigma = self.get_magnitudes(ascc)        
        
        trans, nobs = self.get_transmission(ra, dec, lst)
        flag = np.where(nobs < 25, flag+2, flag)
        
        ipx, nobs = self.get_intrapix(ra, dec, lst, x, y)
        flag = np.where(nobs < 25, flag+4, flag)
        
        clouds, nobs, sigma = self.get_clouds(ra, dec, lstseq)
        flag = np.where(nobs < 25, flag+8, flag)
        flag = np.where(sigma > .05, flag+16, flag)
        
        systematics = trans + ipx + clouds
        flag = np.where(np.isnan(systematics), flag + 1, flag)
        
        return mag, trans, ipx, clouds, flag     
    
    
def idxstats(indices, values, statistic='mean', keeplength=False):
    """ Compute a statistic for all values with the same index.
    
    Args:
        indices (int): An array of indices.
        values (float): An array of values.
        statistic (string or function): The statistic to compute on values
            that have the same index may be any of 'mean', 'std', 'count',
            'sum', 'median' or a function. Default is 'mean'.
        keeplength (bool): If True the return will have the same shape as
            values, otherwise it will be the same length as the number of 
            unique indices. Default is False.
            
    Returns:
        result: The statistic computed on values with the same index.
    
    """
    
    # Check that the statistic is valid.
    known_stats = ['mean', 'std', 'count', 'sum', 'median']
    if not callable(statistic) and statistic not in known_stats:
        raise ValueError('invalid statistic {}'.format(statistic))
    
    # Make sure the indices are integers.
    indices = indices.astype('int') # It would be better to return an error if they are not integer, but how...
    
    # If only one scalar index is given bincount will not work.
    if np.isscalar(indices):
        
        if statistic == 'std':
            return 0.
        if statistic == 'count':
            return 1
        elif statistic in ['mean', 'sum', 'median']:
            return values
        else:
            return statistic(values)
    
    # Count the number of datapoints at each index.
    flatcount = np.bincount(indices)
    
    # Obtain the unique indices
    a = flatcount.nonzero()
    
    # Compute the desired statistic.
    if statistic == 'mean':
        flatsum = np.bincount(indices, values)
        result = flatsum / flatcount
        
    elif statistic == 'std':
        flatsum = np.bincount(indices, values)
        flatsum2 = np.bincount(indices, values ** 2)
        result = np.sqrt(flatsum2 / flatcount - (flatsum / flatcount) ** 2)

    elif statistic == 'count':
        result = flatcount
        
    elif statistic == 'sum':
        result = np.bincount(indices, values)
        
    elif statistic == 'median':
        result = np.zeros(len(flatcount))
        for i in np.unique(indices):
            result[i] = np.median(values[indices == i])
        
    else:
        result = np.zeros(len(flatcount))
        for i in np.unique(indices):
            result[i] = statistic(values[indices == i])
    
    # Return the statistic in the desired format.
    if not keeplength: 
        return result[a]
    else:
        return result[indices]
    
    
def apply_calibration(photfile, aper=0, sysfile=None, redfile=None):

    log = logging.getLogger('bringreduce').getChild('calibration')

    # File and aperture to work on.
    if not os.path.isfile(photfile):
        raise IOError('Photometry file not found: {}'.format(photfile))
    else:
        log.info('Applying corrections to aperture {} of file: {}'.format(aper, photfile))
    
    # The systematics file.
    if sysfile is None:
        head, tail = os.path.split(photfile)
        prefix = 'sys{}_vmag_'.format(aper)
        tail = prefix + tail.rsplit('_')[-1]
        sysfile = os.path.join(head, tail)
    
    if not os.path.isfile(sysfile):
        raise IOError('Systematics file not found: {}'.format(sysfile))
    else:
        log.info('Reading corrections from: {}'.format(sysfile))    
    
    # The output file.
    if redfile is None:
        head, tail = os.path.split(photfile)
        prefix = 'red{}_vmag_'.format(aper)
        tail = prefix + tail.rsplit('_')[-1]
        redfile = os.path.join(head, tail)
    
    if os.path.isfile(redfile):
        raise IOError('Output file already exists: {}'.format(redfile))
    else:
        log.info('Writing results to: {}'.format(redfile))

    # Read the stars and station fields from the photometry.
    f = io.PhotFile(photfile)

    stars = f.read_stars()
    station = f.read_station()
    
    # Open the systematics file.
    sys = GetSystematics(sysfile)

    # Fields and datatypes for the binned lightcurves.
    lightcurves = dict()
    names = ['lstseq', 'nobs', 'lst', 'jd', 'exptime', 'x', 'y', 'sky', 'esky', 'mag{}'.format(aper), 'emag{}'.format(aper), 'trans{}'.format(aper), 'etrans{}'.format(aper), 'clouds{}'.format(aper), 'eclouds{}'.format(aper)]
    formats = ['uint32', 'uint8', 'float64', 'float64', 'float32', 'float32', 'float32', 'float64', 'float64', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32']
       
    for i in range(len(stars['ascc'])):        
        
        # Read the lightcurve.
        lc = f.read_lightcurves(ascc=stars['ascc'][i])

        # Remove flagged datapoints.
        mask = (lc['aflag'] == 0) & (lc['pflag'] == 0)
        lc = lc[mask]

        if (len(lc) == 0):
            stars['nobs'][i] = 0
            continue 

        # Perform look-up in station field.
        idx = np.searchsorted(station['lstseq'], lc['lstseq']) 

        # Compute the corrected lightcurve.
        mag, emag = misc.flux2mag(lc['flux{}'.format(aper)]/station['exptime'][idx], lc['eflux{}'.format(aper)]/station['exptime'][idx])
        mag0, trans, ipx, clouds, cflag = sys.get_systematics(stars['ascc'][i], stars['ra'][i], stars['dec'][i], lc['lstseq'], station['lst'][idx], lc['x'], lc['y'])
        mag = mag - trans - ipx - clouds
        
        # Remove flagged datapoints.
        mask = (cflag == 0) 
        lc = lc[mask]        
        
        if (len(lc) == 0):
            stars['nobs'][i] = 0
            continue
        
        jd = station['jd'][idx][mask]
        lst = station['lst'][idx][mask]
        exptime = station['exptime'][idx][mask]
        
        mag = mag[mask]
        emag = emag[mask]
        trans = trans[mask]
        ipx = ipx[mask]
        clouds = clouds[mask]            
        
        # Compute the final binned lightcurve.
        lstseq, binidx, nobs = np.unique(lc['lstseq']//50, return_inverse=True, return_counts=True) 
        
        lc_bin = np.recarray(len(lstseq), names=names, formats=formats)
        
        lc_bin['lstseq'] = lstseq
        lc_bin['nobs'] = nobs # Number of raw points used for each binned point.
 
        lc_bin['jd'] = idxstats(binidx, jd, statistic='mean')
        lc_bin['lst'] = idxstats(binidx, lst, statistic='mean')
        lc_bin['exptime'] = idxstats(binidx, exptime, statistic='sum')            
        
        lc_bin['x'] = idxstats(binidx, lc['x'], statistic='mean')
        lc_bin['y'] = idxstats(binidx, lc['y'], statistic='mean')            
        
        lc_bin['mag{}'.format(aper)] = idxstats(binidx, mag, statistic='mean')
        lc_bin['emag{}'.format(aper)] = idxstats(binidx, mag, statistic='std')/np.sqrt(nobs)
        lc_bin['sky'] = idxstats(binidx, lc['sky'], statistic='mean')
        lc_bin['esky'] = idxstats(binidx, lc['sky'], statistic='std')/np.sqrt(nobs)
        
        lc_bin['trans{}'.format(aper)] = idxstats(binidx, trans, statistic='mean')
        lc_bin['etrans{}'.format(aper)] = idxstats(binidx, trans, statistic='std')/np.sqrt(nobs)
        lc_bin['clouds{}'.format(aper)] = idxstats(binidx, clouds, statistic='mean')
        lc_bin['eclouds{}'.format(aper)] = idxstats(binidx, clouds, statistic='std')/np.sqrt(nobs)

        lightcurves[stars['ascc'][i]] = lc_bin
    
        stars['nobs'][i] = len(lstseq)
        stars['lstsqmin'][i] = lstseq[0]
        stars['lstsqmax'][i] = lstseq[-1]
                
    with h5py.File(redfile) as g:
        
        idx, = np.where(stars['nobs'] > 0)            
        
        grp = g.create_group('stars')
        for key in stars.keys():
            grp.create_dataset(key, data=stars[key][idx])
            
        grp = g.create_group('lightcurves')
        for key in lightcurves.keys():
            grp.create_dataset(key, data=lightcurves[key])
        
    return
    
