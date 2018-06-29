#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import h5py
import numpy as np

import logging
import multiprocessing as mp

from . import misc
from . import io
from . import grids
from . import sigmas
from . import astrometry

def cdecor_spatial(idx2, idx3, value, error, x, y, maxiter=100, dtol=1e-3, verbose=False):
    """ Perform a coarse decorrelation with intrapixel variations.
    
    Args:
        idx1 (int): Indices along which to calculate the first parameter.
        idx2 (int): Indices along which to calculate the second parameter.
        idx3 (int): Indices along which to calculate the intrapixel variations.
        value (float): Values to fit.
        error (float): Measurement errors corresponding to the values.
        x (float): The x position corresponding to the values.
        y (float): The y position corresponding to the values. 
        maxiter (int): The maximum number of iterations to perform. Defualt
            is 100.
        dtol (float): Maximum allowed change in the parameters, iteration
            terminates if the cahnge falls below this value. Default is 1e-3.
        verbose (bool): Output the current iteration. Default is True.
        
    Returns:
        par1 (float): The parameters corresponding to idx1.
        par2 (float): The parameters corresponding to idx2.
        par3 (float): The amplitudes of the intrapixel variations corresponding
            to idx3.
        quality: A named tuple with fields niter, chisq, npoints and npars
            describing the number of iterations, the chi-square value, the
            number of datapoints and the number of parameters of the fit.
    
    """
    
    sort = np.argsort(idx3)
    idx2 = idx2[sort]
    idx3 = idx3[sort]
    value = value[sort]
    error = error[sort]
    x = x[sort]
    y = y[sort]
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(value)
    npars2 = np.amax(idx2) + 1
    npars3 = 4*(np.amax(idx3) + 1)
    npars = npars2 + npars3
    
    # Create arrays.
    weights = 1./error**2
    par2 = np.zeros(npars2)
    par3 = np.zeros((npars3/4, 4))
    
    strides = np.cumsum(np.bincount(idx3))
    strides = np.append(0, strides)
    mat = np.vstack([np.sin(2*np.pi*x), np.cos(2*np.pi*x), np.sin(2*np.pi*y), np.cos(2*np.pi*y)]).T
    ipx = np.zeros(len(value))
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = {}'.format(niter)
        
        # Compute the parameters.
        par2 = np.bincount(idx2, weights*(value - ipx))/np.bincount(idx2, weights)
        
        res = value - par2[idx2]
        wsqrt = np.sqrt(weights)
        for i in range(npars3/4):
            par3[i] = np.linalg.lstsq(mat[strides[i]:strides[i+1],:]*wsqrt[strides[i]:strides[i+1]:,None], res[strides[i]:strides[i+1]]*wsqrt[strides[i]:strides[i+1]])[0]
            ipx[strides[i]:strides[i+1]] = np.sum(mat[strides[i]:strides[i+1],:]*par3[i], axis=1)
        
        # Check if the solution has converged.
        if (niter > 0):
            
            dcrit2 = np.nanmax(np.abs(par2 - par2_old))
            dcrit3 = np.nanmax(np.abs(par3 - par3_old))
            
            if (dcrit2 < dtol) & (dcrit3 < dtol):
                break
        
        par2_old = np.copy(par2)
        par3_old = np.copy(par3)
    
    # Compute the chi-square of the fit.
    chisq = weights*(value - par2[idx2] - ipx)**2        
    chisq = np.sum(chisq)
    
    quality = {'niter':niter, 'chisq':chisq, 'npoints':npoints, 'npars':npars}    
    
    return par2, par3, quality
    
def cdecor_temporal(idx1, idx2, value, error, sigma1, sigma2, maxiter=100, dtol=1e-3, verbose=False):
    """ Perform a coarse decorrelation with extra error terms.
    
    Args:
        idx1 (int): Indices along which to calculate the first parameter.
        idx2 (int): Indices along which to calculate the second parameter.
        error (float): Measurement errors corresponding to the values.
        sigma1 (float): Initial value for the extra error corresponding to
            idx1.
        sigma2 (float): Initial value for the extra error corresponding to
            idx2.
        maxiter (int): The maximum number of iterations to perform. Defualt
            is 100.
        dtol (float): Maximum allowed change in the parameters, iteration
            terminates if the cahnge falls below this value. Default is 1e-3.
        verbose (bool): Output the current iteration. Default is True.
        
    Returns:
        par1 (float): The parameters corresponding to idx1.
        par2 (float): The parameters corresponding to idx2.
        sigma1 (float): The extra error corresponding to idx1.
        sigma2 (float): The extra error corresponding to idx2.
        quality: A named tuple with fields niter, chisq, npoints and npars
            describing the number of iterations, the chi-square value, the
            number of datapoints and the number of parameters of the fit.
    
    """
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(value)
    npars1 = np.amax(idx1) + 1
    npars2 = np.amax(idx2) + 1
    npars = npars2
    
    # Create arrays.
    par1 = np.zeros(npars1)
    par2 = np.zeros(npars2)
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = {}'.format(niter)
            
        # Compute the parameters.
        sigma1 = sigmas.find_sigma(idx1, value - par2[idx2], error**2 + (sigma2**2)[idx2])
        par2, sigma2 = sigmas.find_par_sigma(idx2, value, error**2 + (sigma1**2)[idx1])
        
        # Check if the solution has converged.
        if (niter > 0):
            
            dcrit2 = np.nanmax(np.abs(par2 - par2_old))
            
            if (dcrit2 < dtol):
                break
        
        # Check if the solution is oscillating?
        if (niter > 1):
            
            dcrit2 = np.nanmax(np.abs(par2 - par2_older))
            
            if (dcrit2 < dtol):
                break
        
        if (niter > 0):
            par2_older = np.copy(par2_old)
        
        par2_old = np.copy(par2)
        
    # Compute the chi-square of the fit.
    chisq = (value - par2[idx2])**2/(error**2 + (sigma1**2)[idx1] + (sigma2**2)[idx2])     
    chisq = np.sum(chisq)
    
    quality = {'niter':niter, 'chisq':chisq, 'npoints':npoints, 'npars':npars}
    
    return par2, sigma1, sigma2, quality

def spatial_worker(in_queue, out_queue):
    
    while True:
        
        item = in_queue.get()
    
        if (item == 'DONE'):
            break
        else:
            idx, camtransidx, intrapixidx, idx2, idx4, mag, emag, x, y = item
            trans, amp, quality = cdecor_spatial(idx2, idx4, mag, emag, x, y)
            out_queue.put((idx, camtransidx, intrapixidx, trans, amp, quality))

    return
    
def temporal_worker(in_queue, out_queue):
    
    while True:
        
        item = in_queue.get()
        
        if (item == 'DONE'):
            break
        else:
            idx, staridx, lstseq, idx1, idx3, mag, emag, sig1, sig2 = item
            clouds, sig1, sig2, quality = cdecor_temporal(idx1, idx3, mag, emag, sig1, sig2)
            out_queue.put((idx, staridx, lstseq, clouds, sig1, sig2, quality))

    return

class CoarseDecorVmag(object):
    
    def __init__(self, photfile, aperture, sysfile=None, **kwargs):
        """ Perform a coarse decorrelation on all data in a given file."""
        
        self.log = logging.getLogger('bringreduce').getChild('calibration')        
        
        # fLC file and aperture to work on.
        self.photfile = photfile
        self.aper = aperture
        
        if not os.path.isfile(self.photfile):
            self.log.warning('File {} not found, terminating.'.format(self.photfile))
            return
        else:
            self.log.info('Calculating corrections for aperture {} of file {}'.format(self.aper, self.photfile))
        
        # The systematics file.
        if sysfile is None:
            head, tail = os.path.split(self.photfile)
            prefix = 'sys%i_vmag_'%self.aper
            tail = prefix + tail.rsplit('_')[-1]
            sysfile = os.path.join(head, tail)
        
        self.sysfile = sysfile
        
        if os.path.isfile(self.sysfile):
            self.log.warning('Systematics file {} already exists, terminating'.format(self.sysfile))
            return
        else:
            self.log.info('Writing results to {}'.format(self.sysfile))
        
        # Initialize with default parameters unless arguments were given.
        self.verbose = kwargs.pop('verbose', False)
        self.dtol = kwargs.pop('dtol', 1e-3)
        self.outer_maxiter = kwargs.pop('outer_maxiter', 5)
        self.inner_maxiter = kwargs.pop('inner_maxiter', 100)
        self.nprocs = kwargs.pop('nprocs', 4)
        
        self.camgrid = 'polar'
        self.camnx = kwargs.pop('camnx', 13500)
        self.camny = kwargs.pop('camny', 720)
        
        self.ipxgrid = 'polar'
        self.ipxnx = kwargs.pop('ipxnx', 270)
        self.ipxny = self.camny
    
        # Perform the coarse decorrelation.
        self._calculate()
    
        return
    
    def _read_data(self, index, mode):
        """ Read a portion of the data, create indices and remove flagged
        datapoints."""
        
        if (mode == 'spatial'):
            mask = (self.decidx == index)  
        elif (mode == 'temporal'):
            mask = (self.skyidx == index)
        else:
            self.log.warning('Unknown mode {}, how did you even get here?').format(mode)
            return
        
        ascc = self.stars['ascc'][mask]
        ra = self.stars['ra'][mask]
        dec = self.stars['dec'][mask]
        nobs = self.stars['nobs'][mask]
        
        staridx = self.staridx[mask]
        skyidx = self.skyidx[mask]
        
        # Read data.
        fields = ['lstseq', 'flux%i'%self.aper, 'eflux%i'%self.aper, 'x', 'y', 'pflag', 'aflag']
        lightcurves = self.f.read_lightcurves(ascc, fields, perstar=False)
        
        lstseq = lightcurves['lstseq']
        flux = lightcurves['flux%i'%self.aper]
        eflux = lightcurves['eflux%i'%self.aper]        
        x = lightcurves['x']
        y = lightcurves['y']
        aflag = lightcurves['aflag']
        pflag = lightcurves['pflag']
        
        station = self.f.read_station(['lst', 'exptime'], lstseq)        
        lst = station['lst']
        exptime = station['exptime']
        
        # Normalize fluxes with exptime.
        flux, eflux = flux/exptime, eflux/exptime
        
        # Compute the hour angle.
        ra = np.repeat(ra, nobs)
        dec = np.repeat(dec, nobs)
        ha = astrometry.ra2ha(ra, lst)
        
        # Create indices.
        staridx = np.repeat(staridx, nobs)
        camtransidx, decidx = self.camgrid.radec2idx(ha, dec)
        intrapixidx, decidx = self.ipxgrid.radec2idx(ha, dec)
        skyidx = np.repeat(skyidx, nobs)
        lstseq = lstseq - self.lstmin        
        
        # Remove bad data.
        mask = (aflag == 0) & (pflag == 0)
        
        flux = flux[mask]
        eflux = eflux[mask]
        x = x[mask]
        y = y[mask]
        
        staridx = staridx[mask]
        decidx = decidx[mask]
        camtransidx = camtransidx[mask]
        intrapixidx = intrapixidx[mask]
        skyidx = skyidx[mask]
        lstseq = lstseq[mask]
        
        # Convert flux to magnitudes.
        mag, emag = misc.flux2mag(flux, eflux)
        mag = mag - self.stars['vmag'][staridx]
        
        return mag, emag, x, y, staridx, decidx, camtransidx, intrapixidx, skyidx, lstseq
    
    def _spatial(self):
        """ Solve for the time-independent camera transmission and intrapixel 
        variations.
        """
        
        mngr = mp.Manager()
        in_queue = mp.Queue(2*self.nprocs)
        out_queue = mngr.Queue()
        the_pool = mp.Pool(self.nprocs, spatial_worker, (in_queue, out_queue))
        
        decidx = np.unique(self.decidx)
    
        for idx in decidx:
            
            # Read data.
            mag, emag, x, y, staridx, _, camtransidx, intrapixidx, skyidx, lstseq = self._read_data(idx, mode='spatial')
            
            if (len(mag) == 0): continue
            
            if not np.all(np.isfinite(mag)):
                self.log.warning('Bad mag in spatial solver on ring {}'.format(idx))
                
            if not np.all(np.isfinite(emag)):
                self.log.warning('Bad emag in spatial solver on ring {}'.format(idx))
                
            if not np.all(np.isfinite(x)):
                self.log.warning('Bad x in spatial solver on ring {}'.format(idx))
                
            if not np.all(np.isfinite(y)):
                self.log.warning('Bad y in spatial solver on ring {}'.format(idx))            
            
            # Apply temporal correction if known.
            if self.got_sky:
                mag = mag - self.clouds['clouds'][skyidx, lstseq]
                
                if not np.all(np.isfinite(mag)):
                    self.log.warning('Bad cloud corrected mag in spatial solver on ring {}'.format(idx))                
                
                emag = np.sqrt(emag**2 + self.magnitudes['sigma'][staridx]**2 + self.clouds['sigma'][skyidx, lstseq]**2)

                if not np.all(np.isfinite(emag)):
                    self.log.warning('Bad enhanced emag in spatial solver on ring {}'.format(idx))            
            
            # Create unique indices.
            camtransidx, idx2 = np.unique(camtransidx, return_inverse=True)
            intrapixidx, idx4 = np.unique(intrapixidx, return_inverse=True)
            
            self.trans['nobs'][camtransidx, idx] = np.bincount(idx2)
            self.intrapix['nobs'][intrapixidx, idx] = np.bincount(idx4)
            
            in_queue.put((idx, camtransidx, intrapixidx, idx2, idx4, mag, emag, x, y))
            
        for i in range(self.nprocs):
            in_queue.put('DONE')
            
        the_pool.close()
        the_pool.join()
                        
        out_queue.put('DONE')
        
        for item in iter(out_queue.get, 'DONE'):
            
            idx, camtransidx, intrapixidx, trans, amplitudes, quality = item
            
            if np.isnan(quality['chisq']):
                self.log.warning('Bad result in spatial solver on ring {}.'.format(idx))
                
            # Store results.
            self.trans['trans'][camtransidx, idx] = trans
            self.intrapix['amplitudes'][intrapixidx, idx] = amplitudes
            
            self.spatial['niter'][idx] = quality['niter']
            self.spatial['chisq'][idx] = quality['chisq']
            self.spatial['npoints'][idx] = quality['npoints']
            self.spatial['npars'][idx] = quality['npars'] 
            
        return
        
    def _temporal(self):
        """ Solve for the time-dependent sky transmission."""

        mngr = mp.Manager()
        in_queue = mp.Queue(2*self.nprocs)
        out_queue = mngr.Queue()
        the_pool = mp.Pool(self.nprocs, temporal_worker, (in_queue, out_queue))

        skyidx = np.unique(self.skyidx)
        
        for idx in skyidx:           
            
            # Read data.
            mag, emag, x, y, staridx, decidx, camtransidx, intrapixidx, _, lstseq = self._read_data(idx, mode='temporal')
            
            if (len(mag) == 0): continue
                
            if not np.all(np.isfinite(mag)):
                self.log.warning('Bad mag in temporal solver on pixel {}'.format(idx))
                
            if not np.all(np.isfinite(emag)):
                self.log.warning('Bad emag in temporal solver on pixel {}'.format(idx))
                
            if not np.all(np.isfinite(x)):
                self.log.warning('Bad x in temporal solver on pixel {}'.format(idx))
                
            if not np.all(np.isfinite(y)):
                self.log.warning('Bad y in temporal solver on pixel {}'.format(idx))
            
            # Apply known spatial correction.
            mag = mag - self.trans['trans'][camtransidx, decidx]
            
            if not np.all(np.isfinite(mag)):
                self.log.warning('Bad trans corrected mag in temporal solver on pixel {}'.format(idx))            
            
            mag = mag - np.sum(self.intrapix['amplitudes'][intrapixidx, decidx]*np.array([np.sin(2*np.pi*x), np.cos(2*np.pi*x), np.sin(2*np.pi*y), np.cos(2*np.pi*y)]).T, axis=1)
            
            if not np.all(np.isfinite(mag)):
                self.log.warning('Bad ipx corrected mag in temporal solver on pixel {}'.format(idx))            
            
            # Create unique indices.
            staridx, idx1 = np.unique(staridx, return_inverse=True)
            lstseq, idx3 = np.unique(lstseq, return_inverse=True)
            
            self.magnitudes['nobs'][staridx] = np.bincount(idx1)
            self.clouds['nobs'][idx, lstseq] = np.bincount(idx3)
            
            in_queue.put((idx, staridx, lstseq, idx1, idx3, mag, emag, self.magnitudes['sigma'][staridx], self.clouds['sigma'][idx, lstseq]))
            
        for i in range(self.nprocs):
            in_queue.put('DONE')
            
        the_pool.close()
        the_pool.join()
                        
        out_queue.put('DONE')
        
        for item in iter(out_queue.get, 'DONE'):
            
            idx, staridx, lstseq, clouds, sigma1, sigma2, quality = item
        
            if np.isnan(quality['chisq']):
                self.log.warning('Bad result in temporal solver on pixel {}.'.format(idx))     
        
            # Store results.
            self.magnitudes['sigma'][staridx] = sigma1
            self.clouds['clouds'][idx, lstseq] = clouds
            self.clouds['sigma'][idx, lstseq] = sigma2
            
            self.temporal['niter'][idx] = quality['niter']
            self.temporal['chisq'][idx] = quality['chisq']
            self.temporal['npoints'][idx] = quality['npoints']
            self.temporal['npars'][idx] = quality['npars']
            
        # Set got_sky to True.
        self.got_sky = True
            
        return
    
    def _calculate(self):
        """ Perform the coarse decorrelation."""

        # Set up the IO and coordinate grids.
        self.f = io.PhotFile(self.photfile)
        self.camgrid = grids.PolarGrid(self.camnx, self.camny)
        self.ipxgrid = grids.PolarGrid(self.ipxnx, self.ipxny)
        
        # Read data.
        self.stars = self.f.read_stars(['ascc', 'ra', 'dec', 'vmag', 'nobs'])
        self.station = self.f.read_station(['lstseq', 'lst'])        
        
        mask = self.stars['vmag'] <= 8.4
        for key in self.stars.keys():
            self.stars[key] = self.stars[key][mask]
        
        # Create indices.
        self.staridx = np.arange(len(self.stars['ascc']))
        _, self.decidx = self.camgrid.radec2idx(self.stars['ra'], self.stars['dec'])
        
        try:
            hg = grids.HealpixGrid(8)
        except:
            self.log.info('Failed to create instance of HealpixGrid, using catalogue look-up instead.')
            self.skyidx = io.read_skyidx(self.stars['ascc'])
        else:
            self.log.info('Using healpy for sky-index generation.')
            self.skyidx = hg.radec2idx(self.stars['ra'], self.stars['dec'])
        
        self.lstmin = np.amin(self.station['lstseq'])
        self.lstmax = np.amax(self.station['lstseq'])
        self.lstlen = self.lstmax - self.lstmin + 1 
        
        # The spatial calculation statistics.
        self.spatial = dict()
        self.spatial['niter'] = np.zeros(self.camgrid.ny+2, dtype='uint')
        self.spatial['chisq'] = np.full(self.camgrid.ny+2, fill_value=np.nan)
        self.spatial['npoints'] = np.zeros(self.camgrid.ny+2, dtype='uint')
        self.spatial['npars'] = np.zeros(self.camgrid.ny+2, dtype='uint')
        
        # The temporal calculation statistics.
        self.temporal = dict()
        self.temporal['niter'] = np.zeros(768, dtype='uint')
        self.temporal['chisq'] = np.full(768, fill_value=np.nan)
        self.temporal['npoints'] = np.zeros(768, dtype='uint')
        self.temporal['npars'] = np.zeros(768, dtype='uint')
        
        # The magnitudes.
        self.magnitudes = dict()
        self.magnitudes['ascc'] = self.stars['ascc']        
        self.magnitudes['vmag'] = self.stars['vmag']
        self.magnitudes['nobs'] = np.zeros(len(self.stars['ascc']), dtype='uint32')
        self.magnitudes['mag'] = self.stars['vmag']
        self.magnitudes['sigma'] = np.zeros(len(self.stars['ascc']), dtype='float32')
        
        # The transmission map.
        self.trans = dict()
        self.trans['nobs'] = np.zeros((self.camgrid.nx+2, self.camgrid.ny+2), dtype='uint32')
        self.trans['trans'] = np.full((self.camgrid.nx+2, self.camgrid.ny+2), fill_value=np.nan, dtype='float32')
        
        # The intrapixel amplitudes.
        self.intrapix = dict()
        self.intrapix['nobs'] = np.zeros((self.ipxgrid.nx+2, self.ipxgrid.ny+2), dtype='uint32')
        self.intrapix['amplitudes'] = np.full((self.ipxgrid.nx+2, self.ipxgrid.ny+2, 4), fill_value=np.nan, dtype='float32')
        
        # The cloud corrections.
        self.clouds = dict()
        self.clouds['nobs'] = np.zeros((768, self.lstlen), dtype='uint32')
        self.clouds['clouds'] = np.full((768, self.lstlen), fill_value=np.nan, dtype='float32')
        self.clouds['sigma'] = np.zeros((768, self.lstlen), dtype='float32')

        self.got_sky = False
        
        # Perform the coarse decorrelation.
        self.log.info('Performing coarse decorrelation for: {}'.format(self.photfile))
        for niter in range(self.outer_maxiter):
        
            self.log.info('Iteration {} out of {}:'.format(niter + 1, self.outer_maxiter))
            
            self.log.info('    Calculating camera systematics...')
            self._spatial()
            
            self.log.info('    Calculating atmospheric systematics...')
            self._temporal()
        
        # Write the results to file.
        with h5py.File(self.sysfile) as f:
    
            # Write the header.
            hdr = f.create_group('header')
            
            hdr.attrs['dtol'] = self.dtol
            hdr.attrs['outer_maxiter'] = self.outer_maxiter
            hdr.attrs['inner_maxiter'] = self.inner_maxiter
            
            for key in self.spatial.keys():
                hdr.create_dataset('spatial/' + key, data=self.spatial[key])
            
            for key in self.temporal.keys():
                hdr.create_dataset('temporal/' + key, data=self.temporal[key])
            
            # Write the data.
            grp = f.create_group('data')
            
            # Write the magnitudes.
            for key in self.magnitudes.keys():
                grp.create_dataset('magnitudes/' + key, data=self.magnitudes[key])
            
            # Write the camera transmission.
            idx1, idx2 = np.where(self.trans['nobs'] > 0)            
            grp.create_dataset('transmission/idx1', data=idx1, dtype='uint32')
            grp.create_dataset('transmission/idx2', data=idx2, dtype='uint32')          
            
            for key in self.trans.keys():
                grp.create_dataset('transmission/' + key, data=self.trans[key][idx1,idx2])

            grp['transmission'].attrs['grid'] = 'polar'
            grp['transmission'].attrs['nx'] = self.camnx
            grp['transmission'].attrs['ny'] = self.camny
            
            # Write the intrapixel variations.
            idx1, idx2 = np.where(self.intrapix['nobs'] > 0)            
            grp.create_dataset('intrapix/idx1', data=idx1, dtype='uint32')
            grp.create_dataset('intrapix/idx2', data=idx2, dtype='uint32')
            
            for key in self.intrapix.keys():
                grp.create_dataset('intrapix/' + key, data=self.intrapix[key][idx1,idx2])
            
            grp['intrapix'].attrs['grid'] = 'polar'
            grp['intrapix'].attrs['nx'] = self.ipxnx
            grp['intrapix'].attrs['ny'] = self.ipxny
            
            # Write the sky transmission.
            idx, lstseq = np.where(self.clouds['nobs'] > 0)            
            grp.create_dataset('clouds/idx', data=idx, dtype='uint32')
            grp.create_dataset('clouds/lstseq', data=lstseq+self.lstmin, dtype='uint32')
            
            for key in self.clouds.keys():
                grp.create_dataset('clouds/' + key, data=self.clouds[key][idx, lstseq])
                
            grp['clouds'].attrs['grid'] = 'healpix'
            grp['clouds'].attrs['nx'] = 8
            grp['clouds'].attrs['lstmin'] = self.lstmin
            grp['clouds'].attrs['lstmax'] = self.lstmax
            grp['clouds'].attrs['lstlen'] = self.lstlen
        
        return

def main():
    
    return

if __name__ == '__main__':
    
    main()
