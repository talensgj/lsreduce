# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 17:21:45 2016

@author: talens
"""

import os
import logging

from multiprocessing import Pool
import time
import numpy as np
from threading import Thread
import Queue
import atexit

import h5py
from astropy.io import fits
    
from . import configuration as cfg
from . import grids

log = logging.getLogger('lsreduce')

_archivequeue = Queue.Queue()
_archiverclosed = True


def stop_archiver():
    if _archivequeue is not None:
        _archivequeue.put(None)

atexit.register(stop_archiver)


def _archive_worker(file, dest, delete):
    import os
    import traceback
    from astropy.io import fits
    try:
        out = os.path.join(dest, os.path.basename(file)) + '.gz'
        with fits.open(file) as f:
            f.writeto(out, clobber=True)
        if delete:
            os.remove(file)
    except Exception as e:
        pass #return e, traceback.format_exc()


def _archiver(queue, nworkers=5):
    """queue contains (file, destfolder) tuples. None is sentinal value and
    shuts  whole thing down"""
    global _archiverclosed

    pool = Pool(nworkers)

    try:
        while True:
            try:
                item = queue.get_nowait()
                if item is None:
                    try:
                        item = queue.get_nowait()
                    except Queue.Empty:
                        _archiverclosed = True
                        break

                pool.apply_async(_archive_worker, item)
            except Queue.Empty:
                pass

            time.sleep(.1)

    finally:
        pool.close()
        pool.join()


def read_siteinfo(filename, sitename):

    # Read the siteinfo file.
    dtype = [('sitename', '|S20'), ('lat', 'float32'), ('lon', 'float32'), ('height', 'float32'), ('ID', '|S2')]
    siteinfo = np.genfromtxt(filename, dtype=dtype)

    # Select the desired site.    
    mask = siteinfo['sitename'] == sitename
    siteinfo = siteinfo[mask]

    return siteinfo  
    
def read_targets(filename):
    
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Remove end of line and blank lines.
    lines = [line.rstrip('\n') for line in lines] 
    lines = [line for line in lines if line]
    
    # Split names and targets.
    targets = {line.rsplit(':')[0]:line.rsplit(':')[1] for line in lines}
    
    # Split targets.
    for key in targets.keys():
        line = targets[key]
        line = line.rsplit(',')
        line = [i.strip() for i in line] 
        targets[key] = line    
    
    return targets
  
def write_masterframe(filename, image, header, overwrite=True):
    
    image = image.astype('float32')    
    
    hdu = fits.PrimaryHDU(data=image, header=header)
    hdu.writeto(filename, clobber=overwrite)
    
    return
    
def update_darktable(filename, darktable, nmax=5):  
    
    data = np.loadtxt(darktable, dtype='S')    
    data = np.atleast_1d(data)
    
    if (len(data) == nmax):
        
        data = np.roll(data, 1)
        data[0] = filename
        
    else:
        
        data = np.append(filename, data)
        
    np.savetxt(darktable, data, fmt='%s')
    
    return
  
def read_masterdark(darktable):
    
    data = np.loadtxt(darktable, dtype='S')  
    
    try:
        filename = data[0]
    except:
        filename = data[None][0]
    
    return fits.getdata(filename, header=True)

def archive_files(filelist, destination, delete=True, nworkers=None):

    global _archiverclosed

    for f in filelist:
        _archivequeue.put((f, destination, delete))
    _archivequeue.put(None)

    if _archiverclosed:
        nworkers = 5 if nworkers is None else nworkers

        archive_thread = Thread(target=_archiver,
                                name='Raw Archiver Daemon',
                                args=(_archivequeue,),
                                kwargs={'nworkers': nworkers})
        archive_thread.start()
        _archiverclosed = False

    elif nworkers is not None:
        raise ValueError('nworkers may only be specified on the first call')

def read_astromaster(filename):

    wcspars = dict()
    polpars = dict()

    with h5py.File(filename, 'r') as f:
        
        grp = f['wcspars']
        wcspars['lst'] = grp.attrs['lst']
        wcspars['crval'] = grp['crval'].value
        wcspars['crpix'] = grp['crpix'].value
        wcspars['cdelt'] = grp['cdelt'].value
        wcspars['pc'] = grp['pc'].value
        
        grp = f['polpars']
        polpars['nx'] = grp.attrs['nx']
        polpars['ny'] = grp.attrs['ny']
        polpars['order'] = grp.attrs['order']
        polpars['x_wcs2pix'] = grp['x_wcs2pix'].value
        polpars['y_wcs2pix'] = grp['y_wcs2pix'].value
        polpars['x_pix2wcs'] = grp['x_pix2wcs'].value
        polpars['y_pix2wcs'] = grp['y_pix2wcs'].value        
        
    return wcspars, polpars
    
def write_astromaster(filename, wcspars, polpars):

    with h5py.File(filename) as f:
        
        grp = f.create_group('wcspars')
        grp.attrs['lst'] = wcspars['lst']
        grp.create_dataset('crval', data=wcspars['crval'])
        grp.create_dataset('crpix', data=wcspars['crpix'])
        grp.create_dataset('cdelt', data=wcspars['cdelt'])
        grp.create_dataset('pc', data=wcspars['pc'])
        
        grp = f.create_group('polpars')
        grp.attrs['nx'] = polpars['nx']
        grp.attrs['ny'] = polpars['ny']
        grp.attrs['order'] = polpars['order']
        grp.create_dataset('x_wcs2pix', data=polpars['x_wcs2pix'])
        grp.create_dataset('y_wcs2pix', data=polpars['y_wcs2pix'])
        grp.create_dataset('x_pix2wcs', data=polpars['x_pix2wcs'])
        grp.create_dataset('y_pix2wcs', data=polpars['y_pix2wcs'])  
        
    return      
    
def read_catalogue(starcat=cfg.starcat):
    
    catalogue = dict()

    data = fits.getdata(starcat) 
    
    catalogue['ascc'] = data['ASCC']
    catalogue['ra'] = data['_RAJ2000']
    catalogue['dec'] = data['_DEJ2000']
    catalogue['vmag'] = data['Vmag']
    catalogue['bmag'] = data['Bmag']
    catalogue['sptype'] = data['SpType']
    catalogue['blend'] = data['Blend']
    catalogue['var'] = data['Var']
    
    catalogue = np.rec.fromarrays([catalogue[key] for key in catalogue.keys()], names=[key for key in catalogue.keys()])    
    
    sort = np.argsort(catalogue['ascc'])
    catalogue = catalogue[sort]
    
    return catalogue    
    
def read_skyidx(ascc0, starcat=cfg.starcat):

    data = fits.getdata(starcat)

    ascc = data['ASCC']
    skyidx = data['SKYIDX']    
    
    sort = np.argsort(ascc)
    ascc = ascc[sort]
    skyidx = skyidx[sort]

    args = np.searchsorted(ascc, ascc0)

    return skyidx[args]    
    
def read_stack(filelist):

    nimages = len(filelist)

    stack = []
    headers = []

    for i in range(nimages):

        filename = filelist[i]

        try:
            image, header = fits.getdata(filename, header=True)
        except:
            log.warn('Failed to read image {}'.format(filename))
            continue
        else:
            stack.append(image)
            headers.append(header)
        
    stack = np.array(stack)

    return stack, headers 
    
def write_binimage(filename, image, header, overwrite=True):

    image = image.astype('float32')
    
    hdu = fits.PrimaryHDU(data=image, header=header)
    hdu.writeto(filename, clobber=overwrite)
    
    return
    
class TmpFile(object):
    
    def __init__(self, filename, overwrite=True):
        
        if overwrite:
            try:
                os.remove(filename)
            except:
                pass
        
        self.filename = filename
        
        return
        
    def add_header(self, header, siteinfo, aper, skyrad):

         with h5py.File(self.filename) as f:
            
            grp = f.create_group('header')
            grp.attrs['site-obs'] = header['SITE-OBS']
            grp.attrs['cam-obs'] = header['CAMERA']
            grp.attrs['lat'] = siteinfo['lat']
            grp.attrs['lon'] = siteinfo['lon']
            grp.attrs['alt'] = siteinfo['height']
            grp.attrs['aversion'] = header['APIVER']
            grp.attrs['cversion'] = header['SWVER']
            grp.attrs['rversion'] = '2017.00' # TODO
            grp.attrs['exptime'] = [6.4, 3.2] # TODO
            grp.attrs['ccdtemp'] = -15 # TODO
            grp.attrs['aper'] = aper
            grp.attrs['skyrad'] = skyrad
        
    def add_station(self, station):
        
        # Sort the array by lstseq.
        sort = np.argsort(station['lstseq'])
        station = station[sort]        
        
        with h5py.File(self.filename) as f:
            
            grp = f.create_group('station')
            
            for field in station.dtype.names:
                grp.create_dataset(field, data=station[field])
                
        return
    
    def add_lightcurves(self, curves, cat):

        # Sort the array be ascc first and lstseq second.        
        sort = np.lexsort((curves['lstseq'], curves['ascc']))
        curves = curves[sort]         
        
        # Get the unique stars and the number of points for each star.
        ascc, nobs = np.unique(curves['ascc'], return_counts=True)
               
        strides = np.append(0, np.cumsum(nobs))        
        names = list(curves.dtype.names)
        names.remove('ascc')
        
        # Create arrays.
        lstsqmin = np.zeros(len(ascc), dtype='uint32')
        lstsqmax = np.zeros(len(ascc), dtype='uint32')
        
        with h5py.File(self.filename) as f:
            
            grp = f.create_group('lightcurves')            
            
            for i in range(len(ascc)):
                
                lc = curves[strides[i]:strides[i+1]][names]
                
                lstsqmin[i] = lc['lstseq'][0]
                lstsqmax[i] = lc['lstseq'][-1]                
                
                grp.create_dataset(ascc[i], data=lc)
                
            grp = f.create_group('stars')
            grp.create_dataset('ascc', data=ascc)
            grp.create_dataset('nobs', data=nobs.astype('uint16'))
            grp.create_dataset('lstsqmin', data=lstsqmin)
            grp.create_dataset('lstsqmax', data=lstsqmax)
                
            idx = np.searchsorted(cat['ascc'], ascc)
            grp.create_dataset('ra', data=cat[idx]['ra'])
            grp.create_dataset('dec', data=cat[idx]['dec'])
            grp.create_dataset('vmag', data=cat[idx]['vmag'])
            grp.create_dataset('bmag', data=cat[idx]['bmag'])
            grp.create_dataset('sptype', data=cat[idx]['sptype'])
            grp.create_dataset('blend', data=cat[idx]['blend'])
            grp.create_dataset('var', data=cat[idx]['var'])
            
        return
        
    def add_astrometry(self, astrosol):
        
        if astrosol:        
        
            with h5py.File(self.filename) as f:
                
                grp = f.create_group('astrometry')
                
                for field in astrosol.keys():
                    grp.create_dataset(field, data=astrosol[field])
                
        return
        
def write_stamps(filename, stamps, curves, station, cat, overwrite=True):
    
    if overwrite:
        try:
            os.remove(filename)
        except:
            pass
     
    # Sort the array by ascc first and lstseq second.        
    sort = np.lexsort((curves['lstseq'], curves['ascc']))
    curves = curves[sort] 
    stamps = stamps[sort]        
    
    # Get the unique stars and the number of points for each star.
    ascc, nobs = np.unique(curves['ascc'], return_counts=True)
           
    strides = np.append(0, np.cumsum(nobs))        
    names = list(curves.dtype.names)
    names.remove('ascc')
    
    lstsqmin = np.zeros(len(ascc), dtype='uint32')
    lstsqmax = np.zeros(len(ascc), dtype='uint32')    
    
    with h5py.File(filename) as f:
           
        grp = f.create_group('stamps')
        for i in range(len(ascc)):
            
            grp.create_dataset(ascc[i]+'/lstseq', data=curves['lstseq'][strides[i]:strides[i+1]])
            grp.create_dataset(ascc[i]+'/stamps', data=stamps[strides[i]:strides[i+1]])
        
        grp = f.create_group('lightcurves')
        for i in range(len(ascc)):
            
            lc = curves[strides[i]:strides[i+1]][names]
            
            lstsqmin[i] = lc['lstseq'][0]
            lstsqmax[i] = lc['lstseq'][-1]                
            
            grp.create_dataset(ascc[i], data=lc)
        
        grp = f.create_group('stars')
        grp.create_dataset('ascc', data=ascc)
        grp.create_dataset('nobs', data=nobs.astype('uint16'))
        grp.create_dataset('lstsqmin', data=lstsqmin)
        grp.create_dataset('lstsqmax', data=lstsqmax)
            
        idx = np.searchsorted(cat['ascc'], ascc)
        grp.create_dataset('ra', data=cat[idx]['ra'])
        grp.create_dataset('dec', data=cat[idx]['dec'])
        grp.create_dataset('vmag', data=cat[idx]['vmag'])
        grp.create_dataset('bmag', data=cat[idx]['bmag'])
        grp.create_dataset('sptype', data=cat[idx]['sptype'])
        grp.create_dataset('blend', data=cat[idx]['blend'])
        grp.create_dataset('var', data=cat[idx]['var'])
            
        grp = f.create_group('station')
        for key in station.dtype.names:
            grp.create_dataset(key, data=station[key])
            
    return

###############################################################################
### Keep track of processed files so we can restart.
###############################################################################

def write_in_queue(filename, filelist):
    
    import json    
    
    filelist = list(filelist)    
    
    with open(filename, mode='w') as f:
        json.dump(filelist, f)
        
    return
    
def read_in_queue(filename):
    
    import json
    
    with open(filename, mode='r') as f:
        filelist = json.load(f)
    
    filelist = set(filelist)
    
    return filelist

###############################################################################
### Functions for combining the photometry files.
###############################################################################

def _index_files(filelist):
    
    nfiles = len(filelist)
    
    idx1 = np.array([], dtype='uint16')
    
    stars = dict()
    for i in range(nfiles):
        
        filename = filelist[i]        
        
        with h5py.File(filename, 'r') as f:
            
            grp = f['stars']
            
            for key in grp.keys():
                
                ascc_ = grp['ascc'].value
                
                if key not in stars.keys():
                    stars[key] = grp[key].value
                else:
                    stars[key] = np.append(stars[key], grp[key].value)

            grp = f['lightcurves']

            dtype = grp[ascc_[0]].dtype

        idx1 = np.append(idx1, np.repeat(i, len(ascc_)))
    
    ascc, args, idx2 = np.unique(stars['ascc'], return_index=True, return_inverse=True)
    nstars = len(ascc)    
    nobs = np.zeros((nfiles, nstars), dtype='uint32')
    stars['nobs'] = stars['nobs'].astype('uint32')
    nobs[idx1, idx2] = stars['nobs']
    
    for key in stars.keys():
        stars[key] = stars[key][args]
    
    return stars, nobs, dtype

def _read_header(filelist):
    
    header = dict()
    with h5py.File(filelist[0], 'r') as f:
        
        grp = f['header']
        for key in grp.attrs.keys():
            header[key] = grp.attrs[key]
            
    return header

def _read_curves(filelist, ascc, nobs, dtype):
    
    nfiles = len(filelist)
    nstars = len(ascc)
    
    strides = np.row_stack([nstars*[0], np.cumsum(nobs, axis=0)])    
    curves = {ascc[i]:np.recarray(strides[-1,i], dtype=dtype) for i in range(nstars)}
    
    for i in range(nfiles):
        
        filename = filelist[i]
        
        with h5py.File(filename, 'r') as f:
            
            grp = f['lightcurves']
            
            for j in range(nstars):
                
                if (nobs[i,j] > 0):
                    
                    curves[ascc[j]][strides[i,j]:strides[i+1,j]] = grp[ascc[j]].value
                    
    return curves

def _read_station(filelist):
    
    station = dict()
    for filename in filelist:
        
        with h5py.File(filename, 'r') as f:
            
            grp = f['station']
            
            for key in grp.keys():
                
                if key not in station.keys():
                    station[key] = grp[key].value
                else:
                    station[key] = np.append(station[key], grp[key].value)
                
    return station
    
def _read_astrometry(filelist):
    
    astrometry = dict()
    for filename in filelist:
        
        with h5py.File(filename, 'r') as f:
            
            try:
                grp = f['astrometry']
            except:
                continue
            else:
                
                for key in grp.keys():
                    
                    if key not in astrometry.keys():
                        astrometry[key] = [grp[key].value]
                    else:
                        astrometry[key].append(grp[key].value)
                    
    for key in astrometry.keys():
        astrometry[key] = np.array(astrometry[key])

    return astrometry

def combine_photometry(filename, filelist, overwrite=True, astrometry=True, nsteps=1000):
    
    if overwrite:
        try:
            os.remove(filename)
        except:
            pass    
    
    filelist = np.sort(filelist)    
    
    # Read the combined stars field and index the files.
    stars, nobs, dtype = _index_files(filelist)
    
    nstars = len(stars['ascc'])
    for i in range(0, nstars, nsteps):

        # Read the combined lightcurves for a group of stars.
        curves = _read_curves(filelist, stars['ascc'][i:i+nsteps], nobs[:,i:i+nsteps], dtype)
             
        # Write the combined lightcurves for a group of stars.
        with h5py.File(filename) as f:
            
            for j in range(i, i+len(stars['ascc'][i:i+nsteps])):
             
                tmp = curves[stars['ascc'][j]]
                
                stars['nobs'][j] = len(tmp)
                stars['lstsqmin'][j] = tmp['lstseq'][0]
                stars['lstsqmax'][j] = tmp['lstseq'][-1]
                
                f.create_dataset('lightcurves/{}'.format(stars['ascc'][j]), data=tmp)
             
    # Write the combined "stars" field.
    with h5py.File(filename) as f:
        
        grp = f.create_group('stars')
        for key in stars.keys():
            grp.create_dataset(key, data=stars[key])
        
    # Read the "header" field.
    header = _read_header(filelist)        
    
    # Write the "header" field..
    with h5py.File(filename) as f:
        
        grp = f.create_group('header')
        for key in header.keys():
            grp.attrs[key] = header[key]
        
    # Read the combined "station" field.
    station = _read_station(filelist)
    
    # Write the combined "station" field.
    with h5py.File(filename) as f:
        
        grp = f.create_group('station')
        for key in station.keys():
            grp.create_dataset(key, data=station[key])
    
    # Read the combined "astrometry" field.  
    if astrometry:
        astrometry = _read_astrometry(filelist)
         
    if astrometry:
        
        # Write the combined "astrometry" field.
        with h5py.File(filename) as f:
        
            grp = f.create_group('astrometry')
            for key in astrometry.keys():
                grp.create_dataset(key, data=astrometry[key])
            
    return

###############################################################################
### Function for reading the photometry.
###############################################################################

class PhotFile(object):
    
    def __init__(self, filename):
        
        self.filename = filename
        
        return
        
    def read_stars(self, fields=None):
        
        stars = dict()

        with h5py.File(self.filename, 'r') as f:
            
            grp = f['stars']
            
            if fields is None:
                fields = grp.keys()
            
            for field in fields:
                
                if field in grp.keys():
                    stars[field] = grp[field].value
                else:
                    print 'Warning: skipping field {}, field not found.'.format(field)
        
        return stars
        
    def read_station(self, fields=None, lstseq=None):
        
        station = dict()
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['station']
            
            lstseq_station = grp['lstseq'].value            
            
            if fields is None:
                fields = grp.keys()
            
            for field in fields:
                
                if field in grp.keys():
                    station[field] = grp[field].value
                else:
                    print 'Warning: skipping field {}, field not found.'.format(field)
        
        
        if lstseq is not None:
            
            idx = np.searchsorted(lstseq_station, lstseq)
            for key in station.keys():
                station[key] = station[key][idx]
                
        return station
            
    def read_lightcurves(self, ascc=None, fields=None, perstar=True):
        
        onestar = False        
        
        if ascc is None:
            stars = self.read_stars(['ascc'])
            ascc = stars['ascc']
            
        elif isinstance(ascc, basestring):
            onestar = True
            ascc = [ascc]
        
        nstars = len(ascc)        
        curves = dict()
        nobs = np.zeros(nstars, dtype='int')        
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['lightcurves']
            
            for i in range(nstars):
                
                if ascc[i] in grp.keys():
                    curves[ascc[i]] = grp[ascc[i]].value
                    nobs[i] = len(curves[ascc[i]])
                else:
                    print 'Warning: skipping star {}, star not found.'.format(ascc[i])
        
        if not curves:
            return curves
        
        # Select specified fields.
        if fields is not None:
            
            for i in range(nstars):
                curves[ascc[i]] = curves[ascc[i]][fields]
                    
        # Combine lightcurves.
        if not perstar:
            
            strides = np.append(0, np.cumsum(nobs))
            tmp = np.recarray(strides[-1], dtype=curves[ascc[0]].dtype)
            
            for i in range(nstars):
                tmp[strides[i]:strides[i+1]] = curves[ascc[i]]
                        
            curves = tmp                        
        
        if onestar:
            return curves[ascc[0]]             
            
        return curves

###############################################################################
### Function for reading the systematics.
###############################################################################

def update_systable(filename, systable, nmax=5):  
    
    data = np.loadtxt(systable, dtype='S')    
    data = np.atleast_1d(data)
    
    if (len(data) == nmax):
        
        data = np.roll(data, 1)
        data[0] = filename
        
    else:
        
        data = np.append(filename, data)
        
    np.savetxt(systable, data, fmt='%s')
    
    return
    
def get_sysfile(systable):

    data = np.loadtxt(systable, dtype='S')  
    
    try:
        sysfile = data[0]
    except:
        sysfile = data[None][0]
        
    return sysfile
    
def tmp_clouds(filename, nobs, clouds, sigma, lstmin, lstmax, lstlen, overwrite=True):
    
    if overwrite:
        try:
            os.remove(filename)
        except:
            pass    
        
    with h5py.File(filename) as f:
        
        f.create_dataset('nobs', data=nobs)
        f.create_dataset('clouds', data=clouds)
        f.create_dataset('sigma', data=sigma)
        f.attrs['grid'] = 'healpix'
        f.attrs['nx'] = 8
        f.attrs['lstmin'] = lstmin
        f.attrs['lstmax'] = lstmax
        f.attrs['lstlen'] = lstlen
        
    return

class SysFile(object):    
    
    def __init__(self, filename):
        
        self.filename = filename
        
        return
        
    def read_magnitudes(self):
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['data/magnitudes']
            
            ascc = grp['ascc'].value
            nobs = grp['nobs'].value
            mag = grp['mag'].value
            sigma = grp['sigma'].value
        
        return ascc, nobs, mag, sigma
        
    def read_transmission(self):
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['data/transmission']
            
            idx1 = grp['idx1'].value
            idx2 = grp['idx2'].value
            nobs = grp['nobs'].value
            trans = grp['trans'].value
            
            nx = grp.attrs['nx']
            ny = grp.attrs['ny']
            
        pg = grids.PolarGrid(nx, ny)
        
        nobs = pg.values2grid(idx1, idx2, nobs, np.nan)
        trans = pg.values2grid(idx1, idx2, trans, np.nan)
        
        return pg, nobs, trans
        
    def read_intrapix(self):
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['data/intrapix']
            
            idx1 = grp['idx1'].value
            idx2 = grp['idx2'].value
            nobs = grp['nobs'].value
            amplitudes = grp['amplitudes'].value
            
            nx = grp.attrs['nx']
            ny = grp.attrs['ny']
            
        pg = grids.PolarGrid(nx, ny)
        
        nobs = pg.values2grid(idx1, idx2, nobs, np.nan) 
        amplitudes = np.stack([pg.values2grid(idx1, idx2, amplitudes[:,i], np.nan) for i in range(4)], axis=-1)

        return pg, nobs, amplitudes
        
    def read_clouds(self):
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['data/clouds']
            
            idx = grp['idx'].value
            lstseq = grp['lstseq'].value
            nobs = grp['nobs'].value
            clouds = grp['clouds'].value
            sigma = grp['sigma'].value
            
            nx = grp.attrs['nx']
            lstmin = grp.attrs['lstmin']
            lstmax = grp.attrs['lstmax']
            lstlen = grp.attrs['lstlen']
    
        lstseq = lstseq - lstmin
    
        try:
            hg = grids.HealpixGrid(nx)
        except:
            hg = None
        
        tmp = np.full((768, lstlen), fill_value = np.nan)
        tmp[idx, lstseq] = nobs
        nobs = tmp
        
        tmp = np.full((768, lstlen), fill_value = np.nan)
        tmp[idx, lstseq] = clouds
        clouds = tmp
            
        tmp = np.full((768, lstlen), fill_value = np.nan)
        tmp[idx, lstseq] = sigma
        sigma = tmp        
        
        return hg, nobs, clouds, sigma, lstmin, lstmax
