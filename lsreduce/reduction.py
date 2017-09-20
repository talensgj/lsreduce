# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 14:10:15 2016

@author: talens
"""

import os
import logging

import numpy as np

from astropy.io import fits

from . import io
from . import stacking
from . import astrometry
from . import photometry
from . import configuration as cfg
from . import misc
from . import sigmas

log = logging.getLogger('lsreduce')

###############################################################################
### Helper functions.
###############################################################################

def headers2recarray(headers, fields):
    
    nimages = len(headers)

    dtype = [(key, fields[key]) for key in fields.keys()]
    array = np.recarray(nimages, dtype=dtype)    
    
    for i in range(nimages):
        
        header = headers[i]
        
        for key in fields.keys():
            array[i][key] = header[key]
    
    return array
   
def expand_recarray(array, fields):
    
    dtype = array.dtype.descr + [(key, fields[key]) for key in fields.keys()]
    new_array = np.recarray(array.shape, dtype=dtype)
    
    for key in array.dtype.names:
        new_array[key] = array[key] 
    
    return new_array

###############################################################################
### Functions for reducing dark frames.
############################################################################### 

def combine_bias(stack, headers):
    """ Create a masterbias. """
    
    nbias = len(stack)
    
    # Create the master bias.
    masterbias = np.nanmean(stack, axis=0)
    
    # Create the header of the masterdark.
    header = fits.Header()    
    header['NOBS'] = (nbias, 'Number of images used')
    header['METHOD'] = ('mean', 'Method used to combine the images')   
    header['SITE-OBS'] = (headers[0]['SITE-OBS'], 'Observation site')
    header['CAMERA'] = (headers[0]['CAMERA'], 'Camera name on site')
    header['IMAGETYP'] = ('MASTERBIAS', 'Type of image')
    header['EXPTIME'] = (np.mean(headers['EXPTIME']), 'Exposure time in seconds')
    header['CCDTEMP'] = (np.mean(headers['CCDTEMP']), 'Average CCD temperature (C)')
    header['LSTSEQ'] = (headers[0]['LSTSEQ'], 'Exposure num of first frame used')  
    header['X0'] = headers[0]['X0']
    header['XSIZE'] = headers[0]['XSIZE']
    header['Y0'] = headers[0]['Y0']
    header['YSIZE'] = headers[0]['YSIZE']
    
    return masterbias, header

def combine_darks(stack, headers):
    """ Create a masterdark. """    
    
    ndark = len(stack)
        
    # Create the master dark.
    masterdark = np.nanmean(stack, axis=0)
    
    # Create the header of the masterdark.
    header = fits.Header()    
    header['NOBS'] = (ndark, 'Number of images used')
    header['METHOD'] = ('mean', 'Method used to combine the images')   
    header['SITE-OBS'] = (headers[0]['SITE-OBS'], 'Observation site')
    header['CAMERA'] = (headers[0]['CAMERA'], 'Camera name on site')
    header['IMAGETYP'] = ('MASTERDARK', 'Type of image')
    header['EXPTIME'] = (np.mean(headers['EXPTIME']), 'Exposure time in seconds')
    header['CCDTEMP'] = (np.mean(headers['CCDTEMP']), 'Average CCD temperature (C)')
    header['LSTSEQ'] = (headers[0]['LSTSEQ'], 'Exposure num of first frame used')  
    header['X0'] = headers[0]['X0']
    header['XSIZE'] = headers[0]['XSIZE']
    header['Y0'] = headers[0]['Y0']
    header['YSIZE'] = headers[0]['YSIZE']
    
    return masterdark, header
    
def combine_flats(stack, headers, dark):
    """ Create a masterflat. """
    
    nflat = len(stack)
    
    # Create the master flat.
    stack = stack - dark
    stack = stack/np.nanmean(stack, axis=(1,2), keepdims=True)
    masterflat = np.nanmean(stack, axis=0)
    
    # Create the header of the masterflat.
    header = fits.Header()    
    header['NOBS'] = (nflat, 'Number of images used')
    header['METHOD'] = ('mean', 'Method used to combine the images')   
    header['SITE-OBS'] = (headers[0]['SITE-OBS'], 'Observation site')
    header['CAMERA'] = (headers[0]['CAMERA'], 'Camera name on site')
    header['IMAGETYP'] = ('MASTERFLAT', 'Type of image')
    header['EXPTIME'] = (np.mean(headers['EXPTIME']), 'Exposure time in seconds')
    header['CCDTEMP'] = (np.mean(headers['CCDTEMP']), 'Average CCD temperature (C)')
    header['LSTSEQ'] = (headers[0]['LSTSEQ'], 'Exposure num of first frame used')  
    header['X0'] = headers[0]['X0']
    header['XSIZE'] = headers[0]['XSIZE']
    header['Y0'] = headers[0]['Y0']
    header['YSIZE'] = headers[0]['YSIZE']
    
    return masterflat, header
    
def reduce_bias_frames(camid, filelist, dirtree, nmin=cfg.minbias):
    """ Create a masterbias. """    
    
    log.info('Received {} bias frames.'.format(len(filelist)))
    
    # Read the files.
    stack, headers = io.read_stack(filelist)

    if (len(stack) == 0):
        log.warn('No images were succesfully read.')
        return

    # Create a recarray containing the required header fields.
    fields = {'LSTSEQ':'uint32', 'JD':'float64', 'LST':'float64',
              'EXPTIME':'float32', 'CCDTEMP':'float32',
              'SITE-OBS':'|S20', 'CAMERA':'|S1',
              'X0':'uint16', 'XSIZE':'uint16',
              'Y0':'uint16', 'YSIZE':'uint16'}
    headers = headers2recarray(headers, fields)
   
    nbias = len(stack)
    if (nbias >= nmin):   
        
        log.info('Creating masterbias from {} exposures.'.format(nbias))
    
        masterbias, header = combine_bias(stack, headers)    
    
        filename = '{:08d}{}masterbias.fits'.format(header['LSTSEQ'], camid)
        filename = os.path.join(dirtree['binned'], filename)
    
        log.info('Saving masterbias to {}'.format(filename))    
    
        io.write_masterframe(filename, masterbias, header)

    else:
        log.info('Not enough valid bias frames.')
        
    return    
    
def reduce_dark_frames(camid, filelist, dirtree, darktable, nmin=cfg.mindark):
    """ Create a masterdark. """    
    
    log.info('Received {} dark frames.'.format(len(filelist)))
    
    # Read the files.
    stack, headers = io.read_stack(filelist)

    if (len(stack) == 0):
        log.warn('No images were succesfully read.')
        return

    # Create a recarray containing the required header fields.
    fields = {'LSTSEQ':'uint32', 'JD':'float64', 'LST':'float64',
              'EXPTIME':'float32', 'CCDTEMP':'float32',
              'SITE-OBS':'|S20', 'CAMERA':'|S1',
              'X0':'uint16', 'XSIZE':'uint16',
              'Y0':'uint16', 'YSIZE':'uint16'}
    headers = headers2recarray(headers, fields)
   
    ndarks = len(stack)
    if (ndarks >= nmin):   
        
        log.info('Creating masterdark from {} exposures.'.format(ndarks))
    
        masterdark, header = combine_darks(stack, headers)    
    
        filename = '{:08d}{}masterdark.fits'.format(header['LSTSEQ'], camid)
        filename = os.path.join(dirtree['binned'], filename)
    
        log.info('Saving masterdark to {}'.format(filename))    
    
        io.write_masterframe(filename, masterdark, header)
        io.update_darktable(filename, darktable)
        #log.warn('Running with standard dark, not updating darktable.')
    else:
        log.info('Not enough valid darks.')
        
    return
    
def reduce_flat_frames(camid, filelist, dirtree, darktable, nmin=cfg.minflat):
    """ Create a masterflat. """    
    
    log.info('Received {} flat frames.'.format(len(filelist)))
    
    # Read the files.
    stack, headers = io.read_stack(filelist)
    dark, header = io.read_masterdark(darktable)

    if (len(stack) == 0):
        log.warn('No images were succesfully read.')
        return

    # Create a recarray containing the required header fields.
    fields = {'LSTSEQ':'uint32', 'JD':'float64', 'LST':'float64',
              'EXPTIME':'float32', 'CCDTEMP':'float32',
              'SITE-OBS':'|S20', 'CAMERA':'|S1',
              'X0':'uint16', 'XSIZE':'uint16',
              'Y0':'uint16', 'YSIZE':'uint16'}
    headers = headers2recarray(headers, fields)
   
    nflats = len(stack)
    if (nflats >= nmin):   
        
        log.info('Creating masterflat from {} exposures.'.format(nflats))
    
        masterflat, header = combine_flats(stack, headers, dark)    
    
        filename = '{:08d}{}masterflat.fits'.format(header['LSTSEQ'], camid)
        filename = os.path.join(dirtree['binned'], filename)
    
        log.info('Saving masterflat to {}'.format(filename))    
    
        io.write_masterframe(filename, masterflat, header)

    else:
        log.info('Not enough valid flats.')
        
    return    
    
###############################################################################
### Functions for reducing science frames.
###############################################################################    
   
def preprocess(stack, headers, darktable):
    
    import datetime    
    
    nimages = len(stack)
    
    # Create basic station field from the headers.
    fields = {'lstseq':'uint32', 'jd':'float64', 'lst':'float64',
              'exptime':'float32', 'ccdtemp':'float32',
              'dwldtime':'float32', 'margtime':'float32'}
    station = headers2recarray(headers, fields)    

    # Add fields for derived information.
    fields = {'imnzero':'uint32', 'imnsat':'uint32', 'imstd':'float32', 'immed':'uint16',
              'oversca0':'uint16', 'oversca1':'uint16', 'oversca2':'uint16', 'oversca3':'uint16',
              'lstday':'uint32', 'lstidx':'uint32', 'darkfile':'uint32',
              'year':'uint16', 'month':'uint8', 'day':'uint8', 'hour':'uint8', 'min':'uint8', 'sec':'uint8', 'fsec':'float32'}
    station = expand_recarray(station, fields)     
    
    # Add derivatives of the header fields.
    station['lstday'] = station['lstseq']//13500
    station['lstidx'] = station['lstseq']%13500
    
    for i in range(nimages):
        
        tmp = datetime.datetime.strptime(headers[i]['UTC-OBS'], '%Y-%m-%d %H:%M:%S.%f')
    
        station[i]['year'] = tmp.year
        station[i]['month'] = tmp.month
        station[i]['day'] = tmp.day
        station[i]['hour'] = tmp.hour
        station[i]['min'] = tmp.minute
        station[i]['sec'] = tmp.second
        station[i]['fsec'] = tmp.microsecond      
    
    # Boundaries of the overscan region.
    header = headers[0]
    lx = header['X0']
    ux = lx + header['XSIZE']
    ly = header['Y0']
    uy = ly + header['YSIZE']     
    
    # Compute statistics on the raw images.
    station['imnzero'] = np.sum((stack == 0), axis=(1,2))
    station['imnsat'] = np.sum((stack > 64000), axis=(1,2))
    station['imstd'] = np.nanstd(stack, axis=(1,2))
    station['immed'] = np.nanmedian(stack, axis=(1,2))
    
    station['oversca0'] = np.nanmedian(stack[:,uy:,lx:ux], axis=(1,2)) # Top
    station['oversca1'] = np.nanmedian(stack[:,ly:uy,:lx], axis=(1,2)) # Left
    station['oversca2'] = np.nanmedian(stack[:,:ly,lx:ux], axis=(1,2)) # Bottom
    station['oversca3'] = np.nanmedian(stack[:,ly:uy,ux:], axis=(1,2)) # Right         
        
    stack = stack.astype('float32')        
        
    # Subtract the masterdark.
    masterdark, darkheader = io.read_masterdark(darktable)
    station['darkfile'] = darkheader['LSTSEQ'] 
    stack = stack - masterdark
    
    # Remove the overscan region.
    stack = stack[:,ly:uy,lx:ux]
        
    return stack, station  
   
def sunmoon2station(station, siteinfo, astro):
    
    import ephem    
    
    nimages = len(station)    
    
    fields = {'sunra':'float32', 'sundec':'float32', 'sunalt':'float32',
              'moonra':'float32', 'moondec':'float32', 'moonalt':'float32',
              'moonx':'uint16', 'moony':'uint16',
              'moonph':'float32', 'moonmag':'float32'}
    dtype = [(key, fields[key]) for key in fields.keys()]
    dtype = station.dtype.descr + dtype
    nstation = np.recarray(nimages, dtype=dtype)
    
    # Copy the existing fields.
    for key in station.dtype.names:
        nstation[key] = station[key]
    station = nstation    
    
    # Compute the sky position of the sun.
    ra, dec, alt = astrometry.sun_position(siteinfo, station['jd'])
    
    station['sunra'] = ra
    station['sundec'] = dec
    station['sunalt'] = alt
    
    # Compute the sky position of the moon.
    ra, dec, alt = astrometry.moon_position(siteinfo, station['jd'])    
    
    station['moonra'] = ra
    station['moondec'] = dec
    station['moonalt'] = alt

    # Compute the CCD position of the moon.
    moon = ephem.Moon()

    for i in range(nimages):
        
        x, y, mask = astro.world2pix(station[i]['lst'], ra[i], dec[i])
       
        if mask:
            station[i]['moonx'] = np.around(x)
            station[i]['moony'] = np.around(y)
        else:
            station[i]['moonx'] = 9999
            station[i]['moony'] = 9999
            
        moon.compute(station[i]['jd'] - 2415020) # ephem uses Dublin JD  
        station[i]['moonph'] = moon.moon_phase
        station[i]['moonmag'] = moon.mag
    
    return station
    
def lightcurves(stack, station, astro, cat, aper, skyrad, maglim):
    """ Perform aperture photometry on the images. """    
    
    nimages = len(stack)    
    naper = len(aper)    
    
    # Initialize the photometry.
    phot = photometry.Photometry(aper, skyrad)    
    
    # Select stars brighter than the magnitude limit.
    select = (cat['vmag'] <= maglim)
    ascc = cat['ascc'][select]
    ra, dec = cat['ra'][select], cat['dec'][select]    
    
    # Create a recarray to hold the results.
    names = ['ascc', 'lstseq', 'sky', 'esky', 'peak', 'x', 'y', 'pflag', 'aflag']
    names = names + ['flux{}'.format(i) for i in range(naper)]
    names = names + ['eflux{}'.format(i) for i in range(naper)]    
    
    formats = ['|S32', 'uint32', 'float64', 'float32', 'float32', 'float32', 'float32', 'uint8', 'uint8']
    formats = formats + naper*['float64'] + naper*['float32']    
    
    curves = np.recarray((0,), names=names, formats=formats)    
    
    for i in range(nimages):    

        image = stack[i]
        lst = station[i]['lst']
        lstseq = station[i]['lstseq']
            
        # Compute positions and photometry.    
        x, y, mask = astro.world2pix(lst, ra, dec)  
        flux, eflux, sky, esky, peak, pflag = phot.get_phot(image, x, y)
    
        nstars = len(x)    
        
        # Put the results in the recarray.
        curves_ = np.recarray(nstars, names=names, formats=formats)
        curves_['ascc'] = ascc[mask]
        curves_['lstseq'] = np.repeat(lstseq, nstars)
        curves_['sky'] = sky
        curves_['esky'] = esky
        curves_['peak'] = peak
        curves_['x'] = x
        curves_['y'] = y
        curves_['pflag'] = pflag
        curves_['aflag'] = 0 # Will be set later.
        for i in range(naper):
            curves_['flux{}'.format(i)] = flux[:,i]
            curves_['eflux{}'.format(i)] = eflux[:,i]
            
        curves = np.append(curves, curves_)

    # Sort the array be ascc first and lstseq second.        
    sort = np.lexsort((curves['lstseq'], curves['ascc']))
    curves = curves[sort] 

    return curves

def calibrate_clouds(lstidx, skyidx, res, errsq, lstlen):
      
    nobs = np.zeros((lstlen, 768), dtype='uint32') 
    clouds = np.full((lstlen, 768), fill_value=np.nan, dtype='float32')
    sigma = np.full((lstlen, 768), fill_value=np.nan, dtype='float32')      
      
    if (len(lstidx) == 0):
        return nobs, clouds, sigma
      
    newidx = np.ravel_multi_index((lstidx, skyidx), (lstlen, 768))
    newidx, idx = np.unique(newidx, return_inverse=True) 
    lstidx, skyidx = np.unravel_index(newidx, (lstlen, 768)) 
    
    # Compute the clouds calibration.
    nobs[lstidx,skyidx] = np.bincount(idx)
    clouds[lstidx,skyidx], sigma[lstidx,skyidx] = sigmas.find_par_sigma(idx, res, errsq)      
    
    return nobs, clouds, sigma

def live_calibration(station, curves, cat, systable): 
    
    # Create empty arrays.
    sys = np.zeros(len(curves), dtype='float32')
    cflag = np.zeros(len(curves), dtype='uint8')
    
    # Perform look-up in station field.
    args = np.searchsorted(station['lstseq'], curves['lstseq'])    
    lst = station['lst'][args]   
    exptime = station['exptime'][args]    
    
    # Perform look-up in catalogue.
    args = np.searchsorted(cat['ascc'], curves['ascc'])
    ra = cat['ra'][args]
    dec = cat['dec'][args]    
    vmag = cat['vmag'][args]        
    
    # Convert fluxes to magnitudes.
    mag, emag = misc.flux2mag(curves['flux0']/exptime, curves['eflux0']/exptime)
    
    # Convert RA to HA.
    ha = astrometry.ra2ha(ra, lst)
    
    # Open the latest systematics file.
    sysfile = io.get_sysfile(systable)
    f = io.SysFile(sysfile)  
    
    # Compute the magnitude calibration.
    ascc, nobs, _, sigma = f.read_magnitudes() 
    
    tmp = np.zeros(len(curves))
    for i in range(len(curves)):
        
        if curves['ascc'][i] in ascc:
            tmp[i] = sigma[ascc == curves['ascc'][i]]

    sigma = tmp
    
    # Compute the transmission calibration.
    pg, nobs, trans = f.read_transmission()
    idx1, idx2 = pg.radec2idx(ha, dec)
    sys = sys + trans[idx1,idx2]    
    cflag = np.where(nobs[idx1,idx2] < 25, cflag+2, cflag)    
    
    # Compute the intrapixel calibration. 
    pg, nobs, amp = f.read_intrapix()
    idx1, idx2 = pg.radec2idx(ha, dec)
    ipx_x = amp[idx1,idx2,0]*np.sin(2*np.pi*curves['x']) + amp[idx1,idx2,1]*np.cos(2*np.pi*curves['x'])
    ipx_y = amp[idx1,idx2,2]*np.sin(2*np.pi*curves['y']) + amp[idx1,idx2,3]*np.cos(2*np.pi*curves['y'])
    sys = sys + ipx_x + ipx_y
    cflag = np.where(nobs[idx1,idx2] < 25, cflag+4, cflag)     
    
    # Get indices for computing the cloud calibration.
    lstmin = np.amin(station['lstseq'])
    lstmax = np.amax(station['lstseq'])
    lstlen = lstmax - lstmin + 1
    lstidx = curves['lstseq'] - lstmin   
    skyidx = io.read_skyidx(curves['ascc'])  
    
    # Compute the cloud calibrations.
    mask = (curves['pflag'] == 0) & (curves['aflag'] == 0) & np.isfinite(sys)
    nobs, clouds, sigma = calibrate_clouds(lstidx[mask], skyidx[mask], (mag - vmag - sys)[mask], (emag**2 + sigma**2)[mask], lstlen)
        
    sys = sys + clouds[lstidx,skyidx]
    cflag = np.where(nobs[lstidx,skyidx] < 25, cflag+8, cflag) 
    cflag = np.where(sigma[lstidx,skyidx] > .05, cflag+16, cflag) 
    cflag = np.where(np.isnan(sys), cflag+1, cflag)    
    
    # Add the temporary calibration to the lightcurves.
    fields = {'tmpmag0':'float32', 'tmpemag0':'float32', 'cflag':'uint8'}
    curves = expand_recarray(curves, fields)     
    
    curves['tmpmag0'] = mag - sys
    curves['tmpemag0'] = emag
    curves['cflag'] = cflag    
    
    return curves, nobs, clouds, sigma, lstmin, lstmax, lstlen

def binned_image(stack, station, header, astro):
    """ Create a binned image. """
    
    nimages = len(stack)

    # Choose a reference image.    
    mid_idx = nimages//2 # Central image for uneven number of images.
    
    # Get the WCS parameters and initialize the stacker.
    wcspars = astro.get_wcspars()
    stacker = stacking.Stacker(wcspars)    

    # Stack the images together.
    for i in range(nimages):
        stacker.add_image(station[i]['lst'], stack[i], station[i]['exptime'], station[mid_idx]['lst'])

    # Get the resulting image and a header containing the WCS transformations.
    binimage, binheader = stacker.get_image(station[mid_idx]['lst'])
    
    # Add additional information to the header.
    binheader['NIMAGES'] = (nimages, 'Number of images combined.')
    binheader['SITE-OBS'] = (header['site-obs'], 'Observation site')
    binheader['CAMERA'] = (header['camera'], 'Camera name on site')
    binheader['IMAGETYP'] = ('BINIMAGE', 'Type of image')
    binheader['EXPTIME'] = (np.mean(station['exptime']), 'Exposure time in seconds')    
    binheader['CCDTEMP'] = (np.mean(station['ccdtemp']), 'Average CCD temperature (C)')
    binheader['LSTSEQ'] = (station[mid_idx]['lstseq'], 'Exposure num of the reference frame')    
    binheader['LST'] = (station[mid_idx]['lst'], 'Sidereal time of the reference frame')
    binheader['JD'] = (station[mid_idx]['jd'], 'JD of the reference frame')
    binheader['LSTMID'] = (np.mean(station['lst']), 'Average LST in hours')
    binheader['JDMID'] = (np.mean(station['jd']), 'Average JD')    
    
    return binimage, binheader

def postage_stamps(curves, stack, station, nhalf=cfg.nhalf):
    
    npoints = len(curves)  
    
    idx = np.searchsorted(station['lstseq'], curves['lstseq'])
    stamps = np.zeros((npoints, 2*nhalf+1, 2*nhalf+1), dtype='float32')    
    
    for i in range(npoints):      
        
        x, y = curves['x'][i], curves['y'][i]
        xi, yi = np.around(x), np.around(y) 
        
        lx = int(xi - nhalf)
        ux = int(xi + nhalf + 1)
        ly = int(yi - nhalf)
        uy = int(yi + nhalf + 1)
        
        stamps[i] = stack[idx[i],ly:uy,lx:ux]
    
    return stamps

def reduce_science_frames(camid, filelist, siteinfo, dirtree, darktable, astromaster, systable):
    """ Process the raw data to produce lightcurves and binned images. """    
    
    log.info('Received {} science frames.'.format(len(filelist)))  
    
    # Read the catalogue.
    cat = io.read_catalogue()     
    
    # Initialize the astrometry.
    wcspars, polpars = io.read_astromaster(astromaster)
    astro = astrometry.Astrometry(wcspars, polpars)         
    
    # Read the files.
    stack, headers = io.read_stack(filelist)

    if (len(stack) == 0):
        log.warn('No images were succesfully read.')
        return

    # Preprocess the files.
    stack, station = preprocess(stack, headers, darktable) 

    # Perform astrometry.
    select = (cat['vmag'] <= cfg.maglim_astro) & (cat['vmag'] > 4.)           
    ra, dec = cat['ra'][select], cat['dec'][select]   
         
    aflag, astrosol = astro.solve(stack[0], station[0], ra, dec) 
    
    # Add sun and moon information.
    station = sunmoon2station(station, siteinfo, astro)    
    
    # Perform the fast photometry.
    log.info('Performing fast photometry.')
    fast_curves = lightcurves(stack, station, astro, cat, cfg.aper_fast, cfg.skyrad, cfg.maglim_fast)
    fast_curves['aflag'] = aflag

    # Perform the live calibration.
    fast_curves, nobs, clouds, sigma, lstmin, lstmax, lstlen = live_calibration(station, fast_curves, cat, systable)

    # Save the live calibration.
    filename = 'tmp_clouds{:08d}.hdf5'.format(station[0]['lstseq'])
    filename = os.path.join(dirtree['sys'], filename)
    
    log.info('Saving the live calibration to {}'.format(filename))
    
    io.tmp_clouds(filename, nobs, clouds, sigma, lstmin, lstmax, lstlen)

    # Save the fast lightcurves.
    filename = 'tmp_fast{:08d}.hdf5'.format(station[0]['lstseq'])        
    filename = os.path.join(dirtree['tmp'], filename)
    
    log.info('Saving fast photometry to {}'.format(filename))    
    
    f = io.TmpFile(filename) 
    f.add_header(headers[0], siteinfo, cfg.aper_fast, cfg.skyrad)      
    f.add_station(station)    
    f.add_lightcurves(fast_curves, cat)  
    f.add_astrometry(astrosol)
    
    # Create postage stamps.
    targets = io.read_targets(cfg.targets)
    for key in targets.keys():
        
        log.info('Trying to make postage stamps for target group {}'.format(key))     
        
        mask = np.in1d(fast_curves['ascc'], targets[key])
        
        if np.any(mask):
            
            stamp_curves = fast_curves[mask]
            stamps = postage_stamps(stamp_curves, stack, station)            
            
            filename = '{}_{:08d}{}.hdf5'.format(key, station[0]['lstseq'], camid)
            filename = os.path.join(dirtree['targets'], filename)
            
            log.info('Saving postage stamps for target group {} to {}'.format(key, filename))             
            
            io.write_stamps(filename, stamps, stamp_curves, station, cat)
            
        else:
            log.info('No stars present for target group {}'.format(key))
        
    # Create binned images.
    binstack = [] 
    binheaders = []

    log.info('Creating binned image from {} exposures.'.format(len(stack)))        
    
    binimage, binheader = binned_image(stack, station, headers[0], astro)
        
    filename = 'bin_{:08d}{}.fits.gz'.format(binheader['lstseq'], camid)
    filename = os.path.join(dirtree['binned'], filename)
    
    log.info('Saving binned image to {}'.format(filename))        
    
    io.write_binimage(filename, binimage, binheader) 
    
    binheaders = [binheader]
    binstack = np.array([binimage])
    
    fields = {'lstseq':'uint32', 'nimages':'uint8', 'exptime':'float32',
              'lst':'float64', 'jd':'float64', 'lstmid':'float64', 'jdmid':'float64'}
    binstation = headers2recarray(binheaders, fields)
    binstation = sunmoon2station(binstation, siteinfo, astro)

    # Perform the slow photometry.
    log.info('Performing slow photometry.')
    slow_curves = lightcurves(binstack, binstation, astro, cat, cfg.aper_slow, cfg.skyrad, cfg.maglim_slow)    
    slow_curves['aflag'] = aflag    
    
    # Save the slow lightcurves.
    filename = 'tmp_slow{:08d}.hdf5'.format(binstation[0]['lstseq'])        
    filename = os.path.join(dirtree['tmp'], filename)
    
    log.info('Saving slow photometry to {}'.format(filename))        
    
    f = io.TmpFile(filename) 
    f.add_header(headers[0], siteinfo, cfg.aper_slow, cfg.skyrad)       
    f.add_station(binstation)    
    f.add_lightcurves(slow_curves, cat)       
         
    return

def main():
    
    return
    
if __name__ == '__main__':
    
    main()
