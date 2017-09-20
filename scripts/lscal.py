# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 14:58:39 2017

@author: talens
"""

import os
import glob
import datetime

import numpy as np

from lsreduce import io
from lsreduce.cdecor_vmag import CoarseDecorVmag
from lsreduce.apply_calibration import apply_calibration

import logging.config

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('lsreduce')

def baseline(date, part, camera, source, dest):
    
    # Check that the date is valid.
    try:
        datetime.datetime.strptime(date, '%Y%m')
    except ValueError:
        raise ValueError("Incorrect date format, should be YYYYMM")
     
    # Get the photometry files matching the date, part and camera.
    if (part == 'A'): 
     
        filelist = glob.glob(os.path.join(source, '*/lightcurves/fast_{}0?{}.hdf5'.format(date, camera)))
        filelist = filelist + glob.glob(os.path.join(source, '*/lightcurves/fast_{}1[0-5]{}.hdf5'.format(date, camera)))
    
    elif (part == 'B'):
    
        filelist = glob.glob(os.path.join(source, '*/lightcurves/fast_{}1[6-9]{}.hdf5'.format(date, camera)))
        filelist = filelist + glob.glob(os.path.join(source, '*/lightcurves/fast_{}[23]?{}.hdf5'.format(date, camera)))
        
    elif (part == 'C'):
    
        filelist = glob.glob(os.path.join(source, '*/lightcurves/fast_{}*{}.hdf5'.format(date, camera)))
       
    else:
        raise ValueError("part should be either 'A', 'B' or 'C'.")
    
    # Check that there are valid files.
    if len(filelist) == 0:
        log.info('No valid data found.')
        return None

    # Sort the filelist.
    filelist = np.sort(filelist)    
    
    log.info('Combining files:')
    for filename in filelist:
        log.info('    {}'.format(filename)) 

    # Create the destination directory.    
    outpath = os.path.join(dest, camera)    
    io.ensure_dir(outpath)
    
    # Check that the destination file does not exist.
    photfile = os.path.join(outpath, 'fast_{}{}{}.hdf5'.format(date, part, camera)) 
    if os.path.isfile(photfile):
        raise IOError('Output file already exists: {}'.format(photfile))
    else:
        log.info('Writing results to: {}'.format(photfile))
    
    # Combine the files.
    io.combine_photometry(photfile, filelist, astrometry=False, overwrite=False)
    
    return photfile

def lscal(date, part, cameras, aper, source, dest):
    
    for cam in cameras:
        
        photfile = baseline(date, part, cam, source, dest)

        if photfile is not None:        
        
            sys = CoarseDecorVmag(photfile, aper)
            
            apply_calibration(photfile, aper)
    
    return
    
if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Perform the coarse decorrelation on a baseline.')
    parser.add_argument('date', type=str,
                        help='a date of the form YYYYMM')
    parser.add_argument('part', type=str, choices=['A','B','C'],
                        help='letter indicating which baseline to create') 
    parser.add_argument('dest', type=str,
                        help='Location where the products will be written, e.g. /data3/talens/2017Q4. If the path does not exist it will be created. Subdirectories are generated automatically.')
    parser.add_argument('-c', '--cam', type=str, nargs='+', default=['LSN', 'LSE', 'LSS', 'LSW', 'LSC'],
                        help ='the camera(s) to perform the combination for', dest='cameras')                        
    parser.add_argument('-a', '--aper', type=int, choices=[0,1], default=0,
                        help ='the aperture to perform the coarse decorrelation on', dest='aper')
    parser.add_argument('-d', '--data', type=str, default='/data5/mascara/LaSilla',
                        help='Location of the raw data.', dest='source')
    args = parser.parse_args()

    lscal(args.date, args.part, args.cameras, args.aper, args.source, args.dest)    
    