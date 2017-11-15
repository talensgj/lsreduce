#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import time
import logging

import numpy as np

from . import io
from . import reduction    
from . import astrometry
from . import configuration as cfg
from . import summarize

###############################################################################
### Helper functions.
###############################################################################

def get_datestring():

    import datetime
    
    t = datetime.datetime.utcnow()

    return t.strftime('%Y%m%d')

def listdir_fullpath(d):
    
    return [os.path.join(d, f) for f in os.listdir(d)]

###############################################################################
### Placeholders.
###############################################################################

def combine_temporary_files(date, camid, dirtree):
    
    log = logging.getLogger('lsreduce')
    
    # Combine the fast lightcurve files.
    filelist = glob.glob(os.path.join(dirtree['tmp'], 'tmp_fast*.hdf5'))

    if (len(filelist) > 0):

        filename = 'fast_{}{}.hdf5'.format(date, camid)
        filename = os.path.join(dirtree['lightcurves'], filename)

        log.info('Combining {} temporary fast files, and writig results to {}'.format(len(filelist), filename))

        io.combine_photometry(filename, filelist)

        for filename in filelist:
            os.remove(filename)
    
    # Combine the slow lightcurve files.
    filelist = glob.glob(os.path.join(dirtree['tmp'], 'tmp_slow*.hdf5'))

    if (len(filelist) > 0):

        filename = 'slow_{}{}.hdf5'.format(date, camid)
        filename = os.path.join(dirtree['lightcurves'], filename)

        log.info('Combining {} temporary slow files, and writig results to {}'.format(len(filelist), filename))

        io.combine_photometry(filename, filelist)

        for filename in filelist:
            os.remove(filename)

    return    

def calibrate_photometry(date, camid, dirtree, systable):
    
    import datetime
    from . import cdecor_vmag    
    
    log = logging.getLogger('lsreduce')
    
    # Get all fast lightcurves in products directory.
    log.info('Getting list of fast lightcurve files.')
    
    filelist = glob.glob(os.path.join(cfg.reddir, r'*\lightcurves\fast_*.hdf5'))
    filelist = np.sort(filelist)

    if (len(filelist) == 0):
        log.warn('No lightcurves found.')
        return

    # Find fast lightcurves from past 15 nights.
    log.info('Obtaining fast lightcurves from the past 15 nights.')
    
    date0 = datetime.datetime.strptime(date, '%Y%m%d')
    date0 = date0 - datetime.timedelta(days=15)
    date0 = date0.strftime('%Y%m%d')
    dates = np.array([os.path.split(filename)[-1].rsplit('.')[0].rsplit('_')[-1] for filename in filelist])
    
    mask = (dates > (date0 + camid)) & (dates > ('20170320' + camid))
    filelist = filelist[mask]

    if (len(filelist) == 0):
        log.warn('No lightcurves found in the past 15 nights.')
        return

    log.info('Found {} good nights.'.format(len(filelist)))
    
    # Save the combined lightcurves.
    filename = 'cfast_{}{}.hdf5'.format(date, camid)
    filename = os.path.join(dirtree['sys'], filename)
    
    log.info('Saving combined photometry to {}'.format(filename))
    
    io.combine_photometry(filename, filelist, astrometry=False)
    
    # Perform the coarse decorrelation.
    log.info('Running the coarse decorrelation.')    
    
    sys = cdecor_vmag.CoarseDecorVmag(filename, 0)
    
    # Update the systable.
    log.info('Updating the systematics table.')
    io.update_systable(sys.sysfile, systable)
    
    # Remove the combined lightcurves.
    os.remove(filename)
    
    return


###############################################################################
### Main reduction loop.
###############################################################################

def build_dirtree(date, camid):

    log = logging.getLogger('lsreduce')    
    log.info('Building directory tree for {}{}.'.format(date, camid))
    
    dirtree = dict()
    dirtree['rawarchive'] = os.path.join(cfg.arcdir, date + camid)
    dirtree['tmp'] = cfg.tmpdir
    dirtree['binned'] = os.path.join(cfg.reddir, date + camid, 'binned')
    dirtree['lightcurves'] = os.path.join(cfg.reddir, date + camid, 'lightcurves')
    dirtree['thumbs'] = os.path.join(cfg.reddir, date + camid, 'thumbs')
    dirtree['sys'] = os.path.join(cfg.reddir, date + camid, 'sys')    
    dirtree['targets'] = os.path.join(cfg.reddir, date + camid, 'targets')    
    
    for key in dirtree.keys():
        
        try:
            os.makedirs(dirtree[key])
        except:
            log.warn('Directory exists, {}'.format(dirtree[key]))        
            
    return dirtree

def parse_files(filelist, nbias, ndark, nflat):
    
    filelist = np.sort(filelist)

    # Separate calibration and science frames.
    bias_frames = [filename for filename in filelist if 'bias' in filename]
    dark_frames = [filename for filename in filelist if 'dark' in filename]
    flat_frames = [filename for filename in filelist if 'flat' in filename]
    
    filelist = [filename for filename in filelist if 'bias' not in filename]
    filelist = [filename for filename in filelist if 'dark' not in filename]
    filelist = [filename for filename in filelist if 'flat' not in filename]
    
    science_frames = filelist
    
    if (len(science_frames) > 0):    
    
        # Extract the lst sequence number from the filename.
        lstseq = np.array([int(os.path.split(filename)[-1][:8]) for filename in science_frames])
        
        # Select the oldest set of 50 science frames present in the directory.
        setidx = lstseq//50
        nim = sum(setidx == np.amin(setidx))
        
    else:
        
        nim = 0

    bias_frames = bias_frames[:nbias]
    dark_frames = dark_frames[:ndark]
    flat_frames = flat_frames[:nflat]
    science_frames = science_frames[:nim]  
    remainder = bias_frames[nbias:] + dark_frames[ndark:] + flat_frames[nflat:] + science_frames[nim:]
    
    return science_frames, bias_frames, dark_frames, flat_frames, remainder

#def rereduce(rawdir, outdir, nscience=50, ndark=50):
#
#    import shutil
#
#    log = logging.getLogger('rereduce')
#    log.info('Starting rereduction.')
#
#    # Check if the directories exist.
#    if not os.path.exists(rawdir):
#        raise IOError('Directory {} does not exist.'.format(rawdir))
#        
#    if not os.path.exists(outdir):
#        raise IOError('Directory {} does not exist.'.format(outdir))
#
#    # Parse the raw location to obtain the date and camid.
#    head, tail = os.path.split(rawdir)
#    date, camid = tail[:8], tail[8:]
#
#    # Clean up the directory where the results will be written.
#    contents = listdir_fullpath(outdir)
#    olddir = os.path.join(outdir, 'old')
#    os.makedirs(olddir)
#    for name in contents:
#        shutil.move(name, olddir)
#        
#    # Create the dirtree for rereduction.
#    dirtree = dict()
#    dirtree['tmp'] = outdir
#    dirtree['binned'] = os.path.join(outdir, 'binned')
#    dirtree['lightcurves'] = os.path.join(outdir, 'lightcurves')
#    dirtree['thumbs'] = os.path.join(outdir, 'thumbs')
#    dirtree['sys'] = os.path.join(outdir, 'sys')    
#    dirtree['targets'] = os.path.join(outdir, 'targets')  
#    
#    for key in dirtree.keys():
#        
#        try:
#            os.makedirs(dirtree[key])
#        except:
#            log.warn('Directory exists, {}'.format(dirtree[key])) 
#
#    # Get the siteinfo, darktables and astrometry.
#    siteinfo = io.read_siteinfo(cfg.siteinfo, cfg.sitename)
#    darktables = cfg.darktables[camid] 
#    astromaster = cfg.astromaster[camid]
#    systable = cfg.systable[camid]
#
#    # Get the files to be processed.
#    filelist = listdir_fullpath(rawdir)
#
#    while (len(filelist) > 0):            
#
#        # Parse files.
#        science_frames, dark_frames, filelist = parse_files(filelist, nscience, ndark)
#
#        if (len(dark_frames) > 0):
#
#            reduction.reduce_dark_frames(camid, dark_frames, dirtree, darktables)
#
#        if (len(science_frames) > 0):
#
#            reduction.reduce_science_frames(camid, science_frames, siteinfo, dirtree, darktables, astromaster, systable)
#            
#    # Combine temporary lightcurve files.
#    combine_temporary_files(date, camid, dirtree)
#
#    # Send summary email.
#    summarize.reduction_summary(dirtree['lightcurves']) 
#
#    return

def rawdir_monitor(camid, twilight=5, nscience=49, nbias=20, ndark=20, nflat=20, timeout=10, step=6.4):    

    log = logging.getLogger('lsreduce')
    log.info('Initializing main reduction loop for camera {}.'.format(camid))    
    
    # Set up the listener.
    siteinfo = io.read_siteinfo(cfg.siteinfo, cfg.sitename)
    darktable = cfg.darktable[camid] 
    astromaster = cfg.astromaster[camid]
    systable = cfg.systable[camid]
    rawdir = cfg.rawdir
    
    in_queue = set()
    day_tasks_finished = True
    bias_time = 0
    dark_time = 0
    flat_time = 0
    science_time = 0
    
    # See if the reduction is being restarted.
    test = glob.glob(os.path.join(cfg.tmpdir, '*.txt'))
    
    if (len(test)  == 1):

        log.info('Restarting night.')        
        
        queuefile = test[0]
        
        # Get directory tree for this night.
        head, tail = os.path.split(queuefile)
        date = tail.rstrip('.txt')
        dirtree = build_dirtree(date, camid)
        
        # Get processed files for this night.
        in_queue = io.read_in_queue(queuefile)
        
        # Get contents of rawdir and add any processed files to the archive.
        filelist = listdir_fullpath(rawdir)
        filelist = in_queue.intersection(set(filelist))
        if (len(filelist) > 0):        
            io.archive_files(list(filelist), dirtree['rawarchive'])
                
        day_tasks_finished = False
        
    elif (len(test) == 0):
        queuefile = None
        
    else:
        log.warning('Found multiple .txt files in tmp directory, unknown situation.')
        exit()
   
    while True:
        
        # Get the sun altitude.
        log.info('Computing current sun position.')
        ra, dec, sun_alt = astrometry.sun_position(siteinfo)

        # Obtain the current list of raw files.
        log.info('Getting filelist.')
        filelist = listdir_fullpath(rawdir)
        filelist = set(filelist)        
        
        # See which files are not in the processing queue yet.
        log.info('Getting unprocesed files.')
        new_files = filelist.difference(in_queue)
        new_files = list(new_files)

        time.sleep(step)

        if (sun_alt < twilight):
        
            if day_tasks_finished:
                
                log.info('Beginning of new night.')
                
                # Build directory tree for the current night.
                date = get_datestring()                
                dirtree = build_dirtree(date, camid)                               
                
                # Write a queuefile to indicate the night has begun.
                queuefile = os.path.join(dirtree['tmp'], date + '.txt')                
                io.write_in_queue(queuefile, set())                
                
                day_tasks_finished = False
        
            # Parse files.
            log.info('Parsing files.')
            science_frames, bias_frames, dark_frames, flat_frames, remainder = parse_files(new_files, nscience, nbias, ndark, nflat)        
            log.debug('{}, {}, {}, {}'.format(bias_time, dark_time, flat_time, science_time))
            
            if (len(bias_frames) == nbias) | (bias_time >= (nbias + timeout)):
                
                # Reduce the frames.
                reduction.reduce_bias_frames(camid, bias_frames, dirtree)
                
                # Archive the files.
                in_queue.update(bias_frames)
                io.write_in_queue(queuefile, in_queue)
                io.archive_files(bias_frames, dirtree['rawarchive'])

                bias_time = 0
                
            elif (len(bias_frames) > 0):
                
                if (bias_time == 0):
                    bias_time = len(bias_frames)
                else:
                    bias_time += 1
            
            
            if (len(dark_frames) == ndark) | (dark_time >= (ndark + timeout)):
                
                # Reduce the frames.
                reduction.reduce_dark_frames(camid, dark_frames, dirtree, darktable)
                
                # Archive the files.
                in_queue.update(dark_frames)
                io.write_in_queue(queuefile, in_queue)
                io.archive_files(dark_frames, dirtree['rawarchive'])

                dark_time = 0
                
            elif (len(dark_frames) > 0):
                
                if (dark_time == 0):
                    dark_time = len(dark_frames)
                else:
                    dark_time += 1
                    
            if (len(flat_frames) == nflat) | (flat_time >= (nflat + timeout)):
                
                # Reduce the frames.
                reduction.reduce_flat_frames(camid, flat_frames, dirtree, darktable)
                
                # Archive the files.
                in_queue.update(flat_frames)
                io.write_in_queue(queuefile, in_queue)
                io.archive_files(flat_frames, dirtree['rawarchive'])

                flat_time = 0
                
            elif (len(flat_frames) > 0):
                
                if (flat_time == 0):
                    flat_time = len(flat_frames)
                else:
                    flat_time += 1
        
            if (len(science_frames) == nscience) | (science_time >= (nscience + timeout)):
                
                # Reduce the frames.
                reduction.reduce_science_frames(camid, science_frames, siteinfo, dirtree, darktable, astromaster, systable)

                # Archive the files.
                in_queue.update(science_frames)
                io.write_in_queue(queuefile, in_queue)
                io.archive_files(science_frames, dirtree['rawarchive'])
                         
                science_time = 0            
                
            elif (len(science_frames) > 0):
                
                if (science_time == 0):
                    science_time = len(science_frames)
                else:
                    science_time += 1
            
        elif (sun_alt >= twilight) & (~day_tasks_finished):
            
            log.info('End of night, processing all remaining raw files.')            
            
            # Process all remaining files.
            while (len(new_files) > 0):            
            
                # Parse files.
                science_frames, bias_frames, dark_frames, flat_frames, remainder = parse_files(new_files, nscience, nbias, ndark, nflat)
            
                if (len(bias_frames) > 0):
                    
                    reduction.reduce_bias_frames(camid, bias_frames, dirtree)
                    
                    # Archive the files.
                    in_queue.update(bias_frames)  
                    io.write_in_queue(queuefile, in_queue)
                    io.archive_files(bias_frames, dirtree['rawarchive'])            
            
                if (len(dark_frames) > 0):
                    
                    reduction.reduce_dark_frames(camid, dark_frames, dirtree, darktable)
                    
                    # Archive the files.
                    in_queue.update(dark_frames)  
                    io.write_in_queue(queuefile, in_queue)
                    io.archive_files(dark_frames, dirtree['rawarchive'])
                    
                if (len(flat_frames) > 0):
                    
                    reduction.reduce_flat_frames(camid, flat_frames, dirtree, darktable)
                    
                    # Archive the files.
                    in_queue.update(flat_frames)  
                    io.write_in_queue(queuefile, in_queue)
                    io.archive_files(flat_frames, dirtree['rawarchive']) 
            
                if (len(science_frames) > 0):
                    
                    reduction.reduce_science_frames(camid, science_frames, siteinfo, dirtree, darktable, astromaster, systable)
                    
                    # Archive the files.
                    in_queue.update(science_frames)
                    io.write_in_queue(queuefile, in_queue)
                    io.archive_files(science_frames, dirtree['rawarchive'])
                  
            log.info('Finished processing raw files, starting daytime tasks.')                  
                  
            # Combine temporary lightcurve files.
            log.info('Combing temporary lightcurves.')
            combine_temporary_files(date, camid, dirtree)
            os.remove(queuefile)            
            
            # Send summary email.
            summarize.reduction_summary(dirtree['lightcurves'])            
            
            # Run daytime calibration. 
            log.info('Running daytime calibration.')
            calibrate_photometry(date, camid, dirtree, systable)
            
            # Send summary email.
            summarize.calibration_summary(dirtree['sys'], astromaster[1])            
            
            log.info('Finished daytime taks.')            
            
            # Reset.
            day_tasks_finished = True
            science_time = 0
            bias_time = 0
            dark_time = 0
            flat_time = 0
            in_queue = set()               
            
        elif day_tasks_finished:
            log.info('Waiting for sunset.')
                    
        else:
            log.warn('Encountered unknown situation.')
          
    return


def main():
    return


if __name__ == '__main__':
    logging.basicConfig(filename='/home/talens/MASCARA/lsreduce/example.log', level=logging.DEBUG)
    log = logging.getLogger('lsreduce')

    main()
