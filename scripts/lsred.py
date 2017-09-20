#!/usr/bin/env python

import sys, os
sys.path.insert(0, os.getcwd())
import traceback

from astropy.utils.iers import conf as iersconf
iersconf.iers_auto_url = 'http://www.strw.leidenuniv.nl/~stuik/finals2000A.all'

import logging.config
logging.config.fileConfig('./lsreduce/logging.conf')
log = logging.getLogger('lsreduce')

from lsreduce import rawdir_monitor

if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Initialize the reduction loop.')
    parser.add_argument('camera', type=str, choices=['LSN', 'LSE', 'LSS', 'LSW', 'LSC'],
                        help='the camera to initialize the reduction loop for')
    args = parser.parse_args()

    try:
        rawdir_monitor.rawdir_monitor(args.camera)
    except Exception as e:
        log.error('Reduction Fault:' + traceback.format_exc())
        raise e

