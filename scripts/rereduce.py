#!/usr/bin/env python

import sys, os
sys.path.insert(0, os.getcwd())
from astropy.utils.iers import conf as iersconf
iersconf.iers_auto_url = 'http://www.strw.leidenuniv.nl/~stuik/finals2000A.all'
import logging.config

logging.basicConfig(level=logging.DEBUG)

from bringreduce import rawdir_monitor

if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Initialize the reduction loop.')
    parser.add_argument('rawdir', type=str,
                        help='the location of the files to be rereduced')
    parser.add_argument('outdir', type=str,
                        help='the location the results should be written to')
    args = parser.parse_args()

    rawdir_monitor.rereduce(args.rawdir, args.outdir, ndark=50)

