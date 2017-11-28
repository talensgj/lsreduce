#!/usr/bin/env python

import numpy as np

from lsreduce import astromaster

# Fiducial pointings: ha0, dec0, x0, y0, theta0
fiducial = dict()

fiducial['LSC'] = np.array([-0.1, -32.5 ,2004., 1336., 183.2])
fiducial['LSN'] = np.array([1.4, 12.1 ,2004., 1336., 1.7])
fiducial['LSE'] = np.array([-45.1, -19.8 ,2004., 1336., 291.6])
fiducial['LSS'] = np.array([-3.6, -69.7 ,2004., 1336., 184.7])
fiducial['LSW'] = np.array([45.0, -23.9, 2004., 1336., 72.6])

# Masked regions.
astromasks = dict()
dtype = [('lx', 'uint16'), ('ux', 'uint16'), ('ly', 'uint16'), ('uy', 'uint16')]

astromasks['LSC'] = None
astromasks['LSN'] = None
astromasks['LSE'] = None
astromasks['LSS'] = np.array([(2850, 4008, 2275, 2672)], dtype=dtype)
astromasks['LSW'] = np.array([(0, 400, 2275, 2672)], dtype=dtype)

if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Create astrometric solution.')
    parser.add_argument('image', type=str, 
                        help='The path to the image to use for creating the solution.')
    parser.add_argument('camera', type=str, choices=['LSN', 'LSE', 'LSS', 'LSW', 'LSC'],
                        help='The camera for which we are creating a solution.')
    parser.add_argument('catalogue', type=str, 
                        help='The path to the stellar catalogue.')
    parser.add_argument('--dark', type=str, default=None,
                        help='The path to the dark to use for creating the solution.')
    parser.add_argument('--outpath', type=str, default='.',
                        help='The directory where the solution should be written.')
    args = parser.parse_args()    
       
    astromaster.astromaster(args.image, fiducial[args.camera], args.catalogue, astromasks[args.camera], args.dark, args.outpath)