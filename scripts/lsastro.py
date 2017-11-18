#!/usr/bin/env python

import numpy as np

from lsreduce import astromaster

# Fiducial pointings: ha0, dec0, x0, y0, theta0
fiducial = dict()

# MASCARA La Silla.
fiducial['LSC'] = [np.array([-0.6, -32.2 ,2004., 1336., 180.]), None]
fiducial['LSN'] = np.array([6.7, 11.2 ,2004., 1336., 0.])
fiducial['LSE'] = np.array([-42.1, -15.6 ,2004., 1336., 290.])
fiducial['LSW'] = np.array([47.0, -27.9, 2004., 1336., 70.])

astromask = np.recarray((1,), dtype=[('lx', 'uint16'), ('ux', 'uint16'), ('ly', 'uint16'), ('uy', 'uint16')])
astromask['lx'] = 2850 
astromask['ux'] = 4008  
astromask['ly'] = 2275 
astromask['uy'] = 2672 
fiducial['LSS'] = [np.array([-6., -70. ,2004., 1336., 185.]), astromask]

astromask = np.recarray((1,), dtype=[('lx', 'uint16'), ('ux', 'uint16'), ('ly', 'uint16'), ('uy', 'uint16')])
astromask['lx'] = 0 
astromask['ux'] = 400  
astromask['ly'] = 2275 
astromask['uy'] = 2672 
fiducial['LSW'] = [np.array([47.0, -27.9, 2004., 1336., 70.]), astromask]

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
       
    astromaster.astromaster(args.image, fiducial[args.camera], args.catalogue, args.dark, args.outpath)