# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 17:18:01 2016

@author: talens
"""
import os
# Name of the observing site.
sitename = 'Sutherland' 

confdir = r'y:\conf'
rawdir = r'd:\raw'
reddir = r'd:\products' # Location were data products are written to.
arcdir = r'd:\rawarchive'
tmpdir = r'd:\tmp'

###############################################################################
### Locations of various files.
###############################################################################

starcat = os.path.join(confdir, 'bringcat20170327.fits') # Location of the stellar catalogue.
siteinfo = os.path.join(confdir, 'siteinfo.dat') # Location of the siteinfo file.
targets = os.path.join(confdir, 'targets.dat') # Location of the targets file.

# Locations of the darktables.
darktables = {'SAE': [os.path.join(confdir, 'darktableSAElong.dat'),
                      os.path.join(confdir, 'darktableSAEshort.dat')],
              'SAW': [os.path.join(confdir, 'darktableSAWlong.dat'),
                      os.path.join(confdir, 'darktableSAWshort.dat')]}
              
# Locations of the astrometry master solutions.
astromaster = {'SAE': os.path.join(confdir, 'astromasterSAE.hdf5'),
               'SAW': os.path.join(confdir, 'astromasterSAW.hdf5')}

# Locations of the systables.
systable = {'SAE': os.path.join(confdir, 'systableSAE.dat'),
            'SAW': os.path.join(confdir, 'systableSAW.dat')}

###############################################################################
### Parameters of the reduction.
###############################################################################

mindarks = 10 # Minimum number of darks needed to create a new masterdark.

maglim_astro = 7.5 # Magnitude limit for solving the astrometry.
maglim_fast = 8.4 # Magnitude limit for the fast photometry.
maglim_slow = 10. # Magnitude limit for the slow photometry.

aper_fast = [2.5, 4.5] # Apertures for the fast photometry.
aper_slow = [2.5, 3.5, 4.5, 5.5] # Apertures for the slow photometry.
skyrad = [6, 21] # Inner and outer radius of the sky annulus.

nhalf = 25 # Size of postage stamps.