# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 17:18:01 2016

@author: talens
"""

import os

# Name of the observing site.
sitename = 'LaSilla' 

confdir = r'd:\conf'
rawdir = r'd:\raw'
reddir = r'd:\products' # Location were data products are written to.
arcdir = r'd:\rawarchive'
tmpdir = r'd:\tmp'

###############################################################################
### Locations of various files.
###############################################################################

starcat = os.path.join(confdir, 'bringcat20180428.fits') # Location of the stellar catalogue.
siteinfo = os.path.join(confdir, 'siteinfo.dat') # Location of the siteinfo file.
targets = os.path.join(confdir, 'targets.dat') # Location of the targets file.

# Locations of the darktables.
darktable = {'LSN': os.path.join(confdir, 'darktableLSN.dat'),
             'LSE': os.path.join(confdir, 'darktableLSE.dat'),
             'LSS': os.path.join(confdir, 'darktableLSS.dat'),
             'LSW': os.path.join(confdir, 'darktableLSW.dat'),
             'LSC': os.path.join(confdir, 'darktableLSC.dat')}
              
# Locations of the astrometry master solutions.
astromaster = {'LSN': [os.path.join(confdir, 'astromasterLSNdusk.hdf5'),
                       os.path.join(confdir, 'astromasterLSNnight.hdf5'),
                       os.path.join(confdir, 'astromasterLSNdawn.hdf5')],
               'LSE': [os.path.join(confdir, 'astromasterLSEdusk.hdf5'),
                       os.path.join(confdir, 'astromasterLSEnight.hdf5'),
                       os.path.join(confdir, 'astromasterLSEdawn.hdf5')],
               'LSS': [os.path.join(confdir, 'astromasterLSSdusk.hdf5'),
                       os.path.join(confdir, 'astromasterLSSnight.hdf5'),
                       os.path.join(confdir, 'astromasterLSSdawn.hdf5')],
               'LSW': [os.path.join(confdir, 'astromasterLSWdusk.hdf5'),
                       os.path.join(confdir, 'astromasterLSWnight.hdf5'),
                       os.path.join(confdir, 'astromasterLSWdawn.hdf5')],
               'LSC': [os.path.join(confdir, 'astromasterLSCdusk.hdf5'),
                       os.path.join(confdir, 'astromasterLSCnight.hdf5'),
                       os.path.join(confdir, 'astromasterLSCdawn.hdf5')]}

# Locations of the systables.
systable = {'LSN': os.path.join(confdir, 'systableLSN.dat'),
            'LSE': os.path.join(confdir, 'systableLSE.dat'),
            'LSS': os.path.join(confdir, 'systableLSS.dat'),
            'LSW': os.path.join(confdir, 'systableLSW.dat'),
            'LSC': os.path.join(confdir, 'systableLSC.dat')}

###############################################################################
### Parameters of the reduction.
###############################################################################

minbias = 10 # Minimum number of bias frames needed to create a new masterbias.
mindark = 10 # Minimum number of dark frames needed to create a new masterdark.
minflat = 10 # Minimum number of flat frames needed to create a new masterflat.

maglim_astro = 7.5 # Magnitude limit for solving the astrometry.
maglim_fast = 8.4 # Magnitude limit for the fast photometry.
maglim_slow = 10. # Magnitude limit for the slow photometry.

aper_fast = [2.5, 4.5] # Apertures for the fast photometry.
aper_slow = [2.5, 3.5, 4.5, 5.5] # Apertures for the slow photometry.
skyrad = [6, 21] # Inner and outer radius of the sky annulus.

nhalf = 25 # Size of postage stamps.