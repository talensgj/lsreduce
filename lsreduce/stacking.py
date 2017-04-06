# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 16:32:11 2016

@author: talens
"""

import numpy as np

from . import astrometry

class Stacker(object):
    
    def __init__(self, wcspars, **kwargs):
        
        self.wcspars = wcspars
        
        self.nx = kwargs.get('nx', 4008)
        self.ny = kwargs.get('ny', 2672)

        self.sidex = kwargs.get('sidex', 32)
        self.sidey = kwargs.get('sidey', 32)
        self.margin = kwargs.get('margin', 50)
               
        self._make_grid()
        self._make_images()

        return
           
    def _make_grid(self):
        
        # Number of tiles.
        self.ntile_x = self.nx//self.sidex
        self.ntile_y = self.ny//self.sidey   
     
        # Number of discarded pixels.
        resx = self.nx%self.sidex
        resy = self.ny%self.sidey
        
        if (resx%2 > 0) | (resy%2 > 0):
            print 'Uneven number of discarded pixels, please adjust patch size.'
            exit()
            
        # Centers of the tiles.
        gx = np.arange(self.ntile_x)*self.sidex + self.sidex/2. + resx/2.
        gy = np.arange(self.ntile_y)*self.sidey + self.sidey/2. + resy/2.
        gx, gy = np.meshgrid(gx, gy)        

        self.gx = gx.ravel()
        self.gy = gy.ravel() 

        # Size of binned image.       
        self.bin_nx = self.nx + 2*self.margin
        self.bin_ny = self.ny + 2*self.margin  
       
        return        
    
    def _make_images(self):
        
        self.bin_image = np.zeros((self.bin_ny, self.bin_nx), dtype='float64')
        self.bin_count = np.zeros((self.bin_ny, self.bin_nx), dtype='int') 
        self.exptime = 0. 
        self.nimages = 0       

        return
         
    def add_image(self, lst, image, exptime, lst0):

        from itertools import izip

        # Shift tile coordinates to reference time.  
        ra, dec = astrometry.wcs2world(self.wcspars, self.gx-.5, self.gy-.5, lst)
        cx, cy = astrometry.world2wcs(self.wcspars, ra, dec, lst0)        
        
        cx, cy = cx+.5, cy+.5        
        cx, cy = np.around(cx), np.around(cy)
        cx, cy = cx.astype('int'), cy.astype('int') # Speed gain!

        gx, gy = np.around(self.gx), np.around(self.gy)
        gx, gy = gx.astype('int'), gy.astype('int')

        for i, j, k, l in izip(gy, gx, cy, cx):

            try:
                self.bin_image[self.margin+k-self.sidey/2:self.margin+k+self.sidey/2, self.margin+l-self.sidex/2:self.margin+l+self.sidex/2] += image[i-self.sidey/2:i+self.sidey/2,j-self.sidex/2:j+self.sidex/2]
                self.bin_count[self.margin+k-self.sidey/2:self.margin+k+self.sidey/2, self.margin+l-self.sidex/2:self.margin+l+self.sidex/2] += 1
            except ValueError:
                continue

        self.exptime += exptime        
        self.nimages += 1

        return
        
    def get_image(self, lst0, reset=True):

        # Get the final binned image.
        bin_image = self.bin_image[self.margin:-self.margin, \
                                   self.margin:-self.margin]
        bin_count = self.bin_count[self.margin:-self.margin, \
                                   self.margin:-self.margin]
        
        with np.errstate(invalid='ignore'):
            bin_image = bin_image/bin_count
        
        bin_image = np.nan_to_num(bin_image) 

        # Get the exposure time and number of images.
        exptime = self.exptime
        nimages = self.nimages
          
        if (nimages > 0):     
            exptime = exptime/nimages

        # Reset the images.
        if reset:
            self._make_images()

        w = astrometry.create_wcs(self.wcspars, lst0)
        bin_header = w.to_header()

        return bin_image, bin_header

