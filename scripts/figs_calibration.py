# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 16:24:07 2017

@author: talens
"""

from lsreduce import summarize

def main(filelist, astromaster):
    
    for filename in filelist:
        
        summarize.fig_transmission(filename, astromaster)
        summarize.fig_intrapix(filename, astromaster)
        summarize.fig_clouds(filename)
    
    return
    
if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Make figures of calibration terms.')
    parser.add_argument('files', type=str, nargs='+',
                        help='the file(s) to create the figures for')
    parser.add_argument('astromaster', type=str,
                        help='the astrometric solution to use when creating the figures')
    args = parser.parse_args()

    main(args.files, args.astromaster)
