#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  6 14:50:12 2025

@author: mrutyunjaya
"""

import pyCallisto as pyc
import pyCallistoUtils as utils
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from astropy.visualization import time_support
import numpy as np

plt.close('all')

fits_paths = ['data/GREENLAND_20240716_130442_62.fit.gz',
              'data/GREENLAND_20240716_131942_62.fit.gz',
              'data/GREENLAND_20240716_133442_62.fit.gz', 
              'data/GREENLAND_20240716_134942_62.fit.gz',
              'data/GREENLAND_20240716_140443_62.fit.gz', 
              'data/GREENLAND_20240716_141942_62.fit.gz']
#Initialise the joined obejct with the first FITS file
joined=pyc.pyCallisto.fromFile(fits_paths[0])

for fits_path in fits_paths:
    fits=pyc.pyCallisto.fromFile(fits_path)
    
    #Generate spectrogram
    #plt.figure()
    #fits.spectrogram()
    #plt.show()
for fits_path in fits_paths[1:]:
    joined = joined.appendTimeAxis(fits_path)
#joined.spectrogram() 
#plt1 = joined.subtractBackground()
#plt1.spectrogram(xtick=15,fontsize=20)
#plt0=joined.sliceTimeAxis("01:00:00","04:15:00")

#plt2= plt0.sliceFrequencyAxis(12,80)
#plt2.spectrogram()
#plt2= plt0.sliceFrequencyAxis(11.8,85)
plt3= joined.subtractBackground()
plt3.spectrogram(xtick=2)
joined_bg_subtracted = joined.subtractBackground()
joined_bg_subtracted.spectrogramWithGOES(
    tstart="2024-07-16 13:05",
    tend="2024-07-16 14:34",
    satellite_number=16,    
    xtick=5,
    save_path='spectrogram_with_goes_overlay_2024-07-16.png'
)
plt.show()