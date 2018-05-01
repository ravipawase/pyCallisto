
"""
Purpose of this program is to test different functions together

"""

import sys
sys.path.append('../src/')	#this and above line is added because 'pyCallisto.py' and 
							#'pyCallistoUtils.py' are not in the same folder as this script.
import pyCallisto as pyc
import pyCallistoUtils as utils
#import pyfits
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

fits1_path = '../data/IISERP_20151104_031152_59.fit'
fits2_path = '../data/IISERP_20151104_032652_59.fit'


#plot multiple files
fits1 = pyc.pyCallisto.fromFile(fits1_path)
plt = fits1.spectrogram() #this will show in imshow thing
plt.savefig("fits1.png")


#join time axis
joined1 = fits1.appendTimeAxis(fits2_path)
plt = joined1.spectrogram() #this will show in imshow thing
plt.savefig("joined.png")

#slice in frequency axis
freq_sliced = joined1.sliceFrequencyAxis("200", "600")
plt = freq_sliced.spectrogram() #this will show in imshow thing
plt.savefig("freq_sliced.png")


#slice in time axis
time_sliced = freq_sliced.sliceTimeAxis("03:23:00", "03:29:00")
plt = time_sliced.spectrogram() #this will show in imshow thing
plt.savefig("time_sliced.png")


#do background subtraction
background_subtracted = freq_sliced.subtractBackground()
plt = background_subtracted.spectrogram()
plt.savefig("background_subtracted.png")


#get meanlightcurve
background_subtracted.meanLightCurve(outImage ="mean_Light_Curve.png", grid=True)

#get meanSpectrum
background_subtracted.meanSpectrum(outImage ="mean_spectrum.png", grid=True)

#get light curve at one frequency
background_subtracted.lightCurve(300, outImage ="Lightcurve.png", grid=True)


#get spectrum
background_subtracted.spectrum( '2015/11/04','03:26:00', outImage ="singletimespectrum.png", grid=True)


