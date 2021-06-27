
"""
Purpose of this program is to test different functions together

"""

import sys
sys.path.append('../src/')	#this and above line is added because 'PyCallisto.py' and
							#'pyCallistoUtils.py' are not in the same folder as this script.
import pyCallisto as pyc
import pyCallistoUtils as utils
#import pyfits
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

fits1_path = '../data/IISERP_20151104_031152_59.fit'
fits2_path = '../data/IISERP_20151104_032652_59.fit'

#fits1_path = '../data/GAURI_20151104_042959_59.fit'
#fits2_path =  '../data/GAURI_20151104_044459_59.fit'

#fits1_path = '../data/GAURI_20151104_033000_59.fit'			#####
#fits2_path =  '../data/GAURI_20151104_034500_59.fit'		#####

#fits1_path = '../data/GAURI_20151104_035959_59.fit'
#fits2_path =  '../data/GAURI_20151104_041459_59.fit'

#fits1_path = '../data/data.fits'
#fits2_path =  '../data/data_1.fits'


#plot multiple files
fits1 = pyc.PyCallisto.from_file(fits1_path)
plt = fits1.spectrogram() #this will show in imshow thing
plt.savefig("fits1.png")


#join time axis
joined1 = fits1.append_time_axis(fits2_path)
plt = joined1.spectrogram() #this will show in imshow thing
plt.savefig("joined.png")

#slice in frequency axis
freq_sliced = joined1.slice_frequency_axis("200", "400")
plt = freq_sliced.spectrogram() #this will show in imshow thing
plt.savefig("freq_sliced.png")


#slice in time axis
time_sliced = freq_sliced.slice_time_axis("03:21:00", "03:31:00")
plt = time_sliced.spectrogram() #this will show in imshow thing
plt.savefig("time_sliced.png")


#do background subtraction
background_subtracted = time_sliced.subtract_background()
plt = background_subtracted.spectrogram()
plt.savefig("background_subtracted.png")


#get meanlightcurve
background_subtracted.mean_light_curve(out_image="mean_Light_Curve.png", grid=True)

#get meanSpectrum
background_subtracted.mean_spectrum(out_image="mean_spectrum.png", grid=True)

#get light curve at one frequency
background_subtracted.light_curve(300, out_image="Lightcurve.png", grid=True)


#get spectrum
background_subtracted.spectrum( '2015/11/04','03:25:00', out_image="singletimespectrum.png", grid=True)


