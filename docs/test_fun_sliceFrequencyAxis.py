import sys
sys.path.append('../src/')	#this and above line is added because 'pyCallisto.py' and 
							#'pyCallistoUtils.py' are not in the same folder as this script.
import pyCallisto as pyc
import pyCallistoUtils as utils
#import pyfits
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

fits1_path = '../data/IISERP_20151104_031152_59.fit'

fits1 = pyc.pyCallisto.fromFile(fits1_path)
plt = fits1.spectrogram()
plt.savefig('test.png')
plt.clf()

freq_sliced = fits1.sliceFrequencyAxis(300, 400)

plt = freq_sliced.spectrogram()
plt.savefig('freq_sliced_300_400.png')


