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



fits1 = pyc.pyCallisto.fromFile(fits1_path)
plt = fits1.spectrogram()
plt.savefig('test.png')
fits1.spectrum("2015/11/04", "03:12:53")	#both python objects

