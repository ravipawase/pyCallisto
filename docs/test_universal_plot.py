
"""
Purpose of this program is to test different functions together

"""

"""import sys
sys.path.append('../src/')	#this and above line is added because 'PyCallisto.py' and
							#'pyCallistoUtils.py' are not in the same folder as this script."""
import pyCallisto as pyc
import pyCallistoUtils as utils
#import pyfits
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

#fits1_path = '../data/IISERP_20151104_031152_59.fit'
#fits2_path = '../data/IISERP_20151104_032652_59.fit'

#fits1_path = '../data/GAURI_20151104_042959_59.fit'
#fits2_path =  '../data/GAURI_20151104_044459_59.fit'

fits1_path = 'data/GAURI_20151104_033000_59.fit'			#####
fits2_path =  'data/GAURI_20151104_034500_59.fit'		#####

#plot multiple files
fits1 = pyc.pyCallisto.fromFile(fits1_path)
plt = fits1.spectrogram() #this will show in imshow thing
plt.savefig("fits1.png")


fits1.universalPlot()


