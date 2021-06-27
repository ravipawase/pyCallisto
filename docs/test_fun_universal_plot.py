import sys
sys.path.append('../src/')	#this and above line is added because 'PyCallisto.py' and
							#'pyCallistoUtils.py' are not in the same folder as this script.
import pycallisto as pyc
import pycallisto_utils as utils
#import pyfits
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

fits1_path = '../data/IISERP_20151104_031152_59.fit'
fits2_path = '../data/IISERP_20151104_032652_59.fit'



fits1 = pyc.PyCallisto.from_file(fits1_path)
plt = fits1.spectrogram(option=3)
plt.savefig('spectrogram_option1.png')

fits1.universal_plot()

