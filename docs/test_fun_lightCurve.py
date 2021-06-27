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
plt = fits1.spectrogram()
plt.savefig('test.png')
fits1.light_curve(300, out_image="Lightcurve_300.png")
fits1.light_curve(400, out_image="Lightcurve_400.png")
fits1.light_curve(500, out_image="Lightcurve_500.png")


#def lightCurve(self, frequency, plot=True, outimage ="Lightcurve.png", returndata=False):
