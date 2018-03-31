import pyCallisto as pyc
import pyCallistoUtils as utils
#import pyfits
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

fits1_path = '../Data/IISERP_20151104_031152_59.fit'
fits2_path = '../Data/IISERP_20151104_032652_59.fit'


#plot multiple files
fits1 = pyc.pyCallisto.fromFile(fits1_path)
#plot = fits1.showSpectrum() #this will show in imshow thing
#plot.savefig("testtttttttttttttttttttt.png")


#join time axis
joined1 = fits1.joinTimeaxis(fits2_path)
sliced = joined1.sliceTimeaxis("03:23:00", "03:29:00")
freq_sliced = sliced.sliceFrequencyaxis("200", "600")
background_subtracted = freq_sliced.subtractBackground()
plot = background_subtracted.plotSpectrum()
#utils.visualise(plot)
utils.visualise(plot, show=False, outpath= 'test.png')


#background_subtracted.gettimeseries()
#background_subtracted.getfreqseries()




