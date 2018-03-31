import datetime as dt
#import pyfits
import astropy.io.fits as pyfits
#import matplotlib.pyplot as plt
#import matplotlib.dates as mdates
#from matplotlib import cm
#import os
#import numpy as np
#from matplotlib.dates import  DateFormatter
#import sys



def checkFitsCallisto(fitsfile):
	"""Check whether fits file has two HDUs or not
	"""
	hdus = pyfits.open(fitsfile)
	if len(hdus) == 2:
		hdus.close()
		#del hdus
		return True
	else:
		hdus.close()
		#del hdus
		return False

def checkBInTable(fitsfile):
	"""Check whether fits file bintable is valid
	"""
	hdus = pyfits.open(fitsfile)
	bintablehdu = hdus[1]

	try:
		bintablehdu.data
		hdus.close()
		#del hdus
		return True
	except:
		hdus.close()
		#del hdus
		return False


def tosec(td):
	"""Calculate the total seconds in a timedate.timedate object
	"""
	return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10 ** 6) / 10 ** 6


def toDate(string):
	"""Break the string with "/" separator
	return the equivalent datetime.date object
	"""
	yr, mnth, day = string.split('/')
	yr, mnth, day = int(yr), int(mnth), int(day) 
	return dt.date(yr, mnth, day)


def toTime(string):
	"""
	break the string with "/" separator
	return the datetime.time object
	"""
	hr, mn, sec = string.split(':')
	sec = sec.split('.')[0]  # if sec has a fractional value
	hr, mn, sec = int(hr), int(mn), int(sec)
	#if second's value is 60, replace it to 59
	#some old data from some observatories has this "leap second" problem
	if sec == 60:
		sec = 59
	return dt.time(hr, mn, sec)


def visualise(plt_object,  show =True, outpath= 'tmp.png'):
	"""
	input 
		1) matplotlib.pyplot object (plt)
		2) outpath : path of the image to be stored
		3) show (keyword) = if True, just show, if False save instead plotting
	
	returns
		returns nothing , saves or shows plot
	"""
	
	if show:
		plt_object.show()
	else:
		plt_object.savefig(outpath)

