import datetime as dt
#import pyfits
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import cm
import os
import numpy as np
from matplotlib.dates import  DateFormatter
import sys
import matplotlib
import io
import pyCallistoUtils as utils		#local utility file
import copy

class pyCallisto:

	def __init__(self, HDUList):
		"""
		Create pyCallisto object from HDUList
		"""
		self.hdus = HDUList
		self.imagehdu = self.hdus[0]
		self.bintablehdu = self.hdus[1]
		self.imageheader = self.imagehdu.header
		# get the datamin amd maximum for colorbar
		self.datamin = int(self.imagehdu.data.min())
		self.datamax = int(self.imagehdu.data.max())
		self.datamid = (self.datamin + self.datamax) / 2

	def __del__(self):
		#print("Clearing the memory by deleting the HDUlist object nothandled correctly by pyfits")
		self.hdus.close()


	@classmethod
	def fromFile(cls, infits):
		"""
		Create pyCallisto object from fits file on disc
		"""
		
		# checking number of HDUs
		if not utils.checkFitsCallisto(infits):
			#through an error here
			raise ValueError("No. of HDUs are wrong in input fits file, not a proper callisto file , cannot proceed")
			#print("No. of HDUs are wrong, not a proper callisto file ")
			#return -1
		
		#check bintable data accessible
		if not utils.checkBInTable(infits):
			#through an error here
			raise ValueError("Bintable data may be corrupted in input fits file, cannot proceed")
			#print "Problem with bintabledata"
			#return -1
		HDUList = pyfits.open(infits)
		return cls(HDUList)






	def plotSpectrum(self, option=3, xtick= 2, blevel=0, endpts= [False, False], cmap=cm.jet, cbar=True, cbar_ori='vertical', fontsize=14):

		#plot the fig 
		fig, ax = plt.subplots(figsize=(7,7))
		
		if option == 1:
			y, x = self.imagehdu.data.shape
			cax = ax.imshow(self.imagehdu.data, extent=[0, x, 0, y], aspect='auto', cmap=cmap, vmin=blevel)
			if cbar == True:
				ticks =list( np.linspace(blevel,self.datamax, 10).astype('int')) #calculate 10 ticks positins 
				if cbar_ori == 'horizontal':
					cbar = fig.colorbar(cax, ticks = ticks, orientation='horizontal')
				else:
					cbar = fig.colorbar(cax, ticks = ticks)
				cbar.set_label('Intensity', rotation= 90)		
		
			plt.xlabel('Row Count')
			plt.ylabel('Column Count')


		if option == 2:
			xstart = int(self.imageheader['CRVAL1'])
			xstep = float(self.imageheader['CDELT1'])  
			xlength = int(self.imageheader['NAXIS1'])
			xend = xstart + xstep * xlength
			freqs = self.bintablehdu.data['frequency'][0]  # these are true frequencies
			cax = ax.imshow(self.imagehdu.data, extent=[xstart, xend, freqs[-1], freqs[0]], aspect='auto', cmap=cmap, vmin=blevel)
		
			if cbar == True:
				ticks =list( np.linspace(blevel,self.datamax, 10).astype('int')) #calculate 10 ticks positins 
				if cbar_ori == 'horizontal':
					cbar = fig.colorbar(cax, ticks = ticks, orientation='horizontal')
				else:
					cbar = fig.colorbar(cax, ticks = ticks)
				cbar.set_label('Intensity', rotation= 90)
			
			if endpts[1]:
				#set y minorticks for first and last frequency
				y_lims = [freqs[-1], freqs[0]]
				#print(y_lims)
				plt.yticks(list(plt.yticks()[0]) + y_lims)
		
			plt.xlabel('Time (sec of day)')
			plt.ylabel('Frequency [MHz]')

		if option == 3:
			# get the start time and end time of observation
			start_date = utils.toDate(self.imageheader['DATE-OBS'])
			starttime = utils.toTime(self.imageheader['TIME-OBS'])
	
	
			starttime = dt.datetime.combine(start_date, starttime)  # make a datetime object
			endtime = utils.toTime(self.imageheader['TIME-END'])
			endtime = dt.datetime.combine(start_date, endtime)
	
			# get the frequencies
			freqs = self.bintablehdu.data['frequency'][0]  # these are true frequencies
	
			# set the limits for plotting 
			x_lims = [starttime, endtime]
			x_lims = mdates.date2num(x_lims)  # dates to numeric values
			y_lims = [freqs[-1], freqs[0]]

	
	
			cax = ax.imshow(self.imagehdu.data, extent=[x_lims[0], x_lims[1], y_lims[0], y_lims[1]], aspect='auto', cmap=cmap, vmin=blevel)
			if cbar == True:
				ticks =list( np.linspace(blevel,self.datamax, 10).astype('int')) #calculate 10 ticks positins 
				if cbar_ori == 'horizontal':
					cbar = fig.colorbar(cax, ticks = ticks, orientation='horizontal')
				else:
					cbar = fig.colorbar(cax, ticks = ticks)
				cbar.set_label('Intensity', rotation= 90)
		
			ax.xaxis_date()  # x axis has a date data
	
			ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=range(60), interval=xtick, tz=None))
			ax.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))

			#set ytick labels to be intigers only, in case we have a small frequency range, pyplot tends to show fraction 
			ax.get_yaxis().get_major_formatter().set_useOffset(False)
	
			if endpts[0]:
				#get the start time and end time ans set ticks at those position son x-axis using minor_locator
				total_sec = utils.tosec(endtime - starttime)
				ax.xaxis.set_minor_locator(mdates.SecondLocator(bysecond=None, interval=total_sec, tz=None))
				ax.xaxis.set_minor_formatter(DateFormatter('%H:%M:%S'))
			fig.autofmt_xdate()
	
			if endpts[1]:
				#set y minorticks for first and last frequency
				plt.yticks(list(plt.yticks()[0]) + y_lims)
				
			#ax.xaxis.set_minor_formatter(DateFormatter('%H:%M:%S'))
			plt.xlabel('Universal Time')
			plt.ylabel('Frequency [MHz]')

	
		#title = fitsfile + "\n" + imageheader['DATE-OBS'] + "::" + imageheader['TIME-OBS']
		title = self.imageheader['CONTENT']
		plt.title(title)	
		#plt.show()
		return plt



	def joinTimeaxis(self, fits2):
		"""Take second radiohelliograph observation fits file and join two in time axis (x axis) and return new pyCallisto object.

		Args:
			fits2 (string or pyCallisto object): Second input fits file
			
		Returns:
			Returns pyCallisto object
		"""


		#get imagedata of first fits file
		imagedata1 = self.imagehdu.data



		#get the image data of second fits file 
		#(which is the argument here, could be a string/path of fits file or pyCallisto object)

		#if argument is a string (representing path of fits file)
		if isinstance(fits2, str):
			#check if it is proper callisto file or not
			if not (utils.checkFitsCallisto(fits2)):
				raise Exception("Fits file is not proper callisto files")
				#print("Fits file is not proper callisto files")
				#return -1
			#get the data
			hdus = pyfits.open(fits2)
			imagehdu2 = hdus[0]
			bintablehdu2 = hdus[1]
			imageheader2 = imagehdu2.header
			binheader2 = bintablehdu2.header
			imagedata2 = imagehdu2.data

		else:#if argument is a pyCallisto object
			
			imagehdu2 = fits2.hdus[0]
			bintablehdu2 = fits2.hdus[1]
			imageheader2 = imagehdu2.header
			binheader2 = bintablehdu2.header
			imagedata2 = imagehdu2.data


		# check if both imagedats they have same y dimensions (frequency in out case)
		if not imagedata1.shape[0] == imagedata2.shape[0]:
			raise Exception("Frequency dimensions do not match, cannot concatinate files")
			#print("Frequency dimensions do not match, cannot concatinate files")
			#hdus.close()
			#return -1		

		# check the order of fits files in time axis, if not correct, flip it
		
		#get start time of first file
		startdate = utils.toDate(self.imageheader['DATE-OBS'])
		starttime = utils.toTime(self.imageheader['TIME-OBS'])
		starttime1 = dt.datetime.combine(startdate, starttime)  # make a datetime object

		#get start time of second
		startdate = utils.toDate(imageheader2['DATE-OBS'])
		starttime = utils.toTime(imageheader2['TIME-OBS'])
		starttime2 = dt.datetime.combine(startdate, starttime)  # make a datetime object

		#if input files in not proper  order in time swap them
		if not starttime1 < starttime2:
			imagehdu2, bintablehdu2, imageheader2, binheader2, imagedata2, self.imagehdu, self.bintablehdu, self.imageheader, self.binheader, imagedata1 = self.imagehdu, self.bintablehdu, self.imageheader, self.binheader, imagedata1, imagehdu2, bintablehdu2, imageheader2, binheader2, imagedata2 
			starttime1, starttime2 = starttime2, starttime1 
		
		# endtime of first file
		enddate = utils.toDate(self.imageheader['DATE-END'])
		endtime = utils.toTime(self.imageheader['TIME-END'])
		endtime1 = dt.datetime.combine(enddate, endtime)  # make a datetime object
	
		# check that two fits files are continuous in time axis
		td_sec = dt.timedelta(seconds=1)
		if starttime2 - endtime1 > td_sec:
			raise Exception("Fits  files are not continuous in time axis, cannot join them")
			#print("Fits  files are not continuous in time axis, cannot join them")
			#hdus.close()
			#return -1		
			#sys.exit(1)
	
		# check if both the files have same sampling /
		if not self.imageheader['CDELT1'] == imageheader2['CDELT1']:
			raise Exception("Two fits files do not have the same sampling in time axis, cannot join them")
			#print("Two fits files do not have the same sampling in time axis, cannot join them")
			#hdus.close()
			#return -1		
			#sys.exit(1)
	
		# join the numpy Array 
		imagedata = np.concatenate((imagedata1, imagedata2), axis=1)
	
		# compose a new image header
		imageheader = copy.deepcopy(self.imageheader)
		imageheader['NAXIS1'] = imagedata.shape[1]
		imageheader['DATE-END'] = imageheader2['DATE-END']
		imageheader['TIME-END'] = imageheader2['TIME-END']
		imageheader['DATAMIN'] = imagedata.min()
		imageheader['DATAMAX'] = imagedata.max()
		imageheader['COMMENT'] = "created on " + str(dt.datetime.now()) + " by joining fits files " # + str(self.infits) + " " + str(fits2)

		# create a primary hdu containing this image and header data
		imagehdu = pyfits.PrimaryHDU(imagedata, header=imageheader)
		newhdulist = pyfits.HDUList([imagehdu]) 
	
		# bintabledata
		xlength = int(self.imageheader['NAXIS1']) + int(imageheader2['NAXIS1']) 
		xlength = int(xlength)
		 
		rangelist = [x * 1.0 for x in range(xlength)]
		bintabledatatime = np.array([ rangelist ])
		bintabledatafreqs = list(self.bintablehdu.data[0][1].copy())
		bintabledatafreqs = np.array([bintabledatafreqs])
	
		# these two arrays needs to be two dimensional, though the content is one dimensional,
		# to be in sync with original data format, which is 1RX2C 
		#print bintabledatatime.shape
		#print bintabledatafreqs.shape
	
		# generalize the "format" 
		format1 = str(bintabledatatime.shape[1]) + "D8.3"
		format2 = str(bintabledatafreqs.shape[1]) + "D8.3"
		col1 = pyfits.Column(name='TIME', format=format1, array=bintabledatatime)
		col2 = pyfits.Column(name='FREQUENCY', format=format2, array=bintabledatafreqs)
		cols = pyfits.ColDefs(np.asarray([col1, col2]))
	
		# create bintable and  add it to hdulist
		tbhdu = pyfits.BinTableHDU.from_columns(cols)
		newhdulist.append(tbhdu)

		#return new pycallisto object
		return pyCallisto(newhdulist)





	def sliceTimeaxis(self, time1, time2):
		"""Make a slice of input radiohelliograph observation fits file along a time axis and return a new object
	
		Args:
			time1 (string): start of a slice
							time in HH:MM:SSformat 
			time1 (string): end of a slice
							time in HH:MM:SS format
		Returns
			pyCallisto object
		"""

		# assuming that the input file has same start date and end date  (i.e. observed on the same day)
		# time1 and time 2 is in HH:MM:SS
		if not (len(time1.split(':')) == 3):
			raise Exception("Time format not proper, please provide time in HH:MM:SS format")
			#print "Time format not proper, please provide time in HH:MM:SS format"
			#return -1		
			#sys.exit(1)
	
		if not (len(time2.split(':')) == 3):
			raise Exception("Time format not proper, please provide time in HH:MM:SS format")
			#print "Time format not proper, please provide time in HH:MM:SS format"
			#return -1
			#sys.exit(1)
	
#		# open input fits file
#		hdus = pyfits.open(infits)
#		imagehdu = hdus[0]
#		bintablehdu = hdus[1]
#		imageheader = imagehdu.header
#		imagedata = imagehdu.data
	
		# check the dates
		startdate = utils.toDate(self.imageheader['DATE-OBS'])
		endtdate = utils.toDate(self.imageheader['DATE-END'])
		if not startdate == endtdate:
			raise Exception("Startdate and enddate differ, right now we do not support this")
			#print "Startdate and enddate differ, right now we do not support this"
			#hdus.close()
			#return -1
			#sys.exit(1)
		
		# get the times in datetime.date format for easy manipulation
		starttime = utils.toTime(self.imageheader['TIME-OBS'])  # datetime.date object
		starttime = dt.datetime.combine(startdate, starttime)  # datetime.datetime object
		endtime = utils.toTime(self.imageheader['TIME-END'])
		endtime = dt.datetime.combine(startdate, endtime)
	
		time1 = utils.toTime(time1)
		time1 = dt.datetime.combine(startdate, time1)
		time2 = utils.toTime(time2)
		time2 = dt.datetime.combine(startdate, time2)
	
		# check the "time constraints"
		if not (time1 < time2):
			time1, time2 = time2, time1
		if not (starttime < time1):
			#print("Time1 out of bound, can't slice!")
			print("Start time of input file : ", starttime)
			print("End time of input file : ", endtime)
			raise Exception("Time1 out of bound, can't slice!")
			#hdus.close()
			#return -1
			#sys.exit(1)
		if not (endtime > time2):
			#print("Time2 out of bound, can't slice!")
			print("Start time of input file : ", starttime)
			print("End time of input file : ", endtime)
			raise Exception("Time2 out of bound, can't slice!")
			#hdus.close()
			#return -1
			#sys.exit(1)
		
		# get the startpixel and endpixel in piselterms from input of time1 and time2 
		startpixel = time1 - starttime
		startpixel = utils.tosec(startpixel)  # in seconds
		startoffset = startpixel
		startpixel = int(startpixel / float(self.imageheader['CDELT1']))
	
		endpixel = time2 - starttime
		endpixel = utils.tosec(endpixel)
		endpixel = int(endpixel / float(self.imageheader['CDELT1']))
	
		# do the actual slicing
		imagedata = copy.deepcopy(self.imagehdu.data)
		imagedata = imagedata[:, startpixel:endpixel]
	
		# update the imageheader accordingly
		imageheader = copy.deepcopy(self.imageheader)
		imageheader['NAXIS1'] = imagedata.shape[1]
		imageheader['TIME-OBS'] = str(time1.time())
		imageheader['TIME-END'] = str(time2.time())
		imageheader['DATAMIN'] = imagedata.min()
		imageheader['DATAMAX'] = imagedata.max()
		imageheader['CRVAL1'] = int(self.imageheader['CRVAL1']) + startoffset 
	
		# create a primary hdu containing this image and header data
		imagehdu = pyfits.PrimaryHDU(imagedata, header=imageheader)
		newhdulist = pyfits.HDUList([imagehdu])
	
		# create a new bintable and update it in data structure
		xlength = endpixel - startpixel 
		rangelist = [x * 1.0 for x in range(xlength)]
		bintabledatatime = np.array([ rangelist ])
		bintabledatafreqs = list(self.bintablehdu.data[0][1].copy())
		bintabledatafreqs = np.array([bintabledatafreqs]) 
	
		format1 = str(bintabledatatime.shape[1]) + "D8.3"
		format2 = str(bintabledatafreqs.shape[1]) + "D8.3" 
		col1 = pyfits.Column(name='TIME', format=format1, array=bintabledatatime)
		col2 = pyfits.Column(name='FREQUENCY', format=format2, array=bintabledatafreqs)
		cols = pyfits.ColDefs(np.asarray([col1, col2]))
		tbhdu = pyfits.BinTableHDU.from_columns(cols)
	
		# add it to new hdulist we have created
		newhdulist.append(tbhdu)

		return pyCallisto(newhdulist)




	def sliceFrequencyaxis(self, freq1, freq2):
		"""Make a slice of input radiohelliograph observation fits file along a frequency axis
	
		Args:
			freq1 (int): start of a slice
			freq2 (int): end of a slice

		Returns:
			new pyCallisto object
		"""

		# convert  the input frequency from string to int
		freq1 = int(freq1)
		freq2 = int(freq2)
		#swap frequencies if not in proper orde
		if freq1 > freq2:
			freq1, freq2 = freq2, freq1


	#	# open input fits file
		hdus = self.hdus
		imagehdu = hdus[0]
		bintablehdu = hdus[1]
		imageheader = imagehdu.header
		imagedata = imagehdu.data
		bintblfreqdata = bintablehdu.data[0][1]

		# check the frequencies
		startfreq = int(bintblfreqdata[-1])
		endfreq = int(bintblfreqdata[0])

	
		if (freq1 < startfreq or freq1 > endfreq):
			#print("Frequency out of bound, cannot slice")
			print("Start Frequency of input file : ", startfreq)
			print("End Frequency of input file ", endfreq)
			raise Exception("Frequency out of bound, cannot slice")
			#hdus.close()
			#return -1
			#raise SystemExit, 0
			#sys.exit()
		
		if freq2 < startfreq or freq2 > endfreq :
			#print "Frequency out of bound, cannot slice"
			print("Start Frequency of input file : ", startfreq)
			print("End Frequency of input file ", endfreq)
			raise Exception("Frequency out of bound, cannot slice")
			#hdus.close()
			#return -1
			#raise SystemExit, 0
			#sys.exit()
	
		if  (freq2 - freq1 < 1):
			#print "Too thin slice demanded, cannot slice thinner than 1 unit"
			print("Start Frequency of input file : ", startfreq)
			print("End Frequency of input file ", endfreq)
			raise Exception("Too thin slice demanded, cannot slice thinner than 1 unit")
			#hdus.close()
			#return -1
			#raise SystemExit, 0
			#sys.exit()
	
		# update the bintabledata according to freq1 and freq2
		try:
			slicept1 = np.argwhere((bintblfreqdata > freq1) & (bintblfreqdata < freq2))[0][0]
			slicept2 = np.argwhere((bintblfreqdata > freq1) & (bintblfreqdata < freq2))[-1][-1]
	
		except:
			#print("frequency limits given are smaller than single channel, please increase it and retry")
			print("Start Frequency of input file : ", startfreq)
			print("End Frequency of input file : ", endfreq)
			raise Exception("frequency limits given are smaller than single channel, please increase it and retry")
			#hdus.close()
			#return -1	
		
	
		if (slicept2 - slicept1) < 1:
			#print "frequency limits given are smaller than single channel, please increase it and retry"
			print("Start Frequency of input file : ", startfreq)
			print("End Frequency of input file : ", endfreq)
			raise Exception("frequency limits given are smaller than single channel, please increase it and retry")
			#hdus.close()
			#return -1
			#raise SystemExit, 0
			#sys.exit(1)
	
		bintabledatafreqs = bintblfreqdata[np.argwhere((bintblfreqdata > freq1) & (bintblfreqdata < freq2))]
		bintabledatafreqs = np.array([bintabledatafreqs])

		xlength = int(imagedata.shape[1])
		rangelist = [x * 1.0 for x in range(xlength)]
		bintabledatatime = np.array([ rangelist ])
	
		format1 = str(bintabledatatime.shape[1]) + "D8.3"
		format2 = str(bintabledatafreqs.shape[1]) + "D8.3" 
		col1 = pyfits.Column(name='TIME', format=format1, array=bintabledatatime)
		col2 = pyfits.Column(name='FREQUENCY', format=format2, array=bintabledatafreqs)
		cols = pyfits.ColDefs(np.asarray([col1, col2]))
		tbhdu = pyfits.BinTableHDU.from_columns(cols)
	
		# update the imageheader accordingly
		imageheader = copy.deepcopy(self.imageheader)
		imagedata = copy.deepcopy(imagedata)
		imagedata = imagedata[slicept1:slicept2, :]
		imageheader['NAXIS2'] = imagedata.shape[0]
		imageheader['DATAMIN'] = imagedata.min()
		imageheader['DATAMAX'] = imagedata.max()
		imageheader['CRVAL2'] = imagedata.shape[0]

		# create a primary hdu containing this image and header data
		imagehdu = pyfits.PrimaryHDU(imagedata, header=imageheader)
		newhdulist = pyfits.HDUList([imagehdu])
	 
		# add it to hdulist, and write hdulist  to  fits file
		newhdulist.append(tbhdu)
	
		return pyCallisto(newhdulist)

	def subtractBackground(self):
		"""estimate and subtract background from a fits file
	
		Args:
			input " No input
			output:
				returns a pyCallisto object
		"""
	
		# open input fits file
		hdus = self.hdus
		imagehdu = hdus[0]
		bintablehdu = hdus[1]
		imageheader = imagehdu.header
		imagedata = imagehdu.data
	
		print("size of the bintblfreqdata", type(bintablehdu.data), len(bintablehdu.data))
		bintblfreqdata = bintablehdu.data[0][1]
	
	
	
		#estimate background and subtract from original data
		med_col = np.median(imagedata, axis=1, keepdims=1)
		imagedata = copy.deepcopy(imagedata)
		imagedata = imagedata - med_col


		#we need to update the header accordingly  ##############
		# compose a new image header
		imageheader = copy.deepcopy(imageheader)
		imageheader['DATAMIN'] = imagedata.min()
		imageheader['DATAMAX'] = imagedata.max()
	
	
		# create a primary hdu containing this updated image data
		imagehdu = pyfits.PrimaryHDU(imagedata, header=imageheader)
		newhdulist = pyfits.HDUList([imagehdu])
	 
		# add the unchanged bintablehdu to hdulist
		newhdulist.append(bintablehdu)

		return pyCallisto(newhdulist)






	def getTimeSeries(self, plot=True, outpng ="timeseries.png", returndata=False):
		"""Collapse the 2D fits image along time axis and return 1D array
	
		Args:
			plot (Boolean) : plot it or not, Default True
			returndata (Boolean) :   default False
	
		returns:if returndata is set to True 
				list of collapsed 1d array (1d numpy array), 
				 respective time in sec of a day(1d numpy array),
				 respective time in datetime.datetime object(1d list)			 
		"""
		# open input fits file
		hdus = self.hdus
		imagehdu = hdus[0]
		bintablehdu = hdus[1]
		imageheader = imagehdu.header
		imagedata = imagehdu.data
		bintblfreqdata = bintablehdu.data[0][1]
	
		# Sum the data along oth axis
		sumimage = np.sum(imagedata, axis=0)
	
		start_date = utils.toDate(imageheader['DATE-OBS'])
		starttime = utils.toTime(imageheader['TIME-OBS'])
		starttime = dt.datetime.combine(start_date, starttime)  # make a datetime object
		endtime = utils.toTime(imageheader['TIME-END'])
		endtime = dt.datetime.combine(start_date, endtime)

		def gettimeaxis(start, end, delta):
			curr = start
			while curr < end:
				yield curr
				curr += delta
	
		timeaxis = [time for time in gettimeaxis(starttime, endtime, (endtime - starttime) / sumimage.shape[0])]
	
		if plot:
			plt.clf()
			fig, ax = plt.subplots()
			ax.plot_date(timeaxis, sumimage, 'b-', xdate=True)
			ax.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
			plt.xlabel('Time')
			plt.ylabel('Total count')
			#title = infits + "\n" + imageheader['DATE-OBS'] + "::" + imageheader['TIME-OBS']
			title = self.imageheader['CONTENT']
			plt.title(title)
			plt.grid()
			#plt.show()
			plt.savefig(outpng)
		
		if returndata:	
			# generate numpy array of sec. of a day
			xstart = int(imageheader['CRVAL1'])
			xstep = float(imageheader['CDELT1'])  
			xlength = int(imageheader['NAXIS1'])
			xend = xstart + (xstep * xlength)
			timesec = np.arange(xstart, xend, xstep)
		
			# compose data into a list
			data = [sumimage, timesec, timeaxis]
			#hdus.close()
			return data


def showSpectrum(self, filepath, filetype):


	def getFrequencySeries(self, plot=True, outpng ="frequency_series.png", returndata=False):
		"""Collapse the 2D fits image along frequency axis and return 1D array
	
		Args:
			plot (Boolean): plot or not ?, Default True
			outpng (String):  name of the png file to plot
			returndata (Boolean) : return data or not, Default False
		
		returns : list of  collapsed 1d array (1d numpy array),
				respective frequency channels (1d numpy array)
		"""
	
		# open input fits file
		hdus = self.hdus
		imagehdu = hdus[0]
		bintablehdu = hdus[1]
		imageheader = imagehdu.header
		imagedata = imagehdu.data
		bintblfreqdata = bintablehdu.data[0][1]
	
		# Sum the data along oth axis
		sumimage = np.sum(imagedata, axis=1)
		 
		if plot:
			# fig, ax = plt.subplots()
			plt.clf()
			plt.plot(bintblfreqdata[:-1], sumimage, 'b-',)	#check why we need to add -1 here
			plt.axis([bintblfreqdata[-1], bintblfreqdata[0], sumimage.min(), sumimage.max()])
			plt.xlabel(imageheader['CTYPE2'])
			plt.ylabel('Total count')
			#title = infits + "\n" + imageheader['DATE-OBS'] + "::" + imageheader['TIME-OBS'] + " To " + imageheader['DATE-END'] + "::" + imageheader['TIME-END']
			title = self.imageheader['CONTENT']
			plt.title(title)
			plt.grid()
			plt.savefig(outpng)
			#hdus.close()
		
		if returndata: 
			# compose data into a list
			data = [sumimage, bintblfreqdata]
			hdus.close()
			return data

