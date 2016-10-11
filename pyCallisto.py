'''
This module is a collection of utility functions for radio Hellio spectrometer data like callisto.

Functions avaialble here can plot the data, join or slice the data along timeaxis as well as frequency axis.And all data can be summed along axis to give a "timesries" or "frequencyseries".
Input to functions is a single or multiple fits file/s, output can be a fits file in joining and slicing operation.In case of plotting an image is saved at given path. And in case of "timeseries" and "frequencyseries"  you  can ask for plot or you can return the data and do different task with it like plotting, fitting a function etc. 
'''

import datetime as dt
import pyfits
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import cm
import os
import numpy as np
from matplotlib.dates import  DateFormatter
import sys
import matplotlib
#from matplotlib.ticker import MaxNLocator

def checkFitsCallisto(fitsfile):
    """Check whether fits file has two HDUs or not
    """
    hdus = pyfits.open(fitsfile)
    if len(hdus) == 2:
        return True
    else:
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
    return dt.time(hr, mn, sec)


def plot(fitsfile, imagepath, option=3, xtick= 5, blevel=0, endpts= [False, False], cmap=cm.jet, cbar=True, cbar_ori='vertical', fontsize=12):
    """Plot a radiohelliograph observation file and save at given path
    
    Args:
        fitsfile (string): path of a fitsfile to be plotted
        imagepath (string): path to store imagefile 
        option (int): option 1/2/3, default 3
                options :   1 raw : plot column no. and row no. on x and y axis 
                			2 time and frequncy : time in sec like 1-900, and frequency
                			3 absolute time and frequency 
        xtick = x axis tick interval, default is every 5 minutes
        blevel (int): background level of your data, default 0
        endpts = show end points on x and y axis, deafuslt false
        cmap (matplotlib colormap)= default cm.jet, 
        cbar (Boolean) : Plot colorbar or not ? default True, 
        cbar_ori (string )= colorbar orientation, 'vertical' / 'horizontal' , default 'vertical'  
        fontsize (int) = fontsize for the labels and tttles etc
    """
	#set the default fontsize
    matplotlib.rcParams.update({'font.size': fontsize})
	
    # checking number of HDUs
    if not checkFitsCallisto(fitsfile):
        print "No. of HDUs are wrong, not a proper callisto file "
        return -1
        #sys.exit(1)
  
    #open input fits image 
    hdus = pyfits.open(fitsfile)
    imagehdu = hdus[0]      
    bintablehdu = hdus[1]
    imageheader = imagehdu.header
    
    # get the datamin amd maximum for colorbar
    datamin = int(imagehdu.data.min())
    datamax = int(imagehdu.data.max())
    datamid = (datamin + datamax) / 2
    
    #plot the fig 
    fig, ax = plt.subplots(figsize=(20,10))
    if option == 1:
        y, x = imagehdu.data.shape
        cax = ax.imshow(imagehdu.data, extent=[0, x, 0, y], aspect='auto', cmap=cmap, vmin=blevel)
        if cbar == True:
            if cbar_ori == 'horizontal':
                cbar = fig.colorbar(cax, ticks=[blevel, (datamax+blevel)/2, datamax], orientation='horizontal')
            else:
                cbar = fig.colorbar(cax, ticks=[blevel, (datamax+blevel)/2, datamax])
        
        
        plt.xlabel('Row Count')
        plt.ylabel('Column Count')
             
    if option == 2:
        xstart = int(imageheader['CRVAL1'])
        xstep = float(imageheader['CDELT1'])  
        xlength = int(imageheader['NAXIS1'])
        xend = xstart + xstep * xlength
        freqs = bintablehdu.data['frequency'][0]  # these are true frequencies
        cax = ax.imshow(imagehdu.data, extent=[xstart, xend, freqs[-1], freqs[0]], aspect='auto', cmap=cmap, vmin=blevel)
        if cbar == True:
            if cbar_ori == 'horizontal':
                cbar = fig.colorbar(cax, ticks=[blevel, (datamax+blevel)/2, datamax], orientation='horizontal')
            else:
                cbar = fig.colorbar(cax, ticks=[blevel, (datamax+blevel)/2, datamax])
        
        if endpts[1]:
            #set y minorticks for first and last frequency
            y_lims = [freqs[-1], freqs[0]]
            print y_lims
            plt.yticks(list(plt.yticks()[0]) + y_lims)
        
        plt.xlabel('Time (sec of day)')
        plt.ylabel('Frequency [MHz]')
         
    if option == 3:
        # get the start time and end time of observation
        start_date = toDate(imageheader['DATE-OBS'])
        starttime = toTime(imageheader['TIME-OBS'])
        
        
        starttime = dt.datetime.combine(start_date, starttime)  # make a datetime object
        endtime = toTime(imageheader['TIME-END'])
        endtime = dt.datetime.combine(start_date, endtime)
        
        # get the frequencies
        freqs = bintablehdu.data['frequency'][0]  # these are true frequencies
        
        # set the limits for plotting 
        x_lims = [starttime, endtime]
        x_lims = mdates.date2num(x_lims)  # dates to numeric values
        y_lims = [freqs[-1], freqs[0]]

		
        
        cax = ax.imshow(imagehdu.data, extent=[x_lims[0], x_lims[1], y_lims[0], y_lims[1]], aspect='auto', cmap=cmap, vmin=blevel)
        if cbar == True:
            if cbar_ori == 'horizontal':
                cbar = fig.colorbar(cax, ticks=[blevel, (datamax+blevel)/2, datamax], orientation='horizontal')
            else:
                cbar = fig.colorbar(cax, ticks=[blevel, (datamax+blevel)/2, datamax])
        ax.xaxis_date()  # x axis has a date data
        
        ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=range(60), interval=xtick, tz=None))
        ax.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))

        #set ytick labels to be intigers only, in case we have a small frequency range, pyplot tends to show fraction 
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        
        if endpts[0]:
		    #get the start time and end time ans set ticks at those position son x-axis using minor_locator
		    total_sec = tosec(endtime - starttime)
		    ax.xaxis.set_minor_locator(mdates.SecondLocator(bysecond=None, interval=total_sec, tz=None))
		    ax.xaxis.set_minor_formatter(DateFormatter('%H:%M:%S'))
        fig.autofmt_xdate()
        
        if endpts[1]:
            #set y minorticks for first and last frequency
            plt.yticks(list(plt.yticks()[0]) + y_lims)
        
        
        
        #ax.xaxis.set_minor_formatter(DateFormatter('%H:%M:%S'))
        
        plt.xlabel('Universal Time')
        plt.ylabel('Frequency [MHz]')
    
    title = fitsfile + "\n" + imageheader['DATE-OBS'] + "::" + imageheader['TIME-OBS']
    plt.title(title)    
    plt.savefig(imagepath)



def jointimeaxis(fits1, fits2, ofile):
    """Take two radiohelliograph observation files and join them in time axis (x axis)
    
    Args:
        fits1 (string): First input fits file  
        fits2 (string): Second input fits file
        ofile (string): path to store the output fits file 
    """
    
    # check fits files
    # open fits1, fits 2
    # get the data and headers, bintables and headers
    if not (checkFitsCallisto(fits1) and checkFitsCallisto(fits2)):
        print "Either one or both fits files are not proper callisto files"
        return -1        
        #sys.exit(1)
    if not ((ofile.split('.')[-1] == "fits") or (ofile.split('.')[-1] == "fit")):
        print "output file should have a fit/fits extension, please provide it"
        return -1        
        #sys.exit(1)
        
    hdus = pyfits.open(fits1)
    imagehdu1 = hdus[0]      
    bintablehdu1 = hdus[1]
    imageheader1 = imagehdu1.header
    binheader1 = bintablehdu1.header 
    imagedata1 = imagehdu1.data
    
    hdus = pyfits.open(fits2)
    imagehdu2 = hdus[0]      
    bintablehdu2 = hdus[1]
    imageheader2 = imagehdu2.header
    binheader2 = bintablehdu2.header
    imagedata2 = imagehdu2.data
    
    # check they have same y dimensions
    if not imagedata1.shape[0] == imagedata2.shape[0]:
        print "Frequency dimensions do not match, cannot concatinate files"
        return -1        
        #sys.exit(1)
         
    # check the order of fits files in time axis, if not correct flip it
    startdate = toDate(imageheader1['DATE-OBS'])
    starttime = toTime(imageheader1['TIME-OBS'])
    starttime1 = dt.datetime.combine(startdate, starttime)  # make a datetime object

    startdate = toDate(imageheader1['DATE-OBS'])
    starttime = toTime(imageheader2['TIME-OBS'])
    starttime2 = dt.datetime.combine(startdate, starttime)  # make a datetime object
    
    #if input files in not proper  order in time swap them
    if not starttime1 < starttime2:
        imagehdu2, bintablehdu2, imageheader2, binheader2, imagedata2, imagehdu1, bintablehdu1, imageheader1, binheader1, imagedata1 = imagehdu1, bintablehdu1, imageheader1, binheader1, imagedata1, imagehdu2, bintablehdu2, imageheader2, binheader2, imagedata2 
        starttime1, starttime2 = starttime2, starttime1 
    
    # endtime of first file
    enddate = toDate(imageheader1['DATE-END'])
    endtime = toTime(imageheader1['TIME-END'])
    endtime1 = dt.datetime.combine(enddate, endtime)  # make a datetime object
    
    # check that two fits files are continuous in time axis
    td_sec = dt.timedelta(seconds=1)
    if starttime2 - endtime1 > td_sec:
        print "Fits  files are not continuous in time axis, cannot join them"
        return -1        
        #sys.exit(1)
    
    # check if both the files have same sampling /
    if not imageheader1['CDELT1'] == imageheader2['CDELT1']:
        print "Two fits files do not have the same sampling in time axis, cannot join them" 
        return -1        
        #sys.exit(1)
    
    # join the numpy Array 
    imagedata = np.concatenate((imagedata1, imagedata2), axis=1)
    
    # compose a new image header
    imageheader = imageheader1.copy()
    imageheader['NAXIS1'] = imagedata.shape[1]
    imageheader['DATE-END'] = imageheader2['DATE-END']
    imageheader['TIME-END'] = imageheader2['TIME-END']
    imageheader['DATAMIN'] = imagedata.min()
    imageheader['DATAMAX'] = imagedata.max()
    imageheader['COMMENT'] = "created on " + str(dt.datetime.now()) + " by joining fits files " + str(fits1) + " " + str(fits2)
    
    # create a primary hdu containing this image and header data
    imagehdu = pyfits.PrimaryHDU(imagedata, header=imageheader)
    newhdulist = pyfits.HDUList([imagehdu]) 
    
    # bintabledata
    xlength = int(imageheader1['NAXIS1']) + int(imageheader2['NAXIS1']) 
    xlength = int(xlength)
     
    rangelist = [x * 1.0 for x in range(xlength)]
    bintabledatatime = np.array([ rangelist ])
    bintabledatafreqs = list(bintablehdu1.data[0][1].copy())
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
    
    # create bintable
    tbhdu = pyfits.BinTableHDU.from_columns(cols)
    # add it to hdulist
    newhdulist.append(tbhdu)
    # write hdulist  to  fits file
    if os.path.isfile(ofile):
        print "file %s already exists, cannot overwrite, please choose a different output file name"
        return -1
    newhdulist.writeto(ofile)



def slicetimeaxis(infits, time1, time2, ofile):
    """Make a slice of input radiohelliograph observation fits file along a time axis
    
    Args:
        infits (string): path of a fits file
        time1 (string): start of a slice
                        time in HH:MM:SSformat 
        time1 (string): end of a slice
                        time in HH:MM:SS format
        ofile (string): path of the output fits file
    """
    # check for the inputs
    if not os.path.isfile(infits):
        print "input fits file %s does not exists"
        return -1
        #sys.exit(1)
    
    # checking number of HDUs
    if not checkFitsCallisto(infits):
        print "No. of HDUs are wrong, not a proper callisto file "
        return -1
        #sys.exit(1)

    # assuming that the input file has same start date and end date 
    # time1 and time 2 is in HH:MM:SS
    if not (len(time1.split(':')) == 3):
        print "Time format not proper, please provide time in HH:MM:SS format"
        return -1        
        #sys.exit(1)
    
    if not (len(time2.split(':')) == 3):
        print "Time format not proper, please provide time in HH:MM:SS format"
        return -1
        #sys.exit(1)
    
    # open input fits file
    hdus = pyfits.open(infits)
    imagehdu = hdus[0]
    bintablehdu = hdus[1]
    imageheader = imagehdu.header
    imagedata = imagehdu.data
    
    # check the dates
    startdate = toDate(imageheader['DATE-OBS'])
    endtdate = toDate(imageheader['DATE-END'])
    if not startdate == endtdate:
        print "Startdate and enddate differ, right now we do not support this"
        return -1
        #sys.exit(1)
        
    # get the times in datetime.date format for easy manipulation
    starttime = toTime(imageheader['TIME-OBS'])  # datetime.date object
    starttime = dt.datetime.combine(startdate, starttime)  # datetime.datetime object
    endtime = toTime(imageheader['TIME-END'])
    endtime = dt.datetime.combine(startdate, endtime)
    
    time1 = toTime(time1)
    time1 = dt.datetime.combine(startdate, time1)
    time2 = toTime(time2)
    time2 = dt.datetime.combine(startdate, time2)
    
    # check the "time constraints"
    if not (time1 < time2):
        time1, time2 = time2, time1
    if not (starttime < time1):
        print "Time1 out of bound, can't slice!"
        print "Start time of input file : %s"%(starttime)
        print "End time of input file : %s"%(endtime)
        return -1
        #sys.exit(1)
    if not (endtime > time2):
        print "Time2 out of bound, can't slice!"
        print "Start time of input file : %s"%(starttime)
        print "End time of input file : %s"%(endtime)
        return -1
        #sys.exit(1)
        
    # get the startpixel and endpixel in piselterms from input of time1 and time2 
    startpixel = time1 - starttime
    startpixel = tosec(startpixel)  # in seconds
    startoffset = startpixel
    startpixel = int(startpixel / float(imageheader['CDELT1']))
    
    endpixel = time2 - starttime
    endpixel = tosec(endpixel)
    endpixel = int(endpixel / float(imageheader['CDELT1']))
    
    # do the actual slicing
    imagedata = imagedata[:, startpixel:endpixel]
    
    # update the imageheader accordingly
    imageheader['NAXIS1'] = imagedata.shape[1]
    imageheader['TIME-OBS'] = str(time1.time())
    imageheader['TIME-END'] = str(time2.time())
    imageheader['DATAMIN'] = imagedata.min()
    imageheader['DATAMAX'] = imagedata.max()
    imageheader['CRVAL1'] = int(imageheader['CRVAL1']) + startoffset 
    
    # create a primary hdu containing this image and header data
    imagehdu = pyfits.PrimaryHDU(imagedata, header=imageheader)
    newhdulist = pyfits.HDUList([imagehdu])
    
    # update the bintable accordingly
    xlength = endpixel - startpixel 
    rangelist = [x * 1.0 for x in range(xlength)]
    bintabledatatime = np.array([ rangelist ])
    bintabledatafreqs = list(bintablehdu.data[0][1].copy())
    bintabledatafreqs = np.array([bintabledatafreqs]) 
    
    format1 = str(bintabledatatime.shape[1]) + "D8.3"
    format2 = str(bintabledatafreqs.shape[1]) + "D8.3" 
    col1 = pyfits.Column(name='TIME', format=format1, array=bintabledatatime)
    col2 = pyfits.Column(name='FREQUENCY', format=format2, array=bintabledatafreqs)
    cols = pyfits.ColDefs(np.asarray([col1, col2]))
    tbhdu = pyfits.BinTableHDU.from_columns(cols)
    
    # add it to hdulist, a nd write hdulist  to  fits file
    newhdulist.append(tbhdu)
    if os.path.isfile(ofile):
        print "file %s already exists, cannot overwrite, please choose a different output file name"
        return -1
    newhdulist.writeto(ofile)


def slicefreqaxis(infits, freq1, freq2, ofile):
    """Make a slice of input radiohelliograph observation fits file along a frequency axis
    
    Args:
        infits (string): path of a fits file
        freq1 (int): start of a slice            
        time2 (int): end of a slice
        ofile (string): path of the output fits file
    """
    
    # check for the inputs
    if not os.path.isfile(infits):
        print "input file %s does not exists"
        return -1        
        #sys.exit(1)
    
    # checking number of HDUs
    if not checkFitsCallisto(infits):
        print "No. of HDUs are wrong, not a proper callisto file "
        return -1        
        #sys.exit(1)
    
    # convert  the input frequency from string to int
    freq1 = int(freq1)
    freq2 = int(freq2)
    #swap frequencies if not in proper orde
    if freq1 > freq2:
        freq1, freq2 = freq2, freq1
    
    # open input fits file
    hdus = pyfits.open(infits)
    imagehdu = hdus[0]
    bintablehdu = hdus[1]
    imageheader = imagehdu.header
    imagedata = imagehdu.data
    bintblfreqdata = bintablehdu.data[0][1]
        
    # check the frequencies
    startfreq = int(bintblfreqdata[-1])
    endfreq = int(bintblfreqdata[0])

    
    if (freq1 < startfreq or freq1 > endfreq):
        print "Frequency out of bound, cannot slice"
        print "Start Frequency of input file : %d"%(startfreq)
        print "End Frequency of input file : %d"%(endfreq)
        return -1
        #raise SystemExit, 0
        #sys.exit()
        
    if freq2 < startfreq or freq2 > endfreq :
        print "Frequency out of bound, cannot slice"
        print "Start Frequency of input file : %d"%(startfreq)
        print "End Frequency of input file : %d"%(endfreq)
        return -1
        #raise SystemExit, 0
        #sys.exit()
    
    if  (freq2 - freq1 < 1):
        print "Too thin slice demanded, cannot slice thinner than 1 unit"
        print "Start Frequency : %d"%(startfreq)
        print "End Frequency : %d"%(endfreq)
        return -1
        #raise SystemExit, 0
        #sys.exit()
    
    # update the bintabledata according to freq1 and freq2
    try:
        slicept1 = np.argwhere((bintblfreqdata > freq1) & (bintblfreqdata < freq2))[0][0]
        slicept2 = np.argwhere((bintblfreqdata > freq1) & (bintblfreqdata < freq2))[-1][-1]
    
    except:
        print "frequency limits given are smaller than single channel, please increase it and retry"
        print "Start Frequency : %d"%(startfreq)
        print "End Frequency : %d"%(endfreq)        
        return -1	
        
    
    if (slicept2 - slicept1) < 1:
        print "frequency limits given are smaller than single channel, please increase it and retry"
        print "Start Frequency : %d"%(startfreq)
        print "End Frequency : %d"%(endfreq)        
        return -1
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
    if os.path.isfile(ofile):
        print "file %s already exists, cannot overwrite, please choose a different output file name"
        return -1

    #print bintabledatafreqs.shape
    print "start frequency of o/p slice is "+str(bintabledatafreqs[0][-1][0])
    print "end frequency of o/p slice is "+str(bintabledatafreqs[0][0][0])
    newhdulist.writeto(ofile)



def gettimeseries(infits, plot=True, returndata=False):
    """Collapse the 2D fits image along time axis and return 1D array
    
    Args:
        infits (string) : input fits file
        plot (Boolean) : plot it or not, Default True
        returndata (Boolean) :   default False
    
    returns: list of collapsed 1d array (1d numpy array), 
             respective time in sec of a day(1d numpy array),
             respective time in datetime.datetime object(1d list)             
    """
    
    # check for the inputs
    if not os.path.isfile(infits):
        print "file %s does not exists"
        return -1
    
    # checking number of HDUs
    if not checkFitsCallisto(infits):
        print "No. of HDUs are wrong, not a proper callisto file "
        return -1
    
    # open input fits file
    hdus = pyfits.open(infits)
    imagehdu = hdus[0]
    bintablehdu = hdus[1]
    imageheader = imagehdu.header
    imagedata = imagehdu.data
    bintblfreqdata = bintablehdu.data[0][1]
    
    # Sum the data along oth axis
    sumimage = np.sum(imagedata, axis=0)
    
    start_date = toDate(imageheader['DATE-OBS'])
    starttime = toTime(imageheader['TIME-OBS'])
    starttime = dt.datetime.combine(start_date, starttime)  # make a datetime object
    endtime = toTime(imageheader['TIME-END'])
    endtime = dt.datetime.combine(start_date, endtime)

    def gettimeaxis(start, end, delta):
        curr = start
        while curr < end:
            yield curr
            curr += delta
    
    timeaxis = [time for time in gettimeaxis(starttime, endtime, (endtime - starttime) / sumimage.shape[0])]
    
    if plot:
        fig, ax = plt.subplots()
        ax.plot_date(timeaxis, sumimage, 'b-', xdate=True)
        ax.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
        plt.xlabel('Time')
        plt.ylabel('Total count')
        title = infits + "\n" + imageheader['DATE-OBS'] + "::" + imageheader['TIME-OBS']
        plt.title(title)
        plt.grid()
        #plt.show()
        plt.savefig(infits.split('.')[0] + "_timeseries.png")
        
    if returndata:    
        # generate numpy array of sec. of a day
        xstart = int(imageheader['CRVAL1'])
        xstep = float(imageheader['CDELT1'])  
        xlength = int(imageheader['NAXIS1'])
        xend = xstart + (xstep * xlength)
        timesec = np.arange(xstart, xend, xstep)
        
        # compose data into a list
        data = [sumimage, timesec, timeaxis]
        
        return data



def getfreqseries(infits, plot=True, returndata=False):
    """Collapse the 2D fits image along frequency axis and return 1D array
    
    Args:
        infits (string): input fits file
        plot (Boolean): plot or not ?, Default True
        returndata (Boolean) : return data or not, Default False
        
    returns : list of  collapsed 1d array (1d numpy array),
                        respective frequency channels (1d numpy array)
    
    """
    
    # check for the inputs
    if not os.path.isfile(infits):
        print "file %s does not exists"
        return -1
    
    # checking number of HDUs
    if not checkFitsCallisto(infits):
        print "No. of HDUs are wrong, not a proper callisto file "
        return -1
    
    # open input fits file
    hdus = pyfits.open(infits)
    imagehdu = hdus[0]
    bintablehdu = hdus[1]
    imageheader = imagehdu.header
    imagedata = imagehdu.data
    bintblfreqdata = bintablehdu.data[0][1]
    
    # Sum the data along oth axis
    sumimage = np.sum(imagedata, axis=1)
     
    if plot:
        # fig, ax = plt.subplots()
        plt.plot(bintblfreqdata, sumimage, 'b-',)
        plt.axis([bintblfreqdata[-1], bintblfreqdata[0], sumimage.min(), sumimage.max()])
        plt.xlabel(imageheader['CTYPE2'])
        plt.ylabel('Total count')
        title = infits + "\n" + imageheader['DATE-OBS'] + "::" + imageheader['TIME-OBS'] + " To " + imageheader['DATE-END'] + "::" + imageheader['TIME-END']
        plt.title(title)
        plt.grid()
        plt.savefig(infits.split('.')[0] + "_freqseries.png")
    
    if returndata: 
        # compose data into a list
        data = [sumimage, bintblfreqdata]
        return data
        
        
######################

def subkg(infits, outfits):
    """estimate and subtract background from a fits file
    
    Args:
        infits (string): input fits file
        outfits (string): output fits file
    """
    
    # check for the inputs
    if not os.path.isfile(infits):
        print "file %s does not exists"
        return -1
    
    # checking number of HDUs
    if not checkFitsCallisto(infits):
        print "No. of HDUs are wrong, not a proper callisto file "
        return -1
    
    # open input fits file
    hdus = pyfits.open(infits)
    imagehdu = hdus[0]
    bintablehdu = hdus[1]
    imageheader = imagehdu.header
    imagedata = imagehdu.data
    bintblfreqdata = bintablehdu.data[0][1]
    
    
    
    #estimate background and subtract from original data
    med_col = np.median(imagedata, axis=1, keepdims=1)
    imagedata = imagedata - med_col


    #we need to update the header accordingly  ##############
    # compose a new image header
    imageheader['DATAMIN'] = imagedata.min()
    imageheader['DATAMAX'] = imagedata.max()
    
    
    # create a primary hdu containing this updated image data
    imagehdu = pyfits.PrimaryHDU(imagedata, header=imageheader)
    newhdulist = pyfits.HDUList([imagehdu])
 
    # add the unchanged bintablehdu to hdulist
    newhdulist.append(bintablehdu)
    if os.path.isfile(outfits):
        print "file %s already exists, cannot overwrite, please choose a different output file name"
        return -1
    newhdulist.writeto(outfits)
    
