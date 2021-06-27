import datetime as dt
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import cm
import numpy as np
from matplotlib.dates import DateFormatter
import pycallisto_utils as utils  # local utility file
import copy
import sys

sys.path.append('.')


class PyCallisto:
    def __init__(self, hdu_list):
        """
        Create PyCallisto object from hdu_list
        """
        self.hdus = hdu_list
        self.imageHdu = self.hdus[0]
        self.binTableHdu = self.hdus[1]
        self.imageHeader = self.imageHdu.header
        # get the dataMin amd maximum for colorbar
        self.dataMin = int(self.imageHdu.data.min())
        self.dataMax = int(self.imageHdu.data.max())
        self.dataMid = (self.dataMin + self.dataMax) / 2

    def __del__(self):
        # print("Clearing the memory by deleting the HDUlist object not handled correctly by pyfits")
        self.hdus.close()

    @classmethod
    def from_file(cls, infits):
        """
        Create PyCallisto object from fits file on disc

        input arguments:
            infits: input fits file

        """

        # checking number of HDUs
        if not utils.check_fits_callisto(infits):
            # throw an error here
            raise ValueError("No. of HDUs are wrong in input fits file, not a proper callisto file , cannot proceed")

        # check bintable data accessible
        if not utils.check_bin_table(infits):
            # throw an error here
            raise ValueError("Bintable data may be corrupted in input fits file, cannot proceed")
        hdu_list = pyfits.open(infits)
        return cls(hdu_list)

    def spectrogram(self, option=3, xtick=2, blevel=0, end_pts=[False, False], fig_size=(8, 6), cmap=cm.jet,
                    color_bar=True, color_bar_ori='vertical', font_size=14):

        """
        Return a matplotlib.pyplot.plot object representing the plot of the input fits file.
        input arguments:
            xticks: frequency of xticks in mins (default =2 )
            blevel: background level (default = 0)
            figsize: a tuple representing size of image (default = (8,6))
            cmap: a matplotlib colormap object (default = cm.jet)
            cbar: to plot a colorbar or not (default = true)
            color_bar_ori: colorbar orientation (default = 'vertical')
            font_size: font_size used in plot (default = 14)

        """

        fig, ax = plt.subplots(figsize=(fig_size[0], fig_size[1]))

        if option == 1:
            y, x = self.imageHdu.data.shape
            cax = ax.imshow(self.imageHdu.data, extent=[0, x, 0, y], aspect='auto', cmap=cmap, vmin=blevel)
            if color_bar:
                ticks = list(np.linspace(blevel, self.dataMax, 10).astype('int'))  # calculate 10 ticks positions
                if color_bar_ori == 'horizontal':
                    color_bar = fig.colorbar(cax, ticks=ticks, orientation='horizontal')
                else:
                    color_bar = fig.colorbar(cax, ticks=ticks)
                color_bar.set_label('Intensity', rotation=90)

            plt.xlabel('Row Count')
            plt.ylabel('Column Count')

        if option == 2:
            x_start = int(self.imageHeader['CRVAL1'])
            x_step = float(self.imageHeader['CDELT1'])
            x_length = int(self.imageHeader['NAXIS1'])
            x_end = x_start + x_step * x_length
            freqs = self.binTableHdu.data['frequency'][0]  # these are true frequencies
            cax = ax.imshow(self.imageHdu.data, extent=[x_start, x_end, freqs[-1], freqs[0]], aspect='auto', cmap=cmap,
                            vmin=blevel)

            if color_bar:
                ticks = list(np.linspace(blevel, self.dataMax, 10).astype('int'))  # calculate 10 ticks positins
                if color_bar_ori == 'horizontal':
                    color_bar = fig.colorbar(cax, ticks=ticks, orientation='horizontal')
                else:
                    color_bar = fig.colorbar(cax, ticks=ticks)
                color_bar.set_label('Intensity', rotation=90)

            if end_pts[1]:
                # set y minorticks for first and last frequency
                y_lims = [freqs[-1], freqs[0]]
                plt.yticks(list(plt.yticks()[0]) + y_lims)

            plt.xlabel('Time (sec of day)')
            plt.ylabel('Frequency (MHz)')

        if option == 3:
            # get the start time and end time of observation
            start_date = utils.to_date(self.imageHeader['DATE-OBS'])
            start_time = utils.to_time(self.imageHeader['TIME-OBS'])

            start_time = dt.datetime.combine(start_date, start_time)  # make a datetime object
            end_time = utils.to_time(self.imageHeader['TIME-END'])
            end_time = dt.datetime.combine(start_date, end_time)

            # get the frequencies
            freqs = self.binTableHdu.data['frequency'][0]  # these are true frequencies

            # set the limits for plotting
            x_lims = [start_time, end_time]
            x_lims = mdates.date2num(x_lims)  # dates to numeric values
            y_lims = [freqs[-1], freqs[0]]

            cax = ax.imshow(self.imageHdu.data, extent=[x_lims[0], x_lims[1], y_lims[0], y_lims[1]], aspect='auto',
                            cmap=cmap, vmin=blevel)
            if color_bar:
                ticks = list(np.linspace(blevel, self.dataMax, 10).astype('int'))  # calculate 10 ticks positions
                if color_bar_ori == 'horizontal':
                    color_bar = fig.colorbar(cax, ticks=ticks, orientation='horizontal')
                else:
                    color_bar = fig.colorbar(cax, ticks=ticks)
                color_bar.set_label('Intensity', rotation=90)

            ax.xaxis_date()  # x axis has a date data

            ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=range(60), interval=xtick, tz=None))
            ax.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))

            # set ytick labels to be integers only, in case we have a small frequency range, pyplot tends to show
            # fraction
            ax.get_yaxis().get_major_formatter().set_useOffset(False)

            if end_pts[0]:
                # get the start time and end time ans set ticks at those positions on x-axis using minor_locator
                total_sec = utils.tosec(end_time - start_time)
                ax.xaxis.set_minor_locator(mdates.SecondLocator(bysecond=None, interval=total_sec, tz=None))
                ax.xaxis.set_minor_formatter(DateFormatter('%H:%M:%S'))
            fig.autofmt_xdate()

            if end_pts[1]:
                # set y minorticks for first and last frequency
                plt.yticks(list(plt.yticks()[0]) + y_lims)

            # ax.xaxis.set_minor_formatter(DateFormatter('%H:%M:%S'))
            plt.xlabel('Universal Time')
            plt.ylabel('Frequency (MHz)')

        title = self.imageHeader['CONTENT']
        plt.title(title)
        return plt

    def append_time_axis(self, fits2):
        """
        Take second radiohelliograph observation fits file and join two in time axis (x axis) and return new
        PyCallisto object.


        input arguments:
            fits2: path of a Second input fits file or PyCallisto object

        """

        # get image_data of first fits file
        image_data1 = self.imageHdu.data

        # get the image data of second fits file
        # (which is the argument here, could be a string/path of fits file or PyCallisto object)

        # if argument is a string (representing path of fits file)
        if isinstance(fits2, str):
            # check if it is proper callisto file or not
            if not (utils.check_fits_callisto(fits2)):
                raise Exception("Fits file is not proper callisto files")

            # get the data
            hdus = pyfits.open(fits2)
            image_hdu2 = hdus[0]
            bin_table_hdu2 = hdus[1]
            image_header2 = image_hdu2.header
            bin_header2 = bin_table_hdu2.header
            image_data2 = image_hdu2.data

        else:  # if argument is a PyCallisto object

            image_hdu2 = fits2.hdus[0]
            bin_table_hdu2 = fits2.hdus[1]
            image_header2 = image_hdu2.header
            bin_header2 = bin_table_hdu2.header
            image_data2 = image_hdu2.data

        # check if both imagedats they have same y dimensions (frequency in out case)
        if not image_data1.shape[0] == image_data2.shape[0]:
            raise Exception("Frequency dimensions do not match, cannot concatinate files")

        # check the order of fits files in time axis, if not correct, flip it

        # get start time of first file
        start_date = utils.to_date(self.imageHeader['DATE-OBS'])
        start_time = utils.to_time(self.imageHeader['TIME-OBS'])
        start_time1 = dt.datetime.combine(start_date, start_time)  # make a datetime object

        # get start time of second
        start_date = utils.to_date(image_header2['DATE-OBS'])
        start_time = utils.to_time(image_header2['TIME-OBS'])
        start_time2 = dt.datetime.combine(start_date, start_time)  # make a datetime object

        # if input files in not proper  order in time swap them
        if not start_time1 < start_time2:
            image_hdu2, bin_table_hdu2, image_header2, bin_header2, image_data2, self.imageHdu,  = \
                self.imageHdu, self.binTableHdu, self.imageHeader, self.binheader, image_data1, image_hdu2,

            self.binTableHdu, self.imageHeader, self.binheader, image_data1 = \
                image_header2, bin_header2, image_data2, bin_table_hdu2

            start_time1, start_time2 = start_time2, start_time1

        # end_time of first file
        end_date = utils.to_date(self.imageHeader['DATE-END'])
        end_time = utils.to_time(self.imageHeader['TIME-END'])
        end_time1 = dt.datetime.combine(end_date, end_time)  # make a datetime object

        # check that two fits files are continuous in time axis
        td_sec = dt.timedelta(seconds=1)
        if start_time2 - end_time1 > td_sec:
            raise Exception("Fits  files are not continuous in time axis, cannot join them")

        # check if both the files have same sampling /
        if not self.imageHeader['CDELT1'] == image_header2['CDELT1']:
            raise Exception("Two fits files do not have the same sampling in time axis, cannot join them")

        # join the numpy Array
        image_data = np.concatenate((image_data1, image_data2), axis=1)

        # compose a new image header
        image_header = copy.deepcopy(self.imageHeader)
        image_header['NAXIS1'] = image_data.shape[1]
        image_header['DATE-END'] = image_header2['DATE-END']
        image_header['TIME-END'] = image_header2['TIME-END']
        image_header['DATAMIN'] = image_data.min()
        image_header['DATAMAX'] = image_data.max()
        image_header['COMMENT'] = "created on " + str(
            dt.datetime.now()) + " by joining fits files "  # + str(self.infits) + " " + str(fits2)

        # create a primary hdu containing this image and header data
        image_hdu = pyfits.PrimaryHDU(image_data, header=image_header)
        new_hdulist = pyfits.HDUList([image_hdu])

        # bintabledata
        x_length = int(self.imageHeader['NAXIS1']) + int(image_header2['NAXIS1'])
        x_length = int(x_length)

        range_list = [x * 1.0 for x in range(x_length)]
        bin_table_data_time = np.array([range_list])
        bintable_data_freqs = list(self.binTableHdu.data[0][1].copy())
        bintable_data_freqs = np.array([bintable_data_freqs])

        # these two arrays needs to be two dimensional, though the content is one dimensional,
        # to be in sync with original data format, which is 1RX2C

        # generalize the "format"
        format1 = str(bin_table_data_time.shape[1]) + "D8.3"
        format2 = str(bintable_data_freqs.shape[1]) + "D8.3"
        col1 = pyfits.Column(name='TIME', format=format1, array=bin_table_data_time)
        col2 = pyfits.Column(name='FREQUENCY', format=format2, array=bintable_data_freqs)
        cols = pyfits.ColDefs(np.asarray([col1, col2]))

        # create bintable and  add it to hdulist
        tbhdu = pyfits.BinTableHDU.from_columns(cols)
        new_hdulist.append(tbhdu)

        # return new pycallisto object
        return PyCallisto(new_hdulist)

    def slice_time_axis(self, time1, time2):
        """
        Make a slice of input radiohelliograph observation fits file along a time axis and return a new object

        input arguments:
            time1 (string): start of a slice, time in HH:MM:SSformat
            time2 (string): end of a slice, time in HH:MM:SS format

        Returns
            PyCallisto object
        """

        # assuming that the input file has same start date and end date  (i.e. observed on the same day)
        # time1 and time 2 is in HH:MM:SS
        if not (len(time1.split(':')) == 3):
            raise Exception("Time format not proper, please provide time in HH:MM:SS format")

        if not (len(time2.split(':')) == 3):
            raise Exception("Time format not proper, please provide time in HH:MM:SS format")

        # check the dates
        start_date = utils.to_date(self.imageHeader['DATE-OBS'])
        end_date = utils.to_date(self.imageHeader['DATE-END'])
        if not start_date == end_date:
            raise Exception("start_date and end date differ, right now we do not support this")

        # get the times in datetime.date format for easy manipulation
        start_time = utils.to_time(self.imageHeader['TIME-OBS'])  # datetime.date object
        start_time = dt.datetime.combine(start_date, start_time)  # datetime.datetime object
        end_time = utils.to_time(self.imageHeader['TIME-END'])
        end_time = dt.datetime.combine(start_date, end_time)

        time1 = utils.to_time(time1)
        time1 = dt.datetime.combine(start_date, time1)
        time2 = utils.to_time(time2)
        time2 = dt.datetime.combine(start_date, time2)

        # check the "time constraints"
        if not (time1 < time2):
            time1, time2 = time2, time1
        if not (start_time < time1):
            # print("Time1 out of bound, can't slice!")
            print("Start time of input file : ", start_time)
            print("End time of input file : ", end_time)
            raise Exception("Time1 out of bound, can't slice!")

        if not (end_time > time2):
            # print("Time2 out of bound, can't slice!")
            print("Start time of input file : ", start_time)
            print("End time of input file : ", end_time)
            raise Exception("Time2 out of bound, can't slice!")

        # get the start_pixel and end_pixel in piselterms from input of time1 and time2
        start_pixel = time1 - start_time
        start_pixel = utils.tosec(start_pixel)  # in seconds
        start_offset = start_pixel
        start_pixel = int(start_pixel / float(self.imageHeader['CDELT1']))

        end_pixel = time2 - start_time
        end_pixel = utils.tosec(end_pixel)
        end_pixel = int(end_pixel / float(self.imageHeader['CDELT1']))

        # do the actual slicing
        image_data = copy.deepcopy(self.imageHdu.data)
        image_data = image_data[:, start_pixel:end_pixel]

        # update the image_header accordingly
        image_header = copy.deepcopy(self.imageHeader)
        image_header['NAXIS1'] = image_data.shape[1]
        image_header['TIME-OBS'] = str(time1.time())
        image_header['TIME-END'] = str(time2.time())
        image_header['DATAMIN'] = image_data.min()
        image_header['DATAMAX'] = image_data.max()
        image_header['CRVAL1'] = int(self.imageHeader['CRVAL1']) + start_offset

        # create a primary hdu containing this image and header data
        image_hdu = pyfits.PrimaryHDU(image_data, header=image_header)
        new_hdu_list = pyfits.HDUList([image_hdu])

        # create a new bintable and update it in data structure
        x_length = end_pixel - start_pixel
        range_list = [x * 1.0 for x in range(x_length)]
        bin_table_data_time = np.array([range_list])
        bin_table_data_freqs = list(self.binTableHdu.data[0][1].copy())
        bin_table_data_freqs = np.array([bin_table_data_freqs])

        format1 = str(bin_table_data_time.shape[1]) + "D8.3"
        format2 = str(bin_table_data_freqs.shape[1]) + "D8.3"
        col1 = pyfits.Column(name='TIME', format=format1, array=bin_table_data_time)
        col2 = pyfits.Column(name='FREQUENCY', format=format2, array=bin_table_data_freqs)
        cols = pyfits.ColDefs(np.asarray([col1, col2]))
        tb_hdu = pyfits.BinTableHDU.from_columns(cols)

        # add it to new hdulist we have created
        new_hdu_list.append(tb_hdu)

        return PyCallisto(new_hdu_list)

    def slice_frequency_axis(self, freq1, freq2):
        """
        Make a slice of input radiohelliograph observation fits file along a frequency axis

        Input arguments:
            freq1 (int): start of a slice
            freq2 (int): end of a slice

        Returns:
            new PyCallisto object
        """

        # convert  the input frequency from string to int
        freq1 = int(freq1)
        freq2 = int(freq2)
        # swap frequencies if not in proper orde
        if freq1 > freq2:
            freq1, freq2 = freq2, freq1

        # open input fits file
        hdus = self.hdus
        image_hdu = hdus[0]
        bin_table_hdu = hdus[1]
        image_header = image_hdu.header
        image_data = image_hdu.data
        bintbl_freq_data = bin_table_hdu.data[0][1]

        # check the frequencies
        start_freq = int(bintbl_freq_data[-1])
        end_freq = int(bintbl_freq_data[0])

        if freq1 < start_freq or freq1 > end_freq:
            # print("Frequency out of bound, cannot slice")
            print("Start Frequency of input file : ", start_freq)
            print("End Frequency of input file ", end_freq)
            raise Exception("Frequency out of bound, cannot slice")

        if freq2 < start_freq or freq2 > end_freq:
            # print "Frequency out of bound, cannot slice"
            print("Start Frequency of input file : ", start_freq)
            print("End Frequency of input file ", end_freq)
            raise Exception("Frequency out of bound, cannot slice")

        if freq2 - freq1 < 1:
            # print "Too thin slice demanded, cannot slice thinner than 1 unit"
            print("Start Frequency of input file : ", start_freq)
            print("End Frequency of input file ", end_freq)
            raise Exception("Too thin slice demanded, cannot slice thinner than 1 unit")

        # update the bintabledata according to freq1 and freq2
        try:
            slice_pt1 = np.argwhere((bintbl_freq_data >= freq1) & (bintbl_freq_data <= freq2))[0][0]
            slice_pt2 = np.argwhere((bintbl_freq_data > freq1) & (bintbl_freq_data < freq2))[-1][-1] + 1
        # +1 at slicept2 is needed to accomodate the slicing style of python to get the slice of exact size

        except:
            # print("frequency limits given are smaller than single channel, please increase it and retry")
            print("Start Frequency of input file : ", start_freq)
            print("End Frequency of input file : ", end_freq)
            raise Exception("frequency limits given are smaller than single channel, please increase it and retry")

        if (slice_pt2 - slice_pt1) < 1:
            print("Start Frequency of input file : ", start_freq)
            print("End Frequency of input file : ", end_freq)
            raise Exception("frequency limits given are smaller than single channel, please increase it and retry")

        bin_table_data_freqs = bintbl_freq_data[np.argwhere((bintbl_freq_data >= freq1) & (bintbl_freq_data <= freq2))]
        bin_table_data_freqs = np.array([bin_table_data_freqs])

        x_length = int(image_data.shape[1])
        range_list = [x * 1.0 for x in range(x_length)]
        bintable_data_time = np.array([range_list])

        format1 = str(bintable_data_time.shape[1]) + "D8.3"
        format2 = str(bin_table_data_freqs.shape[1]) + "D8.3"
        col1 = pyfits.Column(name='TIME', format=format1, array=bintable_data_time)
        col2 = pyfits.Column(name='FREQUENCY', format=format2, array=bin_table_data_freqs)

        cols = pyfits.ColDefs(np.asarray([col1, col2]))
        tb_hdu = pyfits.BinTableHDU.from_columns(cols)

        # update the image_header accordingly
        image_header = copy.deepcopy(self.imageHeader)
        image_data = copy.deepcopy(image_data)
        image_data = image_data[slice_pt1:slice_pt2, :]
        image_header['NAXIS2'] = image_data.shape[0]
        image_header['DATAMIN'] = image_data.min()
        image_header['DATAMAX'] = image_data.max()
        image_header['CRVAL2'] = image_data.shape[0]

        # create a primary hdu containing this image and header data
        image_hdu = pyfits.PrimaryHDU(image_data, header=image_header)
        new_hdu_list = pyfits.HDUList([image_hdu])

        # add it to hdulist, and write hdulist  to  fits file
        new_hdu_list.append(tb_hdu)

        return PyCallisto(new_hdu_list)

    def subtract_background(self):
        """
        Estimate and subtract background from a fits file

        returns a PyCallisto object
        """

        # open input fits file
        hdus = self.hdus
        image_hdu = hdus[0]
        bin_table_hdu = hdus[1]
        image_header = image_hdu.header
        image_data = image_hdu.data

        # bintblfreqdata = bin_table_hdu.data[0][1]

        # estimate background and subtract from original data
        med_col = np.median(image_data, axis=1, keepdims=True)
        image_data = copy.deepcopy(image_data)
        image_data = image_data - med_col

        # we need to update the header accordingly
        # compose a new image header
        image_header = copy.deepcopy(image_header)
        image_header['DATAMIN'] = image_data.min()
        image_header['DATAMAX'] = image_data.max()

        # create a primary hdu containing this updated image data
        image_hdu = pyfits.PrimaryHDU(image_data, header=image_header)
        new_hdu_list = pyfits.HDUList([image_hdu])

        # add the unchanged bin_table_hdu to hdulist
        new_hdu_list.append(bin_table_hdu)

        return PyCallisto(new_hdu_list)

    def mean_light_curve(self, plot=True, out_image="timeseries.png", return_data=False, fig_size=(8, 8), font_size=14,
                         grid=True):
        """
        Create mean light curve by Collapse the 2D fits image along time axis.
        Eithers saves the image or returns the data

        Input arguments:
            plot (Boolean): plot it or not, Default True
            out_image(string): name of the image to be saved
            return_data (Boolean):  return data or not (default False)

        returns:if return_data is set to True
                list of collapsed 1d array (1d numpy array),
                 respective time in sec of a day(1d numpy array),
                 respective time in datetime.datetime object(1d list)
        """
        # open input fits file
        hdus = self.hdus
        image_hdu = hdus[0]
        image_header = image_hdu.header
        image_data = image_hdu.data
        bin_table_hdu = hdus[1]
        # bintbl_freq_data = bin_table_hdu.data[0][1]

        # Sum the data along oth axis
        sum_image = np.sum(image_data, axis=0)
        sum_image = sum_image / image_data.shape[0]  # divide by 200 to normalize

        start_date = utils.to_date(image_header['DATE-OBS'])
        start_time = utils.to_time(image_header['TIME-OBS'])
        start_time = dt.datetime.combine(start_date, start_time)  # make a datetime object
        end_time = utils.to_time(image_header['TIME-END'])
        end_time = dt.datetime.combine(start_date, end_time)
        time_delta = dt.timedelta(days=0, hours=0, minutes=0, seconds=image_header['CDELT1'])

        def get_time_axis(start, end, delta):
            curr = start
            while curr < end:
                yield curr
                curr += delta

        time_axis = [time for time in get_time_axis(start_time, end_time, time_delta)]

        if plot:
            fig, ax = plt.subplots(figsize=(fig_size[0], fig_size[1]))
            ax.plot_date(time_axis, sum_image, 'b-', xdate=True)
            ax.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
            plt.xticks(rotation=45)
            plt.xlabel('Universal Time', fontSize=font_size)
            plt.ylabel('Total count', fontSize=font_size)
            # title = infits + "\n" + image_header['DATE-OBS'] + "::" + image_header['TIME-OBS']
            # title = self.image_header['CONTENT']

            plt.title("Mean light curve", fontSize=font_size)
            if grid:
                plt.grid()
            plt.savefig(out_image)

        else:
            if return_data:
                # generate numpy array of sec. of a day
                x_start = int(image_header['CRVAL1'])
                x_step = float(image_header['CDELT1'])
                x_length = int(image_header['NAXIS1'])
                x_end = x_start + (x_step * x_length)
                time_in_sec = np.arange(x_start, x_end, x_step)

                # compose data into a list
                data = [sum_image, time_in_sec, time_axis]
                # hdus.close()
                return data

    def mean_spectrum(self, plot=True, out_image="frequency_series.png", return_data=False, font_size=14, grid=True):
        """
        Create mean spectrum by collapsing the 2D fits image along frequency axis.
        Either plots image or returns data.

        INput arguments:
            plot (Boolean): plot or not ?, Default True
            out_image (String):  name of the png file to plot
            return_data (Boolean) : return data or not (Default False)

        returns : list of  collapsed 1d array (1d numpy array),
                respective frequency channels (1d numpy array)
        """

        # open input fits file
        hdus = self.hdus
        image_hdu = hdus[0]
        bin_table_hdu = hdus[1]
        image_header = image_hdu.header
        image_data = image_hdu.data
        bintbl_freq_data = bin_table_hdu.data[0][1]

        # Sum the data along oth axis
        sum_image = np.sum(image_data, axis=1)
        sum_image = sum_image / image_data.shape[1]  # divide by xsize to normalize

        if plot:
            plt.clf()
            plt.plot(bintbl_freq_data, sum_image, 'b-', )  # check why we need to add -1 here
            plt.axis([bintbl_freq_data[-1], bintbl_freq_data[0], sum_image.min(), sum_image.max()])
            plt.xlabel("Frequency (MHz)", fontSize=font_size)
            plt.ylabel('Total count', fontSize=font_size)
            # title = infits + "\n" + image_header['DATE-OBS'] + "::" + image_header['TIME-OBS'] + " To " +
            # image_header['DATE-END'] + "::" + image_header['TIME-END'] title = self.image_header['CONTENT']
            plt.title("Mean Spectrum", fontSize=font_size)
            if grid:
                plt.grid()
            plt.savefig(out_image)
        else:
            if return_data:
                # compose data into a list
                data = [sum_image, bintbl_freq_data]
                hdus.close()
                return data

    def print_frequencies(self):
        """
        Print the list of frequencies.
        """
        # open input fits file
        hdus = self.hdus
        image_hdu = hdus[0]
        bin_table_hdu = hdus[1]
        # image_header = image_hdu.header
        # image_data = image_hdu.data
        # bintbl_freq_data = bin_table_hdu.data[0][1]

        print(np.array2string(bin_table_hdu.data[0][1]))

    def light_curve(self, frequency, plot=True, out_image="Lightcurve.png", return_data=False, fig_size=(8, 8),
                    font_size=14, grid=True):
        """
        Plot the lightcurve for given frequency, i.e. time vs amplitude, or return the data

        Input arguments:
            frequency: frequency to plot lightcurve
            plot (Boolean): plot it or not (Default = True)
            out_image: Name of the image to save, default is "Lightcurve.png"
            return_data (Boolean) : return the data or not (default = False)
            figsize: size of the fig (default = (8, 8))
            font_size: font_size of the used in plots (default =14)
            grid: to plot agrid or not (default = True)
        returns:if return_data is set to True
#				return a tuple of (time_axis, light_curve)
                where
                    time_axis is array of python dataetime object
                    light_curve is numpy array
        """

        # open input fits file
        hdus = self.hdus
        image_hdu = hdus[0]
        bin_table_hdu = hdus[1]
        image_header = image_hdu.header
        image_data = image_hdu.data
        # bintbl_freq_data = bin_table_hdu.data[0][1]
        bintbl_freq_data = bin_table_hdu.data[0][1]

        # check input frequency
        start_freq = int(bintbl_freq_data[-1])
        end_freq = int(bintbl_freq_data[0])
        if frequency < start_freq or frequency > end_freq:
            raise Exception("Input frequency is out of limit for this data, aborting the operation")

        # get the nearest value in the bintable to the user input frequency
        min_diff = (np.abs(bintbl_freq_data - frequency)).min()
        idx = (np.abs(bintbl_freq_data - frequency)).argmin()
        nearest_frequency = bintbl_freq_data[idx]

        # give warning if the difference between input frequency and the nearest frequency is greater than 5
        if min_diff > 5:
            print(
                "Please note that the difference between demanded frequency and the nearest one in data is %g"
                % min_diff)

        # get the corresponding row from image_data
        light_curve = image_data[idx, :]

        # get startdate-time and enddate-time
        start_date = utils.to_date(image_header['DATE-OBS'])
        start_time = utils.to_time(image_header['TIME-OBS'])
        start_time = dt.datetime.combine(start_date, start_time)  # make a datetime object
        end_time = utils.to_time(image_header['TIME-END'])
        end_time = dt.datetime.combine(start_date, end_time)

        def get_time_axis(start, end, delta):
            curr = start
            while curr < end:
                yield curr
                curr += delta

        time_axis = [time for time in get_time_axis(start_time, end_time,
                                                    (end_time - start_time) / light_curve.shape[0])]

        if plot:
            plt.clf()
            fig, ax = plt.subplots(figsize=(fig_size[0], fig_size[1]))
            ax.plot_date(time_axis, light_curve, 'b-', xdate=True)
            ax.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
            plt.xticks(rotation=45)

            plt.xlabel('Universal Time', fontSize=font_size)
            plt.ylabel('Amplitude', fontSize=font_size)
            # title = infits + "\n" + image_header['DATE-OBS'] + "::" + image_header['TIME-OBS']
            # title = self.image_header['CONTENT']
            title = "Light curve - " + str(round(nearest_frequency, 2)) + " MHz"
            plt.title(title, fontSize=font_size)
            if grid:
                plt.grid()
            plt.savefig(out_image)
        else:
            if return_data:
                return time_axis, light_curve

    def spectrum(self, in_date, in_time, binning=2, binning_method='avg', plot=True, out_image="singletimespectrum.png",
                 return_data=False, fig_size=(8, 6), font_size=14, grid=True):
        """
        Plot the spectrum for a given time, i.e. amplitude at all frequencies at given time

        Input arguments:
            indate :date to plot spectrum of,
                    should be a python datetime object or
                    string of format 'YYYY/MM/DD'
            intime : time to plot spectrum of,
                    should be a python datetime.time object or
                    string of format 'HH:MM:SS'
            binning: pixel-level binning while selecting a row (default = 2)
                caution, only advanced users should alter this parameter
            binnigmethod : either of avg, sum, med (default = 'avg')
            plot (Boolean) : plot it or not, Default True
            outimage : Name of the image to save, default is "singletimespectrum.png"
            returndata (Boolean) :   default False
            figsize: size of the fig (default = (8, 8))
            font_size: font_size of the used in plots (default =14)
            grid: to plot agrid or not (default = True)

        returns:if returndata is set to True
                return a tuple of (bintbl_freq_data, spectrum)
                where
                    bintbl_freq_data :  list of frequencies
                    spectrum :  intensity count at these respective ffrequencies

        """

        # open input fits file
        hdus = self.hdus
        image_hdu = hdus[0]
        bin_table_hdu = hdus[1]
        image_header = image_hdu.header
        image_data = image_hdu.data
        bintbl_freq_data = bin_table_hdu.data[0][1]
        bintbl_freq_data = bin_table_hdu.data[0][1]

        # get startdate-time and enddate-time of the file
        start_date = utils.to_date(image_header['DATE-OBS'])
        start_time = utils.to_time(image_header['TIME-OBS'])
        start_time = dt.datetime.combine(start_date, start_time)  # make a datetime object
        end_time = utils.to_time(image_header['TIME-END'])
        end_time = dt.datetime.combine(start_date, end_time)
        # check the format of in date and intime

        if isinstance(in_date, dt.date) and isinstance(in_time, dt.time):
            in_date_time = dt.datetime.combine(in_date, in_time)
        else:
            if isinstance(in_date, str) and isinstance(in_time, str):
                if len(in_date.split('/')) == 3:
                    # year, month, day = indate.split('/')
                    in_date = utils.to_date(in_date)
                else:
                    raise Exception("Date string not in proper format")

                if len(in_time.split(':')) == 3:
                    # hr, mint, sec = intime.split(':')
                    in_time = utils.to_time(in_time)
                else:
                    raise Exception("Time string not in proper format")
                in_date_time = dt.datetime.combine(in_date, in_time)
            else:
                raise Exception("Date and/or time not in proper format")

        # check input time
        if in_date_time < start_time or in_date_time > end_time:
            print("Provided datetime is " + str(in_date_time))
            print("Start datetime for this file is " + str(start_time))
            print("End datetime for this file is" + str(end_time))
            raise Exception("Input time is out of limit for this data, aborting the operation")

        # get the columnnumber of given datetime in a image_data
        def get_time_axis(start, end, delta):
            curr = start
            while curr < end:
                yield curr
                curr += delta

        time_axis = [time for time in get_time_axis(start_time, end_time,
                                                    (end_time - start_time) / image_data.shape[1])]

        time_index = time_axis.index(in_date_time)
        spectrum = image_data[:, time_index - binning:time_index + binning + 1]

        # print(spectrum.shape)
        if binning_method == 'med':
            spectrum = np.median(spectrum, axis=1)
        if binning_method == 'sum':
            spectrum = spectrum.sum(axis=1)
        if binning_method == 'avg':
            spectrum = np.average(spectrum, axis=1)

        if plot:
            plt.clf()
            fig, ax = plt.subplots(figsize=(fig_size[0], fig_size[1]))
            ax.plot_date(bintbl_freq_data, spectrum, 'b-', xdate=True)
            ticks = list(
                np.linspace(bintbl_freq_data[-1], bintbl_freq_data[0], 20).astype('int'))  # calcult 200 ticks positions
            plt.xticks(ticks, ticks, rotation=45)
            # ax.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
            plt.xlabel('Frequency (MHz)', fontSize=font_size)
            plt.ylabel('Amplitude', fontSize=font_size)
            # title = infits + "\n" + image_header['DATE-OBS'] + "::" + image_header['TIME-OBS']
            # title = self.image_header['CONTENT']
            title = "Spectrum - " + str(in_date_time)
            plt.title(title)
            if grid:
                plt.grid()
            plt.savefig(out_image)

        else:
            if return_data:
                return bintbl_freq_data, spectrum

    def universal_plot(self, plot_name="universal_plot_with_add_processing.png", title='Universal Plot',
                       return_plot=False, xtick=3, ytick=3, end_pts=[False, False], blevel=0, fig_size=(10, 8),
                       cmap=cm.jet, label_font_size=10, title_font_size=14, color_bar=True, color_bar_ori='horizontal'):
        """
        plot universal plot


        input arguments:
            plot_name: name of the image to be saved
            title: title of plot
            return_plot: to return matplotlib plot object or not (defayult = False)
            xtick: frequency of xticks in minutes (default = 3)
            ytick: frequency of yticks in minutes (default = 3)
            end_pts:	plot endpoints on x and y axis irrespective of ticks (fefault= (False, false))
            blevel: background level (default = 0)
            fig_size:  tuple representing size of image (default = (8,6))
            cmap: a matplotlib colormap object (default = cm.jet)
            label_font_size: font_size used in plot for labels (default = 10)
            title_font_size: font_size used in plot for title (default = 14)
            color_bar: to plot a colorbar or not (default = true)
            color_bar_ori: colorbar orientation (default = 'vertical')
        """

        # create global plot object with subplot
        fig = plt.figure(figsize=(fig_size[0], fig_size[1]))

        # plot the spectrogram in the one subplot
        ax1 = plt.subplot2grid(shape=(8, 8), loc=(0, 0), rowspan=5, colspan=5)

        # get the start time and end time of observation
        start_date = utils.to_date(self.imageHeader['DATE-OBS'])
        start_time = utils.to_time(self.imageHeader['TIME-OBS'])

        start_time = dt.datetime.combine(start_date, start_time)  # make a datetime object
        end_time = utils.to_time(self.imageHeader['TIME-END'])
        end_time = dt.datetime.combine(start_date, end_time)

        # get the frequencies
        freqs = self.binTableHdu.data['frequency'][0]  # these are true frequencies

        # set the limits for plotting
        x_lims = [start_time, end_time]
        x_lims = mdates.date2num(x_lims)  # dates to numeric values
        y_lims = [freqs[-1], freqs[0]]

        im1 = ax1.imshow(self.imageHdu.data, extent=[x_lims[0], x_lims[1], y_lims[0], y_lims[1]], aspect='auto',
                         cmap=cmap, vmin=blevel)

        # set ytick labels to be intigers only, in case we have a small frequency range, pyplot tends to show fraction
        ax1.get_yaxis().get_major_formatter().set_useOffset(False)
        plt.minorticks_on()

        if end_pts[0]:
            # get the start time and end time ans set ticks at those position son x-axis using minor_locator
            total_sec = utils.tosec(end_time - start_time)
            ax1.xaxis.set_minor_locator(mdates.SecondLocator(bysecond=None, interval=total_sec, tz=None))
            ax1.xaxis.set_minor_formatter(DateFormatter('%H:%M:%S'))
        fig.autofmt_xdate()

        if end_pts[1]:
            # set y minorticks for first and last frequency
            plt.yticks(list(plt.yticks()[0]) + y_lims)

        # plt.xlabel('Universal Time')
        plt.ylabel('Frequency (MHz)', fontSize=label_font_size)

        ax1.tick_params(direction='in', axis='both', which='both')
        ax1.yaxis.set_ticks_position('both')
        ax1.xaxis.set_ticks_position('both')
        ax1.tick_params(labelbottom=False)

        ax2 = plt.subplot2grid(shape=(8, 8), loc=(5, 0), rowspan=3, colspan=5, sharex=ax1)

        # use meanlightcurve to get the lightcurve data
        data = self.mean_light_curve(plot=False, return_data=True)
        sum_image, time_in_sec, time_axis = data

        plt.minorticks_on()
        plt.xticks(rotation=45)
        ax2.plot(time_axis, sum_image, 'k')
        ax2.tick_params(direction='in', axis='both', which='both')
        ax2.yaxis.set_ticks_position('both')
        ax2.xaxis.set_ticks_position('both')
        ax2.margins(0.0)
        ax2.set_xlabel('Universal Time')

        # define ax3
        ax3 = plt.subplot2grid(shape=(8, 8), loc=(0, 5), rowspan=5, colspan=3, sharey=ax1)
        # use meanspectrum to get spectrum data
        data = self.mean_spectrum(plot=False, return_data=True)
        sum_image, bintblfreqdata = data
        # use spectrum data to plot spectrum in third subplot
        ax3.plot(sum_image, bintblfreqdata, 'k')
        plt.minorticks_on()
        ax3.tick_params(direction='in', axis='both', which='both')
        ax3.yaxis.set_ticks_position('both')
        ax3.xaxis.set_ticks_position('both')
        ax3.tick_params(labelleft=False)
        ax3.tick_params(labelright=True)
        ax3.invert_xaxis()
        plt.xticks(rotation=45)
        ax3.margins(0.0)

        plt.suptitle(title, fontsize=title_font_size + 4)
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
        fig.colorbar(im1, cax=cbar_ax)
        plt.savefig(plot_name)
