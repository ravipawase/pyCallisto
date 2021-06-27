import datetime as dt
import astropy.io.fits as pyfits


def check_fits_callisto(fitsfile):
    """Check whether fits file has two HDUs or not
    """
    hdus = pyfits.open(fitsfile)
    if len(hdus) == 2:
        hdus.close()
        return True
    else:
        hdus.close()
        return False


def check_bin_table(fitsfile):
    """Check whether fits file bintable is valid
    """
    hdus = pyfits.open(fitsfile)
    bintablehdu = hdus[1]

    try:
        bin_data = bintablehdu.data  # try accessing the data
        hdus.close()
        return True
    except:
        hdus.close()
        return False


def tosec(td):
    """Calculate the total seconds in a timedate.timedate object
    """
    return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10 ** 6) / 10 ** 6


def to_date(string):
    """Break the string with "/" separator
    return the equivalent datetime.date object
    """
    yr, mnth, day = string.split('/')
    yr, mnth, day = int(yr), int(mnth), int(day)
    return dt.date(yr, mnth, day)


def to_time(string):
    """
    break the string with "/" separator
    return the datetime.time object
    """
    hr, mn, sec = string.split(':')
    sec = sec.split('.')[0]  # if sec has a fractional value
    hr, mn, sec = int(hr), int(mn), int(sec)
    # if second's value is 60, replace it to 00 and update and manage propogative changes in the minute and hour values
    # some old data from some observatories has this "leap second" problem
    if sec == 60:
        sec = 00
        mn = mn + 1
        if mn > 59:
            hr = hr + 1
            mn = mn - 60
            if hr > 23:
                hr = hr - 24
    return dt.time(hr, mn, sec)


def visualise(plt_object, show=True, outpath='tmp.png'):
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
