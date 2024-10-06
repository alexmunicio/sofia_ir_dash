import astropy.io.fits as fits
import pandas as pd
import numpy as np
from astropy.convolution import convolve
from astropy.convolution import Box1DKernel
import datetime,pytz

MONTH_CONVERT = {1: "January", 2: "February", 3: "March", 4: "April", 5: "May", 6: "June", 7: "July", 8: "August", 9: "September", 10: "October", 11: "November", 12: "Decmeber"}

# converts altitude from feet to kilometers
FEET_TO_KILOMETERS = 1/3280.84

# function to return normalization value
def get_fluxnorm(flux,trans):

    fti=np.where((np.isfinite(trans))&(np.isfinite(flux)))[0]
    ftil=len(fti)
    ftil4=int(ftil/4.)
    fti=fti[ftil4:ftil-ftil4] #take only the middle half of non-nan flux
    fintrans=trans[fti]
    finflux=flux[fti]

    #tmp 21/5/18 changed 0.99 to 0.95 to accomodate texes data; make sure it still works for exes
    transmax=np.max(fintrans)
    finfluxbase=finflux[fintrans>=0.95*transmax] #take only flux where trans is near max (to filter out trans lines)
    ffb_mean=np.mean(finfluxbase)
    ffb_std=np.std(finfluxbase)
    ffb_useful=(finfluxbase>=ffb_mean-ffb_std)&(finfluxbase<=ffb_mean+ffb_std) #take only parts that are sig away from mean

    norm=np.mean(finfluxbase[ffb_useful])/transmax #divide by trans max to get ~1

    return norm

# function to pull EXES data from a FITS file and returns data as a dictionary
def get_exes_file_data(file_name, smooth_width = 9, dir = "EXES_Files"):

    path = dir + "/" + file_name
    # print("\n", path, "\n")
    with fits.open(path) as hdu:
        primary_hdu = hdu[0]
        primary_header = primary_hdu.header
        primary_data = primary_hdu.data
    
    # label the lists within the data
    wavenumber = primary_data[0]
    flux = primary_data[1]
    uncertainty = primary_data[2]
    atran = primary_data[3]

    # calculate the wavelengths for each wavenumber
    wavelength = 10000/wavenumber

    # calculate the values for normalized flux
    norm = get_fluxnorm(flux, atran)
    norm_flux = flux/norm

    if(smooth_width is not None):
        # convolve flux if a smooth width is given
        smooth_flux = convolve(norm_flux, Box1DKernel(smooth_width), preserve_nan = True)
    else: smooth_flux = None

    # initialize useful header info
    obj = primary_header["OBJECT"]
    TELEL = primary_header["TELEL"]
    lat = primary_header['LAT_STA']
    lon = primary_header['LON_STA']
    ALTI_END = primary_header['ALTI_END']
    ALTI_STA = primary_header['ALTI_STA']
    avg_ALTI = (ALTI_END + ALTI_STA)/2
    avg_ALTI_kilometers = avg_ALTI * FEET_TO_KILOMETERS
    Tout = primary_header['TEMP_OUT']

    # initialize experiment date
    year = int(primary_header['DATE-OBS'][0:4])
    month = int(primary_header['DATE-OBS'][5:7])
    day = int(primary_header['DATE-OBS'][8:10])
    hour = int(primary_header['DATE-OBS'][11:13])
    minute = int(primary_header['DATE-OBS'][14:16])
    date = datetime.datetime(year, month, day, hour, minute, tzinfo = pytz.utc)

    # create a summarative title string for the spectrum plot
    spectrum_title = obj + " | " + str(TELEL) + "$^\circ$" + " | " + str(avg_ALTI) + " | " + MONTH_CONVERT[month] + "," + str(year)

    # create a neatly organized list for the map table's "cellText" argument
    map_table_array = [obj, str(TELEL) + "\N{DEGREE SIGN}", str(avg_ALTI) + " ft", str(round(-1 * lon, 3)) + "\N{DEGREE SIGN}W, " + str(round(lat, 3)) + "\N{DEGREE SIGN}N"]

    # return information as a dictionary
    return {
        "hdu": primary_hdu,
        "header": primary_header,
        "data": primary_data,
        "wavenumber": wavenumber,
        "flux": flux,
        "atran": atran,
        "norm_flux": norm_flux,
        "smooth_flux": smooth_flux,
        "uncertainty": uncertainty,
        "wavelength": wavelength,
        "wavelength_range": [min(wavelength), max(wavelength)],
        "object": obj,
        "latitude": lat,
        "longitude": lon,
        "start_altitude": ALTI_STA, # altitude is in units of feet
        "end_altitude": ALTI_END,
        "avg_altitude": avg_ALTI,
        "avg_altitude_km": avg_ALTI_kilometers, # altitude is in units of kilometers
        "temperature": Tout,
        "telescope_elevation_angle": TELEL,
        "map_table_array": map_table_array,
        "spectrum_title": spectrum_title,
        "date": date
    }

# this function was made to do exactly the same thing as "get_exes_file_data()" with one distinct difference:
# this function will switch any big endian arrays stored in the FITS file into little endian arrays
# this was done because the pandas library uses little endian arrays instead of big endians, and so FITS files can sometimes be incompitble with pandas
# Consequently, I made this function to deal with that
# if you don't know what a little or big endian is, don't worry, you're not missing out on anything cool
# little and big endians are two distinct memory allocation conventions used to store arrays in low level code like Assembly.
# it's usually not an issue in high level code, like python, but it came up while I was working on this
# so if you don't know which function to use, both will probably work fine, but this one is probably a safer bet than get_exes_file_data()

# TLDR: this does the same thing as "get_exes_file_data()", but it prevents a bug from occuring when the FITS data is stored in a pandas dataframe
def get_endian_modified_exes_file_data(file_name, smooth_width = 9, dir = "EXES_Files"):

    # initialize full path name
    path = dir + "/" + file_name
    # print("\n", path, "\n")

    hdulist = fits.open(path)
    primary_hdu = hdulist[0]
    primary_header = primary_hdu.header
    primary_data = primary_hdu.data

    #check to see if its big endian
    if primary_data.dtype.byteorder == '>' : 
        primary_data = primary_data.byteswap().newbyteorder()

    wavenumber = primary_data[0]
    flux = primary_data[1]
    uncertainty = primary_data[2]
    atran = primary_data[3]

    # calculate the wavelengths for each wavenumber
    wavelength = 10000.0 / wavenumber

    # calculate the values for normalized flux
    norm = get_fluxnorm(flux, atran)
    norm_flux = flux/norm

    if(smooth_width is not None):
        # convolve flux if a smooth width is given
        smooth_flux = convolve(norm_flux, Box1DKernel(smooth_width), preserve_nan = True)
    else: smooth_flux = None

    df = pd.DataFrame({
        "wavenumber": wavenumber,
        "atran": atran,
        "flux": flux,
        "uncertainty": uncertainty,
        "norm_flux": norm_flux,
        "smooth_flux": smooth_flux,
    })

    # initialize useful header info
    obj = primary_header["OBJECT"]
    TELEL = primary_header["TELEL"]
    lat = primary_header['LAT_STA']
    lon = primary_header['LON_STA']
    ALTI_END = primary_header['ALTI_END']
    ALTI_STA = primary_header['ALTI_STA']
    avg_ALTI = (ALTI_END + ALTI_STA)/2
    avg_ALTI_kilometers = avg_ALTI * FEET_TO_KILOMETERS
    Tout = primary_header['TEMP_OUT']

    # initialize experiment date
    year = int(primary_header['DATE-OBS'][0:4])
    month = int(primary_header['DATE-OBS'][5:7])
    day = int(primary_header['DATE-OBS'][8:10])
    hour = int(primary_header['DATE-OBS'][11:13])
    minute = int(primary_header['DATE-OBS'][14:16])
    date = datetime.datetime(year, month, day, hour, minute, tzinfo = pytz.utc)

    # create a summarative title string for the spectrum plot
    spectrum_title = obj + " | " + str(TELEL) + "$^\circ$" + " | " + str(avg_ALTI) + " | " + MONTH_CONVERT[month] + "," + str(year)

    # create one for the dashboard that works better with the formatting
    dashboard_spectrum_title = obj + " | " + str(TELEL) + " deg." + " | " + str(avg_ALTI) + " | " + MONTH_CONVERT[month] + "," + str(year)

    # create a neatly organized list for the map table's "cellText" argument
    map_table_array = [obj, str(TELEL) + "\N{DEGREE SIGN} ", str(avg_ALTI) + "ft f", str(round(-1 * lon, 3)) + "\N{DEGREE SIGN}W, " + str(round(lat, 3)) + "\N{DEGREE SIGN}N"]

    hdulist.close()

    # return information as a dictionary
    return {
        "hdu": primary_hdu,
        "header": primary_header,
        "data": primary_data,
        "wavenumber": wavenumber,
        "wavenumber_range":[min(wavenumber), max(wavenumber)],
        "flux": flux,
        "atran": atran,
        "dataframe": df,
        "norm_flux": norm_flux,
        "smooth_flux": smooth_flux,
        "uncertainty": uncertainty,
        "wavelength": wavelength,
        "wavelength_range": [min(wavelength), max(wavelength)],
        "object": obj,
        "latitude": lat,
        "longitude": lon,
        "start_altitude": ALTI_STA, # altitude is in units of feet
        "end_altitude": ALTI_END,
        "avg_altitude": avg_ALTI,
        "avg_altitude_km": avg_ALTI_kilometers, # altitude is in units of kilometers
        "temperature": Tout,
        "telescope_elevation_angle": TELEL,
        "map_table_array": map_table_array,
        "spectrum_title": spectrum_title,
        "dashboard_spectrum_title": dashboard_spectrum_title,
        "date": date
    }

  