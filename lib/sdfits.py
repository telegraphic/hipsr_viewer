#!/usr/bin/env python
"""
sdfits.py
----------

Functions and helpers to convert HIPSR5 files into SD-FITS

"""

import sys, os, re, time
from datetime import datetime
import pyfits as pf, numpy as np, tables as tb
import pylab 
from optparse import OptionParser

import config as config
from printers import LinePrint
from hipsr6 import Hipsr6


__version__  = config.__version__
__author__   = config.__author__
__email__    = config.__email__
__license__  = config.__license__
__modified__ = datetime.fromtimestamp(os.path.getmtime(os.path.abspath( __file__ )))

path = os.getcwd()

def extractMid(x):
    """ Extract the mid part of an array """
    return x[len(x)/4:3*len(x)/4]

def fitLine(x, y, n_chans):
    """ Fit a line to data, only using central channels """
    x_fine = np.linspace(x[0], x[-1], n_chans)
    x = x[len(x)/4:3*len(x)/4]
    y = y[len(y)/4:3*len(y)/4]
    p = np.polyfit(x[::-1], y, 1)  # linear fit
    v = np.polyval(p, x_fine)
    return v
    
def loadDiodeTemp(h6, filename_x, filename_y):
    """ Load a diode temp csv """
    
    f_fine = h6.freqs
    f      = h6.freqs_cal
    
    temps_x = np.fromfile(filename_x).reshape([13,16])
    temps_y = np.fromfile(filename_y).reshape([13,16])
    
    temps_fine_x = np.zeros([13, 8192])
    temps_fine_y = np.zeros([13, 8192])
    
    for i in range(0,13):
        temps_fine_x[i] = fitLine(f, temps_x[i], 8192)
        temps_fine_y[i] = fitLine(f, temps_y[i], 8192)
        
    return temps_fine_x, temps_fine_y

def applyCal(beam, row, freqs, freqs_cal, cf, T_d_x, T_d_y):
    """ Apply basic calibration in Jy
    
    P_sys / (CF*P_on- CF* p_off) * T_d  = T_sys 
    
    P_sys: Total power in channel
    P_on:  Noise diode on
    P_off: Noise diode off
    T_d:   Temperature of the diode
    CF:    Cal factor which relates diode measurement P_on to P_sys
    
    """
    
    P_sys_xx = beam.cols.xx[row].astype('float')
    xx_on    = beam.cols.xx_cal_on[row].astype('float')
    xx_off   = beam.cols.xx_cal_off[row].astype('float')
    P_on_xx  = np.average(extractMid(xx_on))
    P_off_xx = np.average(extractMid(xx_off))
    
    #P_on_xx  = fitLine(freqs_cal, xx_on, len(freqs))
    #P_off_xx = fitLine(freqs_cal, xx_off, len(freqs))

    P_sys_yy = beam.cols.yy[row].astype('float')
    yy_on    = beam.cols.yy_cal_on[row].astype('float')
    yy_off   = beam.cols.yy_cal_off[row].astype('float')
    P_on_yy  = np.average(extractMid(yy_on))
    P_off_yy = np.average(extractMid(yy_off))
    
    #P_on_yy  = fitLine(freqs_cal, yy_on, len(freqs))
    #P_off_yy = fitLine(freqs_cal, yy_off, len(freqs))

    
    T_sys_xx = P_sys_xx / (cf*P_on_xx - cf*P_off_xx) * T_d_x
    T_sys_yy = P_sys_yy / (cf*P_on_yy - cf*P_off_yy) * T_d_y
    
    return T_sys_xx, T_sys_yy

def computeTsys(beam, row, T_d_x, T_d_y):
    """ Compute Tsys from frequency data 
    
    T_sys = Td / (Pon/Poff - 1)
    
    """
    
    xx_on    = beam.cols.xx_cal_on[row].astype('float')
    xx_off   = beam.cols.xx_cal_off[row].astype('float')
    
    yy_on    = beam.cols.yy_cal_on[row].astype('float')
    yy_off   = beam.cols.yy_cal_off[row].astype('float')
    
    T_sys_x = np.average(T_d_x[len(T_d_x)/4:3*len(T_d_x)/4]) / (xx_on/xx_off -1)
    T_sys_y = np.average(T_d_y[len(T_d_x)/4:3*len(T_d_x)/4]) / (yy_on/yy_off -1)
    
    l = len(T_sys_x)
    return np.average(T_sys_x[l/4:3*l/4]), np.average(T_sys_y[l/4:3*l/4])

def generateCards(filename):
  """
  Parses a text file and generates a pyfits card list.
  Do NOT feed this a full FITS file, feed it only a human-readable 
  FITS header template. 
  
  A text file is opened, acard is created from each line, then verified. 
  If the line does not pass verification, no card is appended.
  
  Parameters
  ----------
  filename: str
      name of text file header to open and parse
  """
  infile = open(filename)

  header = pf.Header()

  # Loop through each line, converting to a pyfits card
  for line in infile.readlines():
      line = line.rstrip('\n')
      line = line.strip()
      if(line == 'END'):
        break
      else:
        c = pf.Card().fromstring(line)
        c.verify() # This will attempt to fix issuesx[1]
        header.append(c)
        
  return header.cards

def fitsFormatLookup(x):
    """ Helper function to map FITS format codes into numpy format codes
    
    Notes
    -----
    FITS format code         Description                     8-bit bytes
    L                        logical (Boolean)               1
    X                        bit                             *
    B                        Unsigned byte                   1
    I                        16-bit integer                  2
    J                        32-bit integer                  4
    K                        64-bit integer                  4
    A                        character                       1
    E                        single precision floating point 4
    D                        double precision floating point 8
    C                        single precision complex        8
    M                        double precision complex        16
    P                        array descriptor                8    
    """
    
    # This is a python dictionary to lookup mappings.
    return {
           'L' : 'bool_',
           'X' : 'bool_',
           'B' : 'ubyte',
           'I' : 'int16',
           'J' : 'int32',
           'K' : 'int64',
           'A' : 'str_',
           'E' : 'float32',
           'D' : 'float64',
           'C' : 'complex64',
           'M' : 'complex128',
           'P' : 'float32'
    }.get(x, 'float32')

def formatLookup(format_str):
    """ Look up the format of a FITS string """
    pat   = '(\d+)([A-Z])'
    match = re.search(pat, format_str)
    #print match.group()
    
    data_len = int(match.group(1))
    data_fmt = str(match.group(2))
    np_fmt   = fitsFormatLookup(data_fmt)
    np_dtype = '%i%s'%(data_len, np_fmt)
    
    return np_dtype, data_len, np_fmt 

def generateZeros(num_rows, format, dim=None):
    """ Generate blank data to populate binary table 
    
    Used by generateSDFits() to form column definitions.
    
    Parameters
    ----------
    format: str
        FITS format code, e.g. 16A, 2048E
    dim: str
        dimensions for multidimensional data array, e.g. (1024,2,1,1)
        Defaults to None
    """
    
    np_dtype, data_len, np_fmt   = formatLookup(format)
    
    return np.zeros(num_rows, dtype=np_dtype)

def generatePrimaryHDU(hdu_header='header_primaryHDU.txt'):
    """ Generates the Primary HDU
    
    Parameters
    ----------
    hdu_header: string
        Name of the HDU header file to parse to generate the header.
        Defaults to header_primaryHDU.txt
    """
    
    hdu   = pf.PrimaryHDU()
    cards = generateCards(hdu_header)
    
    for card in cards:
        #print card
        if card.keyword == 'COMMENT':
            pass
            hdu.header.add_comment(card.value)
        elif card.keyword == 'HISTORY':
            pass
            hdu.header.add_history(card.value)
        else:
            hdu.header.set(card.keyword, card.value, card.comment)
    
    return hdu

def generateBlankDataHDU(num_rows=1, header_file='header_dataHDU.txt',
                   coldef_file='coldefs_dataHDU.txt'):
    """ Generate a blank data table with N rows.
    
    Parameters
    ----------
    num_rows: int
        The number of rows in the binary table.
    header_file: str
        Path to the header file. Defaults to 'header_dataHDU.txt'
    coldef_file: str
        Path to the file containing column definitions.
        Defaults to 'coldefs_dataHDU.txt'
    
    """
    
    cols = []
    
    # The column definitions are loaded from an external file, which is
    # parsed line-by-line, using regular experssions.
    
    unit_pat   = "unit\s*\=\s*'([\w/%]+)'"
    name_pat   = "name\s*\=\s*'([\w-]+)'"
    dim_pat    = "dim\s*\=\s*'(\([\d,]+\))'"
    format_pat = "format\s*\=\s*'(\w+)'" 

    # Loop through, matching on each line
    cfile = open(coldef_file)
    for line in cfile.readlines():
        unit = name = dim = format = None
        name_match = re.search(name_pat, line)
        if name_match:
            name = name_match.group(1)
             
            format_match = re.search(format_pat, line)
            dim_match    = re.search(dim_pat, line)
            unit_match   = re.search(unit_pat, line)

            if unit_match:   unit = unit_match.group(1)
            if dim_match:    dim  = dim_match.group(1)
                        
            if format_match: 
                fits_fmt = format_match.group(1)
                zarr     = generateZeros(num_rows, fits_fmt, dim)

            
            # Append the column to the column list
            cols.append(pf.Column(name=name, format=fits_fmt, unit=unit, dim=dim, array=zarr))
    
    # Now we have made a list of columns, we can make a new table
    coldefs = pf.ColDefs(cols)
    #print coldefs
    tbhdu   = pf.new_table(coldefs)
    
    # If that all worked, we can populate with the final header values
    cards = generateCards(header_file)
    
    for card in cards:
        if card.key == 'COMMENT':
            pass
            tbhdu.header.add_comment(card.value)
        elif card.key == 'HISTORY':
            pass
            tbhdu.header.add_history(card.value)
        else:
            tbhdu.header.set(card.key, card.value, card.comment)
    
    return tbhdu

def generateDataHDU(input_file, 
                    header_file='hipsr_core/header_dataHDU.txt',
                    coldef_file='hipsr_core/coldefs_dataHDU.txt'):
    """ Generate a new data table based upon an input file
    
    Parameters
    ----------
    header_file: str
        Path to the header file. Defaults to 'header_dataHDU.txt'
    coldef_file: str
        Path to the file containing column definitions.
        Defaults to 'coldefs_dataHDU.txt'
    input_file: str
        String to the input file to grab data from. Defaults to none.
    
    """
    
    sd_in      = pf.open(input_file)
    sd_data    = sd_in[1].data
    num_rows   = sd_data.shape[0]
    
    cols = []
    
    # The column definitions are loaded from an external file, which is
    # parsed line-by-line, using regular experssions.
    
    unit_pat   = "unit\s*\=\s*'([\w/%]+)'"
    name_pat   = "name\s*\=\s*'([\w-]+)'"
    dim_pat    = "dim\s*\=\s*'(\([\d,]+\))'"
    format_pat = "format\s*\=\s*'(\w+)'" 
    
    # Loop through, matching on each line
    cfile = open(coldef_file)
    for line in cfile.readlines():
        unit = name = dim = format = None
        name_match = re.search(name_pat, line)
        if name_match:
            name = name_match.group(1)
             
            format_match = re.search(format_pat, line)
            dim_match    = re.search(dim_pat, line)
            unit_match   = re.search(unit_pat, line)
    
            if unit_match:   
                unit = unit_match.group(1)
            
            
            if dim_match:    
                dim       = dim_match.group(1)
            
            arr_shape = sd_data[name].shape
                    
            if format_match: 
                fits_fmt = format_match.group(1)
                zarr=None

                try:
                    if name == 'DATA' or name == 'FLAGGED':
                        np_dtype, data_len, data_fmt = formatLookup(fits_fmt)
                        print name, " no data"
                    else:
                        # Data array must be flattened (e.g. (2,2) -> 4)
                        np_dtype, data_len, data_fmt = formatLookup(fits_fmt)
                        if data_len > 1 and data_fmt != 'str_':
                            z_shape = (sd_data[name].shape[0], data_len)
                        else:
                             z_shape = sd_data[name].shape
                        #print name, z_shape, sd_data[name].shape
                        zarr     = sd_data[name].reshape(z_shape)
                        
                except:
                    print "Error with %s"%name
            
                # Append the column to the column list
                cols.append(pf.Column(name=name, format=fits_fmt, unit=unit, dim=dim, array=zarr))
    
    # Now we have made a list of columns, we can make a new table
    #print cols
    coldefs = pf.ColDefs(cols)
    #print coldefs
    tbhdu   = pf.new_table(coldefs)
    
    # If that all worked, we can populate with the final header values
    cards = generateCards(header_file)
    
    for card in cards:
        if card.keyword == 'COMMENT':
            pass
            tbhdu.header.add_comment(card.value)
        elif card.keyword == 'HISTORY':
            pass
            tbhdu.header.add_history(card.value)
        else:
            tbhdu.header.set(card.keyword, card.value, card.comment)
    
    return tbhdu


def timestamp2dt(timestamp):
    """ Convert timestamp to date and time for SD-FITS """
    
    dt = datetime.utcfromtimestamp(timestamp)
    
    date = dt.strftime("%Y-%m-%d")
    # TODO: Check this is correct
    time = dt.hour * 3600 + dt.minute * 60 + dt.second + dt.microsecond * 1e-6
    return (date, time)
    
def findMatchingTimestamps(h5, sd):
    """ Compare HIPSR and SD-FITS timestamps

    Returns an array of indexes corresponding to the best matches.
    These indexes are for the HDF5 file, and correspond to data rows.
    Throws and error if the t_diff is over 2s (one integration).

    """

    sd_data = sd[1].data
    hp_ts = h5.root.raw_data.beam_01.col("timestamp")
    hp_dts = np.array([datetime.utcfromtimestamp(ts) for ts in hp_ts])
    
    t_idx = []
    for row in range(len(sd_data['TIME'])):
            
        utime = sd_data['TIME'][row]
        udate = sd_data['DATE-OBS'][row]
        
        # From string to datetime obj
        d_d  = datetime.strptime(udate, "%Y-%m-%d")
        # from datetime obj to timestamp
        d_ts = time.mktime(d_d.utctimetuple())
        # date + time into timestamp
        dt_ts = d_ts + utime
        # Creating overall timestamp
        dt = datetime.utcfromtimestamp(dt_ts)

        # TODO: Figure out why an 8hr offset is required??!
        t_diffs = hp_ts - dt_ts - 60*60*8
        idx = np.argmin(np.abs(t_diffs))

        if np.abs(t_diffs[idx]) >= 1.1:
            print "Warning: large t_diff: ",
            print idx, t_diffs[idx]
        t_idx.append(idx)

        if np.abs(t_diffs[idx]) >= 2:
            print "ERROR: Time difference between two files is too large. No match found."
            raise

    t_idx = np.array(t_idx)

    return t_idx

def generateBlankSDFits(num_rows, header_primary, header_tbl, coldef_file):
    """ Generate a blank SD-FITS file

    This function returns a blank SD-FITS file, with num_rows in the binary table.
    It generates all the required columns, then fills them with blank data (zeros).

    Parameters
    ----------
    num_rows: int
        The number of rows in the binary table.
    header_primary: str
            Path to the primaryHDU header file. Defaults to 'header_primaryHDU.txt'
    header_data: str
        Path to the binary table header file. Defaults to 'header_dataHDU.txt'
    coldef_file: str
        Path to the file containing column definitions for the binary table.
        Defaults to 'coldefs_dataHDU.txt'

    """
    print header_primary
    prhdu = generatePrimaryHDU(header_primary)
    tbhdu = generateBlankDataHDU(num_rows, header_tbl, coldef_file)
    
    # Insert creation date
    time_tuple = time.gmtime()
    date_str = time.strftime("%Y-%m-%dT%H:%M:%S", time_tuple)
    prhdu.header['DATE'] = date_str
        
    hdulist = pf.HDUList([prhdu, tbhdu])
    
    return hdulist

def generateSDFitsFromMbcorr(input_file, header_primary, header_tbl, coldef_file):
    """ Generate a blank SD-FITS file

    This function returns a blank SD-FITS file, with num_rows in the binary table.
    It generates all the required columns, then fills them with blank data (zeros).

    Parameters
    ----------
    num_rows: int
        The number of rows in the binary table.
    header_primary: str
            Path to the primaryHDU header file. Defaults to 'header_primaryHDU.txt'
    header_data: str
        Path to the binary table header file. Defaults to 'header_dataHDU.txt'
    coldef_file: str
        Path to the file containing column definitions for the binary table.
        Defaults to 'coldefs_dataHDU.txt'

    """

    prhdu = generatePrimaryHDU(header_primary)
    tbhdu = generateDataHDU(input_file, header_tbl, coldef_file)
    hdulist = pf.HDUList([prhdu, tbhdu])

    return hdulist

def generateSDFitsFromHipsr(filename_in, path_in, filename_out, path_out, write_stokes=False):
    """ Generate an SD-FITS file from a hipsr5 file """
    
    # Open h5 file
    print "\nOpening files"
    print "-------------"
    h5file = os.path.join(path_in, filename_in)
    out_file = os.path.join(path_out, filename_out)
    h6 = Hipsr6(h5file)
   
    
    num_acc  = h6.beams[0].shape[0] 
    num_rows = num_acc * 13
    
    if num_acc == 0:
        print "No data in %s. Skipping."%h5file
        return -1
    
    print "Input file: %s"%h6.h5.filename
    print "No accumulations: %s, no rows: %s"%(num_acc, num_rows)
    print h6
    abspath = os.path.abspath( __file__ ).rstrip('sdfits.py')
    diode_cal_file_x  = "%s/diode_jy_x.cal"%abspath
    diode_cal_file_y  = "%s/diode_jy_y.cal"%abspath
    cal_factor_file   = "%s/cal_factor.cal"%abspath
    
    cf = np.fromfile(cal_factor_file)
    diode_temps_x, diode_temps_y = loadDiodeTemp(h6, diode_cal_file_x, diode_cal_file_y)

    # We now need to generate a blank SD-FITS file, with the same number of rows
    print "\nGenerating blank SD-FITS file with %i rows..."%num_rows
    
    if write_stokes:
        print "Stokes flag found - writing I,Q,U,V"
        header_primary='hipsr_core/header_primaryHDU.txt'
        header_tbl='hipsr_core/header_dataHDU_stokes.txt'
        coldef_file='hipsr_core/coldefs_dataHDU_stokes.txt'
    else:
        print "Writing XX, YY"
        header_primary='hipsr_core/header_primaryHDU.txt'
        header_tbl='hipsr_core/header_dataHDU.txt'
        coldef_file='hipsr_core/coldefs_dataHDU.txt'
    
    hdulist = generateBlankSDFits(num_rows, header_primary, header_tbl, coldef_file)
    print hdulist.info()
    
    # Next, we copy over observation data    
    print "Filling new SD-FITS with HIPSR data..."
    
    pointing = h6.tb_pointing.cols
    obs      = h6.tb_observation.cols
    sdtab    = hdulist[1].data
    sdhead   = hdulist[1].header
    
    freqs     = h6.freqs
    freqs_cal = h6.freqs_cal
    
    # Fill in header values
    sdhead["OBSERVER"] = obs.observer[0]
    sdhead["PROJID"]   = obs.project_id[0]
    
    # Fill in common values
    print "Filling in common values... ",
    sdtab["SCAN"][:]     = 1
    sdtab["EXPOSURE"][:] = obs.acc_len[0]
    sdtab["OBJECT"][:]   = pointing.source[0]
    sdtab["OBJ-RA"][:]   = pointing.ra[0]
    sdtab["OBJ-DEC"][:]  = pointing.dec[0]
    sdtab["RESTFRQ"][:]  = obs.frequency[0]    
    sdtab["FREQRES"][:]  = np.abs(obs.bandwidth[0])*1e6 / 8192
    sdtab["BANDWID"][:]  = np.abs(obs.bandwidth[0])
    sdtab["CRPIX1"][:]   = 4095
    sdtab["CRVAL1"][:]   = obs.frequency[0] * 1e6
    sdtab["CDELT1"][:]   = np.abs(obs.bandwidth[0])*1e6 / 8192
    sdtab["FLAGGED"][:]  = 0
    sdtab["SCANRATE"][:] = obs.scan_rate[0]


    # TCS INFO
    sdtab["OBSMODE"][:]  = obs.obs_mode[0] 
    sdtab["IF"][:]       = 1
    print "OK."
    
    row_sd   = 0
    cycle_id = 0
    scaling = 2**22 # Divide through to change 32-bit to 
    
    flipped = False
    if obs.bandwidth[0] < 0:
        flipped = True
    
    print "Filling in unique values... "
    scan_pointing_len = h6.tb_scan_pointing.shape[0]
    
    for row_h5 in range(num_acc):
        cycle_id += 1 # Starts at 1 in SD-FITS file
        for beam in h6.h5.root.raw_data:
            LinePrint("%i of %i"%(row_sd, num_rows))
            
            if cycle_id <= scan_pointing_len:
                raj_id = "mb%s_raj"%beam.name.lstrip('beam_')
                dcj_id = "mb%s_dcj"%beam.name.lstrip('beam_')
                
                sdtab["CYCLE"][row_sd]   = cycle_id
                beam_id = int(beam.name.lstrip('beam_'))
                
                # Fix beam mapping (remove after fixing mapping)
                sdtab["BEAM"][row_sd]     = beam_id
                
                sdtab["CRVAL3"][row_sd]   = h6.tb_scan_pointing.col(raj_id)[cycle_id-1]
                sdtab["CRVAL4"][row_sd]   = h6.tb_scan_pointing.col(dcj_id)[cycle_id-1]
                sdtab["AZIMUTH"][row_sd]  = h6.tb_scan_pointing.col("azimuth")[cycle_id-1]
                sdtab["ELEVATIO"][row_sd] = h6.tb_scan_pointing.col("elevation")[cycle_id-1]
                

                try:
                    timestamp  = beam.cols.timestamp[row_h5]
                    date_obs, time = timestamp2dt(timestamp)
                    sdtab["DATE-OBS"][row_sd] = date_obs
                    sdtab["TIME"][row_sd]     = time
                    
                    # Compute T_sys for each beam
                    T_d_x = diode_temps_x[beam_id-1]
                    T_d_y = diode_temps_y[beam_id-1]
                    T_sys_x, T_sys_y = computeTsys(beam, row_h5, T_d_x, T_d_y)
                
                    #print T_sys_x, T_sys_y
                    sdtab["TSYS"][row_sd] = (T_sys_x, T_sys_y)
                    
                    
                    
                    if write_stokes:
                        # Currently not calibrating!
                        xx = beam.cols.xx[row_h5].astype('float32') 
                        yy = beam.cols.yy[row_h5].astype('float32') 
                        re_xy = beam.cols.re_xy[row_h5].astype('float32') 
                        im_xy = beam.cols.im_xy[row_h5].astype('float32')
                        
                        # Blank DC bin
                        xx[0], yy[0],re_xy[0], im_xy[0] = np.zeros(4)
                        if flipped:
                            xx, yy, re_xy, im_xy = xx[::-1], yy[::-1], re_xy[::-1], im_xy[::-1]
                        
                        xx = xx / np.average(extractMid(xx)) * T_sys_x
                        yy = yy / np.average(extractMid(yy)) * T_sys_y
                        re_xy = re_xy / np.average(extractMid(re_xy)) * (T_sys_x + T_sys_y)/2
                        im_xy = im_xy / np.average(extractMid(im_xy)) * (T_sys_x + T_sys_y)/2
                        
                        # Ettore tells me Parkes uses this definition
                        # i.e. that I is the average of xx + yy
                        ii = (xx + yy) / 2
                        qq = (xx - yy) / 2
                        uu = re_xy
                        vv = im_xy
                        
                        # Form one data vector
                        data1 = np.append(ii, qq)
                        data2 = np.append(uu, vv)
                        data  = np.append(data1, data2)
                        data  = data.reshape([1,1,4,8192])
                    else:
                        
                        xx = beam.cols.xx[row_h5].astype('float32') 
                        yy = beam.cols.yy[row_h5].astype('float32') 
                        
                        # Blank DC bin
                        xx[0], yy[0] = 0,0
                        if flipped:
                            xx, yy = xx[::-1], yy[::-1]                           
                        
                        #print "cal factor: %2.3f"%cf
                        #print "Diode temp: %s"%T_d
                        #xx, yy = applyCal(beam, row_h5, freqs, freqs_cal, cf, T_d_x, T_d_y)
                        
                        xx = xx / np.average(extractMid(xx)) * T_sys_x
                        yy = yy / np.average(extractMid(yy)) * T_sys_y

                        data = np.append(xx, yy)
                        data = data.reshape([1,1,2,8192]) 
                    
                    sdtab["DATA"][row_sd] = data
                    

                    
                except:
                    if beam.name != 'beam_02':
                        print "\nWARNING: missing row in %s"%beam.name
                        print "Current index: %i"%row_h5
                        print "Row length: %i"%beam.shape[0]
                        raise

                row_sd += 1
            else:
                print "WARNING: scan_pointing table is not complete."
                print "%s table length: %i"%(beam.name, beam.shape[0])
                print "scan_pointing table length: %i"%scan_pointing_len
    
    h6.h5.close()
    
    if os.path.exists(out_file):
        print "\nInfo: File exists, deleting..."
        os.remove(out_file)

    print "\nInfo: Saving to file"
    hdulist.writeto(out_file)
    hdulist.close()

