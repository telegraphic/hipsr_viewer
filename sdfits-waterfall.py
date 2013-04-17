#!/usr/bin/env python
"""
sdfits_viewer.py
----------------

Viewer for Parkes multibeam SD-FITS files.

"""

import pyfits, numpy, pylab, sys, os, re
import pylab as plt

path = os.getcwd()

# Regular expression to match SDFITS extension
regex = '([0-9A-Za-z-_]+).sdfits'
filelist = []

for filename in os.listdir(path):
    match = re.search(regex, filename)
    if match:
        filelist.append(filename)
    
for iname in filelist:
    name = os.path.join(path, iname)
    print "File: %s"%iname
    imname = 'tmp'
        
    fits     = pyfits.open(name,mode='update')
    header   = fits[1].header
    data     = fits[1].data
    data_or  = fits[1].data['DATA']
    flagg_or = fits[1].data['FLAGGED']
    
    # Frequency-like axis values
    ref_pix   = fits[1].data['CRPIX1']  # Reference pixel
    ref_val   = fits[1].data['CRVAL1']  # Value at reference pixel (in Hz)
    ref_delt  = fits[1].data['CDELT1']  # Delta between pixels
    num_pix   = data_or.shape[-1]       # Last axis of data array (multidimensional)
    
    data_name = fits[1].header.get('TTYPE25')   # Might be better to look this up than assume it's col 25
    data_unit = fits[1].header.get('TUNIT25')
    
    freq_type = fits[1].header.get('CTYPE1')
    freq_unit = 'MHz'                           # Might be better to check this
    
    freqs = (numpy.arange(0,num_pix,1) * ref_delt[0] + ( ref_val[0] - ref_pix[0] * ref_delt[0] ) ) / 1e6 
    
    data_or_shape = numpy.shape(data_or)
    
    
        
    print "Data array shape: ",
    print data_or_shape
    if len(data_or_shape)>2:
    	print "warning: N-dimensional data array. Attempting to reshape..."
    	try: 
    		data_or = data_or.reshape(data_or_shape[0],data_or_shape[-1]*data_or_shape[-2])
    		
    	except:
    		print "couldn't figure out how to reshape."
    		break
    
        
    #plt.plot(freqs, numpy.average(data_or[:,::-2],axis=0))
    pltdata = data_or[:,::-2] / numpy.average(data_or[:,::-2],axis=0)
    stdev   = numpy.std(pltdata)
    avg     = numpy.average(data_or)
    numpy.putmask(pltdata, pltdata >= 2*stdev, numpy.nan)
    
    fig = plt.figure(figsize=(4*3,3*3))
    ax  = fig.add_subplot(111)
    plt.imshow(pltdata[::13])
    plt.colorbar()
    print pltdata.shape[1]/pltdata.shape[0]
    ax.set_aspect(pltdata.shape[1]/pltdata.shape[0])
    plt.xticks(numpy.arange(0,8192)[::8192/10], ["%2.2f"%fr for fr in freqs[::8192/10]])
    #plt.ylabel('%s [%s]'%(data_name, data_unit))
    #plt.xlabel('%s [%s]'%(freq_type, freq_unit))
    plt.show()