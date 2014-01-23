#!/usr/bin/env python
"""
hipsr-sdview.py
==================

SD-FITS viewer designed for HIPSR data.

"""



# Imports
import sys
from optparse import OptionParser

import lib.config as config
from lib.sdfits import *

# Python metadata
__version__  = config.__version__
__author__   = config.__author__
__email__    = config.__email__
__license__  = config.__license__
__modified__ = datetime.fromtimestamp(os.path.getmtime(os.path.abspath( __file__ )))

USES_PYSIDE = False

import matplotlib

if matplotlib.__version__ == '0.99.3':
    print "Error: your matplotlib version is too old to run this. Please upgrade."
    exit()
else:
    matplotlib.use('Qt4Agg')
    if USES_PYSIDE:
        matplotlib.rcParams['backend.qt4']='PySide'
    else:
        matplotlib.rcParams['backend.qt4']='PyQt4'
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
    from matplotlib.figure import Figure
try:
    import matplotlib.pylab as plt
except:
    print "Error: cannot load Pylab. Check your matplotlib install."
    exit()
    
try:
    #import lib.qt_compat as qt_compat
    #QtGui = qt_compat.import_module("QtGui")
    #QtCore = qt_compat.QtCore
    
    
    
    import PyQt4
    from PyQt4 import QtGui, QtCore
except:
    print "Error: cannot load PySide or PyQt4. Please check your install."
    exit()
    
try:    
    import numpy as np
except:
    print "Error: cannot load Numpy. Please check your install."
    exit()

try:    
    import pyfits as pf
except:
    print "Error: cannot load PyFITS. Please check your install."
    exit()

from lib.printers import LinePrint




class HipsrGui(QtGui.QWidget):
    """ HIPSR GUI class
    
    A Qt4 Widget that uses matplotlib to display antenna configurations
    
    Parameters
    ----------
    filename: str
        Name of file to open. Defaults to blank, in which case no file is opened.
    """
    def __init__(self, filename=''):
        super(HipsrGui, self).__init__()
        
        # Initialize user interface
        self.filename = filename
        self.current_row = 0
        self.initUI(width=800, height=800)
       

    def initUI(self, width=1024, height=768):
        """ Initialize the User Interface 
        
        Parameters
        ----------
        width: int
            width of the UI, in pixels. Defaults to 1024px
        height: int
            height of the UI, in pixels. Defaults to 768px
        """
        
        self.main_frame = QtGui.QWidget()   
        self.setWindowIcon(QtGui.QIcon('lib/icon.jpg'))  
        #self.gen_gui = generateGui()
        
        # Create buttons/widgets
        self.but_open     = QtGui.QPushButton("Open")
        self.but_open.clicked.connect(self.onButOpen)
        self.lab_info     = QtGui.QLabel(" ")
        self.slider       = QtGui.QSlider()
        self.slider.valueChanged.connect(self.updatePlot)
        
        self.spinner     = QtGui.QSpinBox()
        self.spinner.setRange(0, 1)
        self.spinner.setSingleStep(1)
        self.spinner.valueChanged.connect(self.updateSlider)
        
        self.fits_data = np.zeros([8192])
        
        # Create plots
        self.sp_fig, self.sp_ax = self.createSpectrumPlot()
        
        
        # generate the canvas to display the plot
        self.sp_canvas = FigureCanvas(self.sp_fig)
        self.mpl_toolbar = NavigationToolbar(self.sp_canvas, self.main_frame)
        
        # Widget layout
        layout = QtGui.QVBoxLayout()
        
        h_layout = QtGui.QHBoxLayout()
        h_layout.addWidget(self.sp_canvas)
        h_layout.addWidget(self.slider)
        layout.addLayout(h_layout)
        h_layout = QtGui.QHBoxLayout()
        h_layout.addStretch(1)
        h_layout.addWidget(self.spinner)
        layout.addLayout(h_layout)
        layout.addWidget(self.mpl_toolbar)
        
        bbox = QtGui.QHBoxLayout()
        bbox.addWidget(self.lab_info)
        bbox.addStretch(1)
        bbox.addWidget(self.but_open)
        layout.addLayout(bbox)

        self.setLayout(layout)    
        
        self.setGeometry(300, 300, width, height)
        self.setWindowTitle('HIPSR SD-FITS viewer')    
        self.show()
        
        # Load file if command line argument is passed
        if self.filename != '':
            try:
                self.openSdFits(self.filename)
                self.updatePlot()
            except:
                print "Error: cannot open %s"%self.filename
                #raise
        
    
    def openSdFits(self,filename):
        """Opens an SD-FITS file and loads it into memory"""

        fits     = pf.open(filename)
       
        self.fits_filename     = fits.filename()
        self.fits_time    = fits[1].data['TIME']
        self.fits_date    = fits[1].data['DATE-OBS']
        self.fits_data    = fits[1].data['DATA']
        self.fits_beam    = fits[1].data['BEAM']
        self.fits_flagged = fits[1].data['FLAGGED']
        self.fits_tsys    = fits[1].data['TSYS']
        self.fits_header  = fits[1].header
        self.fits_shape   = np.shape(self.fits_data)
        
        # Frequency-like axis values

        try:
            ref_pix   = fits[1].header['CRPIX1']   # Freq Reference pixel
            ref_val   = fits[1].header['CRVAL1']   # Value at reference pixel (in Hz)
            ref_delt  = fits[1].header['CDELT1']   # Delta between pixels
        except KeyError:
            ref_pix   = fits[1].data['CRPIX1'][0]  # Freq Reference pixel
            ref_val   = fits[1].data['CRVAL1'][0]  # Value at reference pixel (in Hz)
            ref_delt  = fits[1].data['CDELT1'][0]  # Delta between pixels
        num_pix   = self.fits_data.shape[-1]# Last axis of data array (multidimensional)
        self.fits_npol = self.fits_data.shape[-2]
        # Might be better to look this up than assume it's col 25
        self.fits_data_name = fits[1].header.get('TTYPE25')   
        self.fits_data_unit = fits[1].header.get('TUNIT25')


        self.fits_freqs = (np.arange(0,num_pix,1) * ref_delt + (ref_val - ref_pix * ref_delt) ) / 1e6
        self.fits_freq_type = fits[1].header.get('CTYPE1')
        self.fits_freq_unit = 'MHz' # TODO: check this
        
        # Stokes codes are the CRVAL2 in the header
        stokes_codes = {
            1  : 'Stokes I',
            2  : 'Stokes Q',
            3  : 'Stokes U',
            4  : 'Stokes V',
            -1 : 'RR',
            -2 : 'LL',
            -3 : 'RL',
            -4 : 'LR',
            -5 : 'XX',
            -6 : 'YY',
            -7 : 'XY',
            -8 : 'YX'
        }
        
        try:
            self.xlabel = stokes_codes[int(fits[1].header['CRVAL2'])]
            self.ylabel = stokes_codes[int(fits[1].header['CRVAL2']) + int(fits[1].header['CDELT2']) ] 
        except:
            print "Warning: unknown stokes code: %s"%fits[1].header['CRVAL2']
            self.xlabel = 'Unknown X'       
            self.ylabel = 'Unknown Y' 
        
        print fits.info()        

        print "File date:  %s"%fits[0].header['DATE']
        print "Telescope:  %s"%self.fits_header['TELESCOP']
        print "Observer:   %s"%self.fits_header['OBSERVER']
        print "Project ID: %s"%self.fits_header['PROJID']               
        print "Data type:  %s %s"%(self.xlabel, self.ylabel)



        
        # Reset slider
        self.current_row = 0
        self.slider.setValue(0)
        self.slider.setMaximum(self.fits_shape[0]-1)
        self.spinner.setValue(0)
        self.spinner.setRange(0, self.fits_shape[0]-1)
        
    def createSpectrumPlot(self):
          """ Creates a single pylab plot for displaying a spectrum """

          fig = plt.figure(figsize=(8,6),dpi=80)
          fig.set_facecolor('#ededed')
          
          # Format plot
          ax = plt.subplot(111)

          ax.set_xlabel("Frequency [MHz]")
          ax.set_ylabel("Power [-]")
          ax.set_xlim(1000,1500)
          ax.set_ylim(0,1000)
        
          fig.canvas.draw()
      
          return fig, ax
    
    def updateSlider(self):
        """ Link spinner to slider """
        self.slider.setValue(self.spinner.value())
        
    
    def updatePlot(self):
        """ Updates the antenna config plot"""
                
        self.sp_ax.clear()
        self.sp_fig.clear()
        
        row = self.slider.value()
        self.spinner.setValue(row)
        
        data_x    = self.fits_data[row,0,0,0,:]
        flagged_x = self.fits_flagged[row,0,0,0,:]

        if self.fits_npol >= 2:
            data_y    = self.fits_data[row,0,0,1,:]
            flagged_y = self.fits_flagged[row,0,0,1,:]
        freqs     = self.fits_freqs
        tsys      = self.fits_tsys[row,:]
        
        #plt.plot(freqs, data_x)
        
        colors = ['#333333','#CC0000',"#5fafc6","#6dae81","#9d81ba","#cd4a4a"]
        if self.xlabel == 'Stokes I':
            c1 = colors[3]
            c2 = colors[2]
        else:
            c1 = colors[0]
            c2 = colors[1]
        
        self.sp_plot =  plt.plot(freqs[flagged_x == 0], data_x[flagged_x == 0], color=c1, 
                        label='%s [%2.1f Jy]'%(self.xlabel, tsys[0]))  
        if self.fits_npol >= 2:
            self.sp_plot2 = plt.plot(freqs[flagged_y == 0], data_y[flagged_y == 0], color=c2,
                        label='%s [%2.1f Jy]'%(self.ylabel, tsys[1]))         
        
        plt.ylabel('%s [%s]'%(self.fits_data_name, self.fits_data_unit))
        plt.xlabel('%s [%s]'%(self.fits_freq_type, self.fits_freq_unit))
        plt.title('Beam %s: %s %s'%(self.fits_beam[row], self.fits_date[row], self.fits_time[row]))
        plt.xlim(np.min(freqs[flagged_x == 0]), np.max(freqs[flagged_x == 0]))
    
        
        plt.legend()
        #plt.show()
        

        self.sp_fig.canvas.draw()
        #self.sp_canvas.draw()
              
        self.lab_info.setText(self.fits_filename)       
                      
    def onButOpen(self):
        """ Button action: Open station file """
        self.file_dialog    = QtGui.QFileDialog()
        fileparts = self.file_dialog.getOpenFileName()
        
        # Compatibility check
        if isinstance(fileparts, str) or isinstance(fileparts, unicode):
            filename = fileparts
        
        elif not USES_PYSIDE:
            if isinstance(fileparts, PyQt4.QtCore.QString):
                filename = str(fileparts)
            else:
                filename = fileparts
        else:
            filename = fileparts[0]
        
        self.openSdFits(filename)
        self.updatePlot()

    def onButStats(self):
        """ Button action: Open statistics """
        pass


def main():
    
    # Basic option parsing 
    p = OptionParser()
    p.set_usage('fits_pattern_viewer.py [filename] [options]')
    p.set_description(__doc__)
    (options, args) = p.parse_args()

    print "Starting HIPSR SD-FITS viewer..."
    global main_gui
    app = QtGui.QApplication(sys.argv)
    
    try:
        filename = args[0]
        main_gui = HipsrGui(filename)
    except:
        main_gui = HipsrGui()
    app.exec_()
    sys.exit()    

if __name__ == '__main__':
    main()
