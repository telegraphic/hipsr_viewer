#!/usr/bin/env python
# encoding: utf-8
"""
config.py
============

This is the main configuration file for the HIPSR wideband spectrometer server script.

"""

###########
# v1 - Analytical Albatross (Oct 2012)
# v2 - Ballistic Bandicoot  (Apr 2013)

__author__    = "Danny Price"
__email__     = "danny.price@astro.ox.ac.uk" 
__license__   = "GNU GPL"
__version__   = "2.0 - Ballistic Bandicoot"

project_id = 'PXXX'

data_dir    = '/data/hipsr'

tcs_port     = 59011             # BPSR is 50910
tcs_server   = '130.155.182.73'  # hipsr-srv0 eth0

plotter_port = 59012
plotter_host = '127.0.0.1'  # Plotting on localhost

katcp_port = 7147

#boffile    = 'hipsr_parkes_400_2012_Mar_15_2209.bof'
boffile     = 'hispec_8192_v2_2013_Apr_18_1730.bof'

reprogram   = True
reconfigure = True

###############
# Roach to beam mappings
# Last checked on 17th April 2013
# Note: Curtis is noise diode master ROACH

roachlist = {
    "drake"    : "beam_01",
    "hendrix"  : "beam_02",
    "reznor"   : "beam_03",
    "keenan"   : "beam_04",
    "mackaye"  : "beam_05", 
    "albarn"   : "beam_06", 
    "willis"   : "beam_07", 
    "waits"    : "beam_08", 
    "yorke"    : "beam_09",
    "cobain"   : "beam_10",
    "patton"   : "beam_11",
    "barrett"  : "beam_12",
    "curtis"   : "beam_13"
    }


############
# FPGA setup
# Last modified 17th April 2013
# FPGA clock:         200 MHz
# Dump rate:          2 s
# Diode switching:    128 Hz
# Vector accumulator: 4096 long (8192/2, as even / odd are parallel streams)
# 200e6 / 4096 = 48828.125

n_sec             = 2
acc_len           = 48828*n_sec
nar_acc_len       = acc_len * 4096 / 8
n_cycles_per_dump = 128*n_sec
sq_wv_period      = 8 * nar_acc_len / n_cycles_per_dump 

fpga_config = {
    "acc_len"               : acc_len,
    "fft_shift"             : -1,
    "quant_xx_gain"         : -1,
    "quant_yy_gain"         : -1,
    "quant_xy_gain"         : -1,
    "mux_sel"               : 0,
    "nar_sq_wave_period"    : sq_wv_period, 
    "nar_quant_yy_gain"     : -1,
    "nar_quant_xx_gain"     : -1,
    "nar_fft_shift"         : -1,
    "nar_acc_len"           : nar_acc_len
}
        
tcs_regex_esc = '\\n' # Escape for message from TCS
#tcs_regex_esc = '' # Escape for dummy TCS data
