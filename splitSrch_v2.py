#!/usr/bin/env python3
import math
import numpy as np

import matplotlib.pyplot as plt
#from matplotlib import rc
#from matplotlib.patches import Circle

#import subprocess
import argparse
#import glob

#import pyfits
from astropy.io import fits

#rc('text', usetex=True)
# read in parameters

parser = argparse.ArgumentParser(description='Read PSRFITS format search mode data')
parser.add_argument('-f',  '--input_file',  metavar='Input file name',  nargs='+', required=True, help='Input file name')
parser.add_argument('-o',  '--output_file', metavar='Output file name', nargs='+', required=True, help='Output file name')
parser.add_argument('-sub0',  '--subband_start', metavar='Starting subband', nargs='+', required=True, help='Starting subband')
parser.add_argument('-sub1',  '--subband_end', metavar='Ending subband', nargs='+', required=True, help='Ending subband')
#parser.add_argument('-dm', '--chn_dm',      metavar='Channel DM',       nargs='+', required=True, help='Channel DM')

args = parser.parse_args()
sub_start = int(args.subband_start[0])
sub_end = int(args.subband_end[0])
#dm = float(args.chn_dm[0])

############### read in data and parameters #################
in_file = args.input_file[0]
hdulist = fits.open(in_file)
tbdata = hdulist['SUBINT'].data
obsbw = hdulist['PRIMARY'].header['OBSBW']
nbits = int(hdulist['SUBINT'].header['NBITS'])
nchan = int(hdulist['SUBINT'].header['NCHAN'])
nsblk = int(hdulist['SUBINT'].header['NSBLK'])
nsub =  int(hdulist['SUBINT'].header['NAXIS2'])
npol =  int(hdulist['SUBINT'].header['NPOL'])
nstot = int(hdulist['SUBINT'].header['NSTOT'])
print (nchan,npol,nsub)

# read in headers from the first file
hdu_primary_header = hdulist['PRIMARY'].header
hdu_subint_header = hdulist['SUBINT'].header

#dat_freq = tbdata['DAT_FREQ']
#dat_wts = tbdata['DAT_WTS']
dat_freq = np.reshape(tbdata['DAT_FREQ'], (nsub, nchan))
#dat_freq = np.reshape(tbdata['DAT_FREQ'], (nsub, npol, nchan))
dat_wts = np.reshape(tbdata['DAT_WTS'], (nsub, nchan))
#dat_wts = np.reshape(tbdata['DAT_WTS'], (nsub, npol, nchan))
#dat_off = tbdata['DAT_OFFS']
#dat_scl = tbdata['DAT_SCL']
dat_off = np.reshape(tbdata['DAT_OFFS'], (nsub, npol, nchan))
dat_scl = np.reshape(tbdata['DAT_SCL'], (nsub, npol, nchan))
#data = tbdata['DATA']
#data = np.reshape(tbdata['DATA'], (nsub*nsblk, nchan/(8/nbits))).astype(int)
#data = np.reshape(tbdata['DATA'], (nsub*nsblk, nchan/(8/nbits))).astype(np.int8)
#data = np.reshape(tbdata['DATA'], (nsub*nsblk, nchan/(8/nbits)))
data = np.reshape(tbdata['DATA'], (nsub*nsblk, npol, int(nchan/(8/nbits))))
print (dat_freq.shape, dat_wts.shape, dat_off.shape, dat_scl.shape, data.shape)
print (nsub, nsblk, nbits)

# copy columns from the first file
indexval = tbdata['INDEXVAL']
tsubint = tbdata['TSUBINT']
offssub = tbdata['OFFS_SUB']
auxdm = tbdata['AUX_DM']
auxrm = tbdata['AUX_RM']

############### determine channels to include #################
sub_nchan = int(nchan/26)
sub_nstot = int(nstot/26)

mask = np.zeros(nchan, dtype=bool)
mask[sub_start*sub_nchan:(sub_end+1)*sub_nchan] = True
#mask = np.zeros((npol,nchan), dtype=bool)
#mask[:,sub_start*sub_nchan:(sub_end+1)*sub_nchan] = True
#print mask.shape
#mask = np.squeeze(mask.flatten())
#print mask.shape
#print mask

nchan_new = (sub_end - sub_start + 1)*sub_nchan
nstot_new = (sub_end - sub_start + 1)*sub_nstot

mask_dat = np.zeros(int(nchan/(8/nbits)), dtype=bool)
mask_dat[sub_start*int(sub_nchan/(8/nbits)):(sub_end+1)*int(sub_nchan/(8/nbits))] = True
#print sub_nchan, sub_nstot
#print nchan, nstot
#print nchan_new, nstot_new

############### create some new Columns #######################
col_index = fits.Column(name='INDEXVAL', format='1D', array=indexval)
col_tsubint = fits.Column(name='TSUBINT', format='1D', array=tsubint, unit='s')
col_offssub = fits.Column(name='OFFS_SUB', format='1D', array=offssub, unit='s')
col_auxdm = fits.Column(name='AUX_DM', format='1D', array=auxdm, unit='CM-3')
col_auxrm = fits.Column(name='AUX_RM', format='1D', array=auxrm, unit='RAD')

#col_freq = fits.Column(name='DAT_FREQ', format='{0}D'.format(nchan_new*npol), array=dat_freq[:,:,mask])
#col_wts = fits.Column(name='DAT_WTS', format='{0}E'.format(nchan_new*npol), array=dat_wts[:,:,mask])
col_freq = fits.Column(name='DAT_FREQ', format='{0}D'.format(nchan_new), array=dat_freq[:,mask])
col_wts = fits.Column(name='DAT_WTS', format='{0}E'.format(nchan_new), array=dat_wts[:,mask])
col_off = fits.Column(name='DAT_OFFS', format='{0}E'.format(nchan_new*npol), array=dat_off[:,:,mask])
col_scl = fits.Column(name='DAT_SCL', format='{0}E'.format(nchan_new*npol), array=dat_scl[:,:,mask])

#data = np.reshape(data, (nsub, nsblk/(8/nbits), npol, nchan_new)).astype(int)
data = data[:,:,mask_dat]
print(data.shape)
data = np.reshape(data, (nsub, int(nsblk/(8/nbits)), npol, nchan_new)).astype(int)
print (data.shape)
col_data = fits.Column(name='DATA', format='{0}B'.format(nchan_new*int(nsblk/(8/nbits))*npol), dim='({0},{1},{2})'.format(nchan_new, npol, int(nsblk/(8/nbits))), array=data)

#obsfreq_new = np.mean(dat_freq[0, mask])
obsfreq_new = np.mean(dat_freq[0, mask])
obsbw_new = (sub_end - sub_start + 1)*128   # MHz
obsnchan_new = nchan_new

################## Output to a new file ###################
out_file = args.output_file[0]

# create new primary HDU and update parameters
primary_hdu = fits.PrimaryHDU(header=hdu_primary_header)
# update OBSFREQ, OBSBW and OBSNCHAN
primary_hdu.header.set('OBSFREQ', obsfreq_new)
primary_hdu.header.set('OBSBW', obsbw_new)
primary_hdu.header.set('OBSNCHAN', obsnchan_new)

# fix other problems
primary_hdu.header.set('OBS_MODE', 'SEARCH')
#primary_hdu.header.set('CHAN_DM', dm)        # required by PRESTO
primary_hdu.header.set('TRK_MODE', 'TRACK')  # required by PRESTO

stt_crd1 = primary_hdu.header['STT_CRD1']
stt_crd2 = primary_hdu.header['STT_CRD2']

# Online cycle time
primary_hdu.header.set('NRCVR', 1)

# Scan length
#hdulist[0].header.set('SCANLEN', 60.)

# Online cycle time
#hdulist[0].header.set('TCYCLE', 10.)

primary_hdu.header.set('STP_CRD1', stt_crd1)
primary_hdu.header.set('STP_CRD2', stt_crd2)
primary_hdu.header.set('CAL_MODE', 'SYNC')
primary_hdu.header.set('CAL_FREQ', 11.123)
primary_hdu.header.set('CAL_DCYC', 0.)
primary_hdu.header.set('CAL_PHS', 0.25)

# Start LST, second
#hdulist[0].header.set('STT_LST', 10.)

# create new subint HDU and update parameters
cols = fits.ColDefs([col_index, col_tsubint, col_offssub, col_auxdm, col_auxrm, col_freq, col_wts, col_scl, col_off, col_data])
hdu_subint = fits.BinTableHDU.from_columns(cols, header=hdu_subint_header)
#print hdu_subint.columns

# update NCHAN, NSTOT, and REFFREQ
hdu_subint.header.set('NCHAN', obsnchan_new)
hdu_subint.header.set('REFFREQ', obsfreq_new)
hdu_subint.header.set('NSTOT', nstot_new)

# fix other problems
# Nr of bins/pulse period (for gated data)
hdu_subint.header.set('NBIN_PRD', 1)
hdu_subint.header.set('NCHNOFFS', 0)   # required by PRESTO

# Phase offset of bin 0 for gated data
hdu_subint.header.set('PHS_OFFS', 0.)
if npol == 4:
    hdu_subint.header.set('POL_TYPE', 'AABBCRCI')

# create new HDU list and write to a new fits
hdul = fits.HDUList([primary_hdu, hdu_subint])
print (hdul.info())
hdul.writeto(out_file, overwrite=True, checksum=True)

hdulist.close()
