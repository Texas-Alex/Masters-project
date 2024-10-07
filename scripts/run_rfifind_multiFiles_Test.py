#!/usr/bin/env python3 
import math
import numpy as np
import subprocess
import argparse
import glob
import os
import re
from collections import OrderedDict
from astropy.io import fits
#import pyfits
#from multiprocessing import Pool
#from itertools import repeat

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

############### PRESTO calls ####################
######## Observation information class ##############

class obs_info_class:
    '''
    class obs_info(args)
        a class containing all the information about the observation and data processing.
    '''
    def __init__(self, args):
        self.infile = args.input_filename
        self.outfile = args.output_filename
        self.rfi_time = args.rfi_time      # Integration time for rfifind
        #self.ncpus = args.num_cpus         # Number of threads to use
        self.smallFile = args.read_chnfreq # Read channel frequencies from a small file?
        self.freqsig = args.freq_sigma     # Frequency domain sigma
        self.timesig = args.time_sigma     # Time domain sigma
        
        print ('Input file name: %s' % " ".join(self.infile))
        #print ('Number of threads to use: %s'%self.ncpus)
        print ('Integration time for rfifind: %s'%self.rfi_time)

        self.zaplist = args.zap_list
        if self.zaplist == '':
            print ('No zapping channels!')
        else:
            print ('Channels to be zapped (channel number): %s'%self.zaplist)

        self.zapBirds = args.zapbirds_list
        if self.zapBirds == '':
            print ('No birdies!')
        else:
            print ('PRESTO birdies list: %s'%self.zapBirds)

        self.zapFreq = args.zapFreq_list
        if self.zapFreq == '':
            print ('No zapping frequency file!')
        else:
            print ('Frequency zapping file: %s'%self.zapFreq)

        ########### read in header information ###############
        print ('\n')
        if self.smallFile == '':
            print ('Reading header information from %s\n'%self.infile[0])
            hdu = fits.open(self.infile[0])
        else:
            print ('Reading header information from %s\n'%self.smallFile)
            hdu = fits.open(self.smallFile)

        self.obsbw = hdu[0].header['OBSBW']             # Total bandwidth (MHz)
        self.obsfreq = hdu[0].header['OBSFREQ']         # Central frequency (MHz)
        self.nchn = hdu['SUBINT'].header['NCHAN']       # Total number of channels
        self.npol = hdu['SUBINT'].header['NPOL']        # Number of polarisation

        #self.nsub = hdu['SUBINT'].header['NAXIS2']      # Number of subint
        self.nsblk = hdu['SUBINT'].header['NSBLK']      # Samples/row
        #self.nout = self.nsub*self.nsblk                # Spectra per file
        self.tsamp = hdu['SUBINT'].header['TBIN']       # Sampling timie (s)

        #tbdata = hdu['SUBINT'].data
        #self.chn_freq = tbdata['DAT_FREQ'][0]
        chn_bw = float(self.obsbw)/float(self.nchn)
        f0 = float(self.obsfreq) - float(self.obsbw)/2. + chn_bw/2.
        self.chn_freq = np.arange(f0, f0+self.obsbw, chn_bw)
        hdu.close()

        print ('Total Bandwidth: {0} MHz'.format(self.obsbw))
        print ('Central frequency: {0} MHz'.format(self.obsfreq))
        print ('Number of channels: {0}'.format(self.nchn))
        print ('Sample time: {0} s'.format(self.nchn))
        #print ('Spectra per file: {0}'.format(self.nout))

        self.cmd = [
                'rfifind -psrfits -noscales -nooffsets -noweights -time {1} -o {2} -freqsig {3} -timesig {4} {0}'.format(
            " ".join(self.infile), 
            self.rfi_time, 
            self.outfile, 
            self.freqsig, 
            self.timesig
            )
                ]

        ########### zapping channels  #############
        if self.zaplist != '' or self.zapFreq != '':
            self.cmd.append('-zapchan')

            chn_list = []
            ########### command line zapping channels  #############
            if self.zaplist != '':
                chn_list.append('%s'%self.zaplist)

            ########### zap frequencies #############
            if self.zapFreq != '':
                zap_freq = np.loadtxt(self.zapFreq)
                for temp in zap_freq:
                    mask = (self.chn_freq > math.floor(temp[0])) & (self.chn_freq < math.ceil(temp[1]))
                    index, = np.where(mask == True)
                    if len(index) != 0:
                        index = np.array(index)
                        chn_list.append('{0}:{1}'.format(np.amin(index), np.amax(index)))
                        print ('zap {0}:{1} ({2}:{3})'.format(np.amin(index), np.amax(index), math.floor(temp[0]), math.ceil(temp[1])))

            chn_list = ','.join(chn_list)
            self.cmd.append(chn_list)

        self.cmd = ' '.join(self.cmd)

#################################################################
####### Main ########
parser = argparse.ArgumentParser(description='Create PRESTO processing database')
# required
parser.add_argument('input_filename', metavar='filenames', nargs="+",  help='Input file name')
# not required
parser.add_argument('-zap',      '--zap_list',          metavar='0:20,30:40',     help='Channels to be zapped (channel number), e.g., "-zap 0:200,250:260"', default = '', type = str)
parser.add_argument('-t_rfi',    '--rfi_time',          metavar='2.0',            help='Integration time for rfifind', default = 2.0, type = float)
parser.add_argument('-zapFreq',  '--zapFreq_list',      metavar='zapFreqFile',    help='File containing all the frequencies to be zapped (MHz)', default = '', type = str)
parser.add_argument('-zapBirds', '--zapbirds_list',     metavar='birdfile',       help='PRESTO birdies list, needs to be prepared separately', default = '', type = str)
parser.add_argument('-readfreq', '--read_chnfreq',      metavar='filename',       help='A small file to read in channel frequencies', default = '', type = str)
parser.add_argument('-freqsig',  '--freq_sigma',        metavar='4',              help='The +/-sigma cutoff to reject time-domain chunks', default = 4, type = int)
parser.add_argument('-timesig',  '--time_sigma',        metavar='10',             help='The +/-sigma cutoff to reject freq-domain chunks', default = 10, type = int)
parser.add_argument('-out',      '--output_filename',   metavar='rfi_root',       help='Output file name', default = 'rfi_root', type = str)
#parser.add_argument('-ncpus',    '--num_cpus',       metavar='4',              help='Number of CPU to use', default = 10, type = int)
args = parser.parse_args()

#################################################################
####### read in observation information ########
obsInfo = obs_info_class(args)

#################################################################
####### run presto  ########

subprocess.call(obsInfo.cmd, shell=True)
