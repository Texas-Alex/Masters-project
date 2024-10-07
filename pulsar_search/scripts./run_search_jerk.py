import math
import numpy as np
import subprocess
import argparse
import glob
import os
import pickle
import re
from collections import OrderedDict
#from astropy.io import fits
from multiprocessing import Pool
from itertools import repeat

############### PRESTO calls ####################

def do_search (cmd):
    print (cmd)
    subprocess.call(cmd, shell=True)

    fftfile = cmd.split()[-1]
    zmax = cmd.split()[4]
    wmax = cmd.split()[6]
    candfile = fftfile[:-4] + '_ACCEL_' + zmax + '_JERK_' + wmax + '.cand'

    if os.path.isfile(candfile):
    	print ('Removing %s'%fftfile)
    	os.remove(fftfile)
    else:	
    	print ('%s accelsearch not finished'%fftfile)

######## Observation information class ##############

class obs_info_class:
    '''
    class obs_info(args)
        a class containing all the information about the observation and data processing.
    '''
    def __init__(self, args):
        #self.infile = args.input_filename        # filename
        self.zmax = args.zmax                    # zmax
        self.wmax = args.wmax                    # wmax

        #self.snr = args.snr_cut           # S/N threashold
        self.numharm = args.num_harm      # number of harmonics to sum
        self.flo = args.freq_low          # lowest spin frequency to search
        self.ncpus = args.num_cpus        # Number of CPU to use
        #self.skipDeDM = args.skipDeDM     # Skip prepfold? Default: False
        
        #print ('Input file name: %s'%self.infile)
        #print ('Signal-to-noise ratio: %s'%self.snr)
        print ('Number of CPU to use: %s'%self.ncpus)
        print ('zmax: %s'%self.zmax)
        print ('Number of harmonic: %s'%self.numharm)
        print ('Low frequency cut of accelsearch: %s'%self.flo)

        ########## plan the searching ################
        print ('\n')
        print ('Planning the searching...\n')
        presto_db = 'search.db'

        if os.path.isfile(presto_db):
            print ('%s alreadly exits!\n'%presto_db)
            self.db = pickle.load(open(presto_db,'rb'))

            for step, value in self.db.items():
                fftfile = step.split()[-1]
                candfile = fftfile[:-4] + '_ACCEL_' + str(self.zmax) + '_JERK_' + str(self.wmax) + '.cand'
                if value == 1:
                    print ('%s done!\n'%step)
                elif value == 0:
                    if os.path.isfile(fftfile):
                        print ('%s not done!\n'%step)
                    elif os.path.isfile(candfile):
                        print ('%s done!\n'%step)
                        self.db.update({step:1})    # 0: unprocessed; 1: processed
                        pickle.dump(self.db, open(presto_db, 'wb'))
                else:
                    print ('Wrong value!\n')
                    exit(1)
        else:
            print ('Creating %s!\n'%presto_db)
            fft_list = sorted(glob.glob('*.fft'))

            if len(fft_list) == 0:
                print ('Please run "run_fft.py"! Or Accelsearch is done!\n')
                exit(1)
            else:
                self.db = OrderedDict()

                num_dat = len(fft_list)
                for i in range(num_dat):
                    #cmd = 'accelsearch -numharm {0} -zmax {1} -flo {2} {3}'.format(self.numharm, self.zmax, self.flo, fft_list[i])
                    cmd = 'accelsearch -numharm {0} -zmax {1} -wmax {2} -flo {3} {4}'.format(self.numharm, self.zmax, self.wmax, self.flo, fft_list[i])
                    self.db.update({cmd:0})    # 0: unprocessed; 1: processed

                pickle.dump(self.db, open(presto_db,'wb'))
                print ('%s created, ready to search!'%presto_db)

#################################################################
####### Main ########
parser = argparse.ArgumentParser(description='Create PRESTO processing database')
# required
#parser.add_argument('-f',        '--input_filename',      metavar='filename', required=True, help='Input file name')
parser.add_argument('-z',        '--zmax',                metavar='0',        required=True, help='zmax')
#parser.add_argument('-nsubfold', '--num_subband_fold',    metavar='64',       required=True, help='Number of subbands for prepfold')
# not required
#parser.add_argument('-snr',      '--snr_cut',        metavar='8.0',            help='Signal-to-noise ratio', default = 8.0, type = float)
parser.add_argument('-nharm',    '--num_harm',       metavar='8',              help='Number of harmonic', default = 8, type = int)
parser.add_argument('-flo',      '--freq_low',       metavar='0.1',            help='Low frequency cut of accelsearch', default = 0.1, type = float)
parser.add_argument('-ncpus',    '--num_cpus',       metavar='10',             help='Number of CPU to use', default = 10, type = int)
# 2020/01/07
parser.add_argument('-w',        '--wmax',           metavar='0',              help='Jerk term', default = 0, type = int)
#parser.add_argument('-skipRFI',  action='store_true', default=False,           help='Skip rfifind?')
#parser.add_argument('-skipFold',  action='store_true', default=False,          help='Skip prepfold?')
#parser.add_argument('-skipDeDM',  action='store_true', default=False,          help='Skip prepsubband?')
args = parser.parse_args()

#################################################################
####### read in observation information ########
obsInfo = obs_info_class(args)

#################################################################
####### run presto  ########
cmd = []
for step, value in obsInfo.db.items():
    if value == 1:
        print ('%s done!'%step)
    else:
        cmd.append(step)
        print ('%s not done!'%step)

# fft
pool = Pool(obsInfo.ncpus)
pool.map(do_search, cmd)
pool.close()
pool.join()
