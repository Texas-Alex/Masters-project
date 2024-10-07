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
#from itertools import repeat
#from functools import partial

############### PRESTO calls ####################

class obs_info_class:
    '''
    class obs_info(args)
        a class containing all the information about the observation and data processing.
    '''
    def __init__(self, args):
        self.ncpus = args.num_cpus        # Number of CPU to use
        print ('Number of CPU to use: %s'%self.ncpus)

        ########## plan the searching ################
        print ('\n')
        print ('Planning the FFT...\n')
        presto_db = 'fft.db'

        if os.path.isfile(presto_db):
            print ('%s alreadly exits!\n'%presto_db)
            self.db = pickle.load(open(presto_db,'rb'))

            for datfile, value in self.db.items():
                redfft_temp = datfile[:-4] + '_red.fft'
                fft_temp = datfile[:-3] + 'fft'
                if value == 1:
                    print ('%s done!'%datfile)
                elif value == 0:
                    if os.path.isfile(datfile) or os.path.isfile(redfft_temp):
                        print ('%s not done!'%datfile)
                    else:
                        print ('%s done!'%datfile)
                        self.db.update({datfile:1})    # 0: unprocessed; 1: processed
                        pickle.dump(self.db, open(presto_db, 'wb'))
                else:
                    print ('Wrong value!\n')
                    exit(1)
        else:
            print ('Creating %s!\n'%presto_db)
            dat_list = sorted(glob.glob('*.dat'))

            if len(dat_list) == 0:
                print ('Please run "run_prepsubband.py"! Or FFT is done!\n')
                exit(1)
            else:
                self.db = OrderedDict()

                num_dat = len(dat_list)
                for i in range(num_dat):
                    self.db.update({dat_list[i]:0})    # 0: unprocessed; 1: processed

                pickle.dump(self.db, open(presto_db,'wb'))
                print ('%s created, ready to search!'%presto_db)

#################################################################
####### Main ########
parser = argparse.ArgumentParser(description='Create PRESTO processing database')
# required
# not required
parser.add_argument('-ncpus',    '--num_cpus',       metavar='10',             help='Number of CPU to use', default = 10, type = int)
args = parser.parse_args()

#################################################################
presto_db = 'fft.db'

def do_fft (datfile):
    cmd = 'realfft -fwd %s'%datfile
    subprocess.call(cmd, shell=True)

    fftfile = datfile[:-3] + 'fft'
    cmd = 'rednoise %s'%fftfile
    subprocess.call(cmd, shell=True)

    print ('Removing %s'%fftfile)
    os.remove(fftfile)

    # rename *_red.fft to *.fft
    redfftfile = fftfile[:-4] + '_red.fft'
    print ("Rename %s to %s"%(redfftfile, fftfile))
    os.rename(redfftfile, fftfile)

    print ('Removing %s'%datfile)
    os.remove(datfile)
            
####### read in observation information ########
obsInfo = obs_info_class(args)

####### processing ########
file2process = []
for dat_name, value in obsInfo.db.items():
    if value == 1:
        print ('%s done!'%dat_name)
    else:
        file2process.append(dat_name)
        print ('%s not done!'%dat_name)

# fft
pool = Pool(obsInfo.ncpus)
pool.map(do_fft, file2process)
pool.close()
pool.join()
