import math
import numpy as np
import subprocess
import argparse
import glob
import os
#import cPickle
import pickle
import re
from collections import OrderedDict
from astropy.io import fits
#import pyfits
from multiprocessing import Pool
from itertools import repeat

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
        self.infile = args.input_filename   # filename list

        self.dmlo = args.DM_low                  # lowest DM
        self.dmhi = args.DM_high                 # highest DM
        self.nsubDM = args.num_subband_dm        # number of subbands for prepsubband
        self.cDM = args.chn_dm                   # Coherent DM in each chan  (default = 0.0)
        self.length = float(args.length)                # second

        self.ncpus = args.num_cpus        # Number of CPU to use
        #self.skipRFI = args.skipRFI       # Skip rfifind? Default: False
        #self.skipFold = args.skipFold     # Skip prepfold? Default: False
        
        print ('Input file name list: %s'%self.infile)
        print ('Number of CPU to use: %s'%self.ncpus)
        print ('Lower bound of DM: %s'%self.dmlo)
        print ('Higher bound of DM: %s'%self.dmhi)
        print ('Number of subbands for prepsubband: %s'%self.nsubDM)

        ########### read in header information ###############
        print ('\n')

        print ('Reading header information from %s\n'%self.infile)
        hdu = fits.open(self.infile)
        self.obsbw = hdu[0].header['OBSBW']             # Total bandwidth (MHz)
        self.obsfreq = hdu[0].header['OBSFREQ']         # Central frequency (MHz)
        self.nchn = hdu['SUBINT'].header['NCHAN']       # Total number of channels
        self.npol = hdu['SUBINT'].header['NPOL']        # Number of polarisation

        self.nsub = hdu['SUBINT'].header['NAXIS2']      # Number of subint
        self.nsblk = hdu['SUBINT'].header['NSBLK']      # Samples/row
        self.nout = self.nsub*self.nsblk                # Spectra per file
        self.tsamp = hdu['SUBINT'].header['TBIN']       # Sampling timie (s)

        hdu.close()

        print ('Total Bandwidth: {0} MHz'.format(self.obsbw))
        print ('Central frequency: {0} MHz'.format(self.obsfreq))
        print ('Number of channels: {0}'.format(self.nchn))
        print ('Sample time: {0} s'.format(self.tsamp))
        print ('Spectra per file: {0}'.format(self.nout))

        #################################
        nseg = int(self.tsamp*self.nout/self.length)
        frac = self.length/(self.tsamp*self.nout)
        s0_list = []
        nout_new = []
        for i in range(nseg):
            s0_list.append(i*frac)
            if (i+1)*frac < 1:
                nout_new.append(int(self.length/self.tsamp))
            else:
                length = self.tsamp*self.nout - i*self.length
                nout_new.append(int(length/self.tsamp))
        print (nseg, nout_new, self.nout, s0_list)

        ########## plan the searching ################
        print ('\n')
        print ('Planning the searching...\n')
        presto_db = 'prepsubband.db'

        if os.path.isfile(presto_db):
            print ('%s alreadly exits!\n'%presto_db)
            self.db = pickle.load(open(presto_db,'rb'))
            #first_dat = min()

            for step, value in self.db.items():
                dm0 = float(step.split()[10])
                dat_test = step.split()[-1] + '_DM' + '%.2f'%dm0 + '.dat'
                inf_test = step.split()[-1] + '_DM' + '%.2f'%dm0 + '.inf'
                if os.path.isfile(dat_test) and os.path.isfile(inf_test):
                    self.db.update({step:1})    # 0: unprocessed; 1: processed
                    pickle.dump(self.db, open(presto_db,'wb'))
                    print ('%s done!'%step)
                else:
                    self.db.update({step:0})    # 0: unprocessed; 1: processed
                    pickle.dump(self.db, open(presto_db,'wb'))
                    print ('%s not done!'%step)
        else:
            print ('Creating %s!\n'%presto_db)
            self.db = OrderedDict()

            # Step 2: plan operations for prepsubbands, fft, rednoise, search...
            ddplan_fig = 'ddplan_fig.eps'
            cmd = 'DDplan.py -l {0} -d {1} -f {2} -b {3} -n {4} -t {5} -s {6} -o {7} -c {8}'.format(self.dmlo, self.dmhi, self.obsfreq, self.obsbw, self.nchn, self.tsamp, self.nsubDM, ddplan_fig, self.cDM)
            print (cmd)
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            lines = proc.stdout.readlines()

            dmlos = []
            ddm = []
            downsamp = []
            dm_call = []
            call = []
            for line in lines:
                numbers = re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", line.decode('utf-8'))
                #numbers = re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", line)
                if len(numbers) == 9:
                    dmlos.append(float(line.split()[0]))
                    ddm.append(float(line.split()[2]))      # DM step
                    downsamp.append(int(line.split()[3]))   # down sample
                    dm_call.append(int(line.split()[6]))    # DMs/call
                    call.append(int(line.split()[7]))       # Total number of call
                        	
            for i in range(len(dmlos)):
                lodm = dmlos[i]
                dmstep = ddm[i]
            	#print lodm, dmstep, nout
                for j in range(call[i]):
                    lodm_temp = lodm + j*dm_call[i]*dmstep
                    #dedispersion = '{0} {1} {2} {3} {4} {5} {6}'.format('Prepsubband', lodm_temp, dmstep, dm_call[i], self.nout, self.nsubDM, downsamp[i])
                    #cmd_dedispersion = 'prepsubband -psrfits -noscales -nooffsets -noweights -ncpus 1 -mask rfi_root_rfifind.mask -lodm {0} -dmstep {1} -numdms {2} -downsamp {3} -numout {4} -nsub {5} -o {6} {7}'.format(lodm_temp, dmstep, dm_call[i], downsamp[i], self.nout, self.nsubDM, self.infile, self.infile)
                    #print ('%s: 0'%cmd_dedispersion)
                    #self.db.update({cmd_dedispersion:0})    # 0: unprocessed; 1: processed
                    for k in range(nseg):
                        cmd_dedispersion = 'prepsubband -psrfits -noscales -nooffsets -noweights -ncpus 1 -mask rfi_root_rfifind.mask -lodm {0} -dmstep {1} -numdms {2} -downsamp {3} -numout {4} -nsub {5} -o {6} -start {7} {8}'.format(lodm_temp, dmstep, dm_call[i], downsamp[i], nout_new[k], self.nsubDM, self.infile + '_s{0:.2f}'.format(s0_list[k]), s0_list[k], self.infile)
                        self.db.update({cmd_dedispersion:0})    # 0: unprocessed; 1: processed
            
            pickle.dump(self.db, open(presto_db,'wb'))
            print ('%s created, ready to search!'%presto_db)
    
#################################################################
####### Main ########
parser = argparse.ArgumentParser(description='Create PRESTO processing database')
# required
parser.add_argument('-f',        '--input_filename',      metavar='filename',      required=True, help='Input file name')
parser.add_argument('-dmlo',     '--DM_low',              metavar='0.0',           required=True, help='Lower bound of DM')
parser.add_argument('-dmhi',     '--DM_high',             metavar='500.0',         required=True, help='Higher bound of DM')
parser.add_argument('-nsubband', '--num_subband_dm',      metavar='64',            required=True, help='Number of subbands for prepsubband')
parser.add_argument('-cDM',      '--chn_dm',              metavar='0.0',           required=True, help='Coherent DM in each chan  (default = 0.0)')
parser.add_argument('-l',        '--length',              metavar='1800',          required=True, help='Length of the block (second)')
# not required
parser.add_argument('-ncpus',    '--num_cpus',       metavar='10',             help='Number of CPU to use', default = 10, type = int)
#parser.add_argument('-skipRFI',  action='store_true', default=False,           help='Skip rfifind?')
#parser.add_argument('-skipFold',  action='store_true', default=False,           help='Skip prepfold?')
args = parser.parse_args()

#################################################################
####### read in observation information ########
obsInfo = obs_info_class(args)

#################################################################
####### run presto  ########
db = pickle.load(open('prepsubband.db','rb'))

def do_job (cmd):
    print (cmd)
    subprocess.call(cmd, shell=True)

cmd = []
for step, value in db.items():
    if value == 1:
        print ('%s done!'%step)
    else:
        cmd.append(step)
        print ('%s not done!'%step)

pool = Pool(obsInfo.ncpus)
pool.map(do_job, cmd)
pool.close()
pool.join()
