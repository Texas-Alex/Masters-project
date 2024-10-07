import math
import numpy as np
import subprocess
import argparse
import glob
import os
import pickle
import re
from collections import OrderedDict
from astropy.io import fits
from multiprocessing import Pool
from itertools import repeat

############### PRESTO calls ####################

def do_fold (cmd):
    print (cmd)
    subprocess.call(cmd, shell=True)

######## Observation information class ##############

class obs_info_class:
    '''
    class obs_info(args)
        a class containing all the information about the observation and data processing.
    '''
    def __init__(self, args):
        self.infile = args.input_filename        # filename

        self.snr = args.snr_cut           # S/N threashold
        self.ncpus = args.num_cpus        # Number of CPU to use
        self.tpart = args.time_part        # Sub-integration time for period search (npart = length/tpart) in seconds
        self.chnbw = args.chan_bw          # Sub-band bandwidth for DM search (nsub = bw/chnbw) in MHz
        #self.skipFold = args.skipFold     # Skip prepfold? Default: False
        self.length = float(args.length)                # second
        self.start = float(args.start)                # second
        
        print ('Input file name: %s'%self.infile)
        print ('Number of CPU to use: %s'%self.ncpus)
        print ('Signal-to-noise ratio: %s'%self.snr)

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
        length = self.tsamp*self.nout
        npart = int(length/self.tpart)        # The number of sub-integrations to use for the period search
        nsub = int(self.obsbw/self.chnbw)
        print ('Number of subbands for prepfold: %s'%nsub)
        print ('Number of sub-int for prepfold: %s'%npart)
        
        self.end = self.start + float(self.length/(self.tsamp*self.nout))
        ########## plan the folding ################
        print ('\n')
        print ('Planning the folding...\n')
        presto_db = 'fold.db'

        if os.path.isfile(presto_db):
            print ('%s alreadly exits!\n'%presto_db)
            self.db = pickle.load(open(presto_db,'rb'))

            for step, value in self.db.items():
                fold_out = step.split()[13]
                num_cand = step.split()[-8]
                pfd = fold_out + '_ACCEL_Cand_' + num_cand + '.pfd'
                print (pfd)
                if os.path.isfile(pfd):
                    self.db.update({step:1})    # 0: unprocessed; 1: processed
                    pickle.dump(self.db, open(presto_db, 'wb'))
                    print ('%s done!\n'%step)
                else:
                    self.db.update({step:0})    # 0: unprocessed; 1: processed
                    pickle.dump(self.db, open(presto_db, 'wb'))
                    print ('%s not done!\n'%step)
        else:
            print ('Creating %s!\n'%presto_db)

            ######### sift ###########
            proc = subprocess.Popen(['python', '/presto/examplescripts/ACCEL_sift.py'], shell=False, stdout=subprocess.PIPE)
            f = open('candlist.txt', 'w')
            #f.write(proc.stdout.read())
            f.write(proc.stdout.read().decode('utf-8'))
            f.close()

            ######## create fold database #######
            self.db = OrderedDict()

            file_cand = []
            num_cand = []
            snr_cand = []
            DM_f = []
            p_f = []
            f = open('candlist.txt', 'r')
            lines = f.readlines()
            f.close()
        
            for line in lines:
                if re.search('ACCEL',line):
                    snr_temp = float(line.split()[2])
                    if snr_temp >= self.snr:
            	        temp = line.split()[0]
            	        file_cand.append(temp.split(":")[0] + '.cand')
            	        num_cand.append(temp.split(":")[1])
            
            	        DM_f.append(float(line.split()[1]))
            	        snr_cand.append(float(line.split()[2]))
            	        p_f.append(float(line.split()[7])*1e-3)   # second
    
            for i in range(len(file_cand)):
                fold_out = self.infile + '_DM' + str(DM_f[i]) + '_P' + str(p_f[i])
                pfd = fold_out + '_ACCEL_Cand_' + str(num_cand[i]) + '.pfd'

                p = p_f[i]
                if p < 0.002:
                    #Mp, Mdm, N = 2, 2, 24
                    #npart = 50
                    #otheropts = "-ndmfact 3"
                    Mp, Mdm, N = 8, 2, 32
                    otheropts = "-pstep 1 -pdstep 1 -dmstep 3"
                elif p < 0.05:
                    Mp, Mdm, N = 8, 2, 64
                    otheropts = "-pstep 1 -pdstep 1 -dmstep 3"
                elif p < 0.5:
                    Mp, Mdm, N = 1, 1, 128
                    otheropts = "-pstep 1 -pdstep 2 -dmstep 1"
                else:
                    Mp, Mdm, N = 1, 1, 128
                    otheropts = "-nopdsearch -pstep 1 -pdstep 2 -dmstep 1"
                
                #cmd = 'prepfold -psrfits -noscales -nooffsets -noweights -noxwin -ncpus 1 -mask rfi_root_rfifind.mask -dm {0} -o {1} -nsub {2} -npart {3} {4} -n {5} -npfact {6} -ndmfact {7} -accelcand {8} -accelfile {9} {10}'.format(DM_f[i], fold_out, nsub, npart, otheropts, N, Mp, Mdm, num_cand[i], file_cand[i], self.infile)
                cmd = 'prepfold -psrfits -noscales -nooffsets -noweights -noxwin -ncpus 1 -mask rfi_root_rfifind.mask -dm {0} -o {1} -nsub {2} -npart {3} {4} -n {5} -npfact {6} -ndmfact {7} -accelcand {8} -accelfile {9} -start {10} -end {11} {12}'.format(DM_f[i], fold_out, nsub, npart, otheropts, N, Mp, Mdm, num_cand[i], file_cand[i], self.start, self.end, self.infile)
                print (cmd)
                self.db.update({cmd:0})    # 0: unprocessed; 1: processed
    
            pickle.dump(self.db, open(presto_db,'wb'))
    
#################################################################
####### Main ########
parser = argparse.ArgumentParser(description='Create PRESTO processing database')
# required
parser.add_argument('-f',        '--input_filename',      metavar='filename', required=True, help='Input file name')
parser.add_argument('-l',        '--length',         metavar='1800',          required=True, help='Length of the block (second)')
parser.add_argument('-s',        '--start',          metavar='0.2',           required=True, help='Starting point')
# not required
parser.add_argument('-snr',      '--snr_cut',        metavar='8.0',            help='Signal-to-noise ratio', default = 8.0, type = float)
parser.add_argument('-ncpus',    '--num_cpus',       metavar='10',             help='Number of CPU to use', default = 10, type = int)
parser.add_argument('-tpart',    '--time_part',      metavar='60.',            help='Sub-integration time for period search (npart = length/tpart) in seconds', default = 60., type = float)
parser.add_argument('-chnbw',    '--chan_bw',        metavar='8.',            help='Sub-band bandwidth for DM search (nsub = bw/chnbw) in MHz', default = 8., type = float)
#parser.add_argument('-skipFold',  action='store_true', default=False,          help='Skip prepfold?')
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

# prepfold
pool = Pool(obsInfo.ncpus)
pool.map(do_fold, cmd)
pool.close()
pool.join()
