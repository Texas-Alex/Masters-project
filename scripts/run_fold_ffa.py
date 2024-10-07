import math
import numpy as np
import subprocess
import argparse
import glob
import os
import pickle
import re
from collections import OrderedDict
import csv
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
        #self.nsubfold = args.num_subband_fold    # number of subbands for prepfold
        self.cand_ffa = args.cand_ffa    # number of subbands for prepfold

        self.snr = args.snr_cut           # S/N threashold
        self.ncpus = args.num_cpus        # Number of CPU to use
        self.tpart = args.time_part        # Sub-integration time for period search (npart = length/tpart) in seconds
        self.chnbw = args.chan_bw          # Sub-band bandwidth for DM search (nsub = bw/chnbw) in MHz
        #self.skipFold = args.skipFold     # Skip prepfold? Default: False

        print ('Input file name: %s'%self.infile)
        print ('Number of CPU to use: %s'%self.ncpus)
        print ('Signal-to-noise ratio: %s'%self.snr)
        #print ('Number of subbands for prepfold: %s'%self.nsubfold)

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
        
        ########## plan the folding ################
        print ('\n')
        print ('Planning the folding...\n')
        presto_db = 'fold.db'

        if os.path.isfile(presto_db):
            print ('%s alreadly exits!\n'%presto_db)
            self.db = pickle.load(open(presto_db,'rb'))

            for step, value in self.db.items():
                fold_out = step.split()[13]
                pfd = fold_out + '.pfd'
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
            ###### read in FFA candidates ###########
            #ffa = np.loadtxt(self.cand_ffa)
            ffa = []
            with open(self.cand_ffa, newline='') as csvfile:
                reader = csv.reader(csvfile, delimiter=',')
                for row in reader:
                    if row[0] != 'period':
                        ffa.append(np.array(row)) 

            ffa = np.array(ffa)
            ######## create fold database #######
            print(ffa.shape)
            self.db = OrderedDict()

            for cand in ffa:
                #DM_f.append(dm)
            	#p_f.append(1./f0)   # second
                p = float(cand[0])
                dm = float(cand[2])
                snr = float(cand[-1])

                if snr > self.snr:
                    fold_out = self.infile + '_DM' + str(dm) + '_P' + str(p)
                    pfd = fold_out + '.pfd'
                    print (fold_out)

                    if p < 0.002:
                        #Mp, Mdm, N = 2, 2, 24
                        #npart = 50
                        #otheropts = "-ndmfact 3"
                        Mp, Mdm, N = 4, 2, 32
                        otheropts = "-pstep 1 -pdstep 1 -dmstep 3"
                    elif p < 0.05:
                        #Mp, Mdm, N = 2, 1, 50
                        #npart = 40
                        #otheropts = "-pstep 1 -pdstep 2 -dmstep 3"
                        Mp, Mdm, N = 4, 2, 64
                        otheropts = "-pstep 1 -pdstep 1 -dmstep 3"
                    elif p < 0.5:
                        #Mp, Mdm, N = 1, 1, 100
                        #npart = 30
                        #otheropts = "-pstep 1 -pdstep 2 -dmstep 1"
                        Mp, Mdm, N = 1, 1, 128
                        otheropts = "-pstep 1 -pdstep 2 -dmstep 1"
                    else:
                        #Mp, Mdm, N = 1, 1, 200
                        #npart = 30
                        #otheropts = "-nopdsearch -pstep 1 -pdstep 2 -dmstep 1"
                        Mp, Mdm, N = 1, 1, 128
                        otheropts = "-nopdsearch -pstep 1 -pdstep 2 -dmstep 1"
    
        	    #cmd = 'prepfold -psrfits -noscales -nooffsets -noweights -noxwin -ncpus 1 -mask rfi_root_rfifind.mask -dm {0} -o {1} -nsub {2} -npart {3} {4} -n {5} -npfact {6} -ndmfact {7} -accelcand {8} -accelfile {9} {10}'.format(DM_f[i], fold_out, self.nsubfold, npart, otheropts, N, Mp, Mdm, num_cand[i], file_cand[i], self.infile)
                    cmd = 'prepfold -psrfits -noscales -nooffsets -noweights -noxwin -ncpus 1 -mask rfi_root_rfifind.mask -dm {0} -o {1} -nsub {2} -npart {3} {4} -n {5} -npfact {6} -ndmfact {7} -p {8} {9}'.format(dm, fold_out, nsub, npart, otheropts, N, Mp, Mdm, p, self.infile)
                    print (cmd)
                    self.db.update({cmd:0})    # 0: unprocessed; 1: processed
    
            pickle.dump(self.db, open(presto_db,'wb'))
    
#################################################################
####### Main ########
parser = argparse.ArgumentParser(description='Create PRESTO processing database')
# required
parser.add_argument('-f',        '--input_filename',      metavar='filename',       required=True, help='Input file name')
parser.add_argument('-c',        '--cand_ffa',            metavar='FFA candidates', required=True, help='FFA candidates from Riptide')
#parser.add_argument('-nsubfold', '--num_subband_fold',    metavar='64',             required=True, help='Number of subbands for prepfold')
# not required
parser.add_argument('-snr',      '--snr_cut',        metavar='8.0',            help='Signal-to-noise ratio', default = 8.0, type = float)
parser.add_argument('-ncpus',    '--num_cpus',       metavar='10',             help='Number of CPU to use', default = 10, type = int)
parser.add_argument('-tpart',    '--time_part',      metavar='60.',            help='Sub-integration time for period search (npart = length/tpart) in seconds', default = 60., type = float)
parser.add_argument('-chnbw',    '--chan_bw',        metavar='10.',            help='Sub-band bandwidth for DM search (nsub = bw/chnbw) in MHz', default = 10., type = float)
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
