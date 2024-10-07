import glob
import numpy as np
import subprocess
import os
from multiprocessing import Pool

def run_job(cmd):
    print(cmd)
    subprocess.call(cmd, shell=True)

##########################

all_sf = glob.glob('*.sf')

cmd_list = []
for filename in all_sf:
    outfile = filename + '.6to9'
    if os.path.isfile(outfile):
        print('%s done!'%outfile)
    else:
        cmd = 'python /DATA/CARINA_6/dai02a/UWL/search/split_search/splitSrch_v2.py -f {0} -o {1} -sub0 6 -sub1 9'.format(filename, outfile)
        cmd_list.append(cmd)

##########################
ncpus = 5
pool = Pool(ncpus)
pool.map(run_job, cmd_list)
pool.close()
pool.join()
