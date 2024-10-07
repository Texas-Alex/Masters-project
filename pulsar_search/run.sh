#!/usr/bin/env bash

scriptsBase=/software/projects/mwavcs/abennett/scripts
dm_max=500.0
dataFile=${@}

python3 ${scriptsBase}/run_rfifind_multiFiles_Test.py -t_rfi 2.0 -zapFreq ${scriptsBase}/uwl_zap_list -freqsig 4 -timesig 10 -out rfi_root ${dataFile}  &&
python3 ${scriptsBase}/run_prepsubband_multiCPU.py -ncpus 8 -dmlo 0.0 -dmhi ${dm_max} -nsubband 128 -cDM 0.0 ${dataFile} &&
python3 ${scriptsBase}/run_fft.py -ncpus 8 &&
rm -f *_red.inf &&
python3 ${scriptsBase}/run_search_accl.py -ncpus 8 -nharm 32 -z 200 &&
python3 ${scriptsBase}/run_fold_accl.py -ncpus 8 -snr 6.0 ${dataFile}
