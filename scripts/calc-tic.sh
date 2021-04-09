#!/bin/bash
# Database searching

source activate openms_36

DIR='../data/mzML-converted/'

python calc_tic.py "pp_normalization_factors.csv" $(ls ../data/mzML-converted/180412_0749_097_S*PP.mzML) &
python calc_tic.py "bf_normalization_factors.csv" $(ls ../data/mzML-converted/180412_0749_097_S*BF.mzML)




