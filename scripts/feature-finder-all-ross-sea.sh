
# pooled database
bash feature-finder-general.sh '../data/pooled-database/' '_pooled' '../mzML-converted/' 'S0[1,4,7]_Rep|S10_Rep' 'S0[2,5,8]_Rep|S11_Rep' 'S0[3,6,9]_Rep|S12_Rep' &
# sample specific databases
bash feature-finder-general.sh '../data/specific-database/' '' '../mzML-converted/' 'S0[1,4,7]_Rep|S10_Rep' 'S0[2,5,8]_Rep|S11_Rep' 'S0[3,6,9]_Rep|S12_Rep' &
# just using sample 927
bash feature-finder-general.sh '../data/one-sample/' '_single' '../mzML-converted/' 'S0[1,4,7]_Rep|S10_Rep' 'S0[2,5,8]_Rep|S11_Rep' 'S0[3,6,9]_Rep|S12_Rep' &
# tfg just t=0 metatranscriptomes
bash feature-finder-general.sh '../data/tfg-t0-database/' '_tfg_t0' '../mzML-converted/' 'S0[1,4,7]_Rep|S10_Rep' 'S0[2,5,8]_Rep|S11_Rep' 'S0[3,6,9]_Rep|S12_Rep' &
# tfg all metaT
bash feature-finder-general.sh '../data/tfg-all-database/' '_tfg_all' '../mzML-converted/' 'S0[1,4,7]_Rep|S10_Rep' 'S0[2,5,8]_Rep|S11_Rep' 'S0[3,6,9]_Rep|S12_Rep' &

