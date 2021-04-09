# FOR EVERY PREDICTION FILE, FORMAT
# use cobia to format the openms modelled retention times
#cobia openms_modelled_rt -f ../data/kleiner_data/Run4_C4_2000ng_FDR_tryptic_peptide_rt_linear.csv -n ../data/kleiner_data/Run4_C4_2000ng_linear
for i in ../data/digested-databases/*oligo.csv
do
core_file=${i::-4}
echo $i
cobia openms_modelled_rt -f $i -n $core_file &
done


