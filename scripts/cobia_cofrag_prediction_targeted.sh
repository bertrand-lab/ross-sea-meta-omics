# FOR EVERY PREDICTION FILE, FORMAT
# use cobia to format the openms modelled retention times
#cobia openms_modelled_rt -f ../data/kleiner_data/Run4_C4_2000ng_FDR_tryptic_peptide_rt_linear.csv -n ../data/kleiner_data/Run4_C4_2000ng_linear

for i in ../data/digested-databases/*oligo_lc-retention-times.csv
do
core_file=${i::-4}

echo $i
cobia cofrag_prediction -f $i -l ../data/cobia-analysis/dda_params_ross_sea_metap.csv --global targeted -n $core_file -t ../data/targeted-proteomics-data/target_peptides_cofrag.csv &

done


