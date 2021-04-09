# RTPredict for every database that is sample specific
## Loop through every sample-specific database to predict retention times

for i in ../data/digested-databases/GOS*txt
do
echo $i

# remove the .txt ending to append to the output file
core_file=${i::-4}
echo $core_file

# predict retention times!
RTPredict -in_text $i -out_text:file $core_file'_tryptic_peptide_rt_oligo.csv' -svm_model ../data/cobia-analysis/subsampled20-pooled-sample-specific-model_oligo.txt -total_gradient_time 7200 &

done

