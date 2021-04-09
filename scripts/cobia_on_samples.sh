#!/bin/bash
### script goal:
# predict cofragmentation scores for every peptide in the database
# examine the distribution of cofragmentation scores for diatoms, is it different from the null distribution?
###

## To Do:
### write the DDA characteristics

# subset the idXML files to train a SVM, one SVM for all the samples because we used the same LC column and MS settings

# IDMerger across all samples

IDMerger -in ../data/specific-database/*FDR.idXML -out ../data/cobia-analysis/pooled-sample-specific.idXML

cobia subsample_idxml -f ../data/cobia-analysis/pooled-sample-specific.idXML -n ../data/cobia-analysis/subsampled100-pooled-sample-specific.txt -s 100

cobia subsample_idxml -f ../data/cobia-analysis/pooled-sample-specific.idXML -n ../data/cobia-analysis/subsampled50-pooled-sample-specific.txt -s 50

cobia subsample_idxml -f ../data/cobia-analysis/pooled-sample-specific.idXML -n ../data/cobia-analysis/subsampled20-pooled-sample-specific.txt -s 20


# RTPredict for every database that is sample specific
#RTPredict -in_text ../data/kleiner_data/Mock_Comm_RefDB_V3_trypsin.txt -out_text:file $temp_string'_tryptic_peptide_rt_oligo.csv' -svm_model $temp_string'model_oligo.txt' -total_gradient_time $retention_time &

# use cobia to format the openms modelled retention times
#cobia openms_modelled_rt -f ../data/kleiner_data/Run4_C4_2000ng_FDR_tryptic_peptide_rt_linear.csv -n ../data/kleiner_data/Run4_C4_2000ng_linear

# predict cofragmentation scores for every dataset
#cobia cofrag_prediction -l ../data/kleiner_data/dda_params_kleiner_460.csv -f ../data/kleiner_data/Run1_C4_832ng_rbf_lc-retention-times.csv -n ../data/kleiner_data/Run1_C4_832ng_rbf_cofrag --global global

