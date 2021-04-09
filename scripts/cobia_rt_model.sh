## training retention time models from observed peptides (using subsampled peptide observations from the pooled idXML file)
# different subsampled peptide numbers refer to different levels of subsampling

retention_time=7200

# fit an RT Model using 'OLIGO'
RTModel -in '../data/cobia-analysis/subsampled100-pooled-sample-specific.txt' -out '../data/cobia-analysis/subsampled100-pooled-sample-specific-model_oligo.txt' -kernel_type 'OLIGO' -total_gradient_time $retention_time &
RTModel -in '../data/cobia-analysis/subsampled50-pooled-sample-specific.txt' -out '../data/cobia-analysis/subsampled50-pooled-sample-specific-model_oligo.txt' -kernel_type 'OLIGO' -total_gradient_time $retention_time &
RTModel -in '../data/cobia-analysis/subsampled20-pooled-sample-specific.txt' -out '../data/cobia-analysis/subsampled20-pooled-sample-specific-model_oligo.txt' -kernel_type 'OLIGO' -total_gradient_time $retention_time &

