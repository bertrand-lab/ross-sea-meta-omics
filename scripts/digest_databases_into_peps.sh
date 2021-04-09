# looping through each database for metaproteomics project and digesting them into peptides to look at theoretical overlap
for i in ../data/digested-databases/*single*;
do

echo $i
sub_db="${i/.fasta/}"

## digesting each database
cobia database_trypsin -f $i -n $sub_db'.txt' -c "no-write"

done
