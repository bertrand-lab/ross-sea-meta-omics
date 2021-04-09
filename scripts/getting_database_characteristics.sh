# getting database characteristics for number of protein sequences

########## POOLED DATABASES
prot_seq_01=$(grep -c '>' ../data/formatted-databases/GOS_0_1_nonredun_crap.fasta)
prot_seq_08=$(grep -c '>' ../data/formatted-databases/GOS_0_8_nonredun_crap.fasta)
prot_seq_30=$(grep -c '>' ../data/formatted-databases/GOS_3_0_nonredun_crap.fasta)

echo 'pooled db 01: '
echo $prot_seq_01

echo 'pooled db 08: '
echo $prot_seq_08

echo 'pooled db 30: '
echo $prot_seq_30


sum_pooled=$(($prot_seq_01 + $prot_seq_08 + $prot_seq_30))
avg_pooled=$(echo $sum_pooled / 3 | bc -l )
echo 'Average number of protein sequences for pooled databases: '
echo $avg_pooled


########### SAMPLE SPECIFIC DATABASES
overall_sample_specific_sum=0
for i in ../data/formatted-databases/GOS_???_?_?_nonredun_crap.fasta;
do
prot_seq_sum_temp=$(grep -c '>' $i)
echo $i
echo $prot_seq_sum_temp
overall_sample_specific_sum=$((overall_sample_specific_sum + prot_seq_sum_temp))
#echo $overall_sample_specific_sum
done

avg_specific=$(echo $overall_sample_specific_sum / 12 | bc -l)
echo 'Average number of protein sequences for sample-specific databases: '
echo $avg_specific

########## ONE-SAMPLE DATABASES
overall_one_sample_specific_sum=0
for i in ../data/formatted-databases/GOS_930_?_?_nonredun_crap.fasta;
do
prot_seq_sum_temp=$(grep -c '>' $i)
echo $i
echo $prot_seq_sum_temp
overall_one_sample_specific_sum=$((overall_one_sample_specific_sum + prot_seq_sum_temp))
#echo $overall_one_sample_specific_sum
done

avg_specific=$(echo $overall_one_sample_specific_sum / 3 | bc -l)
echo 'Average number of protein sequences for one sample databases: '
echo $avg_specific


######### TFG T-0 DATABASE
echo 'Metatranscriptome experiment t0'
grep -c '>' ../data/formatted-databases/tfg_t0_nonredun_crap.fasta

######### TFG T-0 DATABASE
echo 'Metatranscriptome experiment all'
grep -c '>' ../data/formatted-databases/tfg_assembly.orf_nonredun_crap.fasta


