# converting fasta files to single line fasta

for i in ../data/digested-databases/*fasta;
do

echo $i
no_ending_i=${i/.fasta/}'_single.fasta'
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $i > $no_ending_i

done
