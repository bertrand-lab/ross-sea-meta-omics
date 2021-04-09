# bash script for removing duplicate sequences from databases

cd ../data/formatted-databases/
for FILE in *fasta; do
    echo $FILE
    python ~/projects/fasta_utilities/remove_duplicates.py $FILE
done

