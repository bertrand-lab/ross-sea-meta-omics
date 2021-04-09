# bash script for removing duplicate sequences from databases

cd ../data/formatted-databases/
for FILE in GOS_?_?.fasta; do
    echo $FILE
    echo "there are this many sequences in this pooled folder"
    grep -c ">" $FILE
    python ~/projects/fasta_utilities/remove_duplicates.py $FILE
done

