### trypsin digest is required for cobia

# tryptic digest of database file
# for every database that is sample specific
for FILE in  ../data/formatted-databases/GOS_???_?_?_nonredun.fasta; do
    file_no_fasta=${FILE/.fasta}'_trypsin_digest.txt'
    cobia database_trypsin -f $FILE -n $file_no_fasta -c no-write & # this should all be run in parallel
done

mv ../data/formatted-databases/*digest*.txt ../data/cobia-analysis/
