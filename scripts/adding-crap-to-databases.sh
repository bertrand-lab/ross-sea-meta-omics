# bash script to take each of the databases and append the crap repository to it

for FILE in ls ../data/formatted-databases/*_nonredun.fasta; do
    echo $FILE
    temp_file=${FILE/.fasta/}
    echo $temp_file
    cat /var/www/sfolder/general/crap.fasta $FILE > $temp_file'_crap.fasta' 
done
