#!/bin/bash
# Database searching and fdr application

# this loop is just for the pooled fraction specific databases

# loop through each mass spec file, which has the database name encoded that should be used
# for this, we use only part of the mass spec file name to be used, the size fraction part of it.

DIR='../data/mzML-converted/'
for FILE in "$DIR"*PP.mzML
do
	echo "Processing $FILE file..."
        temp_string=${FILE/.mzML/}
# db_string is the name of the database that should be used for a given mass spec file
        db_string=${temp_string:49:11}
# db string for the pooled by fraction
        db_string_fraction='GOS'${db_string:7:4}
        echo $db_string_fraction

# the name of the database has to be adjusted so it's the one with CRAP appended to it
        db_string_adjusted='../data/formatted-databases/'$db_string_fraction'_nonredun_crap.fasta'
# making the name that will be saved as idXML, so first we need to remove the prefix of the directory structure
        prefix="../data/mzML-converted/"
        temp_string2=${temp_string#"$prefix"}
# we also need to remove the '_PP' from the mass spec file name
        idxml_file_name=${temp_string2/_PP/}
# making the final output name, but now saving it to a given directory
        output_file_name_dir='../data/'$db_string_fraction'/'$idxml_file_name'_pooled'
# the PeptideIndexer command needs the reversed database, so a new variable is needed
        db_string_revcat='../data/formatted-databases/'$db_string_fraction'_nonredun_crap.revCat.fasta'
        echo "input file name $FILE"
        echo "input file name without the ending: $temp_string"
        echo "name of the database that should be used for a ms file $db_string" 
        echo " adjusting the name of the database string $db_string_adjusted" 
        echo " output file name with directory that will be an idxml file $output_file_name_dir" 
        echo " name of database with the reverse sequences included $db_string_revcat"
        echo "to write .. $output_file_name_dir.idXML"
# running the database search
	MSGFPlusAdapter -in $FILE -executable /software/MSGFPlus/MSGFPlus_Oct2018.jar -database $db_string_adjusted -out $output_file_name_dir'.idXML' -add_decoys -threads 32
        PeptideIndexer -in $output_file_name_dir'.idXML' -fasta $db_string_revcat -out $output_file_name_dir'_PI.idXML' -threads 32 -decoy_string 'XXX_' -enzyme:specificity 'none'
        FalseDiscoveryRate -in $output_file_name_dir'_PI.idXML' -out $output_file_name_dir'_FDR.idXML' -PSM 'true' -FDR:PSM 0.01 -threads 32

done


