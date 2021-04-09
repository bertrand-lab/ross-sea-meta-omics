#!/bin/bash
# Database searching and fdr application

# this loop is just for the sample specific databases
DIR='../data/mzML-converted/'
for FILE in "$DIR"*PP.mzML
do
	echo "Processing $FILE file..."
        temp_string=${FILE/.mzML/}
# db_string is the name of the database that should be used for a given mass spec file
        db_string=${temp_string:49:11}
# the name of the database has to be adjusted so it's the one with CRAP appended to it
        db_string_adjusted='../data/formatted-databases/'$db_string'_nonredun_crap.fasta'
# making the name that will be saved as idXML, so first we need to remove the prefix of the directory structure
        prefix="../data/mzML-converted/"
        temp_string2=${temp_string#"$prefix"}
# we also need to remove the '_PP' from the mass spec file name
        idxml_file_name=${temp_string2/_PP/}
# making the final output name, but now saving it to a given directory
        output_file_name_dir='../data/'$db_string'/'$idxml_file_name
# the PeptideIndexer command needs the reversed database, so a new variable is needed
        db_string_revcat='../data/formatted-databases/'$db_string'_nonredun_crap.revCat.fasta'
        echo $FILE; echo $temp_string; echo $db_string; echo $db_string_adjusted; echo $output_file_name_dir; echo $db_string_revcat

# running the database search
	MSGFPlusAdapter -in $FILE -executable /software/MSGFPlus/MSGFPlus_Oct2018.jar -database $db_string_adjusted -out $output_file_name_dir'.idXML' -add_decoys -threads 32
        PeptideIndexer -in $output_file_name_dir'.idXML' -fasta $db_string_revcat -out $output_file_name_dir'_PI.idXML' -threads 32 -decoy_string 'XXX_' -enzyme:specificity 'none'
        FalseDiscoveryRate -in $output_file_name_dir'_PI.idXML' -out $output_file_name_dir'_FDR.idXML' -PSM 'true' -FDR:PSM 0.01 -threads 32

done


