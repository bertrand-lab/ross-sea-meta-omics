#!/usr/bin/env bash

#######################################
# SETTING UP INPUT VARIABLES
######################################
# this is the root file folder where the idXML files from your database search are stored, e.g.: '../data/pooled-database/'
database_folder=$1

# this is just a modification of that, so that it lists all the FDR processed idXML files
database_with_fdr_idxml=$database_folder'*FDR.idXML'

# since we are testing different db configurations, we need to get the root file name without the database search appended to it, e.g.: '_pooled'
sub_string_id=$2

# where are the peak picked mzML files stored? this folder below!
# e.g.: '../mzML-converted/'
PP_folder_location=$3
echo $PP_folder_location

# making so it can loop through multiple different types of groupings
# variables number of grep strings, e.g.: 'S0[1,4,7]_Rep|S10_Rep' 'S0[2,5,8]_Rep|S11_Rep' 'S0[3,6,9]_Rep|S12_Rep'

######################################
# LOOPING THROUGH MASS SPEC IDXML FILES: MapAlign, FeatureFinderIdentification
######################################
# loop starts with the fourth argument passed to the bash script, but is flexible to the number of arguments
for i in "${@:4}"
do
subset_string=$i
echo $subset_string

#####################################
# LOOPING THROUGH THE SIMILAR MASS SPEC FILES, DESIGNATED WITH STRING MATCHES ABOVE IN STRINGARRAY
#####################################

# Variable of all files within a given grouping
all_files=$(ls $database_with_fdr_idxml | grep -E $subset_string)

for FILE in $all_files
do

        echo "Processing $FILE file..."

# Variable of all external files (see Weisser 2017 for description of external files)
        all_other_files=$all_files | grep -v $FILE
# Base file name
        file_base=${FILE/_FDR.idXML/}
# Retention times from similar MS runs have to be aligned, but we don't want to permanently make a copy of every aligned file.
# So we make a temporary copy of each idXML file.
       for SUB_FILE in $all_files
       do
               echo $SUB_FILE
               echo ${SUB_FILE/_FDR.idXML/}'-temp.idXML'
               cp $SUB_FILE ${SUB_FILE/_FDR.idXML/}'-temp.idXML'
       done
# Variable with all temporary copies of files
       all_temp_files=$(ls $database_folder | grep -e '-temp.idXML')

# Variable of all temporary copies of all external files
       all_other_temp_files=$(ls $database_folder | grep -e '-temp.idXML' | grep -v ${file_base#"$database_folder"})
# Moving directories
       cd $database_folder
# getting the base file name without the folder prefix
       file_base_no_folder=${file_base#"$database_folder"}

# getting the base file name without the appended database sub string id on the end
       file_base_no_folder_mzml=${file_base_no_folder/$sub_string_id/}
# Align other idXML files with the reference file being the main quantification file being looked at in this loop iteration.
       MapAlignerIdentification -in $all_temp_files -reference:file $file_base_no_folder'-temp.idXML' -out $all_temp_files
# Merge external idXML files
       IDMerger -in $all_other_temp_files -out $file_base_no_folder'-external.idXML'
# Identify features
       FeatureFinderIdentification -in $PP_folder_location$file_base_no_folder_mzml'_PP.mzML' -id $file_base_no_folder'_FDR.idXML' -id_ext $file_base_no_folder'-external.idXML' -out $file_base_no_folder'.featureXML'
# Remove temporary .idXML files.
       rm *-temp.idXML
# Go back to the previous working directory.
       cd -

done

done
