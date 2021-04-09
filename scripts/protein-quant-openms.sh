#!/bin/bash
# Last step of OpenMS pipeline for quantifying peptides.

DIR=$1
for FILE in "$DIR"*.featureXML
do
	echo "Processing $FILE file..."
	temp_file=${FILE/.featureXML/}
	echo $temp_file
	ProteinQuantifier -in $FILE -peptide_out $temp_file'.csv'
done

