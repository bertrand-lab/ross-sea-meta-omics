#!/bin/bash
# Going through raw mass spec data (mzML files) and preprocessing them: noise filtering, baseline filtering, and peak picking.

DIR='../data/mzML-converted/'

for FILE in "$DIR"*Rep??.mzML
do
	echo "Processing $FILE file..."
        temp_string=${FILE/.mzML/}
        echo $temp_string
	NoiseFilterSGolay -in $FILE -out $temp_string'_SG.mzML'
	BaselineFilter -in $temp_string'_SG.mzML' -out $temp_string'_BF.mzML'
	PeakPickerHiRes -in $temp_string'_BF.mzML' -out $temp_string'_PP.mzML'
done

