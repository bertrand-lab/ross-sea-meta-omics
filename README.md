`scripts/` includes all of the scripts for processing the mass spectrometry data. It requires a few installations before using (OpenMS is the biggest thing to install), but there are a few others. Here I detail how to run all the scripts. 

**Note that these data and the output for this pipeline produce results for two distinct papers. One paper uses the metaproteomic data to connect with the proteomic allocation model described in mn-fe-allocation (other github repo). I have indicated where certain scripts are necessary for that repo below.

These scripts assume the following file structure:

```
/scripts/
/data/
    /mzML-converted/
    /tfg-all-database/
    /mcmurdo-metatrans/Bertrand_McCrow_TFG/
        /summaries_combined/
    /new_annotation/
    /metaGT/
        /AAHZ_L6_ANT_pool2_MetaT
        /AAHZ_L5_ANT_pool1_MetaG
    /GOS_927_0_1
    /GOS_927_0_8
    /GOS_927_3_0
    /GOS_930_0_1
    /GOS_930_0_8
    /GOS_930_3_0
    /GOS_933_0_1
    /GOS_933_0_8
    /GOS_933_3_0
    /GOS_935_0_1
    /GOS_935_0_8
    /GOS_935_3_0
    /GOS_0_1
    /GOS_0_8
    /GOS_3_0
    /pooled-database
    /specific-database
    /one-sample
    /tfg-t0-database
    /tfg-all-database
    /formatted-databases
    /go-slim-analysis
    /normalization-data
    /correction-factor-output
    /digested-databases
    /targeted-proteomics-data
    /cobia-analysis
    /pep-sample-sim
    /go-tax-pipeline-output
    /culture-proteomics-data
/figures/
```

### Database Searching and Peptide Quantification

`mass-spectra-preprocessing.sh` is first used to preprocess the mass spec data.

`adding-crap-to-tfg.sh` is used to format the databases used to include the CRAP database.

`remove-duplicate-sequences-tfg.sh` is used to remove duplicate sequences using fasta-utlities from P. Wilmarth.

`database-searching-openms-tfg.sh` runs the MSGF+ database searching.

`feature-finder-general.sh` this is a bash script for running FeatureFinderIdentification (Weisser et al 2017), which is called on by `feature-finder-all-ross-sea.sh`. Note that there are extra databases used for the second script that were not included in the current publication (mn-fe-allocation model paper) but will be included in further.

### Assigning Peptides to Taxa

First we map peptides that are unique to a given taxa using `peptide-assignment-tfg.sh`, which calls on an R script `peptide-assignment-post-mcl.R`. This R script uses functions written in `post_processing_functions.R`.

After peptides are mapped to diatoms, we then need to get the coarse-grained diatom proteome, which can be done using `mapping_diatom_coarse_grains.sh` which calls `mapping_diatom_coarse_grains.R`, which also uses functions from `post_processing_functions.R`.


