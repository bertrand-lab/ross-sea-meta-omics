`scripts/` includes all of the scripts for processing the mass spectrometry data. It requires a few installations before using (OpenMS is the biggest thing to install), but there are a few others. Here I detail how to run all the scripts. **Note that these data and the output for this pipeline are being used to write another paper simultaneously, so there are some components that would need to be tweaked to run this without error. But for the mn-fe-allocation (see other Bertrand Lab github repo with this name) data production, the pipeline will still work.

These scripts assume the following file structure:

```
/scripts
/data
    /mzML-converted/
    /tfg-all-database/
    /mcmurdo-metatrans/Bertrand_McCrow_TFG/
        /summaries_combined/
    /new_annotation/
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


