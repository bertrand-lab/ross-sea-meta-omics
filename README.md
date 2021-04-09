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

`adding-crap-to-tfg.sh` is used to format the databases used to include the CRAP database. THe same was done for the other database configurations, using the following scripts: 'adding-crap-to-databases.sh', 'adding-crap-to-pooled.sh'.

`remove-duplicate-sequences-tfg.sh` is used to remove duplicate sequences using fasta-utlities from P. Wilmarth, the same was done with other databases using 'remove-duplicate-sequences.sh', and 'remove-duplicate-sequences-pooled-fraction.sh'.

`database-searching-openms-tfg.sh` runs the MSGF+ database searching. For other database configurations (described in more depth in paper), the scripts were: 'database-searching-openms-one-sample.sh', 'database-searching-openms-pooled-fraction.sh', 'database-searching-openms-tfg-t0.sh', and 'database-searching-openms.sh' (for the sample-specific databases).

`feature-finder-general.sh` this is a bash script for running FeatureFinderIdentification (Weisser et al 2017), which is called on by `feature-finder-all-ross-sea.sh`. (This second script runs the peptide quant for all database configurations).

### Calculating Total Ion Current

Calculating the total ion current using the python script 'calc_tic.py'.

### Assigning Peptides to Taxa

First we map peptides that are unique to a given taxa using `peptide-assignment-tfg.sh`, which calls on an R script `peptide-assignment-post-mcl.R`. This R script uses functions written in `post_processing_functions.R`.

After peptides are mapped to diatoms, we then need to get the coarse-grained diatom proteome, which can be done using 'mapping_diatom_coarse_grains.sh' which calls 'mapping_diatom_coarse_grains.R', which also uses functions from 'post_processing_functions.R'.

To map coarse grains to all taxa, we run 'mapping_meta_coarse_grains.sh' which calls on 'mapping_meta_coarse_grains.sh', which also uses functions from 'post_processing_functions.R'.

### Metaproteomic simulation

THe main script here is 'metaproteomic_bias_same_pool.py', which runs the generative and sampling model. This model is mostly similar to that published in McCain and Bertrand 2019 (Journal of Proteome Research), except here we included dynamic exclusion.

### Cobia analysis

'digest_databases.sh'

    - digested databases for inputs to cobia

'cobia_on_samples.sh'

    - subsampling the observed peptides for then running the RT model

'cobia_rt_model.sh'

    - Training the RT model

'cobia_rt_predict.sh'

    - Using the RT model to predict RT times

'cobia_cofrag_prediction.sh'

    - Digesting database and assigning RT times etc.

'cobia_cofrag_prediction_targeted.sh'

    - Cofragmentation risk score prediction with a list of peptides

'cobia_cofrag_prediction_targeted_discovery.sh'

    -  Cofragmentation risk score prediction with a another list of peptides

### Plotting and Analyses for Metaproteomic Paper

As mentioned above, the data from this repo are used to make inferences on the parameters in the model discussed in the repo mn-fe-allocation. However, for the paper where these data alone are discussed, the analyses and figures are done with the following:

'MCL-cluster-normalization.R'

    - Script required for summarizing the normalization factors

'normalization-comparison.R'

    - Compares normalization factors and makes supp figures

'plotting_metapep_sim.R'

    - Plots the metaproteomic simulation results

'discoverability_analysis_ribo_photo.R'

    - Plots the number of peptides for each coarse grained protein group and the simulation results. This requires the previous script to run.

'cobia_analysis_summarizing.R' 

    - Summarizing the cobia analysis

'formatting_nunn_proteomic_data.R'

    - Formatting the Nunn et al 2013 proteomic data for the Figure 2 comparison

'schmidt_cv_dist.R'

    - Getting the coefficient of variation distribution from Schmidt et al 2016

'getting_coarse_grained_all_taxa_mcl.R'

    - This is the main script for the figures in the manuscript. It has all main figures, and most of the supplementary figures as well, and a lot of the data munging.

'getting_database_characteristics.sh'

    - This is for making Table 1 in the manuscript.

'database-construction-notebook.ipynb'

    - To make the different database configurations this jupyter notebook was used. Specifically it subsets sequences from assemblies (as reported in the methods).
