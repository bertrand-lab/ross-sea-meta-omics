## script for taking total model output and determining diatom proteome annotations

library(ggplot2)
library(magrittr)
library(dplyr)
library(tidyverse)
library(readxl)
library(cleaver)
library(foreach)
library(doParallel)
library(testthat)

# read in peptide quant data where tax_assigner has been run already

source("post_processing_functions.R")

tax_data <- read.csv("../data/go-tax-pipeline-output/pipeline_outtfg-all_tax_peptides.csv")

mcl_df <- read_excel("../data/mcmurdo-metatrans/Bertrand_McCrow_TFG/annotation_allTFG.grpnorm_mmetsp_fc_pn_reclassified.edgeR.xlsx", 
                     sheet = 1)

# coarse_diatoms <- get_coarse_grains(tax_group = 'Diatom', 
                                    # tax_assign_df = tax_data, 
                                    # filter_type = 3.0)

# write.csv(coarse_diatoms, file = 'coarse_diatoms_tfg_all.csv')

tax_data_30 <- tax_data %>% filter(filter == 3.0, tax_assign == 'Diatom')

coarse_diatoms_specific <- get_coarse_grains(tax_group = 'Diatom', 
                            tax_assign_df = tax_data_30,
                            filter_type = 3.0)

write.csv(coarse_diatoms_specific, file = '../data/correction-factor-output/coarse_diatoms_tfg_all_specific.csv')

