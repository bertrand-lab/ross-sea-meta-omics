# batch script for getting coarse grained assignments for metparoteomic data

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

tax_data_30 <- tax_data %>% dplyr::filter(filter == 3.0)
tax_data_08 <- tax_data %>% dplyr::filter(filter == 0.8)
tax_data_01 <- tax_data %>% dplyr::filter(filter == 0.1)

coarse_all30 <- get_coarse_grains(tax_group = "all", 
                                tax_assign_df = tax_data_30,
                                filter_type = 3.0, all_taxa = TRUE)
coarse_all08 <- get_coarse_grains(tax_group = "all", 
                                tax_assign_df = tax_data_08,
                                filter_type = 0.8, all_taxa = TRUE)
coarse_all01 <- get_coarse_grains(tax_group = "all", 
                                tax_assign_df = tax_data_01,
                                filter_type = 0.1, all_taxa = TRUE)

write.csv(coarse_all30, file = 'coarse_all_correction_factors_30greedy.csv')
write.csv(coarse_all08, file = 'coarse_all_correction_factors_08greedy.csv')
write.csv(coarse_all01, file = 'coarse_all_correction_factors_01greedy.csv')

