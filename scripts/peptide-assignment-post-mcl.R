##########
## Written by Scott McCain, July 3, 2019
##########

# peptide assignment to taxa and processing after MCL assignment and normalization

# loading libraries
library(tidyverse)
library(readxl)
library(cleaver)
library(foreach)
library(doParallel)
library(testthat)

source("post_processing_functions.R")

# reading in command line arguments
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("Two arguments must be supplied, the input directory for csv files, and the id string to write normalization output to", call.=FALSE)
} else if (length(args)==2) {
  input_directory <- args[1]
  id_string <- args[2]
}

# checking arguments
does_dir_exist <- dir.exists(input_directory)
expect_true(does_dir_exist)
expect_type(id_string, "character")

# this is what will be worked on:
print(input_directory)
print(id_string)

# input_directory <- "../data/tfg-all-database/"
# id_string <- '_tfg_all'
# # 

# read in file with contig mappings to clusters
mcl_df <- read_excel("../data/mcmurdo-metatrans/Bertrand_McCrow_TFG/annotation_allTFG.grpnorm_mmetsp_fc_pn_reclassified.edgeR.xlsx", 
                     sheet = 1)
# aggregating centric and pennate diatom category
mcl_df$group <- ifelse(mcl_df$group == 'Centric Diatom' | mcl_df$group == 'Pennate Diatom', 
                       yes = 'Diatom', 
                       no = mcl_df$group)


mcl_annot <- read_excel("../data/mcmurdo-metatrans/Bertrand_McCrow_TFG/summaries_combined/annotation_allTFG.mmetsp_fc_pn_reclassified.summary_cluster.edgeR.xlsx")

lowcount <- read.delim("../data/mcmurdo-metatrans/Bertrand_McCrow_TFG/annotation_allTFG.lowcount.tab",
                       sep = "\t")
nobest <- read.delim("../data/mcmurdo-metatrans/Bertrand_McCrow_TFG/annotation_allTFG.nobesthit.tab",
                     sep = "\t")

annot_with_go_terms <- read.delim("../data/new_annotation/antarctica_metaG_metaT_final_annotated_mcl_clusters/annotation_allTFG.table_with_GO.txt")

# which(annot_with_go_terms$GO_terms != "")

# mcl_df %>% filter(orf_id == 'contig_99308_962_2059_+')


# designate which files are replicates to then normalize the peptide intensities to
input_file_reps <- list(s01 = c("180412_0749_097_S01_Rep_1.csv", 
                                "180412_0749_097_S01_Rep_2.csv"),
                        s04 = c("180412_0749_097_S04_Rep_1.csv", 
                                "180412_0749_097_S04_Rep_2.csv"),
                        s07 = c("180412_0749_097_S07_Rep_1.csv", 
                                "180412_0749_097_S07_Rep_2.csv"),
                        s10 = c("180412_0749_097_S10_Rep_1.csv", 
                                "180412_0749_097_S10_Rep_2.csv"),
                        s02 = c("180412_0749_097_S02_Rep_1.csv", 
                                "180412_0749_097_S02_Rep_2.csv"),
                        s05 = c("180412_0749_097_S05_Rep_1.csv", 
                                "180412_0749_097_S05_Rep_2.csv"),
                        s08 = c("180412_0749_097_S08_Rep_1.csv", 
                                "180412_0749_097_S08_Rep_2.csv"),
                        s11 = c("180412_0749_097_S11_Rep_1.csv", 
                                "180412_0749_097_S11_Rep_2.csv"),
                        s03 = c("180412_0749_097_S03_Rep_1.csv", 
                                "180412_0749_097_S03_Rep_2.csv", 
                                "180412_0749_097_S03_Rep_3.csv"),
                        s06 = c("180412_0749_097_S06_Rep_1.csv", 
                                "180412_0749_097_S06_Rep_2.csv", 
                                "180412_0749_097_S06_Rep_3.csv"),
                        s09 = c("180412_0749_097_S09_Rep_1.csv", 
                                "180412_0749_097_S09_Rep_2.csv", 
                                "180412_0749_097_S09_Rep_3.csv"),
                        s12 = c("180412_0749_097_S12_Rep_1.csv", 
                                "180412_0749_097_S12_Rep_2.csv", 
                                "180412_0749_097_S12_Rep_3.csv"))

# designate which files are replicates to then normalize the peptide intensities to

# all_file_names_formatted_test <- append_dir_to_files(file_list = input_file_reps,
#                                                 input_directory = input_directory)

run_pipeline(file_list = input_file_reps, 
             input_directory = input_directory,
             id_string = id_string)

# s01 <- normalize_out(c("data/tfg-all-database/180412_0749_097_S01_Rep_1_GOS_927_0_1_tfg_all.csv", 
#                        "data/tfg-all-database/180412_0749_097_S01_Rep_2_GOS_927_0_1_tfg_all.csv"))
# t1 <- tax_assigner(s01)
# 
# s04 <- normalize_out(c("data/mzML-converted/180412_0749_097_S04_Rep_1.csv", 
#                        "data/mzML-converted/180412_0749_097_S04_Rep_2.csv"))
# t4 <- tax_assigner(s04)
# 
# s07 <- normalize_out(c("data/mzML-converted/180412_0749_097_S07_Rep_1.csv", 
#                        "data/mzML-converted/180412_0749_097_S07_Rep_2.csv"))
# t7 <- tax_assigner(s07)
# 
# s10 <- normalize_out(c("data/mzML-converted/180412_0749_097_S10_Rep_1.csv", 
#                        "data/mzML-converted/180412_0749_097_S10_Rep_2.csv"))
# t10 <- tax_assigner(s10)
# 
# annot_good01 <- annot_connect_good(normalize_output = s01)
# annot_fuzzy01 <- annot_connect_fuzzy(normalize_output = s01) 
# annot_tax01 <- annot_connect_tax(normalize_output = t1) 
# 
# annot_good04 <- annot_connect_good(normalize_output = s04)
# annot_fuzzy04 <- annot_connect_fuzzy(normalize_output = s04)
# 
# annot_good07 <- annot_connect_good(normalize_output = s07)
# annot_fuzzy07 <- annot_connect_fuzzy(normalize_output = s07)
# 
# annot_good10 <- annot_connect_good(normalize_output = s10)
# annot_fuzzy10 <- annot_connect_fuzzy(normalize_output = s10)
# 
# s02 <- normalize_out(c("data/mzML-converted/180412_0749_097_S02_Rep_1.csv", 
#                        "data/mzML-converted/180412_0749_097_S02_Rep_2.csv"))
# t2 <- tax_assigner(s02)
# 
# s05 <- normalize_out(c("data/mzML-converted/180412_0749_097_S05_Rep_1.csv", 
#                        "data/mzML-converted/180412_0749_097_S05_Rep_2.csv"))
# t5 <- tax_assigner(s05)
# 
# s08 <- normalize_out(c("data/mzML-converted/180412_0749_097_S08_Rep_1.csv", 
#                        "data/mzML-converted/180412_0749_097_S08_Rep_2.csv"))
# t8 <- tax_assigner(s08)
# 
# s11 <- normalize_out(c("data/mzML-converted/180412_0749_097_S11_Rep_1.csv", 
#                        "data/mzML-converted/180412_0749_097_S11_Rep_2.csv"))
# t11 <- tax_assigner(s11)
# 
# annot_good02 <- annot_connect_good(normalize_output = s02)
# annot_fuzzy02 <- annot_connect_fuzzy(normalize_output = s02)
# 
# annot_good05 <- annot_connect_good(normalize_output = s05)
# annot_fuzzy05 <- annot_connect_fuzzy(normalize_output = s05)
# 
# annot_good08 <- annot_connect_good(normalize_output = s08)
# annot_fuzzy08 <- annot_connect_fuzzy(normalize_output = s08)
# 
# annot_good11 <- annot_connect_good(normalize_output = s11)
# annot_fuzzy11 <- annot_connect_fuzzy(normalize_output = s11)
# 
# 
# s03 <- normalize_out(c("data/mzML-converted/180412_0749_097_S03_Rep_1.csv", 
#                        "data/mzML-converted/180412_0749_097_S03_Rep_2.csv", 
#                        "data/mzML-converted/180412_0749_097_S03_Rep_3.csv"))
# t3 <- tax_assigner(s03)
# s06 <- normalize_out(c("data/mzML-converted/180412_0749_097_S06_Rep_1.csv", 
#                        "data/mzML-converted/180412_0749_097_S06_Rep_2.csv", 
#                        "data/mzML-converted/180412_0749_097_S06_Rep_3.csv"))
# t6 <- tax_assigner(s06)
# s09 <- normalize_out(c("data/mzML-converted/180412_0749_097_S09_Rep_1.csv", 
#                        "data/mzML-converted/180412_0749_097_S09_Rep_2.csv", 
#                        "data/mzML-converted/180412_0749_097_S09_Rep_3.csv"))
# t9 <- tax_assigner(s09)
# s12 <- normalize_out(c("data/mzML-converted/180412_0749_097_S12_Rep_1.csv", 
#                        "data/mzML-converted/180412_0749_097_S12_Rep_2.csv", 
#                        "data/mzML-converted/180412_0749_097_S12_Rep_3.csv"))
# t12 <- tax_assigner(s12)
# 
# annot_good03 <- annot_connect_good(normalize_output = s03)
# annot_fuzzy03 <- annot_connect_fuzzy(normalize_output = s03)
# 
# annot_good06 <- annot_connect_good(normalize_output = s06)
# annot_fuzzy06 <- annot_connect_fuzzy(normalize_output = s06)
# 
# annot_good09 <- annot_connect_good(normalize_output = s09)
# annot_fuzzy09 <- annot_connect_fuzzy(normalize_output = s09)
# 
# annot_good12 <- annot_connect_good(normalize_output = s12)
# annot_fuzzy12 <- annot_connect_fuzzy(normalize_output = s12)
# 
# 
# filter01 <- process_results(sample1 = annot_good01, 
#                             sample2 = annot_good04, 
#                             sample3 = annot_good07, 
#                             sample4 = annot_good10)
# filter01_b <- process_results2(sample1 = annot_good01, 
#                                sample2 = annot_good04, 
#                                sample3 = annot_good07, 
#                                sample4 = annot_good10)
# filter01_a <- inner_join(x = filter01, y = mcl_annot, by = "cluster")
# 
# 
# filter08 <- process_results(sample1 = annot_good02, 
#                             sample2 = annot_good05, 
#                             sample3 = annot_good08, 
#                             sample4 = annot_good11)
# filter08_b <- process_results2(sample1 = annot_good02, 
#                                sample2 = annot_good05, 
#                                sample3 = annot_good08, 
#                                sample4 = annot_good11)
# filter08_a <- inner_join(x = filter08, y = mcl_annot, by = "cluster")
# 
# 
# filter30 <- process_results(sample1 = annot_good03, 
#                             sample2 = annot_good06, 
#                             sample3 = annot_good09, 
#                             sample4 = annot_good12)
# 
# filter30_b <- process_results2(sample1 = annot_good03, 
#                                sample2 = annot_good06, 
#                                sample3 = annot_good09, 
#                                sample4 = annot_good12)
# filter30_a <- inner_join(x = filter30, y = mcl_annot, by = "cluster")
# 
# write.csv(filter01_a, '../data/filter01_ross_sea_protein_cluster.csv')
# write.csv(filter08_a, '../data/filter08_ross_sea_protein_cluster.csv')
# write.csv(filter30_a, '../data/filter30_ross_sea_protein_cluster.csv')
