# Script goals:

# Upload mass spec quant data output from OpenMS.
# Upload contig-MCL cluster assignments
# Connect MCL-cluster assignments to quant 
# Normalization: 
#   - Collapse across repeat injections
#   - Total ion current per run
#   - Average # of tryptic peptides per MCL-cluster

library(tidyverse)
library(readxl)
library(cleaver)
library(foreach)
library(doParallel)
library(testthat)

## read in the functions for processing these data:

source("post_processing_functions.R")

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("Two arguments must be supplied, the input directory for csv files, and the id string to write normalization output to", call.=FALSE)
} else if (length(args)==2) {
  input_directory <- args[1]
  id_string <- args[2]
}

does_dir_exist <- dir.exists(input_directory)
expect_true(does_dir_exist)

expect_type(id_string, "character")

print(input_directory)
print(id_string)
# 
## Calculating MCL cluster normalization factors

# Assigning peptides to MCL clusters

# mzmldir <- dir("data/tfg-all-database/")
# csvfiles <- mzmldir[grep(".csv", x = mzmldir, fixed = TRUE)]
# # csvfilesfdr <- mzmldir[grep("_FDR.csv", x = mzmldir, fixed = TRUE)]
# csvdir <- paste0('data/tfg-all_database/', csvfiles)
# # csvfdr <- paste0('../data/mzML-converted/', csvfilesfdr)

# read in file with contig mappings to clusters
# mcl_df <- read_excel("../data/mcmurdo-metatrans/Bertrand_McCrow_TFG/annotation_allTFG.grpnorm_mmetsp_fc_pn_reclassified.edgeR.xlsx", 
                     # sheet = 1)

#mcl_annot <- read_excel("../data/mcmurdo-metatrans/Bertrand_McCrow_TFG/summaries_combined/annotation_allTFG.mmetsp_fc_pn_reclassified.summary_cluster.edgeR.xlsx")

#lowcount <- read.delim("../data/mcmurdo-metatrans/Bertrand_McCrow_TFG/annotation_allTFG.lowcount.tab",
                       # sep = "\t")
#nobest <- read.delim("../data/mcmurdo-metatrans/Bertrand_McCrow_TFG/annotation_allTFG.nobesthit.tab",
                     # sep = "\t")

# id_peps_sub_norm_factors <- id_peps_norm_factors(csvfdr, csvdir)

# assign_mcl_to_dir(target_dir = "../data/tfg-all-database/")
# calc_norm_factors(target_dir = "../data/tfg-all-database/", norm_file_name = "../data/normalization-data/norm-tfg-all-database.csv")

# assign_mcl_to_dir(target_dir = input_directory)
calc_norm_factors(target_dir = input_directory,
                  norm_file_name = paste0("../data/normalization-data/norm-", 
                                          id_string, ".csv"))

# assign_mcl_to_dir(target_dir = "../data/pooled-database/")
# calc_norm_factors(target_dir = "../data/pooled-database/", 
#                   norm_file_name = "../data/normalization-data/norm-pooled-database.csv")
# 
# assign_mcl_to_dir(target_dir = "../data/specific-database/")
# calc_norm_factors(target_dir = "../data/specific-database/", 
#                   norm_file_name = "../data/normalization-data/norm-specific-database.csv")
# 
# assign_mcl_to_dir(target_dir = "../data/one-sample/")
# calc_norm_factors(target_dir = "../data/one-sample/", 
#                   norm_file_name = "../data/normalization-data/norm-one-sample.csv")

# write.csv(id_peps_sub_norm_factors, "../data/db-search-id-only-normalization.csv")
# norm_pipe <- function(rep_file_list, file_out_name){
#   norm_out <- normalize_out(rep_file_list = rep_file_list)
#   annot_good <- annot_connect_good(normalize_output = norm_out)
#   annot_fuzzy <- annot_connect_fuzzy(normalize_output = norm_out)
#   write.csv(x = annot_good, file = file_out_name)
# }
