##########
## Written by Scott McCain, July 3, 2019
##########

# functions for the following:
#### peptide assignment to MCL clusters, 
#### taxonomic assigment, and 
#### normalization factor calculations

# there are lots of functions in here. refer to the readme_post_processing_functions.R
# to get an idea of the overarching modules and the specific functions goals.

library(tidyverse)
library(readxl)
library(cleaver)
library(foreach)
library(doParallel)
library(testthat)
library(GSEABase)
library(GO.db)

intersect_all <- function(a, b, ...){
  Reduce(intersect, list(a, b, ...))
}

# assign mcl clusters to groups
summarise_group <- function(pep_file_name, 
                            peptide_grouping = mcl_df){
  #############
  # goes into peptide file outputs from ProteinQuantifier
  # maps those peptides to the annotation file, assigning peptides to MCL clusters
  #############
  
  
  # pep_file_name <- "../data/pooled-database/180412_0749_097_S01_Rep_2_GOS_927_0_1_pooled.csv"
  # peptide_grouping <- mcl_df
  
  orf_id_in_df <- "orf_id" %in% colnames(peptide_grouping)
  pep_file_exists <- file.exists(pep_file_name)
  
  expect_true(orf_id_in_df)
  expect_true(pep_file_exists)
  
  print('Processing this file:')
  print(pep_file_name)
  # peptide_grouping <- mcl_df
  # pep_file_name <- "../data/mzML-converted/180412_0749_097_S01_Rep_1.csv"
  # read in file with peptide mappings to contigs with peptide intensities
  prot_quant_output <- read.table(pep_file_name, 
                                  col.names = c('peptide', 'protein', 'n_proteins', 'charge', 'abundance'), 
                                  skip = 3)
  
  # assign empty vectors to later append to the file
  mcl_assignments <- vector()
  mcl_numbers <- vector()
  
  patterns <- c("XXX", "|")
  
  for(peptide in 1:nrow(prot_quant_output)){
    # for(peptide in 1:5){
    
    # peptide <- 1
    
    mcl_assignment <- character()
    mcl_number <- numeric()
    
    contigs_prot <- prot_quant_output[peptide,]$protein %>% 
      as.character() %>% 
      strsplit(split = "/", fixed = TRUE) %>% 
      unlist()
    
    # some contigs are decoys, and some are contaminants
    contigs_prot_adj <- contigs_prot[!grepl(pattern = paste(patterns, collapse = "|"), 
                                            fixed = TRUE, 
                                            x = contigs_prot)]
    # mcl cluster assignments
    mcl_clusters <- peptide_grouping[peptide_grouping$orf_id %in% contigs_prot_adj, ]$cluster
    
    # removed NAs, how many remaining unique mcl clusters are there
    unique_mcls <- length(unique(mcl_clusters[!is.na(mcl_clusters)]))
    
    # assigns to no cluster
    if(unique_mcls == 0){
      mcl_assignment <- "no-mcl"  
      mcl_number <- 0
    }
    if(unique_mcls > 1){
      mcl_assignment <- "multi-mcl"
      mcl_number <- unique_mcls
    } else if(unique_mcls == 1){
      mcl_assignment <- unique(mcl_clusters[!is.na(mcl_clusters)])
      mcl_number <- 1
    }
    
    mcl_assignments <- c(mcl_assignments, mcl_assignment)
    mcl_numbers <- c(mcl_numbers, mcl_number)
    
    # print(mcl_assignments)
    # print(mcl_numbers)
  }
  
  prot_quant_output$mcl_assignments <- mcl_assignments
  prot_quant_output$mcl_numbers <- mcl_numbers
  
  prot_quant_output2 <- prot_quant_output
  
  # rewrites the file so that each new peptide file has two new columns
  # write.csv(x = prot_quant_output2, file = pep_file_name)
  
  return(prot_quant_output2)
  
}

get_mcl_clusters <- function(all_contigs_id,
                             peptide_grouping = mcl_df){
  
  # function for getting all mcl clusters associated with a given grouping
  
  # mcl cluster assignments
  mcl_clusters <- peptide_grouping[peptide_grouping$orf_id %in% all_contigs_id, ]$cluster
  
  # removed NAs, how many remaining unique mcl clusters are there
  unique_mcls <- length(unique(mcl_clusters[!is.na(mcl_clusters)]))
  
  # assigns to no cluster
  if(unique_mcls == 0){
    mcl_assignment <- "no-mcl"  
    mcl_number <- 0
  }
  if(unique_mcls > 1){
    mcl_assignment <- "multi-mcl"
    mcl_number <- unique_mcls
  } else if(unique_mcls == 1){
    mcl_assignment <- unique(mcl_clusters[!is.na(mcl_clusters)])
    mcl_number <- 1
  }
  
  return(list(mcl_assignment, mcl_number))
  
}

# assign go slim terms to groups
summarise_group_go_terms <- function(pep_file_name, 
                                     mcl_df_in = mcl_df,
                                     peptide_grouping = annot_with_go_terms,
                                       slim_loc = "../data/go-slim-analysis/goslim_generic.obo"){
  #############
  # goes into peptide file outputs from ProteinQuantifier
  # maps those peptides to the annotation file, assigning peptides to GO slim categories
  # note that the method of mapping GO slim categories to peptides is greedy
  #############
  # 
  # 
  # pep_file_name <- "../data/tfg-all-database/180412_0749_097_S01_Rep_1_GOS_927_0_1_tfg_all.csv"
  # peptide_grouping <- annot_with_go_terms
  
  # checking the file used is appropriate
  orf_id_in_df <- "orf_id" %in% colnames(peptide_grouping)
  expect_true(orf_id_in_df)
  
  # load the generic GO slim associations
  slim <- getOBOCollection(slim_loc)
  
  # peptide_grouping <- mcl_df
  # pep_file_name <- "../data/mzML-converted/180412_0749_097_S01_Rep_1.csv"
  
  prot_quant_output <- pep_file_name
  
  # prot_quant_output <- consensus_file
  # assign empty vectors to later append to the file
  go_assignments <- vector()
  go_numbers <- vector()
  go_descs <- vector()
  
  mcl_assignments <- vector()
  mcl_numbers <- vector()
  
  patterns <- c("XXX")
  
  for(peptide in 1:nrow(prot_quant_output)){

    # making empty variables to later add to the df
    go_assignment <- character()
    go_number <- numeric()
    go_desc <- character()
    
    contigs_prot <- prot_quant_output[peptide,]$protein %>% 
      as.character() %>% 
      strsplit(split = "/", fixed = TRUE) %>% 
      unlist()
    
    # some contigs match to decoys, but we have already passed the FDR, so filter those out
    contigs_prot_adj <- contigs_prot[!grepl(pattern = paste(patterns, collapse = "|"), 
                                            fixed = TRUE, 
                                            x = contigs_prot)]
    
    # get the mcl clusters associated with these contigs
    mcl_clusters_list <- get_mcl_clusters(all_contigs_id = contigs_prot_adj,
                                          peptide_grouping = mcl_df_in)
    # orf id assignments
    go_clusters <- peptide_grouping[peptide_grouping$orf_id %in% contigs_prot_adj, ]$GO_terms %>% as.character()
    
    # remove blank GO terms
    go_clusters_no_blanks <- remove_blank_go_terms(go_clusters)
    
    # check if there are actually go terms 
    if(length(go_clusters_no_blanks) > 0){
      # get the go terms
      go_terms_formatted <- get_go_terms(as.character(go_clusters_no_blanks))
      
      # search what the GO slim terms are
      go_slim_terms <- get_go_slim_terms(go_terms_formatted, slim = slim)
      
      # removed NAs, how many remaining unique mcl clusters are there
      unique_gos <- length(unique(go_slim_terms[!is.na(go_slim_terms)]))
    
    } else if(length(go_clusters_no_blanks) == 0){
      unique_gos <- 0
    }
    
    # assigns to no cluster
    if(unique_gos == 0){
      go_assignment <- "no-go"  
      go_number <- 0
      go_desc <- "none"
    }
    
    if(unique_gos > 0){
      go_assignment <- "multi-go"
      go_number <- unique_gos
      go_desc <- paste(go_slim_terms, collapse = ';')
      
    }
    
    go_assignments <- c(go_assignments, go_assignment)
    go_numbers <- c(go_numbers, go_number)
    go_descs <- c(go_descs, go_desc)
    
    mcl_assignments <- c(mcl_assignments, mcl_clusters_list[[1]])
    mcl_numbers <- c(mcl_numbers, mcl_clusters_list[[2]])
    
  }
  
  prot_quant_output$go_assignments <- go_assignments
  prot_quant_output$go_numbers <- go_numbers
  prot_quant_output$go_descs <- go_descs
  
  prot_quant_output$mcl_assignments <- mcl_assignments
  prot_quant_output$mcl_numbers <- mcl_numbers
  
  return(prot_quant_output)
  
}


# blah <- summarise_group_go_terms(pep_file_name = "../data/tfg-all-database/180412_0749_097_S01_Rep_1_GOS_927_0_1_tfg_all.csv")

parse_contigs <- function(protein_quant_out, row_i){
  
  ## input is the quant output from the FeatureFinderIdentification csv ProteinQUantifier (OpenMS)
  ## output is the list of ORFS associated with a given peptide in row_i
  
  contigs_prot <- protein_quant_out[row_i,]$protein %>% 
  as.character() %>% 
  strsplit(split = "/", fixed = TRUE) %>% 
  unlist()
  
  return(contigs_prot)
  
}
  
read_and_clean_file <- function(pep_file_name){
  # read in file with peptide mappings to contigs with peptide intensities
  prot_quant_output <- read.table(pep_file_name, 
                                  col.names = c('peptide', 
                                                'protein', 
                                                'n_proteins', 
                                                'charge', 
                                                'abundance'), 
                                  skip = 3)
  prot_quant_output_no_contam <- prot_quant_output[!grepl(pattern = "|", 
                                                          x = prot_quant_output$protein, 
                                                          fixed = TRUE), ]
  return(prot_quant_output_no_contam)
}

remove_blank_go_terms <- function(string_vector_go_terms){
  filtered_go_no_blanks <- string_vector_go_terms[which(string_vector_go_terms != "")]
  filtered_go_no_bp <- filtered_go_no_blanks[which(filtered_go_no_blanks != "biological_process")]
  return(filtered_go_no_bp)
}

get_go_terms <- function(string_vector_go_terms){
  ####
  #  takes the go terms and formats them for querying
  ####
  
  # if it's just one GO term ID, then, just add the "GO:" to the front of it
  # we use the sum here to indicate if *one of* the listed hits has multiple go terms associated with it
  if(sum(grepl(pattern = "||", x = string_vector_go_terms, fixed = TRUE)) == 0){
    go_line_i <- paste0("GO:", string_vector_go_terms)
    return(go_line_i)
  }
  
  # if there is more than one GO term
  else if(sum(grepl(pattern = "||", x = string_vector_go_terms, fixed = TRUE)) > 0){
    unlisted_go_terms <- unlist(strsplit(x = string_vector_go_terms,
                                         "||", 
                                         fixed = TRUE))
    
    # putting the "GO:" string before
    for(j in 1:length(unlisted_go_terms)){
      unlisted_go_terms[j] <- paste0("GO:", unlisted_go_terms[j])
    }
    return(unlisted_go_terms)
    }
}

get_go_slim_terms <- function(formatted_go_terms, slim){
  
  #####
  ## inputs are the go terms, outputs are the go slim terms
  ####
  
  myCollection <- GOCollection(formatted_go_terms)
  go_out <- goSlim(myCollection, slim = slim, ontology = 'BP')
  go_slim_terms <- as.character(go_out[go_out$Percent > 0, ]$Term)
  
  return(go_slim_terms)
  
}

get_file_list <- function(target_dir){
  
  expect_is(target_dir, "character")
  expect_true(dir.exists(target_dir))
  
  mzmldir <- dir(target_dir)
  csvfiles <- mzmldir[grep(".csv", x = mzmldir, fixed = TRUE)]
  # csvfilesfdr <- mzmldir[grep("_FDR.csv", x = mzmldir, fixed = TRUE)]
  csvdir <- paste0(target_dir, csvfiles)
  
  return(csvdir)
  
}

assign_mcl_to_dir <- function(target_dir, peptide_grouping = mcl_df){
  
  ###############
  # goes through a directory, gets the csv files from ProteinQUantifier
  # loops through rewriting the peptide summary files and rewrites the files
  ###############
  csvdir <- get_file_list(target_dir = target_dir)
  
  # go through each file and rewrite it, with two extra columns
  for(i in 1:length(csvdir)){
    summarise_group(pep_file_name = csvdir[i], 
                    peptide_grouping = peptide_grouping)
  }
}

#normalization function
normalize_out <- function(rep_file_list){
  
  print('Running normalization...')
  
  #########
  # Normalizes files across replicates
  #########
  
  # read in files with MCL assignments
  file_list <- list()
  for(file in 1:length(rep_file_list)){
    # file_list[[file]] <- read.csv(rep_file_list[file])
    file_list[[file]] <- read_and_clean_file(rep_file_list[file])
  }
  
  # make peptide list per file for determining common set of peptides across runs
  pep_list <- list()
  
  for(file in 1:length(file_list)){
    pep_list[[file]] <- file_list[[file]]$peptide
  }
  
  # find the common set of peptides across runs
  common_peps <- Reduce(intersect, pep_list)
  
  common_file_list <- list()
  
  # subset only common-across-runs peptides (to compare across samples), 
  # calculate db-dependent norm factor, apply norm factor
  for(file in 1:length(rep_file_list)){
    
    temp_df <- file_list[[file]]
    temp_df2 <- temp_df[temp_df$peptide %in% common_peps,]
    temp_df2$peptide <- temp_df2$peptide %>% as.character()
    
    # exclude contaminants from normalization factor calculation
    temp_df_no_contam <- temp_df[!grepl(pattern = "|", 
                                        x = temp_df$protein,
                                        fixed = TRUE), ]
    
    # Normalization factor calculated from all peptides identified (+ not common across MS runs)
    db_dep_norm_factor <- sum(temp_df_no_contam$abundance)
    temp_df2$db_norm_abund <- temp_df2$abundance/db_dep_norm_factor
    common_file_list[[file]] <- temp_df2
  }
  
  # making two matrices that have no normalization, and that has normalization  
  no_normal_abund_matrix <- matrix(nrow = length(common_peps), 
                                   ncol = length(rep_file_list))
  abund_matrix <- matrix(nrow = length(common_peps), 
                         ncol = length(rep_file_list))
  
  for(file in 1:length(rep_file_list)){
    temp_df <- common_file_list[[file]]
    no_normal_abund_matrix[, file] <- temp_df$abundance
    abund_matrix[ ,file] <- temp_df$db_norm_abund
  }
  
  # getting the column of names which will later be used.
  temp_file_no_abund <- dplyr::select(.data = common_file_list[[1]], 
                                      -c(db_norm_abund, abundance))
  
  # this does not change to the sum, it averages the abundance across
  abund_no_norm_avg <- rowMeans(no_normal_abund_matrix)
  db_abund_avg <- rowMeans(abund_matrix)
  
  consensus_file <- cbind(temp_file_no_abund, db_abund_avg, abund_no_norm_avg)
  
  return(consensus_file)
}

# proteotypic taxonomy
tax_assigner <- function(normalize_output){
  
  print('Assigning taxonomic IDs to peptides...')
  
  quant_out <- normalize_output
  # quant_out <- s01
  
  tax_assignment_vec <- vector()
  grpnorm_assignment_vec <- vector()
  
  for(peptide in 1:nrow(quant_out)){
    # print(peptide)
    # print(quant_out[peptide,])
    # peptide <- 1
    # list of all contigs assignment to this peptide
    contigs_prot <- quant_out[peptide,]$protein %>% 
      as.character() %>% 
      strsplit(split = "/", fixed = TRUE) %>% 
      unlist()
    
    # getting all taxonomic strings associated with set of orfs
    # tax_string <- mcl_df[mcl_df$orf_id %in% contigs_prot,]$best_tax_string
    # 
    # if(length(tax_string) == 0){
    #   tax_assignment <- c(tax_assignment, "no-assignment-contig")
    #   next
    # }
    # 
    # # parsing the tax string
    # tax_string_split <- strsplit(x = tax_string, split = ";")
    # # the level of taxonomic assignment
    # class_assignment <- mapply(tax_string_split, FUN = '[', 2) %>% unique()
    
    class_assignment <- mcl_df[mcl_df$orf_id %in% contigs_prot,]$group %>% unique()
    grpnorm_assignment <- mcl_df[mcl_df$orf_id %in% contigs_prot,]$grpnorm_taxgrp %>% unique()
    
    # assigning to class
    if(length(class_assignment) == 0){
      tax_assignment_vec <- c(tax_assignment_vec, "no-assignment-contig")
    }
    else if(length(class_assignment) == 1){
      tax_assignment_vec <- c(tax_assignment_vec, class_assignment)
    }
    else if(length(class_assignment) > 1){
      tax_assignment_vec <- c(tax_assignment_vec, "promiscuous-peptide")
    }
    
    # assigning to grpnorm_taxgrp
    if(length(grpnorm_assignment) == 0){
      grpnorm_assignment_vec <- c(grpnorm_assignment_vec, "no-assignment-contig")
    }
    else if(length(grpnorm_assignment) == 1){
      grpnorm_assignment_vec <- c(grpnorm_assignment_vec, grpnorm_assignment)
    }
    else if(length(grpnorm_assignment) > 1){
      grpnorm_assignment_vec <- c(grpnorm_assignment_vec, "promiscuous-peptide")
    }
  }
  
  quant_out$tax_assign <- tax_assignment_vec
  quant_out$grpnorm_assign <- grpnorm_assignment_vec
  
  return(quant_out)
}

# # function for summarizing mcl clusters and taxonomic groups together
annot_connect_tax <- function(normalize_output, 
                              tax_level = 'tax'){
  
  # mcl_in_df <- "mcl_assignments" %in% colnames(normalize_output)
  # expect_true(mcl_in_df)
  # 
  # mcl_assignments_in_df <- "mcl_assignments" %in% colnames(normalize_output)
  # expect_true(mcl_assignments_in_df)
  
  tax_in_df <- "tax_assign" %in% colnames(normalize_output)
  expect_true(tax_in_df)
  
  if(tax_level == 'tax'){
    tax_summarized_df <- normalize_output %>%
      group_by(tax_assign) %>%
      summarise(db_mcl_sum = sum(db_abund_avg),
                nn_mcl_sum = sum(abund_no_norm_avg),
                number_peps = n(),
                db_mcl_mean = mean(db_abund_avg),
                nn_mcl_mean = mean(abund_no_norm_avg))
  }
  
  if(tax_level == 'grpnorm'){
    tax_summarized_df <- normalize_output %>%
      group_by(grpnorm_assign) %>%
      summarise(db_mcl_sum = sum(db_abund_avg),
                nn_mcl_sum = sum(abund_no_norm_avg),
                number_peps = n(),
                db_mcl_mean = mean(db_abund_avg),
                nn_mcl_mean = mean(abund_no_norm_avg))
  }
  
  return(tax_summarized_df)
}


convert_to_vec_go_descs <- function(string_of_go_descs){
  # connecting go terms to peptide abundances (following two functions)
  vec_of_go_descs <- unlist(strsplit(string_of_go_descs, ";"))
  return(vec_of_go_descs)
}

annot_go_categories <- function(tax_assign_out){
  # get the go annotation categories and add the peptide intensitis to them
  # check that there are no NA's in the go description
  is_na_in_go_descs <- is.na(tax_assign_out$go_descs) %>% sum()
  expect_true(is_na_in_go_descs == 0)
  
  master_df <- data.frame(go_slim_category = character(),
                          db_abund_avg = numeric(),
                          tax_assign = numeric(),
                          grpnorm_assign = numeric())
  
  # go through each peptide
  for(i in 1:nrow(tax_assign_out)){
    # i <- 3
    # collect the go_descs and the peptide_abundance
    
    vec_of_go_descs <- convert_to_vec_go_descs(string_of_go_descs = tax_assign_out[i, ]$go_descs)
    vec_of_pep_abundances <- rep(tax_assign_out[i, ]$db_abund_avg, length(vec_of_go_descs))
    
    vec_of_pep_tax <- rep(tax_assign_out[i, ]$tax_assign, length(vec_of_go_descs))
    vec_of_pep_grpnorm <- rep(tax_assign_out[i, ]$grpnorm_assign, length(vec_of_go_descs))
    
    sub_df <- data.frame(go_slim_category = vec_of_go_descs,
                         db_abund_avg = vec_of_pep_abundances,
                         tax_assign = vec_of_pep_tax,
                         grpnorm_assign = vec_of_pep_grpnorm)
    
    master_df <- rbind(master_df, sub_df)
    
  }
  return(master_df)
}

get_mcl_weights <- function(tax_assign_out){
  # this function gets the mcl/species specific peptides, subsets them, and then applies the correction factors
  
  # subset the tax/mcl only peptides
  tax_mcl_only_peps <- tax_assign_out %>% 
    filter(tax_assign != 'promiscuous-peptide',
           mcl_assignments != 'multi-mcl',
           mcl_assignments != 'no-mcl')
  
  # get the mcl informative peptides
  mcl_only_peps <- tax_assign_out %>% 
    filter(mcl_assignments != 'multi-mcl',
           mcl_assignments != 'no-mcl') %>% 
    group_by(mcl_assignments) %>% 
    summarize(number_of_peps = n(),
              min_pep_intensity_mcl = min(db_abund_avg),
              sum_pep_intensity_mcl = sum(db_abund_avg))
  
  return(list(mcl_only_peps, tax_mcl_only_peps))
}

get_tax_go_weights_apply <- function(go_master_out, tax_assign_out){
  
  ##### get the go data categories of peptides that are species specific
  tax_only_peps <- go_master_out %>% 
    filter(tax_assign != 'promiscuous-peptide')
  
  ##### calculate weighting factors for applying to data
  go_correction_df <- go_master_out %>% 
    group_by(go_slim_category) %>% 
    summarize(number_of_peps_go = n(),
              min_pep_intensity_go = min(db_abund_avg),
              sum_pep_intensity_go = sum(db_abund_avg))
  
  mcl_weights_list <- get_mcl_weights(tax_assign_out = tax_assign_out)
  mcl_weights_df <- mcl_weights_list[[1]]
  mcl_weights_tax <- mcl_weights_list[[2]]  
  
  tax_correction_df <- tax_assign_out %>% 
    group_by(tax_assign) %>% 
    summarize(number_of_peps_tax = n(),
              min_pep_intensity_tax = min(db_abund_avg),
              sum_pep_intensity_tax = sum(db_abund_avg))
  
  # connect the go calcs with the tax calcs
  tax_go_corrected <- tax_only_peps %>% 
    inner_join(go_correction_df, by = 'go_slim_category') %>% 
    inner_join(tax_correction_df, by = 'tax_assign') %>% 
    mutate(db_abund_avg_corrected = db_abund_avg*sum_pep_intensity_tax*sum_pep_intensity_go)
  
  tax_mcl_corrected <- mcl_weights_tax %>% 
    inner_join(mcl_weights_df, by = 'mcl_assignments') %>% 
    inner_join(tax_correction_df, by = 'tax_assign') %>% 
    mutate(db_abund_avg_corrected = db_abund_avg*sum_pep_intensity_tax*sum_pep_intensity_mcl)
    
  # return four objects. the first is the tax/go correction term applied abundances, and 
  # the second is the tax summarized dataframe, the third is the go summarized dataframe
  return(list(tax_go_corrected, tax_correction_df, go_correction_df, tax_mcl_corrected))
  
}

calc_norm_factors <- function(target_dir, norm_file_name, save_file = TRUE){
  
  rep_file_list <- get_file_list(target_dir = target_dir)  
  # rep_file_list <- csvdir
  
  file_list <- list()
  for(file in 1:length(rep_file_list)){
    
    print(rep_file_list[file])
    
    file_list[[file]] <- read_and_clean_file(rep_file_list[file])
  }
  
  # make peptide list per file for determining common set of peptides across runs
  norm_factors <- vector()
  
  for(file in 1:length(file_list)){
    print(names(file_list[[file]]))
    # has to be changed for mean vs. sum abundance calculations
    norm_factor_i <- file_list[[file]]$abundance %>% sum(na.rm = TRUE)
    print(norm_factor_i)
    norm_factors <- c(norm_factors, norm_factor_i)
  }
  
  df1 <- data.frame(rep_file_list, norm_factors)
  
  if(save_file){
    write.csv(x = df1, file = norm_file_name)
    return(df1)

  }
  
  if(!save_file){
    return(df1)
  }
}

id_peps_norm_factors <- function(rep_file_list_id, rep_file_list_extended){
  
  # rep_file_list_id <- csvfdr
  # rep_file_list_extended <- csvdir
  
  file_list <- list()
  
  norm_factors <- vector()
  
  for(file in 1:length(rep_file_list_id)){
    
    print(rep_file_list_id[file])
    
    id_df <- read.table(rep_file_list_id[file], 
                        col.names = c('peptide', 'protein', 'n_proteins', 'charge', 'abundance'), 
                        skip = 3)
    extended_df <- read.csv(rep_file_list_extended[file])
    
    extended_df_sub <- extended_df[extended_df$peptide %in% id_df$peptide, ]
    
    # has to be changed for mean vs. sum
    norm_factor_i <- extended_df_sub$abundance %>% sum()
    norm_factors <- c(norm_factors, norm_factor_i)
  }
  
  temp_df <- data.frame(rep_file_list_id, norm_factors)
  
  return(temp_df)
}

# function to format the file name input, which has pairs of replicate MS injections
append_dir_to_files <- function(file_list, input_directory){
  
  expect_type(file_list, "list")  
  
  for(number_unique_samples in 1:length(file_list)){
    
    sub_list <- file_list[[number_unique_samples]]
    
    for(sample_rep in 1:length(sub_list)){
      new_file_name <- paste0(input_directory, sub_list[sample_rep])
      new_file_name_no_csv <- gsub(x = sub_list[sample_rep], 
                                   pattern = ".csv", 
                                   replacement = "", fixed = TRUE)
      # print(new_file_name_no_csv)
      actual_file_name_index <- grepl(pattern = new_file_name_no_csv, 
                                      x = dir(input_directory), fixed = TRUE)
      # print(sum(actual_file_name_index))
      new_actual_file_name <- dir(input_directory)[which(actual_file_name_index == TRUE)]
      # print(new_actual_file_name)
      
      get_csv_file <- new_actual_file_name[grep(pattern = ".csv", x = new_actual_file_name, fixed = TRUE)]
      print(get_csv_file)
      file_list[[number_unique_samples]][sample_rep] <- paste0(input_directory, get_csv_file)
    }
  }
  return(file_list)
}

make_master_dfs_for_reporting <- function(){
  
  tax_only_master <- data.frame(tax_assign = character(),
                                number_of_peps_tax = integer(),
                                min_pep_intensity_tax = numeric(),
                                sum_pep_intensity_tax = numeric(),
                                day = numeric(),
                                filter = character())
  
  go_only_master <- data.frame(go_slim_category = factor(),
                               number_of_peps_go = integer(),
                               min_pep_intensity_go = numeric(),
                               sum_pep_intensity_go = numeric(),
                               day = numeric(),
                               filter = character())
  
  corrected_only_master <- data.frame(go_slim_category = factor(),
                                      db_abund_avg = numeric(),
                                      tax_assign = character(),
                                      grpnorm_assign = factor(),
                                      number_of_peps_go = integer(),
                                      min_pep_intensity_go = numeric(),
                                      sum_pep_intensity_go = numeric(),
                                      number_of_peps_tax = integer(),
                                      min_pep_intensity_tax = numeric(),
                                      sum_pep_intensity_tax = numeric(),
                                      db_abund_avg_corrected = numeric(),
                                      day = numeric(),
                                      filter = character())
  
  peptides_unique_taxa <- data.frame(peptide = character(),
                                     protein = factor(),
                                     n_proteins = integer(),
                                     charge = integer(),
                                     db_abund_avg = numeric(),
                                     abund_no_norm_avg = numeric(),
                                     go_assignments = character(),
                                     go_numbers = numeric(),
                                     go_descs = character(),
                                     mcl_assignments = character(),
                                     mcl_numbers = character(),
                                     tax_assign = character(),
                                     grpnorm_assign = character(),
                                     day = numeric(),
                                     filter = character())
  
  mcl_summary_df_master <- data.frame(mcl_assignments = factor(),
                                      db_abund_avg = numeric(),
                                      tax_assign = character(),
                                      grpnorm_assign = factor(),
                                      number_of_peps_mcl = integer(),
                                      min_pep_intensity_mcl = numeric(),
                                      sum_pep_intensity_mcl = numeric(),
                                      number_of_peps_tax = integer(),
                                      min_pep_intensity_tax = numeric(),
                                      sum_pep_intensity_tax = numeric(),
                                      db_abund_avg_corrected = numeric(),
                                      day = numeric(),
                                      filter = character())

  return(list(tax_only_master, go_only_master, corrected_only_master, 
              peptides_unique_taxa, mcl_summary_df_master))
}

write_tax_results <- function(formatted_file_name_list, id_string){
  
  # for each input_file_reps, (3 refers to the number of filter sizes)
  counter_i <- 0
  
  sample_grouping_vals <- c('0.1', '0.8', '3.0')
  
  # create 3 master dataframes that will be the output for go terms only, 
  # tax only, and tax/go corrected values
  
  all_dfs <- make_master_dfs_for_reporting()
  
  tax_only_master <- all_dfs[[1]]
  go_only_master <- all_dfs[[2]]
  corrected_only_master <- all_dfs[[3]]
  peptides_unique_taxa_master <- all_dfs[[4]]
  mcl_summary_df_master <- all_dfs[[5]]
  
  # 3 refers to the number of filter sizes
  for(sample_grouping in 1:3){
    
    # sample_grouping <- 2
    
    # 4 refers to the number of sampling days
    for(i in 1:4){
      
      # i <- 1
      # sample_grouping <- 1
      # formatted_file_name_list <- append_dir_to_files(file_list = input_file_reps,
      #                                                 input_directory = "../data/tfg-all-database/")
      # counter_i <- 4
      
      # normalizes each pair of replicate injections by the mean ion intensity
      norm_out_files <- normalize_out(rep_file_list = c(formatted_file_name_list[[1 + counter_i]]))
      print(c(formatted_file_name_list[[1 + counter_i]]))
      
      ## add go slim terms
      go_slim_out <- summarise_group_go_terms(pep_file_name = norm_out_files)
      
      # assigns taxonomy to each peptide
      tax_assign_out <- tax_assigner(go_slim_out)
      
      # connecting the tax assigned results with functional annotations
      # tax_df <- annot_connect_tax(tax_assign_out, tax_level = 'tax')
      # grpnorm_df <- annot_connect_tax(tax_assign_out, tax_level = 'grpnorm')
      
      # annot_go_categories
      go_annotated_df <- annot_go_categories(tax_assign_out = tax_assign_out)
      
      # applying the correction factors, but also returning the taxonomic-only peptides and go-only peptides
      corrected_tax_go_df_and_dfs <- get_tax_go_weights_apply(go_master_out = go_annotated_df,
                                                              tax_assign_out = tax_assign_out)
      
      # get the relevant objects into dataframes
      corrected_tax_go_df <- corrected_tax_go_df_and_dfs[[1]]
      tax_summary_df <- corrected_tax_go_df_and_dfs[[2]]
      go_summary_df <- corrected_tax_go_df_and_dfs[[3]]
      mcl_summary_df <- corrected_tax_go_df_and_dfs[[4]]

      # get a column that refers to the sampling day      
      corrected_tax_go_df$day <- rep(i, nrow(corrected_tax_go_df))
      tax_summary_df$day <- rep(i, nrow(tax_summary_df))
      go_summary_df$day <- rep(i, nrow(go_summary_df))
      tax_assign_out$day <- rep(i, nrow(tax_assign_out))
      mcl_summary_df$day <- rep(i, nrow(mcl_summary_df))
      
      # make a column that refers to the sample filter size
      corrected_tax_go_df$filter <- rep(sample_grouping_vals[sample_grouping], 
                                        nrow(corrected_tax_go_df))
      tax_summary_df$filter <- rep(sample_grouping_vals[sample_grouping], 
                                   nrow(tax_summary_df))
      go_summary_df$filter <- rep(sample_grouping_vals[sample_grouping], 
                                  nrow(go_summary_df))
      tax_assign_out$filter <- rep(sample_grouping_vals[sample_grouping],
                                   nrow(tax_assign_out))
      mcl_summary_df$filter <- rep(sample_grouping_vals[sample_grouping],
                                   nrow(mcl_summary_df))
      
      # making the names informative
      # corrected_tax_go_df <- corrected_tax_go_df %>%
      #   dplyr::rename(number_of_peps_go = number_of_peps.x,
      #          min_pep_intensity_go = min_pep_intensity.x,
      #          sum_pep_intensity_go = sum_pep_intensity.x,
      #          number_of_peps_tax = number_of_peps.y,
      #          min_pep_intensity_tax = min_pep_intensity.y,
      #          sum_pep_intensity_tax = sum_pep_intensity.y)
      # 
      # tax_summary_df <- tax_summary_df %>% 
      #   dplyr::rename(number_of_peps_tax = number_of_peps,
      #         min_pep_intensity_tax = min_pep_intensity,
      #         sum_pep_intensity_tax = sum_pep_intensity)
      # 
      # go_summary_df <- go_summary_df %>% 
      #   dplyr::rename(number_of_peps_go = number_of_peps,
      #          min_pep_intensity_go = min_pep_intensity,
      #          sum_pep_intensity_go = sum_pep_intensity)
      # 
      # mcl_summary_df <- mcl_summary_df %>%
      #   dplyr::rename(number_of_peps_mcl = number_of_peps.x,
      #          min_pep_intensity_mcl = min_pep_intensity.x,
      #          sum_pep_intensity_mcl = sum_pep_intensity.x,
      #          number_of_peps_tax = number_of_peps.y,
      #          min_pep_intensity_tax = min_pep_intensity.y,
      #          sum_pep_intensity_tax = sum_pep_intensity.y)
      
      
      # append tthese processed dataframes to the master dfs
      go_only_master <- rbind(go_only_master, go_summary_df)
      tax_only_master <- rbind(tax_only_master, tax_summary_df)
      corrected_only_master <- rbind(corrected_only_master, corrected_tax_go_df)
      peptides_unique_taxa_master <- rbind(peptides_unique_taxa_master, tax_assign_out)
      mcl_summary_df_master <- rbind(mcl_summary_df_master, mcl_summary_df)
      
      counter_i <- counter_i + 1
    }

  }

  go_only_master$id_string <- rep(id_string, nrow(go_only_master))
  tax_only_master$id_string <- rep(id_string, nrow(tax_only_master))
  corrected_only_master$id_string <- rep(id_string, nrow(corrected_only_master))
  peptides_unique_taxa_master$id_string <- rep(id_string, nrow(peptides_unique_taxa_master))
  mcl_summary_df_master$id_string <- rep(id_string, nrow(mcl_summary_df_master))
  
  return(list(go_only_master, tax_only_master, corrected_only_master, 
              peptides_unique_taxa_master, mcl_summary_df_master))
  
}

save_pipeline_output <- function(tax_results_out,
                                 id_string){
  
  write.csv(tax_results_out[[1]], file = paste0("pipeline_out", id_string, "_go_only.csv"))
  write.csv(tax_results_out[[2]], file = paste0("pipeline_out", id_string, "_tax_only.csv"))
  write.csv(tax_results_out[[3]], file = paste0("pipeline_out", id_string, "_corrected_only.csv"))
  write.csv(tax_results_out[[4]], file = paste0("pipeline_out", id_string, "_tax_peptides.csv"))
  write.csv(tax_results_out[[5]], file = paste0("pipeline_out", id_string, "_mcl_corrected_peptides.csv"))

}

run_pipeline <- function(file_list = input_file_reps,
                         input_directory = input_directory, 
                         id_string = id_string){
  all_file_names_formatted <- append_dir_to_files(file_list = file_list, 
                                                  input_directory = input_directory)
  all_results <- write_tax_results(formatted_file_name_list = all_file_names_formatted,
                    id_string = id_string)
  
  save_pipeline_output(tax_results_out = all_results,
                       id_string = id_string)
  
}

######### mapping data to coarse grained pools from model

get_all_descriptions_for_contigs <- function(parsed_contigs_out){
  
  # collects all the descriptions associated with the contigs
  
  vector_of_desc_master <- c()
  
  vector_kegg <- c()
  vector_ko <- c()
  vector_kogg <- c()
  vector_pfams <- c()
  vector_tigrfams <- c()
  
  # now this function is greedy. as soon as it identifies a coarse grained match,
  # it assigns it to the protein coarse grain.
  for(j in 1:length(parsed_contigs_out)){
    # print(j)
    # j <- 1
    kegg_desc_contig <- mcl_df[mcl_df$orf_id == parsed_contigs_out[j], ]$kegg_desc
    ko_desc_contig <- mcl_df[mcl_df$orf_id == parsed_contigs_out[j], ]$KO_desc
    kogg_desc_contig <- mcl_df[mcl_df$orf_id == parsed_contigs_out[j], ]$KOG_desc
    pfams_desc_contig <- mcl_df[mcl_df$orf_id == parsed_contigs_out[j], ]$PFams_desc
    tigrframs_desc_contig <- mcl_df[mcl_df$orf_id == parsed_contigs_out[j], ]$TIGRFams_desc
    
    vector_of_desc <- c(kegg_desc_contig,
                        ko_desc_contig,
                        kogg_desc_contig,
                        pfams_desc_contig,
                        tigrframs_desc_contig)
    
    vector_kegg <- c(vector_kegg, kegg_desc_contig)
    vector_ko <- c(vector_ko, ko_desc_contig)
    vector_kogg <- c(vector_kogg, kogg_desc_contig)
    vector_pfams <- c(vector_pfams, pfams_desc_contig)
    vector_tigrfams <- c(vector_tigrfams, tigrframs_desc_contig)
    
    # print(temp_match)
    # print(unique(temp_match))
    # print(vector_of_desc)
    vector_of_desc_master <- c(vector_of_desc_master, vector_of_desc)
    # checking to see if there are any non-U values
  }
  
  desc_df_all <- data.frame(kegg_desc = paste(unique(vector_kegg), collapse = '____'),
                            ko_desc = paste(unique(vector_ko), collapse = '____'),
                            kogg_desc = paste(unique(vector_kogg), collapse = '____'),
                            pfams_desc = paste(unique(vector_pfams), collapse = '____'),
                            tigrfams_desc = paste(unique(vector_tigrfams), collapse = '____'))
  
  return(list(vector_of_desc_master, desc_df_all))
  
}

check_length_model_grains <- function(unique_coarse_grain){
  
  # intermediate function that takes the length of coarse grain assignments 
  # and processes them to determine which assignment it should be
  
  # if they are not all identical, assign the peptide to 'U'
  if(length(unique_coarse_grain) > 1){
    coarse_assignment <- "U"
  }
  # if they are all identical, then assign it
  if(length(unique_coarse_grain) == 1){
    coarse_assignment <- unique_coarse_grain
  }
  # if there are none, then assign it to a 0
  if(length(unique_coarse_grain) == 0){
    coarse_assignment <- "U"
  }
  
  return(coarse_assignment)
  
}

match_desc_to_model_grains <- function(description_vec_out){
  
  # matches description of output to model
  
  # if there are no descriptions, then return
  if(length(description_vec_out) == 0){
    coarse_assignment <- 'U'
  }
  
  model_coarse_grain <- c()
  
  if(length(description_vec_out) > 0){
  # matching descriptions to coarse grained models
  # goes through each description in the vector
  for(i in 1:length(description_vec_out)){
      
      antioxidant_string <- c("MnSOD", "superoxide dismutase")
      if(grepl(paste(antioxidant_string, collapse = "|"), 
               description_vec_out[i], ignore.case = TRUE)){
        model_coarse_grain <- c(model_coarse_grain, "A")
      }
      
      ribosome_string <- c("ribosom")
      ribosome_anti_string <- c("synthesis")
      
      if(grepl(paste(ribosome_string, collapse = "|"), 
               description_vec_out[i], ignore.case = TRUE) &
         !grepl(paste(ribosome_anti_string, collapse = "|"), 
                description_vec_out[i], ignore.case = TRUE)){
        model_coarse_grain <- c(model_coarse_grain, "R")
      }
      
      photosynthesis_string <- c("photosynthesis", "photosystem", 
                                 "light harvesting protein", "Chlorophyll A-B binding",
                                 "ATP synthase alpha chain", "ATP synthase",
                                 "plastocyanin", "K03839 flavodoxin I",
                                 "oxygen-evolving enhancer protein")
      if(grepl(paste(photosynthesis_string, collapse = "|"), 
               description_vec_out[i], ignore.case = TRUE)){
        model_coarse_grain <- c(model_coarse_grain, "P")
      }
      
      n_metabolism_string <- c("NRT2", "nitrate transporter",
                               "nitrate reductase", "nitrite reductase", 
                               "glutamine biosynthe", "glutamate biosynthe",
                               "glutamine synth", "glutamate synth")
      if(grepl(paste(n_metabolism_string, collapse = "|"), 
               description_vec_out[i], ignore.case = TRUE)){
        model_coarse_grain <- c(model_coarse_grain, "Tn")
      }
      
      fe_transporter_string <- c("ISIP2a", "ISIP1", "iron permease", "Fe permease", "FTR1",
                                 "iron ion transport", "Fe ion transport", 
                                 "ferriportin", "Fe3+ hydroxamate transport",
                                 "ZIP Zn/Fe", "Ferric reductase", "sideroflexin", 
                                 "Fe transport", "Fe transporter", "Iron transport",
                                 "Iron transporter")
      
      if(grepl(paste(fe_transporter_string, collapse = "|"), 
               description_vec_out[i], ignore.case = TRUE)){
        model_coarse_grain <- c(model_coarse_grain, "Tfe")
      }
      
      mn_transporter_string <- c("NRAMP", "Mn transport", "Mn transporter", 
                                 "Manganese transporter", "Manganese transport")
      if(grepl(paste(mn_transporter_string, collapse = "|"), 
               description_vec_out[i], ignore.case = TRUE)){
        model_coarse_grain <- c(model_coarse_grain, "Tmn")
      }
      # print(model_coarse_grain)
      # model_coarse_grain <- c("U", "R")
    }
    
    if(length(model_coarse_grain) == 0){
      model_coarse_grain <- 'U'
    }
    
    # check that all the coarse grains are identical
    unique_coarse_grain <- unique(model_coarse_grain)
    
    coarse_assignment <- check_length_model_grains(unique_coarse_grain = unique_coarse_grain)
    
  }
  
    return(coarse_assignment)
    
}

format_desc_df <- function(desc_out_intermediate){
  # need to add NA's if there are no annotations for any orfs returned  
  if(nrow(desc_out_intermediate) == 0){
    desc_out_intermediate[nrow(desc_out_intermediate)+1,] <- NA
  }
  return(desc_out_intermediate)
}

get_coarse_grains <- function(tax_group, tax_assign_df, filter_type, all_taxa = FALSE){
  
  tax_assign_in_df <- 'tax_assign' %in% colnames(tax_assign_df)
  expect_true(tax_assign_in_df)
  
  # tax_group <- 'Diatom'
  # all_taxa <- TRUE
  # tax_assign_df <- tax_data
  # filter_type <- 3.0
  # 
  
  if(!all_taxa){
  org_peps <- tax_assign_df %>% 
      dplyr::filter(tax_assign == tax_group, 
                    filter == filter_type)
  }
  if(all_taxa){
  org_peps <- tax_assign_df %>% 
      dplyr::filter(filter == filter_type)
  }
  
  all_coarse_assignments <- c()
  
  all_associated_desc <- data.frame(kegg_desc = character(),
                            ko_desc = character(),
                            kogg_desc = character(),
                            pfams_desc = character(),
                            tigrfams_desc = character())
  
  for(pep_i in 1:nrow(org_peps)){
    
    print(pep_i)
    # print(diatom_peps[i, ]$go_descs)
    # parse the output, getting all associated contigs associated with the peptide
    contigs_parsed <- parse_contigs(protein_quant_out = org_peps,
                                    row_i = pep_i)
    # for each contig, get a vector of descriptions
    all_descs <- get_all_descriptions_for_contigs(parsed_contigs_out = contigs_parsed)
    
    # match descriptions to the model output
    coarse_assignment_val <- match_desc_to_model_grains(all_descs[[1]])
    
    desc_df <- format_desc_df(all_descs[[2]])
    
    # appending the description to the 
    all_associated_desc <- rbind(all_associated_desc, desc_df)
    all_coarse_assignments <- c(all_coarse_assignments, coarse_assignment_val)
    
  }
  
  org_peps$coarse_grains <- all_coarse_assignments
  org_peps <- cbind(org_peps, all_associated_desc)
  # org_peps$coarse_grain_descs <- all_associated_desc
  
  # if(!all_taxa){
  # org_peps$all_coarse_descs <- all_associated_desc
  # }
  
  return(org_peps)
}


# 
# get_all_descriptions_for_contigs_specific <- function(parsed_contigs_out){
#   
#   vector_of_desc_master <- c()
#   
#   for(j in 1:length(parsed_contigs_out)){
#     # print(j)
#     # j <- 1
#     kegg_desc_contig <- mcl_df[mcl_df$orf_id == parsed_contigs_out[j], ]$kegg_desc
#     ko_desc_contig <- mcl_df[mcl_df$orf_id == parsed_contigs_out[j], ]$KO_desc
#     kogg_desc_contig <- mcl_df[mcl_df$orf_id == parsed_contigs_out[j], ]$KOG_desc
#     pfams_desc_contig <- mcl_df[mcl_df$orf_id == parsed_contigs_out[j], ]$PFams_desc
#     tigrframs_desc_contig <- mcl_df[mcl_df$orf_id == parsed_contigs_out[j], ]$TIGRFams_desc
#     
#     vector_of_desc <- c(kegg_desc_contig,
#                         ko_desc_contig,
#                         kogg_desc_contig,
#                         pfams_desc_contig,
#                         tigrframs_desc_contig)
#     vector_of_desc_master <- c(vector_of_desc_master, vector_of_desc)
#     
#   }
#   
#   return(vector_of_desc_master)
# }
# 
# match_desc_to_model_grains_specific <- function(description_vec_out){
#   
#   # print('in match desc func')
#   model_coarse_grain <- c()
#     
#     # matching descriptions to coarse grained models
#     # model_coarse_grain <- c("U")
#     for(i in 1:length(description_vec_out)){
#       
#       antioxidant_string <- c("MnSOD", "superoxide dismutase")
#       if(grepl(paste(antioxidant_string, collapse = "|"), 
#                description_vec_out[i], ignore.case = TRUE)){
#         model_coarse_grain <- c(model_coarse_grain, "A")
#       }
#       
#       ribosome_string <- c("ribosom")
#       ribosome_anti_string <- c("synthesis")
#       if(grepl(paste(ribosome_string, collapse = "|"), 
#                description_vec_out[i], ignore.case = TRUE) &
#          !grepl(paste(ribosome_anti_string, collapse = "|"), 
#                 description_vec_out[i], ignore.case = TRUE)){
#         model_coarse_grain <- c(model_coarse_grain, "R")
#       }
#       
#       photosynthesis_string <- c("photosynthesis", "photosystem", "light harvesting protein", "Chlorophyll A-B binding")
#       if(grepl(paste(photosynthesis_string, collapse = "|"), 
#                description_vec_out[i], ignore.case = TRUE)){
#         model_coarse_grain <- c(model_coarse_grain, "P")
#       }
#       
#       n_metabolism_string <- c("NRT2", "nitrate transporter",
#                                "nitrate reductase", "nitrite reductase", 
#                                "glutamine biosynthe", "glutamate biosynthe",
#                                "glutamine synth", "glutamate synth")
#       if(grepl(paste(n_metabolism_string, collapse = "|"), 
#                description_vec_out[i], ignore.case = TRUE)){
#         model_coarse_grain <- c(model_coarse_grain, "Tn")
#       }
#       
#       fe_transporter_string <- c("ISIP2a", "ISIP1", "iron permease", "Fe permease", "FTR1",
#                                  "iron ion transport", "Fe ion transport", 
#                                  "ferriportin", "Fe3+ hydroxamate transport",
#                                  "ZIP Zn/Fe", "Ferric reductase", "sideroflexin", 
#                                  "Fe transport", "Fe transporter", "Iron transport",
#                                  "Iron transporter")
#       if(grepl(paste(fe_transporter_string, collapse = "|"), 
#                description_vec_out[i], ignore.case = TRUE)){
#         model_coarse_grain <- c(model_coarse_grain, "Tfe")
#       }
#       
#       mn_transporter_string <- c("NRAMP", "Mn transport", "Mn transporter", 
#                                  "Manganese transporter", "Manganese transport")
#       if(grepl(paste(mn_transporter_string, collapse = "|"), 
#                description_vec_out[i], ignore.case = TRUE)){
#         model_coarse_grain <- c(model_coarse_grain, "Tmn")
#       }
#       
#       # print(model_coarse_grain)
#       
#       # model_coarse_grain <- c("U", "R")
#       
#     }
#     
#     # check that all the coarse grains are identical
#     unique_coarse_grain <- unique(model_coarse_grain)
#     # if they are not all identical, assign the peptide to 'U'
#     if(length(unique_coarse_grain) > 1){
#       coarse_assignment <- "U"
#     }
#     
#     if(length(unique_coarse_grain) == 1){
#       coarse_assignment <- unique_coarse_grain
#     }
#     
#     if(length(unique_coarse_grain) == 0){
#       coarse_assignment <- "U"
#     }
#     
#   return(coarse_assignment)
# }
# 
# get_coarse_grains_specific <- function(tax_group, tax_assign_df, 
#                                        filter_type, all_taxa = FALSE){
#   
#   tax_assign_in_df <- 'tax_assign' %in% colnames(tax_assign_df)
#   expect_true(tax_assign_in_df)
#   
#   # tax_group <- 'Diatom'
#   # all_taxa <- TRUE
#   # tax_assign_df <- tax_data
#   # filter_type <- 3.0
#   # 
#   if(!all_taxa){
#     org_peps <- tax_assign_df %>% 
#       dplyr::filter(tax_assign == tax_group, 
#                     filter == filter_type)
#   }
#   if(all_taxa){
#     org_peps <- tax_assign_df %>% 
#       dplyr::filter(filter == filter_type)
#   }
#   
#   all_coarse_assignments <- c()
#   all_associated_desc <- c()
#   
#   for(pep_i in 1:nrow(org_peps)){
#     
#     # pep_i <- 1
#     print(pep_i)
#     # print(diatom_peps[i, ]$go_descs)
#     # parse the output
#     contigs_parsed <- parse_contigs(protein_quant_out = org_peps,
#                                     row_i = pep_i)
#     
#     # for each contig, get a vector of descriptions
#     all_descs <- get_all_descriptions_for_contigs_specific(parsed_contigs_out = contigs_parsed)
#     
#     print(all_descs)
#     # match descriptions to the model output
#     coarse_assignment_val <- match_desc_to_model_grains_specific(all_descs)
#     
#     all_coarse_assignments <- c(all_coarse_assignments, coarse_assignment_val)
#     
#     # if(!all_taxa){
#     
#     # all_associated_desc <- c(all_associated_desc, paste0(all_descs, collapse = ";"))
#     # }
#     
#   }
#   
#   org_peps$coarse_grains <- all_coarse_assignments
#   
#   # if(!all_taxa){
#   # org_peps$all_coarse_descs <- all_associated_desc
#   # }
#   
#   return(org_peps)
# }

# process_results <- function(sample1, sample2, sample3, sample4, tax_level){
#   
#   sample1$day <- rep(1, nrow(sample1))
#   sample2$day <- rep(2, nrow(sample2))
#   sample3$day <- rep(3, nrow(sample3))
#   sample4$day <- rep(4, nrow(sample4))
#   
#   if(tax_level == 'tax'){
#     
#     # adding a column that has group normalized expression levels
#     
#     sample1_t <- sample1 %>% 
#       dplyr::select(cluster, db_mcl_sum, nn_mcl_sum) %>% 
#       dplyr::rename(db_mcl_sum1 = db_mcl_sum, nn_mcl_sum1 = nn_mcl_sum) %>% 
#       group_by(tax_assign, day) %>% 
#       summarise(total_prot_gr = sum(db_mcl_sum)) %>%
#       left_join(sample1, by = c("tax_assign", "day")) %>% 
#       mutate(db_mcl_sum_tax_norm = db_mcl_sum/total_prot_gr) 
#     
#     sample2_t <- sample2 %>% 
#       dplyr::select(cluster, db_mcl_sum, nn_mcl_sum) %>% 
#       dplyr::rename(db_mcl_sum2 = db_mcl_sum, nn_mcl_sum2 = nn_mcl_sum) %>% 
#       group_by(tax_assign, day) %>% 
#       summarise(total_prot_gr = sum(db_mcl_sum)) %>%
#       left_join(sample2, by = c("tax_assign", "day")) %>% 
#       mutate(db_mcl_sum_tax_norm = db_mcl_sum/total_prot_gr)
#     
#     sample3_t <- sample3 %>% 
#       dplyr::select(cluster, db_mcl_sum, nn_mcl_sum) %>% 
#       dplyr::rename(db_mcl_sum3 = db_mcl_sum, nn_mcl_sum3 = nn_mcl_sum) %>% 
#       group_by(tax_assign, day) %>% 
#       summarise(total_prot_gr = sum(db_mcl_sum)) %>%
#       left_join(sample3, by = c("tax_assign", "day")) %>% 
#       mutate(db_mcl_sum_tax_norm = db_mcl_sum/total_prot_gr) 
#     
#     sample4_t <- sample4 %>% 
#       dplyr::select(cluster, db_mcl_sum, nn_mcl_sum) %>% 
#       dplyr::rename(db_mcl_sum4 = db_mcl_sum, nn_mcl_sum4 = nn_mcl_sum) %>% 
#       group_by(tax_assign, day) %>% 
#       summarise(total_prot_gr = sum(db_mcl_sum)) %>%
#       left_join(sample4, by = c("tax_assign", "day")) %>% 
#       mutate(db_mcl_sum_tax_norm = db_mcl_sum/total_prot_gr) 
#   }
#   if(tax_level == 'grpnorm'){
#     # adding a column that has group normalized expression levels
#     
#     sample1_t <- sample1 %>% 
#       dplyr::select(cluster, db_mcl_sum, nn_mcl_sum) %>% 
#       dplyr::rename(db_mcl_sum1 = db_mcl_sum, nn_mcl_sum1 = nn_mcl_sum) %>% 
#       group_by(grpnorm_assign, day) %>% 
#       summarise(total_prot_gr = sum(db_mcl_sum)) %>%
#       left_join(sample1, by = c("grpnorm_assign", "day")) %>% 
#       mutate(db_mcl_sum_grp_norm = db_mcl_sum/total_prot_gr) 
#     
#     sample2_t <- sample2 %>% 
#       dplyr::select(cluster, db_mcl_sum, nn_mcl_sum) %>% 
#       dplyr::rename(db_mcl_sum2 = db_mcl_sum, nn_mcl_sum2 = nn_mcl_sum) %>% 
#       group_by(grpnorm_assign, day) %>% 
#       summarise(total_prot_gr = sum(db_mcl_sum)) %>%
#       left_join(sampl2, by = c("grpnorm_assign", "day")) %>% 
#       mutate(db_mcl_sum_grp_norm = db_mcl_sum/total_prot_gr) 
#     
#     sample3_t <- sample3 %>% 
#       dplyr::select(cluster, db_mcl_sum, nn_mcl_sum) %>% 
#       dplyr::rename(db_mcl_sum3 = db_mcl_sum, nn_mcl_sum3 = nn_mcl_sum) %>% 
#       group_by(grpnorm_assign, day) %>% 
#       summarise(total_prot_gr = sum(db_mcl_sum)) %>%
#       left_join(sample3, by = c("grpnorm_assign", "day")) %>% 
#       mutate(db_mcl_sum_grp_norm = db_mcl_sum/total_prot_gr) 
#     
#     sample4_t <- sample4 %>% 
#       dplyr::select(cluster, db_mcl_sum, nn_mcl_sum) %>% 
#       dplyr::rename(db_mcl_sum4 = db_mcl_sum, nn_mcl_sum4 = nn_mcl_sum) %>% 
#       group_by(grpnorm_assign, day) %>% 
#       summarise(total_prot_gr = sum(db_mcl_sum)) %>%
#       left_join(sample4, by = c("grpnorm_assign", "day")) %>% 
#       mutate(db_mcl_sum_grp_norm = db_mcl_sum/total_prot_gr) 
#   }
#   if(tax_level == 'func'){
#     sample1_t <- sample1 %>% 
#       dplyr::select(cluster, db_mcl_sum, nn_mcl_sum) %>% 
#       dplyr::rename(db_mcl_sum1 = db_mcl_sum, nn_mcl_sum1 = nn_mcl_sum)
#     
#     sample2_t <- sample2 %>% 
#       dplyr::select(cluster, db_mcl_sum, nn_mcl_sum) %>% 
#       dplyr::rename(db_mcl_sum2 = db_mcl_sum, nn_mcl_sum2 = nn_mcl_sum)
#     
#     sample3_t <- sample3 %>% 
#       dplyr::select(cluster, db_mcl_sum, nn_mcl_sum) %>% 
#       dplyr::rename(db_mcl_sum3 = db_mcl_sum, nn_mcl_sum3 = nn_mcl_sum)
#     
#     sample4_t <- sample4 %>% 
#       dplyr::select(cluster, db_mcl_sum, nn_mcl_sum) %>% 
#       dplyr::rename(db_mcl_sum4 = db_mcl_sum, nn_mcl_sum4 = nn_mcl_sum)
#   }
# 
#   
#   s0102 <- full_join(sample1_t, sample2_t)
#   s0304 <- full_join(sample3_t, sample4_t)
#   finale <- full_join(s0102, s0304)
#   return(finale)
# }
# 
# process_results2 <- function(sample1, sample2, sample3, sample4, tax_level){
#   
#   sample1$day <- rep(1, nrow(sample1))
#   sample2$day <- rep(2, nrow(sample2))
#   sample3$day <- rep(3, nrow(sample3))
#   sample4$day <- rep(4, nrow(sample4))
#   
#   finale <- rbind(sample1, sample2, sample3, sample4)
#   
#   if(tax_level == 'tax'){
#     
#     # adding a column that has group normalized expression levels
#     mcl_good2 <- finale %>%
#       group_by(tax_assign, day) %>% 
#       summarise(total_prot_gr = sum(db_mcl_sum)) %>%
#       left_join(finale, by = c("tax_assign", "day")) %>% 
#       mutate(db_mcl_sum_tax_norm = db_mcl_sum/total_prot_gr) 
#   }
#   if(tax_level == 'grpnorm'){
#     # adding a column that has group normalized expression levels
#     mcl_good2 <- finale %>%
#       group_by(grpnorm_assign, day) %>% 
#       summarise(total_prot_gr = sum(db_mcl_sum)) %>%
#       left_join(finale, by = c("grpnorm_assign", "day")) %>% 
#       mutate(db_mcl_sum_grp_norm = db_mcl_sum/total_prot_gr)
#   } 
#   if(tax_level == 'func'){
#     mcl_good2 <- finale
#   }
# 
#   return(mcl_good2)
# }
# function for assigning annotations and summing peptides within mcl clusters
# annot_connect_good <- function(normalize_output, 
#                                annot_df = mcl_annot){
#   print('Assigning annotations to peptides...')
#   
#   mcl_in_df <- "mcl_assignments" %in% colnames(normalize_output)
#   expect_true(mcl_in_df)
#   
#   mcl_assignments_in_df <- "mcl_assignments" %in% colnames(normalize_output)
#   expect_true(mcl_assignments_in_df)
#   
#   mcl_good <- normalize_output %>% filter(mcl_assignments != 'multi-mcl', 
#                                           mcl_assignments != 'no-mcl') %>% 
#     group_by(mcl_assignments) %>%
#     summarise(db_mcl_sum = sum(db_abund_avg), 
#               nn_mcl_sum = sum(abund_no_norm_avg),
#               number_peps = n(),
#               db_mcl_mean = mean(db_abund_avg),
#               nn_mcl_mean = mean(abund_no_norm_avg)) %>%
#     dplyr::rename(cluster = mcl_assignments) %>%
#     inner_join(y = annot_df, by = "cluster")
# 
#   return(mcl_good)
# }
# 

# annot_connect_fuzzy <- function(normalize_output){
#   
#   mcl_fuzzy <- normalize_output %>% 
#     filter(mcl_assignments == 'no-mcl', n_proteins == 1) %>% 
#     dplyr::rename(orf_id = protein) %>% 
#     inner_join(y = mcl_df, by = 'orf_id')
#   
#   return(mcl_fuzzy)
# }




########### this chunk used to be at the end of pipeline
# 
# # aggregated the list of dataframes
# process_func <- process_results(sample1 = df_list_func[[1]],
#                                 sample2 = df_list_func[[2]],
#                                 sample3 = df_list_func[[3]],
#                                 sample4 = df_list_func[[4]], 
#                                 tax_level = 'func')
# process_tax <- process_results(sample1 = df_list_tax[[1]],
#                                sample2 = df_list_tax[[2]],
#                                sample3 = df_list_tax[[3]],
#                                sample4 = df_list_tax[[4]], 
#                                tax_level = 'tax')
# process_grpnorm <- process_results(sample1 = df_list_grpnorm[[1]],
#                                    sample2 = df_list_grpnorm[[2]],
#                                    sample3 = df_list_grpnorm[[3]],
#                                    sample4 = df_list_grpnorm[[4]], 
#                                    tax_level = 'grpnorm')
# 
# process_func2 <- process_results2(sample1 = df_list_func[[1]],
#                                   sample2 = df_list_func[[2]],
#                                   sample3 = df_list_func[[3]],
#                                   sample4 = df_list_func[[4]], 
#                                   tax_level = 'func')
# process_tax2 <- process_results2(sample1 = df_list_tax[[1]],
#                                  sample2 = df_list_tax[[2]],
#                                  sample3 = df_list_tax[[3]],
#                                  sample4 = df_list_tax[[4]],
#                                  tax_level = 'tax')
# process_grpnorm2 <- process_results2(sample1 = df_list_grpnorm[[1]],
#                                      sample2 = df_list_grpnorm[[2]],
#                                      sample3 = df_list_grpnorm[[3]],
#                                      sample4 = df_list_grpnorm[[4]],
#                                      tax_level = 'grpnorm')
# 
# # write the final results! (in the wide format)
# write.csv(process_func, file = paste0('filter_', sample_grouping, 
#                                       "_", 
#                                       id_string, "_func_wide.csv"))
# write.csv(process_tax, file = paste0('filter_', sample_grouping, 
#                                      "_", 
#                                      id_string, "_tax_wide.csv"))
# write.csv(process_grpnorm, file = paste0('filter_', sample_grouping, 
#                                          "_", 
#                                          id_string, "_grpnorm_wide.csv"))
# 
# # write the final results!
# write.csv(process_func2, file = paste0('filter_', sample_grouping, 
#                                        "_", 
#                                        id_string, "_func.csv"))
# write.csv(process_tax2, file = paste0('filter_', sample_grouping, 
#                                       "_", 
#                                       id_string, "_tax.csv"))
# write.csv(process_grpnorm2, file = paste0('filter_', sample_grouping, 
#                                           "_", 
#                                           id_string, "_grpnorm.csv"))
