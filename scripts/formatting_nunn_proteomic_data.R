# script for formatting Nunn 2013 data for diatom proteomics

library(dplyr)
library(ggplot2)
library(magrittr)
library(reshape2)
library(stringr)

# setwd("../../Google Drive/Projects/mn-fe-allocation/")

# nunn_model <- read.csv("data/model_output/nunn2013_paper_comparison.csv")

# nunn_model_out <- read_in_data_loop(data_name_string = 'csv',
                                    # dir_to_look = "data/model_output/nunn_model_out/")
# nunn_model_mean <- nunn_model_out %>% group_by(Mnx, Fex) %>% summarize_all(mean)

nunn_data <- read.csv("../data/culture-proteomics-data/Nunn2013_Table_S1_annotated_formatted_atpsynth_out.csv")

# getting the sum of spectral counts per coarse grained protein group for + Fe
plus_fe <- nunn_data %>%
  filter(Primary.vs.Secondary.IDs.homologous.protein.sequences. == 'primary') %>% 
  rowwise() %>% 
  mutate(tp_fe = mean(c(tp_fe1, tp_fe2, tp_fe3, tp_fe4), na.rm = TRUE),
         tp_fe_sd = sd(c(tp_fe1, tp_fe2, tp_fe3, tp_fe4), na.rm = TRUE)) %>%
  group_by(coares_grained_category) %>% 
  summarize(coarse_sum = sum(tp_fe),
            coarse_sum_sd = mean(tp_fe_sd))

# getting the number of spectral counts per coarse grained protein group
plus_fe_reps <- nunn_data %>%
  filter(Primary.vs.Secondary.IDs.homologous.protein.sequences. == 'primary') %>% 
  dplyr::select(tp_fe1, tp_fe2, tp_fe3, tp_fe4, coares_grained_category) %>% 
  melt(id = 'coares_grained_category', 
       value.name = 'pep_count',
       variable.name = 'replicate')

plus_fe_reps$pep_count <- plus_fe_reps$pep_count %>% as.numeric()

plus_fe_reps2 <- plus_fe_reps %>% 
  group_by(coares_grained_category, replicate) %>% 
  summarize(coarse_sum = sum(pep_count, na.rm = TRUE)) 

# getting the sum of spectral counts per coarse grained protein group for - Fe
minus_fe <- nunn_data %>%
  filter(Primary.vs.Secondary.IDs.homologous.protein.sequences. == 'primary') %>% 
  rowwise() %>% 
  mutate(tp_fe = mean(c(tp_nofe1, tp_nofe2, tp_nofe3, tp_nofe4), na.rm = TRUE),
         tp_fe_sd = sd(c(tp_nofe1, tp_nofe2, tp_nofe3, tp_nofe4), na.rm = TRUE)) %>%
  group_by(coares_grained_category) %>% 
  summarize(coarse_sum = sum(tp_fe),
            coarse_sum_sd = mean(tp_fe_sd))

# getting the number of spectral counts per coarse grained protein group
minus_fe_reps <- nunn_data %>%
  dplyr::filter(Primary.vs.Secondary.IDs.homologous.protein.sequences. == 'primary') %>% 
  dplyr::select(tp_nofe1, tp_nofe2, tp_nofe3, tp_nofe4, coares_grained_category) %>% 
  melt(id = 'coares_grained_category', 
       value.name = 'pep_count',
       variable.name = 'replicate')# %>% 

minus_fe_reps$pep_count <- minus_fe_reps$pep_count %>% as.numeric()

minus_fe_reps2 <- minus_fe_reps %>% 
  group_by(coares_grained_category, replicate) %>% 
  summarize(coarse_sum = sum(pep_count, na.rm = TRUE)) 

# adding the Fe prime concentration
plus_fe_reps2$Fex <- rep(2876, nrow(plus_fe_reps2))

# determining the percentage of spectral counts in each replicate, and making the coarse grained values
# a corresponding percentage
plus_fe_reps2 <- plus_fe_reps2 %>% 
  group_by(replicate) %>% 
  mutate(coarse_sum_per = coarse_sum/sum(coarse_sum))

# adding the Fe prime concentration
minus_fe_reps2$Fex <- rep(87, nrow(minus_fe_reps2))

# determining the percentage of spectral counts in each replicate, and making the coarse grained values
# a corresponding percentage
minus_fe_reps2 <- minus_fe_reps2 %>% 
  group_by(replicate) %>% 
  mutate(coarse_sum_per = coarse_sum/sum(coarse_sum))

# combining the plus and minus Fe treatments
nunn_output <- rbind(minus_fe_reps2, plus_fe_reps2)

names(nunn_output) <- c('coarse_grain', 'replicate', 'sum_spec_counts', 'Fex', 'percentage')

nunn_output$model <- c(rep('Proteomic Data', nrow(nunn_output)))


# loading and formatting phaeo data ---------------------------------------


phaeo_data_with_na <- read.csv(file = '../data/culture-proteomics-data/mn_fe_phaeocystis_natcomms_all_proteins_atpsynth_out.csv')

abundance_value_columns <- phaeo_data_with_na %>% select(contains('Fe'))

phaeo_data <- phaeo_data_with_na[complete.cases(abundance_value_columns), ]

# transforming to long data format
ph_transform <- phaeo_data %>% 
  select(coarse_grain, contains('Fe')) %>% 
  melt(value.name = 'abundance', 
       variable.name = 'treatment')

# getting the replicate number
ph_transform$replicate <- str_sub(ph_transform$treatment %>% as.character(), start = -1)

# getting the total protein per replicate as a normalizing constant
total_prot_per_replicate <- ph_transform %>% 
  group_by(treatment, replicate) %>% 
  summarize(total_prot_abundance = sum(abundance, na.rm = TRUE))

# getting the P and R as mass fraction
ph_culture_mf <- ph_transform %>% 
  inner_join(total_prot_per_replicate, by = c('treatment', 'replicate')) %>% 
  mutate(protein_mass_fraction = abundance/total_prot_abundance) %>% 
  group_by(treatment, replicate, coarse_grain) %>% 
  summarize(coarse_mass_fraction = sum(protein_mass_fraction, na.rm = TRUE)) %>% 
  dcast( replicate + treatment ~ coarse_grain)

# ph_culture_mf %>% ggplot(aes(x = treatment, y = R)) +
#   geom_point()
# ph_culture_mf %>% ggplot(aes(x = treatment, y = P)) +
#   geom_point()

# plotting culture proteomic data together --------------------------------
nunn_output_for_plotting <- nunn_output %>% 
  filter(coarse_grain == 'R' | coarse_grain == 'P') %>% 
  select(-model) %>% 
  dcast(formula = Fex + replicate ~ coarse_grain)




