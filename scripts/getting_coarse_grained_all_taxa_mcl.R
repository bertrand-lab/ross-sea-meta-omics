# getting coarse grained proteomes from the metaproteome

library(dplyr)
library(readxl)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(magrittr)
library(forcats)

# reading in culture proteomics data
source('formatting_nunn_proteomic_data.R')

# taxonomically mapped peptides
tax_data <- read.csv("../data/go-tax-pipeline-output/pipeline_outtfg-all_tax_peptides.csv")

# dataframe with the mcl-to-ORF id mappings
mcl_df <- read_excel("../data/mcmurdo-metatrans/Bertrand_McCrow_TFG/annotation_allTFG.grpnorm_mmetsp_fc_pn_reclassified.edgeR.xlsx", 
                     sheet = 1)

# data frame with the mcl-cluster consensus annotations
mcl_annot <- read_excel("../data/mcmurdo-metatrans/Bertrand_McCrow_TFG/summaries_combined/annotation_allTFG.mmetsp_fc_pn_reclassified.summary_cluster.edgeR.xlsx")

## merge the mcl_cluster ids to the tax data
mcl_annot_annot_subset <- mcl_annot %>% dplyr::select(cluster, 
                                                      cluster_size, 
                                                      ann_type,
                                                      ann_id,
                                                      ann_desc) %>% 
  dplyr::rename(mcl_assignments = cluster)

# make a row that is an empty annotation description. without this row, 
mcl_annot_annot_subset_extra <- rbind(mcl_annot_annot_subset,
                                      data.frame(mcl_assignments = c('no-mcl', 'multi-mcl'),
                                                 cluster_size = c(NA, 
                                                                  NA),
                                                 ann_type = c(NA, 
                                                              NA),
                                                 ann_id = c(NA, 
                                                            NA),
                                                 ann_desc = c(NA, NA)))

tax_data_mcl <- tax_data %>% 
  inner_join(mcl_annot_annot_subset_extra,
             by = c('mcl_assignments'))

# making the NA's for Ann Desc actually an identified cluster with a number
tax_data_mcl[is.na(tax_data_mcl$ann_desc) & 
               tax_data_mcl$mcl_numbers == 1, ]$ann_desc <- tax_data_mcl[is.na(tax_data_mcl$ann_desc) & 
                                        tax_data_mcl$mcl_numbers == 1, ]$mcl_assignments

all_clusters_in_tax <- tax_data_mcl[grepl(x = tax_data_mcl$ann_desc, 
                                          pattern = 'clust_'), ]$ann_desc 

tax_data_mcl[grepl(x = tax_data_mcl$ann_desc, 
                     pattern = 'clust_'), ]$ann_desc <- sub(pattern = 'clust_', 
                                                            replacement = 'Unknown Protein Cluster ', 
                                                            x = all_clusters_in_tax)

# getting all the annotation descriptions for the mcl clusters
ann_desc_manual <- unique(tax_data_mcl['ann_desc'])

write.csv(ann_desc_manual, file = '../data/go-tax-pipeline-output/mcl_desc.csv')
# this file was then manually annotated

# reading in the manually annotated file
ann_desc_coarse <- read.csv('../data/go-tax-pipeline-output/mcl_desc_manual.csv')

# calculate the proportion of the proteome
summed_tax_abundance_per_day <- tax_data_mcl %>% 
  group_by(day, tax_assign, filter) %>% 
  summarize(sum_pep_tax_abundance = sum(db_abund_avg))

coarse_per_day_summarized <- tax_data_mcl %>% 
  left_join(ann_desc_coarse, by = c('ann_desc')) %>%
  inner_join(summed_tax_abundance_per_day, by = c('day', 'tax_assign', 'filter')) %>% 
  mutate(peptide_prot_percentage = db_abund_avg/sum_pep_tax_abundance) %>% 
  group_by(tax_assign, day, filter, coarse_grain) %>% 
  summarize(sum_coarse_prot = sum(peptide_prot_percentage))
    
# checking that the proteomes sum to one for each taxa on each day
coarse_per_day_summarized %>% 
  group_by(day, filter, tax_assign) %>% 
  summarize(sum_test = sum(sum_coarse_prot)) %>% 
  ggplot(aes(x = sum_test)) + 
  geom_histogram()
# 1    0.1 Chlorophyta                  0.407

# getting the list of taxa that should be in each diagram (i.e. the ones that comprimise 75% of total protein biomass)
abundance_rankings <- tax_data_mcl %>% 
  group_by(filter) %>% 
  summarize(total_all_tax = sum(db_abund_avg))

abundance_rankings_temp <- tax_data_mcl %>% 
  inner_join(abundance_rankings, by = c('filter')) %>% 
  group_by(filter, tax_assign) %>% 
  mutate(fraction_biomass = db_abund_avg/total_all_tax) %>% 
  summarize(taxa_sum_fraction = sum(fraction_biomass))
  
selected_taxa <- abundance_rankings_temp %>% 
  filter(taxa_sum_fraction > 0.025)
selected_taxa_unique <- selected_taxa$tax_assign %>% unique() %>% as.character()
selected_taxa_unique_noprom <- selected_taxa_unique[-c(4, 6)]

selected_taxa_big_no_change <- selected_taxa_unique_noprom[c(1, 2, 4, 8)]
selected_taxa_small <- selected_taxa_unique_noprom[c(3, 5, 6, 7, 9)]
selected_taxa_small <- gsub('Other ', '', x = selected_taxa_small)

selected_taxa_big <- c('Ciliates', 'Dinoflagellates', 'Haptophytes', 'Diatoms')

# tax_data_mcl %>% 
#   filter(tax_assign %in% selected_taxa_unique_noprom) %>% 
#   group_by(tax_assign, filter, day) %>% 
#   summarize(tax_abund = sum(db_abund_avg)) %>% 
#   ggplot(aes(x = fct_reorder(tax_assign, 
#                              tax_abund), 
#              y = tax_abund)) +
#   geom_col() +
#   coord_flip() +
#   facet_grid(day~filter)

# looking at the coarse grained proteome
coarse_per_day_summarized %>% 
  filter(coarse_grain != 'U',
         tax_assign %in% selected_taxa_big_no_change) %>% 
  ggplot(aes(x = coarse_grain, y = sum_coarse_prot)) +
  geom_point() +
  facet_grid(day~tax_assign)

# getting coarse grained estimates using greedy approach ------------------

filter_protein_total_amounts <- read.csv('../data/filter_protein_amounts_page77_smccain_labbook.csv')

# number of protein pellets combined with 8M urea
proteins_pellets <- data.frame(filter = c(0.1, 0.8, 3.0),
                               number_pellets = c(1/3, 1/2, 1))

greedy01 <- read.csv('../data/correction-factor-output/coarse_all_correction_factors_01greedy.csv')
greedy08 <- read.csv('../data/correction-factor-output/coarse_all_correction_factors_08greedy.csv')
greedy30 <- read.csv('../data/correction-factor-output/coarse_all_correction_factors_30greedy.csv')

all_greedy <- rbind(greedy01, greedy08, greedy30)

levels(all_greedy$tax_assign)[5] <- 'Ciliates'
levels(all_greedy$tax_assign)[8] <- 'Diatoms'
levels(all_greedy$tax_assign)[9] <- 'Dinoflagellates'
levels(all_greedy$tax_assign)[19] <- 'Haptophytes'
levels(all_greedy$tax_assign)[18] <- 'Gammaproteobacteria'

all_greedy[grepl(pattern = 'ATP', 
                 x = all_greedy$tigrfams_desc) |
             grepl(pattern = 'ATP', 
                   x = all_greedy$kegg_desc) |
             grepl(pattern = 'ATP', 
                   x = all_greedy$ko_desc) |
             grepl(pattern = 'ATP', 
                   x = all_greedy$kogg_desc) |
             grepl(pattern = 'ATP', 
                   x = all_greedy$pfams_desc), ]$coarse_grains <- "U"  
  

# all_greedy <- all_greedy %>% dplyr::filter(!grepl('ATP', kogg_desc, ignore.case = TRUE))
# getting the total peptide abundance per day per filter for each taxa
summed_tax_abundance_per_day_greedy <- all_greedy %>% 
  group_by(day, tax_assign, filter) %>% 
  summarize(sum_pep_tax_abundance = sum(db_abund_avg))

# getting the total peptide abundance per day for filter for each taxa, but using grpnorm taxa
summed_tax_abundance_per_day_greedy_grpnorm <- all_greedy %>% 
  group_by(day, grpnorm_assign, filter) %>% 
  summarize(sum_pep_tax_abundance = sum(db_abund_avg))

# getting the coarse grained abundance by first transforming the peptide abundance
# to a mass fraction, and then summing those mass fractions.
coarse_per_day_summarized_greedy <- all_greedy %>% 
  inner_join(summed_tax_abundance_per_day_greedy, by = c('day', 'tax_assign', 'filter')) %>% 
  mutate(peptide_prot_percentage = db_abund_avg/sum_pep_tax_abundance) %>% 
  group_by(tax_assign, day, filter, coarse_grains) %>% 
  summarize(sum_coarse_prot = sum(peptide_prot_percentage),
            number_peps = n()) 

# calculating the number of peptides observed for a given taxa on a given day
# to calculate weights
sum_of_peptides_greedy <- coarse_per_day_summarized_greedy %>% 
  group_by(tax_assign, day, coarse_grains) %>% 
  summarize(all_peps_per_day = sum(number_peps))

# getting weighting factors for total protein per filter on each d --------
filter_protein_total_amounts_sum_per_day <- filter_protein_total_amounts %>% 
  inner_join(proteins_pellets, by = c('filter')) %>% 
  mutate(amount_aggregated_pellet = amount*number_pellets) %>% 
  group_by(day) %>% 
  summarize(sum_prot_amount = sum(amount_aggregated_pellet))

filter_protein_total_amounts_avg <- filter_protein_total_amounts %>% 
  inner_join(proteins_pellets, by = c('filter')) %>% 
  mutate(amount_aggregated_pellet = amount*number_pellets) %>% 
  inner_join(filter_protein_total_amounts_sum_per_day, by = c('day')) %>% 
  mutate(fraction_per_day = amount_aggregated_pellet/sum_prot_amount) 

# first calculating the weighted mass fraction per filter, then
# summing up the weighted values.
coarse_per_day_summarized_weighted <- coarse_per_day_summarized_greedy %>% 
  inner_join(sum_of_peptides_greedy, 
             by = c('tax_assign', 'day', 'coarse_grains')) %>% 
  # inner_join(filter_protein_total_amounts_avg, by = 'filter') %>% 
  mutate(weighted_sum_coarse_prot = sum_coarse_prot*(number_peps/all_peps_per_day)) %>% 
  group_by(tax_assign, day, coarse_grains) %>% 
  summarize(weighted_summed_coarse_prot_allfilter = sum(weighted_sum_coarse_prot))


# sub grouping taxonomic analysis -----------------------------------------

# getting the total peptide abundance per day per filter for each taxa
summed_grp_abundance_per_day_greedy <- all_greedy %>% 
  group_by(day, grpnorm_assign, filter) %>% 
  summarize(sum_pep_grp_abundance = sum(db_abund_avg))

# getting the coarse grained abundance by first transforming the peptide abundance
# to a mass fraction, and then summing those mass fractions.
coarse_per_day_summarized_greedy_grpnorm <- all_greedy %>% 
  inner_join(summed_grp_abundance_per_day_greedy, by = c('day', 'grpnorm_assign', 'filter')) %>% 
  mutate(peptide_prot_percentage = db_abund_avg/sum_pep_grp_abundance) %>% 
  group_by(grpnorm_assign, day, filter, coarse_grains) %>% 
  summarize(sum_coarse_prot = sum(peptide_prot_percentage),
            number_peps = n()) 

# calculating the number of peptides observed for a given taxa on a given day
# to calculate weights
sum_of_peptides_greedy_grpnorm <- coarse_per_day_summarized_greedy_grpnorm %>% 
  group_by(grpnorm_assign, day, coarse_grains) %>% 
  summarize(all_peps_per_day = sum(number_peps))

# applying the weights
coarse_per_day_summarized_weighted_gprnorm <- coarse_per_day_summarized_greedy_grpnorm %>% 
  inner_join(sum_of_peptides_greedy_grpnorm, 
             by = c('grpnorm_assign', 'day', 'coarse_grains')) %>% 
  # inner_join(filter_protein_total_amounts_avg, by = 'filter') %>% 
  mutate(weighted_sum_coarse_prot = sum_coarse_prot*(number_peps/all_peps_per_day)) %>% 
  group_by(grpnorm_assign, day, coarse_grains) %>% 
  summarize(weighted_summed_coarse_prot_allfilter = sum(weighted_sum_coarse_prot))

levels(coarse_per_day_summarized_weighted_gprnorm$coarse_grains) <- c('A', 'Photosynthesis\nProteins', 'Ribosomal\nProteins', 'Tfe', 'Tn', 'U')

grpnorm_subgrouping_analysis <- coarse_per_day_summarized_weighted_gprnorm %>% 
  filter(coarse_grains == 'Ribosomal\nProteins' | coarse_grains == 'Photosynthesis\nProteins',
         grpnorm_assign == 'Fragilariopsis' | grpnorm_assign == 'Pseudo-nitzschia') %>% 
  ggplot(aes(x = coarse_grains, y = weighted_summed_coarse_prot_allfilter)) +
  geom_point(size = 4, alpha = 0.6) +
  facet_grid(~grpnorm_assign) +
  theme_bw() +
  ylab('Proportion of Proteome') +
  xlab('Proteomic Pool');grpnorm_subgrouping_analysis

where_is_hapto <- grepl(pattern = 'Hapto', 
      x = levels(coarse_per_day_summarized_weighted_gprnorm$grpnorm_assign))
levels(coarse_per_day_summarized_weighted_gprnorm$grpnorm_assign)[where_is_hapto] <- 'Haptophyta'

grpnorm_subgrouping_and_database_analysis <- coarse_per_day_summarized_weighted_gprnorm %>% 
  filter(coarse_grains == 'Ribosomal\nProteins' | coarse_grains == 'Photosynthesis\nProteins',
         grpnorm_assign == 'Fragilariopsis' | 
           grpnorm_assign == 'Pseudo-nitzschia' |
           grpnorm_assign == 'Dinophyta' | 
           grpnorm_assign == 'Haptophyta' |
           grpnorm_assign == 'Ciliophora') %>% 
  ggplot(aes(x = coarse_grains, y = weighted_summed_coarse_prot_allfilter)) +
  geom_point(size = 4, alpha = 0.6) +
  facet_grid(~grpnorm_assign) +
  theme_bw() +
  ylab('Proportion of Proteome') +
  xlab('Proteomic Pool')

ggsave(grpnorm_subgrouping_analysis, 
       filename = '../figures/subgrouping_frag_pseudo.png',
       width = 7.28, height = 4.45)
ggsave(grpnorm_subgrouping_and_database_analysis, 
       filename = '../figures/subgrouping_and_database_analysis.png',
       width = 9.28, height = 4.45)


# plotting coarse grains -------------------------------------------

ribosomes_psus_big_with_all_filters <- coarse_per_day_summarized_weighted %>%
  filter(coarse_grains == 'R' | coarse_grains == 'P', 
         tax_assign %in% selected_taxa_big[selected_taxa_big != 'Ciliates']) %>% 
  ggplot(aes(x = coarse_grains, 
             y = weighted_summed_coarse_prot_allfilter)) +
  geom_point(size = 4, shape = 'square', colour = 'darkblue', alpha = 0.8) +
  facet_grid(day~tax_assign) +
  theme_bw() +
  ylab('Proportion of Proteome') +
  xlab('') +
  geom_point(data = coarse_per_day_summarized_greedy %>%
               filter(coarse_grains == 'R' | coarse_grains == 'P', 
                      tax_assign %in% selected_taxa_big[selected_taxa_big != 'Ciliates']),
             aes(x = coarse_grains, y = sum_coarse_prot,
                 size = number_peps), 
             colour = 'grey30', 
             alpha = 0.2) +
  scale_size_continuous('Number of Peptides');ribosomes_psus_big_with_all_filters

levels(coarse_per_day_summarized_weighted$coarse_grains) <- c('A', 'Photosynthesis\nProteins', 'Ribosomal\nProteins', 'Tfe', 'Tn', 'U')

ribosomes_psus_big <- coarse_per_day_summarized_weighted %>% 
  filter(coarse_grains == 'Ribosomal\nProteins' | coarse_grains == 'Photosynthesis\nProteins', 
         tax_assign %in% selected_taxa_big[selected_taxa_big != 'Ciliates']) %>% 
  ggplot(aes(x = coarse_grains, 
             y = weighted_summed_coarse_prot_allfilter)) +
  # geom_line(aes(group = day, 
                # colour = day), 
            # alpha = 0.9) +
  geom_point(size = 3, aes(colour = day, shape = coarse_grains)) +
  facet_grid(~tax_assign) +
  theme_bw() +
  ylab('Proportion of Proteome') +
  xlab('') +
  guides(shape = FALSE) +
  scale_colour_gradient('Week');ribosomes_psus_big
  
ribosomes_small_hetero <- coarse_per_day_summarized_weighted %>% 
  filter(coarse_grains == 'Ribosomal\nProteins', 
         tax_assign %in% selected_taxa_small) %>% 
  ggplot(aes(x = tax_assign, 
             y = weighted_summed_coarse_prot_allfilter)) +
  geom_point(size = 5, aes(fill = day), 
             shape = 21,
             alpha = 0.8) +
  theme_bw() +
  ylab('Ribosomal Proteins\n(Proportion of Proteome)') +
  xlab('Taxonomic Group') +
  # ylim(0, 0.4) +
  ggtitle('Prokaryotic Taxa') +
  theme(legend.position = c(0.8, 0.8),
        axis.text.x = element_text(angle = 15, hjust = 1)) +
  scale_fill_gradient('Week');ribosomes_small_hetero


ggsave(ribosomes_psus_big_with_all_filters, 
       filename = '../figures/ribosomes_photo_big_all_peptides.png',
       width = 9.71, height = 9.61)

ggsave(ribosomes_small_hetero, 
       filename = '../figures/ribosomes_heterotrophic.png',
       width = 6.93, height = 4.82)

ggsave(ribosomes_psus_big, 
       filename = '../figures/ribosomes_photo_big.png',
       width = 12.6, height = 4.82)

mass_fraction_across_groups_plot <- ggarrange(ribosomes_psus_big,
          ribosomes_small_hetero, nrow = 2,
          common.legend = TRUE, legend = 'right') 

ggsave(mass_fraction_across_groups_plot, 
       filename = '../figures/mass_fraction_across_groups.png',
       width = 8.97, height = 7.52)


# plotting ribosomes vs photosystems --------------------------------------

ribo_photo_weighted_cor <- coarse_per_day_summarized_weighted %>%
  filter(coarse_grains == 'Ribosomal\nProteins' | coarse_grains == 'Photosynthesis\nProteins', 
         tax_assign %in% selected_taxa_big[selected_taxa_big != 'Ciliates']) %>% 
  dcast(formula = tax_assign + day ~ coarse_grains ) %>%
  dplyr::rename(Photosynthesis = 'Photosynthesis\nProteins',
                Ribosomal = 'Ribosomal\nProteins') %>% 
  ggplot(aes_string(x = 'Photosynthesis', y = 'Ribosomal')) +
  # geom_line(aes(group = tax_assign), alpha = 0.2) +
  geom_point(aes_string(shape = "tax_assign", fill = "day"), size = 5) +
  scale_shape_manual(values = c(21, 22, 23, 25),
                     name = 'Taxonomic Group') +
  scale_fill_continuous(name = 'Week') +
  theme_bw() +
  guides(fill = FALSE) +
  # ylim(0, 0.4) +
  ylab('Ribosomal Protein Mass Fraction') +
  xlab('Photosynthesis Protein Mass Fraction') +
  theme(legend.position = c(0.8, 0.8)) +
  ggtitle('Eukaryotic Taxa');ribo_photo_weighted_cor

# +
  # geom_point(data = nunn_output_for_plotting, shape = 22, 
  #            aes(x = P, y = R), size = 2, alpha = 0.8, fill = 'grey30') +
  # geom_point(data = ph_culture_mf, shape = 25, fill = 'grey30', 
  #            aes(x = P, y = R), size = 2, alpha = 0.8);ribo_photo_weighted_cor

nunn_transformed_rp <- nunn_output_for_plotting %>% mutate(r_over_p = R/P) %>% dplyr::select(r_over_p)
nunn_transformed_rp$tax_assign <- rep('Diatoms', nrow(nunn_transformed_rp))

ph_transformed_rp <- ph_culture_mf %>% mutate(r_over_p = R/P) %>% dplyr::select(r_over_p)
ph_transformed_rp$tax_assign <- rep('Haptophytes', nrow(ph_transformed_rp))

# combining culture datasets
culture_proteomes <- rbind(nunn_transformed_rp, 
                           ph_transformed_rp)
culture_proteomes$cultured <- rep('Culture', nrow(culture_proteomes))
culture_proteomes$day <- rep(NA, nrow(culture_proteomes))

# combining the culture and metaproteome estimates
r_over_p_meta <- coarse_per_day_summarized_weighted %>% 
  filter(coarse_grains == 'Ribosomal\nProteins' | coarse_grains == 'Photosynthesis\nProteins', 
         tax_assign %in% selected_taxa_big[selected_taxa_big != 'Ciliates']) %>% 
  dcast(formula = tax_assign + day ~ coarse_grains ) %>% 
  dplyr::rename(PhotosynthesisProteins = 'Photosynthesis\nProteins',
                RibosomalProteins = 'Ribosomal\nProteins') %>% 
  mutate(r_over_p = RibosomalProteins/PhotosynthesisProteins) %>% 
  dplyr::select(tax_assign, r_over_p, day)

# adding the label to the metaproteome dataset
r_over_p_meta$cultured <- rep('Metaproteome', nrow(r_over_p_meta))

meta_with_culture <- rbind(r_over_p_meta, culture_proteomes) %>% droplevels()

meta_culture_r_over_p <- meta_with_culture %>% 
  ggplot(aes(x = tax_assign, group = cultured, y = r_over_p)) +
  # geom_line(aes(group = tax_assign), alpha = 0.2) +
  geom_point(aes(shape = cultured, fill = day), size = 5, alpha = 0.9, 
             position = position_dodge(width = 0.4)) +
  scale_shape_manual(values = c(21, 22),
                     name = '') +
  scale_fill_continuous(name = 'Week') +
  guides(fill = FALSE) +
  theme_bw() +
  scale_y_log10() +
  ylab('Ribosomal Protein /\n Photosynthesis Protein') +
  xlab('Taxonomic Group') +
  ggtitle('Eukaryotic Taxa') +
  theme(legend.position = c(0.8, 0.75));meta_culture_r_over_p

# +
  # geom_point(data = nunn_output_for_plotting %>% mutate(r_over_p = R/P), shape = 22, 
  #            aes(x = 2, y = r_over_p), size = 2, alpha = 0.8, fill = 'grey30') +
  # geom_point(data = ph_culture_mf %>% mutate(r_over_p = R/P), shape = 25, fill = 'grey30', 
  #            aes(x = 4, y = r_over_p), size = 2, alpha = 0.8)

euk_subplot <- ggarrange(ribo_photo_weighted_cor,
                         meta_culture_r_over_p, align = 'hv', nrow = 2,
                         labels = c('a', 'b'))

top_plot <- ggarrange(ribo_photo_weighted_cor, ribosomes_small_hetero, align = 'hv', nrow = 1, labels = c('a', 'b'))
mass_fraction_figure_both <- ggarrange(top_plot, meta_culture_r_over_p, nrow = 2, 
                                       heights = c(2, 1), labels = c('', 'c'))


ggsave(mass_fraction_figure_both, filename = '../figures/mass_fraction_figure_both.png',
       width = 10.3, height = 8.67)

# getting growth rate independent proteomic fraction ----------------------

cv_peptides <- all_greedy %>% 
  inner_join(summed_tax_abundance_per_day_greedy, 
             by = c('tax_assign', 'day', 'filter')) %>% 
  mutate(pep_abundance_percent = db_abund_avg/sum_pep_tax_abundance) %>% 
  group_by(peptide, tax_assign, filter) %>% 
  summarize(peptide_abund_sd = sd(pep_abundance_percent),
            peptide_abund_mean = mean(pep_abundance_percent, na.rm = TRUE),
            peptide_abund_n = n()) %>% # this n is useful to check that peptides
  # are found across days.
  mutate(peptide_abund_cv = peptide_abund_sd/peptide_abund_mean) %>% 
  filter(peptide_abund_n == 4)

cv_peptides$tax_category <- ifelse(cv_peptides$tax_assign %in% selected_taxa_small, 
                                   yes = 'Heterotrophic',
                                   no = 'Phototrophic')

numbers_peps_per_filter <- all_greedy %>%
  group_by(tax_assign, day) %>% 
  summarize(number_peps_per_day = n())

cv_peptides_all_dom_taxa <- cv_peptides %>% 
  dplyr::filter(tax_assign %in% c(selected_taxa_big, selected_taxa_small)) %>% # | tax_assign %in% selected_taxa_big) %>% 
  ggplot(aes(x = peptide_abund_cv)) +
  geom_density(aes(fill = tax_assign), alpha = 0.6) +
  theme_bw() +
  theme(legend.position = 'None') +
  facet_wrap(~tax_assign, nrow = 5) +
  xlab('Peptide Abundance Coefficient of Variation') +
  ylab('Density');cv_peptides_all_dom_taxa

# cv_peptides_big <- cv_peptides %>% 
#   dplyr::filter(tax_assign %in% selected_taxa_big) %>% # | tax_assign %in% selected_taxa_big) %>%
#   ggplot(aes(x = peptide_abund_cv)) +
#   geom_density(aes(fill = tax_assign), alpha = 0.6) +
#   theme_bw() +
#   theme(legend.position = 'None') +
#   facet_wrap(~tax_assign, nrow = 4) +
#   xlab('Peptide Abundance Coefficient of Variation') +
#   ylab('Density');cv_peptides_big

# cv_peptides_small_hist <- cv_peptides %>% 
#   dplyr::filter(tax_assign %in% selected_taxa_unique_noprom) %>% # | tax_assign %in% selected_taxa_big) %>% 
#   ggplot(aes(x = peptide_abund_cv)) +
#   geom_histogram(aes(fill = tax_assign), alpha = 0.6) +
#   theme_bw() +
#   theme(legend.position = 'None') +
#   facet_wrap(~tax_assign, nrow = 5) +
#   xlab('Peptide Abundance Coefficient of Variation') +
#   ylab('Density');cv_peptides_small_hist
# 
# cv_peptides_big_hist <- cv_peptides %>% 
#   dplyr::filter(tax_assign %in% selected_taxa_big) %>% # | tax_assign %in% selected_taxa_big) %>%
#   ggplot(aes(x = peptide_abund_cv)) +
#   geom_histogram(aes(fill = tax_assign), alpha = 0.6) +
#   theme_bw() +
#   theme(legend.position = 'None') +
#   facet_wrap(~tax_assign, nrow = 4) +
#   xlab('Peptide Abundance Coefficient of Variation') +
#   ylab('Density');cv_peptides_big_hist

heterotrophic_taxa_cv_mean <- cv_peptides %>% 
  dplyr::filter(tax_assign %in% selected_taxa_small) %>% 
  ggplot(aes(x = peptide_abund_cv, y = peptide_abund_mean)) +
  geom_point(alpha = 0.4) +
  theme_bw() +
  ylab('Mean Peptide Abundance') +
  xlab('Peptide Abundance Coefficient of Variation') +
  ggtitle('Prokaryotic Taxa') +
  scale_y_log10()

phototrophic_taxa_cv_mean <- cv_peptides %>% 
  dplyr::filter(tax_assign %in% selected_taxa_big) %>% 
  ggplot(aes(x = peptide_abund_cv, y = peptide_abund_mean)) +
  geom_point(alpha = 0.4) +
  theme_bw() +
  ylab('Mean Peptide Abundance') +
  xlab('Peptide Abundance Coefficient of Variation') +
  scale_y_log10() +
  ggtitle('Eukaryotic Taxa')

cv_mean_plots <- ggarrange(heterotrophic_taxa_cv_mean,
                           phototrophic_taxa_cv_mean, 
                           nrow = 1, 
                           ncol = 2, 
                           labels = c('a', 
                                      'b'))


ggsave(cv_mean_plots, filename = '../figures/cv_mean_hetero_photo.png',
       width = 8.59, height = 6.57)

ggsave(cv_peptides_all_dom_taxa, 
       filename = '../figures/cv_peptides_all_dom_taxa.png', width = 5.16, height = 6.52)

# ggsave(cv_peptides_big, filename = '../figures/cv_peptides_phototroph.png', width = 5.16, height = 7.52)

cv_peptides$fixed_proteome <- ifelse(cv_peptides$peptide_abund_cv < 0.23, 
                                     yes = 'Fixed Proteome', 
                                     no = 'Dynamic Proteome')

mapping_df <- data.frame(tax_assign = c(selected_taxa_big, 
                           selected_taxa_small), 
                         tax_group = c(rep('Eukaryotic Taxa', 4), 
                                       rep('Prokaryotic Taxa', 5)))

# calculating the number of peptides observed for a given taxa 
# to calculate weights
sum_of_peptides_greedy_tax_only <- cv_peptides %>% 
  group_by(tax_assign, fixed_proteome) %>% 
  summarize(all_peps_per_taxa = n())

proportion_of_proteome_fixed_plot_euk <- cv_peptides %>% 
  group_by(tax_assign, fixed_proteome, filter) %>% 
  summarize(sum_of_fixed = sum(peptide_abund_mean),
            number_peps_per_filter = n()) %>% 
  inner_join(sum_of_peptides_greedy_tax_only, by = c('tax_assign', 'fixed_proteome')) %>% 
  mutate(weighted_sum_of_fixed = sum_of_fixed*(number_peps_per_filter/all_peps_per_taxa)) %>% 
  group_by(tax_assign, fixed_proteome) %>% 
  summarize(sum_of_weighted_sum_of_fixed = sum(weighted_sum_of_fixed)) %>% 
  dplyr::filter(tax_assign %in% selected_taxa_small | tax_assign %in% selected_taxa_big,
                fixed_proteome == 'Fixed Proteome') %>%
  inner_join(mapping_df, by = 'tax_assign') %>% 
  filter(tax_group == 'Eukaryotic Taxa') %>% 
  ggplot(aes(y = sum_of_weighted_sum_of_fixed, 
             x = tax_assign,
             colour = tax_assign)) +
  geom_point(size = 3, alpha = 0.6) +
  xlab('Taxonomic\nGroup') +
  ggtitle('Eukaryotic Taxa') +
  ylim(0, 0.5) +
  coord_flip() +
  ylab('Environment-Independent\nProteomic Mass Fraction') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none', title = element_text(size = 10)) +
  scale_colour_manual(values = c('black', 'firebrick', 'black', 'black'),
                      name = '');proportion_of_proteome_fixed_plot_euk
  
proportion_of_proteome_fixed_plot_prok <- cv_peptides %>% 
  group_by(tax_assign, fixed_proteome, filter) %>% 
  summarize(sum_of_fixed = sum(peptide_abund_mean),
            number_peps_per_filter = n()) %>% 
  inner_join(sum_of_peptides_greedy_tax_only, by = c('tax_assign', 'fixed_proteome')) %>% 
  mutate(weighted_sum_of_fixed = sum_of_fixed*(number_peps_per_filter/all_peps_per_taxa)) %>% 
  group_by(tax_assign, fixed_proteome) %>% 
  summarize(sum_of_weighted_sum_of_fixed = sum(weighted_sum_of_fixed)) %>% 
  dplyr::filter(tax_assign %in% selected_taxa_small | tax_assign %in% selected_taxa_big,
                fixed_proteome == 'Fixed Proteome') %>%
  inner_join(mapping_df, by = 'tax_assign') %>% 
  filter(tax_group == 'Prokaryotic Taxa') %>% 
  ggplot(aes(y = sum_of_weighted_sum_of_fixed, 
             x = tax_assign,
             colour = tax_assign)) +
  geom_point(size = 3, alpha = 0.6) +
  xlab('Taxonomic\nGroup') +
  ggtitle('Prokaryotic Taxa') +
  ylim(0, 0.5) +
  coord_flip() +
  ylab('Environment-Independent\nProteomic Mass Fraction') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none', title = element_text(size = 10)) +
  scale_colour_manual(values = c('black', 'black', 'black', 'darkblue', 'black'),
                      name = '');proportion_of_proteome_fixed_plot_prok


examples_of_cv <- cv_peptides %>% 
  dplyr::filter(tax_assign %in% c('Diatoms', 'SAR11')) %>% # | tax_assign %in% selected_taxa_big) %>% 
  ggplot(aes(x = peptide_abund_cv)) +
  geom_density(aes(fill = tax_assign), alpha = 0.6) +
  theme_bw() +
  theme(strip.text = element_text(size = 10)) +
  theme(legend.position = 'None') +
  facet_wrap(~tax_assign, nrow = 5) +
  xlab('Peptide Abundance Coefficient of Variation') +
  ylab('Density') +
  scale_fill_manual(values = c('firebrick', 'darkblue')) +
  geom_vline(xintercept = 0.23);examples_of_cv



combined_independent_figure <- ggarrange(examples_of_cv, ggarrange(proportion_of_proteome_fixed_plot_euk, 
                                    proportion_of_proteome_fixed_plot_prok, 
                                    nrow = 2, 
                                    align = 'hv',
                                    labels = c('b', 'c')), 
          nrow = 1,
          labels = c('a'))

ggsave(combined_independent_figure, filename = '../figures/combined_independent_figure.png',
       width = 10.6, height = 7.02)


# examining mcl terms associated with independent proteome -----------------

# joining the GO descriptions
comparing_categories <- cv_peptides %>% 
  inner_join(x = all_greedy %>% 
               dplyr::filter(mcl_numbers == 1), 
             by = c('peptide')) %>%
  group_by(fixed_proteome, mcl_assignments) %>% 
  summarize(number_mcl_total = n())

comparing_categories_just_dyn <-  comparing_categories %>%
  filter(fixed_proteome == 'Dynamic Proteome')
comparing_categories_just_fix <-  comparing_categories %>%
  filter(fixed_proteome == 'Fixed Proteome')


just_ann_desc_mcl <- tax_data_mcl %>% 
  dplyr::select(ann_desc, mcl_assignments)

just_ann_desc_mcl_unique <- unique(just_ann_desc_mcl)

comparison_of_categories_plot <- comparing_categories_just_dyn %>% 
  left_join(comparing_categories_just_fix, 
                                             by = 'mcl_assignments') %>%
  replace_na(replace = list(number_mcl_total.x = 0,
                            number_mcl_total.y = 0)) %>% 
  mutate(difference_in_number_mcl = number_mcl_total.x - number_mcl_total.y,
         total_number_mcl = number_mcl_total.x + number_mcl_total.y) %>% 
  dplyr::filter(difference_in_number_mcl > 300 | difference_in_number_mcl < -0.1) %>%
  inner_join(just_ann_desc_mcl_unique,
             by = c('mcl_assignments')) %>% 
  mutate(negative_or_pos = difference_in_number_mcl < 0) %>% 
  # filter(macl_assignments == 'clust_72') %>% 
  ggplot(aes(y = difference_in_number_mcl,
             x = fct_reorder(ann_desc, 
                             difference_in_number_mcl, 
                             .desc = FALSE))) +
  coord_flip() +
  # geom_col(aes(fill = negative_or_pos), 
  #          width = 0.5,
  #          alpha = 0.4) +
  geom_point(size = 2, shape = 'square', aes(colour = negative_or_pos)) +
  scale_colour_manual(values = c('honeydew3',
                                 'darkolivegreen4')) +
  theme_bw() +
  theme(legend.position = 'none') +
  xlab('MCL Cluster\nConsensus Annotation') +
  ylab('Environment-dependent Peptides - \nEnvironment-independent Peptides');comparison_of_categories_plot
  

ggsave(comparison_of_categories_plot, 
       filename = '../figures/comparison_dep_vs_ind.png',
       width = 7.02, height = 4.48)

combined_independent_with_mcl <- ggarrange(combined_independent_figure, 
                                           comparison_of_categories_plot, 
          ncol = 1, nrow = 2, labels = c('', 'd'),
          heights = c(3, 4))

ggsave(combined_independent_with_mcl, 
       filename = '../figures/comparison_dep_vs_ind_wmcl.png',
       width = 8, height = 8)


# examining mcl terms by taxa with associated with independent proteome -----------------

# joining the GO descriptions
comparing_categories_mcl_tax <- cv_peptides %>% 
  inner_join(x = all_greedy %>% 
               dplyr::filter(mcl_numbers == 1), 
             by = c('peptide', 'tax_assign')) %>%
  group_by(fixed_proteome, 
           mcl_assignments,
           tax_assign) %>%
  summarize(number_mcl_total = n())

comparing_categories_just_dyn_mcl_tax <-  comparing_categories_mcl_tax %>%
  filter(fixed_proteome == 'Dynamic Proteome')
comparing_categories_just_fix_mcl_tax <-  comparing_categories_mcl_tax %>%
  filter(fixed_proteome == 'Fixed Proteome')

just_ann_desc_mcl <- tax_data_mcl %>% 
  dplyr::select(ann_desc, mcl_assignments)

just_ann_desc_mcl_unique <- unique(just_ann_desc_mcl)

comparison_of_categories_plot_mcl_tax <- comparing_categories_just_dyn_mcl_tax %>% 
  inner_join(comparing_categories_just_fix_mcl_tax, 
             by = c('mcl_assignments', 'tax_assign')) %>% 
  mutate(difference_in_number_mcl = number_mcl_total.x - number_mcl_total.y,
         total_number_mcl = number_mcl_total.x + number_mcl_total.y,
         dif_weighted_val = difference_in_number_mcl/total_number_mcl) %>% 
  # filter(dif_weighted_val > 0.93 | dif_weighted_val < -0.1,
  #        tax_assign %in% selected_taxa_big) %>%
  filter(tax_assign %in% selected_taxa_big) %>% 
  inner_join(just_ann_desc_mcl_unique,
             by = c('mcl_assignments')) %>% 
  mutate(negative_or_pos = dif_weighted_val < 0) %>% 
  ggplot(aes(y = dif_weighted_val,
             x = fct_reorder(ann_desc, 
                             dif_weighted_val, 
                             .desc = FALSE))) +
  coord_flip() +
  # geom_col(aes(fill = negative_or_pos), 
  #          width = 0.5,
  #          alpha = 0.4) +
  geom_point(size = 2, shape = 'square', aes(colour = negative_or_pos)) +
  scale_fill_manual(values = c('darkred', 'darkblue')) +
  theme_bw() +
  theme(legend.position = 'none') +
  facet_grid(~tax_assign) +
  xlab('MCL Cluster Consensus Annotation') +
  ylab('(Environment-dependent Peptides - \nEnvironment-independent Peptides)/\n(All peptides identified)');comparison_of_categories_plot_mcl_tax


ggsave(comparison_of_categories_plot_mcl_tax, 
       filename = '../figures/comparison_dep_vs_ind_mcl_tax.png',
       width = 15.02, height = 20.48)

# coarse grained compared with single protein biomarkers ------------------

##### comparing the mass fraction of ribosomes and photosystem proteins with plastocyanin

cofrag_analysis_sum <- read.csv('../data/cobia-analysis/cobia_analysis_summarized.csv')

pep_abundance_greedy <- all_greedy %>% 
  inner_join(summed_tax_abundance_per_day_greedy, 
             by = c('tax_assign', 'day', 'filter')) %>% 
  mutate(pep_abundance_percent = db_abund_avg/sum_pep_tax_abundance)

isip_abund <- pep_abundance_greedy %>% 
  filter(grepl("ISIP", kegg_desc)) %>% 
  group_by(day, filter, tax_assign) %>% 
  summarize(isip_total = sum(pep_abundance_percent))

coarse_per_day_summarized_weighted %>% 
  filter(coarse_grains == 'Photosynthesis\nProteins') %>% 
  inner_join(isip_abund, by = c('tax_assign', 'day')) %>% 
  ggplot(aes(y = isip_total, 
             x = weighted_summed_coarse_prot_allfilter)) +
  geom_point() +
  # facet_grid(~tax_assign) +
  geom_smooth(method = 'lm')

plast_abund <- pep_abundance_greedy %>% 
  filter(grepl("plastocyanin", kegg_desc)) %>% 
  group_by(day, filter, tax_assign) %>% 
  summarize(plast_total = sum(pep_abundance_percent)) %>% 
  inner_join(cofrag_analysis_sum, by = c('tax_assign', 'day'))

plastocyanin_diatoms <- coarse_per_day_summarized_weighted %>% 
  filter(coarse_grains == 'Photosynthesis\nProteins',
         tax_assign == 'Diatoms') %>% 
  inner_join(plast_abund, by = c('tax_assign', 'day')) %>% 
  filter(filter == 3.0) %>% 
  ggplot(aes(y = plast_total, 
             x = weighted_summed_coarse_prot_allfilter,
             shape = coarse_grains)) +
  geom_point(aes(shape = factor(day), 
                 fill = mean_cofrag_score), 
             size = 4) +  theme_bw() +
  scale_fill_gradient('Cofragmentation Score', low = 'gold',
                      high = 'red') +
  scale_shape_manual('Week',values = c(21, 22, 23, 24)) +  # geom_smooth(method = 'lm', 
  #             se = FALSE, 
  #             colour = 'grey30') +
  # geom_smooth(method = 'lm', se = FALSE, colour = 'grey30') +
  theme_bw() +
  ggtitle('Discovery Proteomics: Diatoms') +
  # theme(legend.position = 'None') +
  ylab('Plastocyanin Mass Fraction') +
  xlab('Photosynthetic Protein Mass Fraction');plastocyanin_diatoms

plastocyanin_haptophytes <- coarse_per_day_summarized_weighted %>% 
  filter(coarse_grains == 'Photosynthesis\nProteins',
         tax_assign == 'Haptophytes') %>% 
  inner_join(plast_abund, by = c('tax_assign', 'day')) %>% 
  filter(filter == 3.0) %>% 
  ggplot(aes(y = plast_total, 
             x = weighted_summed_coarse_prot_allfilter,
             shape = coarse_grains)) +
  geom_point(aes(shape = factor(day), fill = mean_cofrag_score), 
             size = 4) +  theme_bw() +
  scale_fill_gradient('Cofragmentation Score', low = 'gold',
                       high = 'red') +
  scale_shape_manual('Week',values = c(21, 22, 23, 24)) +  # geom_smooth(method = 'lm', 
  #             se = FALSE, 
  #             colour = 'grey30') +
  ggtitle('Discovery Proteomics: Haptophytes') +
  # scale_colour_continuous(name = 'Week') +
  ylab('Plastocyanin Mass Fraction') +
  # theme(legend.position = 'None') +
  xlab('Photosynthetic Protein Mass Fraction');plastocyanin_haptophytes



# targeted proteomics data comparison -------------------------------------

targeted_data <- read.csv('../data/targeted-proteomics-data/formatted-target-data.csv') 

phaeo_biomass <- targeted_data %>% 
  filter(Taxa == 'Phaeocystis', Protein == 'SSU1' | Protein == 'SSU2') %>% 
  group_by(week) %>% 
  summarize(total_protein_per_week = mean(fmol_ug)) 

phaeo_coarse_photo <- coarse_per_day_summarized_weighted %>% 
  filter(coarse_grains == 'Photosynthesis\nProteins',
         tax_assign == 'Haptophytes') %>%
  dplyr::rename(week = day)

targeted_phaeo_proteins <- targeted_data %>% 
  filter(Taxa == 'Phaeocystis', Protein == 'PC2') %>% 
  inner_join(phaeo_biomass, by = 'week') %>% 
  mutate(plasto_mass_fraction = fmol_ug/total_protein_per_week) %>%
  inner_join(phaeo_coarse_photo, by = 'week') %>% 
  ggplot(aes(x = weighted_summed_coarse_prot_allfilter, 
             y = plasto_mass_fraction,
             fill = Peptide)) +
  geom_point(aes(shape = factor(week)), size = 4, alpha = 0.7) +
  scale_shape_manual('Week', values = c(21, 22, 23, 24)) +  # geom_smooth(method = 'lm', 
  # theme(legend.position = 'None') +
  scale_fill_manual(values=c(GGPHNVVFVEDAIPK='white',
                             GDSITWINNK='grey10')) +
  theme_bw() +
  guides(fill=guide_legend(override.aes=list(shape=21)),
         shape = 'none') +
  # scale_colour_continuous(name = 'Week') +
  ylab('Plastocyanin / Rubisco Small \nSubunit (1 and 2) ') +
  xlab('Coarse Grained Photosystem Protein Mass Fraction') +
  ggtitle('Targeted Proteomics: Phaeocystis (Haptophytes)');targeted_phaeo_proteins

merged_coarse_comparison <- ggarrange(ggarrange(plastocyanin_diatoms, 
          plastocyanin_haptophytes, nrow = 1, labels = c('a', 'b'),
          common.legend = TRUE),
          targeted_phaeo_proteins, nrow = 2, labels = c('', 'c'))

ggsave(merged_coarse_comparison, 
       filename = '../figures/coarse_grained_single_protein_biomarkers.png',
       width = 8.01, height = 6.47)

# baseline figures for functional characterization  ------------------------------------------------

total_peptides_mcl_weights <- tax_data_mcl %>% 
  group_by(ann_desc, day) %>% 
  summarize(total_peps_per_day_per_ann_desc = n())

weighted_mcl_tax_vals <- tax_data_mcl %>% 
  group_by(ann_desc, day, filter) %>% 
  summarize(total_mcl_mass = sum(db_abund_avg),
            total_mcl_peps_per_day_per_filter = n()) %>% 
  inner_join(total_peptides_mcl_weights, by = c('day', 'ann_desc')) %>% 
  inner_join(filter_protein_total_amounts_avg, by = c('filter', 'day')) %>% 
  mutate(multipled_mcl_vals = fraction_per_day*total_mcl_mass*total_mcl_peps_per_day_per_filter/total_peps_per_day_per_ann_desc) %>% 
  group_by(ann_desc, day) %>% 
  summarize(weighted_mcl_vals = sum(multipled_mcl_vals))

highest_mcl_ann_desc <- weighted_mcl_tax_vals %>%
  group_by(ann_desc) %>% 
  summarize(sum_all_days = sum(weighted_mcl_vals))

highest_mcl_ann_desc_char <- highest_mcl_ann_desc[highest_mcl_ann_desc$sum_all_days > quantile(highest_mcl_ann_desc$sum_all_days, 0.95) %>% as.numeric, ]$ann_desc

weighted_mcl_tax_vals_plot <- weighted_mcl_tax_vals %>% 
  filter(ann_desc %in% highest_mcl_ann_desc_char,
         !is.na(ann_desc)) %>% 
  ggplot() +
  theme_bw() +
  geom_tile(aes(fill = weighted_mcl_vals, 
                y = fct_reorder(ann_desc, weighted_mcl_vals), 
                x = day)) +
  ylab('MCL Cluster Consensus Annotation') +
  xlab('Week') +
  scale_fill_gradient(name = 'Weighted Protein\nGroup Abundance',
                      high = 'darkred', low = 'grey80');weighted_mcl_tax_vals_plot

ggsave(weighted_mcl_tax_vals_plot, 
       filename = "../figures/weighted_mcl_tax_vals_weeks.png",
       width = 10.8, height = 7.52)


# figure for unknown proteins ---------------------------------------------

## making figure showing the abundance of non-annotated MCL clusters

temp_mcl_weights_df <- weighted_mcl_tax_vals %>% 
  filter(is.na(ann_desc) | ann_desc == 'Bacterial extracellular solute-binding proteins, family 3')

temp_mcl_weights_df$ann_desc <- c(rep('Bacterial extracellular\nsolute-binding proteins, family 3', 4), 
                                  rep('Unknown Protein Cluster', 4))
unknown_protein_figure <- temp_mcl_weights_df %>% 
  ggplot() +
  facet_grid(~day) +
  coord_flip() +
  theme_bw() +
  geom_col(aes(y = weighted_mcl_vals, 
               x = ann_desc,
               fill = weighted_mcl_vals)) +
           # fill = 'darkred', 
           # alpha = 0.8) +
  scale_fill_gradient(name = 'Weighted Protein\nGroup Abundance',
                      high = 'darkred', low = 'grey80') +
  ylab('Weighted Protein Group Abundance') +
  theme(legend.position = 'None') +
  xlab('Cluster\nAnnotation');unknown_protein_figure


# baseline figures for taxonomic characterization  ------------------------------------------------
total_peptides_tax_weights <- tax_data_mcl %>% 
  group_by(tax_assign, day) %>% 
  summarize(total_peps_per_day_per_tax = n())

weighted_tax_vals <- tax_data_mcl %>% 
  group_by(tax_assign, day, filter) %>% 
  summarize(total_tax_mass = sum(db_abund_avg),
            total_tax_peps_per_day_per_filter = n()) %>% 
  inner_join(total_peptides_tax_weights, 
             by = c('day', 'tax_assign')) %>% 
  inner_join(filter_protein_total_amounts_avg, 
             by = c('filter', 'day')) %>% 
  mutate(multipled_tax_vals = fraction_per_day*total_tax_mass*total_tax_peps_per_day_per_filter/total_peps_per_day_per_tax) %>%
  # mutate(multipled_tax_vals = total_tax_mass*total_tax_peps_per_day_per_filter/total_peps_per_day_per_tax) %>% 
  group_by(tax_assign, day) %>% 
  summarize(weighted_tax_vals = sum(multipled_tax_vals))

# changing the name of the taxonomic assignments so they are nicer
levels(weighted_tax_vals$tax_assign)[5] <- 'Ciliates'
levels(weighted_tax_vals$tax_assign)[8] <- 'Diatoms'
levels(weighted_tax_vals$tax_assign)[9] <- 'Dinoflagellates'
levels(weighted_tax_vals$tax_assign)[20] <- 'Haptophytes'
levels(weighted_tax_vals$tax_assign)[19] <- 'Gammaproteobacteria'

weighted_tax_summarized <- weighted_tax_vals %>% 
  group_by(tax_assign) %>% 
  summarize(summed_weighted_tax = sum(weighted_tax_vals))

tax_to_subset <- weighted_tax_summarized[weighted_tax_summarized$summed_weighted_tax > quantile(weighted_tax_summarized$summed_weighted_tax, 0.65), ]$tax_assign %>% as.character()

# modifying the date label so it looks nice in ggplot2
weighted_tax_vals$week <- paste0('Week ', weighted_tax_vals$day)

weighted_tax_val_plot <- weighted_tax_vals %>% 
  filter(tax_assign %in% tax_to_subset,
         tax_assign != 'promiscuous-peptide',
         tax_assign != 'no-assignment-contig') %>% 
  ggplot(aes(x = fct_reorder(tax_assign, weighted_tax_vals), 
             y = weighted_tax_vals)) +
  geom_col(aes(fill = tax_assign)) +
  facet_grid(~week) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'None') +
  xlab('Taxonomic Group') +
  ylab('Taxon-Specific Peptide Abundance (weighted across filter sizes)');weighted_tax_val_plot


# merging summary figures -------------------------------------------------

intro_to_samples_figure <- ggarrange(weighted_tax_val_plot, 
          weighted_mcl_tax_vals_plot,
          nrow = 2, heights = c(1, 2.4),
          labels = c('a', 'b'))

intro_to_samples_figure_wrandom <- ggarrange(ggarrange(weighted_tax_val_plot, 
          unknown_protein_figure, nrow = 2, heights = c(2.5, 1), 
          align = 'hv', labels = c('a', 'b')),
          weighted_mcl_tax_vals_plot, nrow = 2, labels = c('', 'c'),
          heights = c(1, 1.25))


ggsave(intro_to_samples_figure, filename = '../figures/intro_to_samples_mcl_tax.png',
       width = 8, height = 9)

ggsave(intro_to_samples_figure_wrandom, filename = '../figures/intro_to_samples_mcl_tax_wrandom.png',
       width = 12.8, height = 9.61)
