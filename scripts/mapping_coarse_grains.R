## script for taking total model output and determining diatom proteome annotations

library(ggplot2)
library(magrittr)
library(dplyr)

# read in peptide quant data where tax_assigner has been run already

source("post_processing_functions.R")

tax_data <- read.csv("../data/go-tax-pipeline-output/pipeline_outtfg-all_tax_peptides.csv")

# reading in the model output
mnfe_meta_day1 <- read.csv("../../mn-fe-allocation/data/model_output/mnfe4_model2_metaprot5_1_day_1_variable_parameters_mn_fe19.csv")
mnfe_meta_day3 <- read.csv("../../mn-fe-allocation/data/model_output/mnfe4_model2_metaprot5_1_day_3_variable_parameters_mn_fe19.csv")

mnfe_meta <- rbind(mnfe_meta_day1, mnfe_meta_day3)
mnfe_meta$day <- c(1, 3)

# transforming the model output to be the same as the metaproteomic output
mnfe3_meta_trans <- mnfe_meta %>% 
  dplyr::select(A_aa, P_aa, R_aa, Tfe_aa, 
                Tn_aa, day) %>% 
  dplyr::rename(A = A_aa, 
         P = P_aa,
         R = R_aa,
         Tfe = Tfe_aa,
         # Tmn = Tmn_aa,
         Tn = Tn_aa) %>% 
  melt(value.name = 'mean_ion_int_coarse_percent', 
       variable.name = 'coarse_grains', id.vars = c("day"))
# 
sum_ints <- mnfe3_meta_trans %>%
  group_by(day) %>%
  summarize(total_fracs = sum(mean_ion_int_coarse_percent))
# 
mnfe3_meta_trans2 <- mnfe3_meta_trans %>%
  inner_join(sum_ints, by = c("day")) %>%
  mutate(trans_coarse_abundance = mean_ion_int_coarse_percent/total_fracs)


### taking just model output

mnfe_meta_trans_just_frac <- mnfe_meta %>% 
  dplyr::select(A_frac, P_frac, R_frac, Tfe_frac, 
                Tn_frac, day) %>% 
  dplyr::rename(A = A_frac, 
                P = P_frac,
                R = R_frac,
                Tfe = Tfe_frac,
                # Tmn = Tmn_aa,
                Tn = Tn_frac) %>% 
  melt(value.name = 'mean_ion_int_coarse_percent', 
       variable.name = 'coarse_grains', id.vars = c("day"))

coarse_diatoms <- read.csv("coarse_diatoms_tfg_all.csv")
coarse_diatoms <- read.csv("../data/go-tax-pipeline-output/coarse_diatoms_tfg_all_specific.csv")

# coarse_diatoms <- get_coarse_grains(tax_group = 'Diatom', 
#                                     tax_assign_df = tax_data, 
#                                     filter_type = 3.0)
# 
# coarse_haptophyta <- get_coarse_grains(tax_group = 'Other Haptophyte',
#                                 tax_assign_df = tax_data,
#                                 filter_type = 3.0)

tax_data_30 <- tax_data %>% filter(filter == 3.0)

# coarse_all <- get_coarse_grains(tax_group = "all", 
#                                 tax_assign_df = tax_data_30[c(1:100),],
#                                 filter_type = 3.0, all_taxa = TRUE)

# write.csv(coarse_all, file = 'coarse_all_correction_factors_greedy.csv')

coarse_diatoms_sum_prot_day <- coarse_diatoms %>% 
  # filter(coarse_grains != "U") %>% 
  group_by(coarse_grains, day, filter) %>% 
  summarize(mean_ion_int_coarse = mean(db_abund_avg)) %>% 
  ungroup() %>% 
  group_by(day, filter) %>% 
  summarize(total_prot_per_day_no_u = sum(mean_ion_int_coarse))

coarse_diatoms_meta_out <- coarse_diatoms %>% 
  # filter(coarse_grains != "U") %>% 
  # filter(day == 1 | day == 3) %>% 
  group_by(coarse_grains, day, filter) %>% 
  summarize(mean_ion_int_coarse = mean(db_abund_avg),
            number_peps = n()) %>% 
  inner_join(coarse_diatoms_sum_prot_day, by = c('day', 'filter')) %>% 
  mutate(mean_ion_int_coarse_percent = mean_ion_int_coarse/total_prot_per_day_no_u) 


coarse_diatoms_meta_out %>% 
  filter(coarse_grains != 'U') %>% 
  ggplot(aes(x = day, 
             y = mean_ion_int_coarse_percent)) +
    geom_point() +
  facet_grid(~coarse_grains) +
  theme_bw() +
  ylab('Fraction of Diatom Proteome') +
  xlab('') +
  scale_size_continuous(name = 'Number of \nMapped \nPeptides') +
  theme(axis.text.x.bottom = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 16)) +
  ggtitle('Coarse Grained Proteins for Diatoms') +
  # geom_point(data = mnfe_meta_trans_just_frac, 
  #            mapping = aes(x = day, y = mean_ion_int_coarse_percent),
  #            colour = 'blue', size = 3) +
  theme(strip.background = element_rect(fill = 'white'))
  

coarse_haptophyta_sum_prot_day <- coarse_haptophyta %>% 
  filter(coarse_grains != "U") %>% 
  group_by(coarse_grains, day, filter) %>% 
  summarize(mean_ion_int_coarse = mean(db_abund_avg)) %>% 
  ungroup() %>% 
  group_by(day, filter) %>% 
  summarize(total_prot_per_day_no_u = sum(mean_ion_int_coarse))

coarse_haptophyta %>% 
  filter(coarse_grains != "U") %>% 
  # filter(day == 1 | day == 3) %>% 
  group_by(coarse_grains, day, filter) %>% 
  summarize(mean_ion_int_coarse = mean(db_abund_avg),
            number_peps = n()) %>% 
  inner_join(coarse_haptophyta_sum_prot_day, by = c('day', 'filter')) %>% 
  mutate(mean_ion_int_coarse_percent = mean_ion_int_coarse/total_prot_per_day_no_u) %>% 
  ggplot(aes(x = day, 
             y = mean_ion_int_coarse_percent,
             size = number_peps)) +
  geom_point() +
  facet_grid(~coarse_grains) +
  theme_bw() +
  ylab('Percent of Mapped Proteome') +
  xlab('') +
  scale_size_continuous(name = 'Number of \nMapped \nPeptides') +
  theme(axis.text.x.bottom = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 16)) +
  ggtitle('Coarse Grained Proteins for Haptophytes')

coarse_haptophyta %>% filter(coarse_grains == "R")





