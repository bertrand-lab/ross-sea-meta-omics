## this script looks at the metaproteomic simulation results. It calculates the proposed correction factors.

library(dplyr)
library(ggplot2)
library(magrittr)
library(gridExtra)
library(ggpubr)
library(forcats)


# reading in data and formatting ------------------------------------------

sampled_data <- read.csv("../data/pep-sample-sim/sampled_proper_observed_orgs30_n_seqs2000_5prot_groups_updated_ambiguous.csv")
actual_data <- read.csv("../data/pep-sample-sim/actual_proper_observed_orgs30_n_seqs2000_5prot_groups_updated_ambiguous.csv")
observed_data <- read.csv("../data/pep-sample-sim/observed_proper_observed_orgs30_n_seqs2000_5prot_groups_updated_ambiguous.csv")

sampled_data$organism <- as.factor(sampled_data$organism)
observed_data$organism <- as.factor(observed_data$organism)
actual_data$organism <- as.factor(actual_data$organism)

sampled_data$protein_pool <- as.factor(sampled_data$protein_pool)
observed_data$protein_pool <- as.factor(observed_data$protein_pool)
actual_data$protein_pool <- as.factor(actual_data$protein_pool)

levels(actual_data$protein_pool) <- c("Super Low", "Lowest Diversity Protein Group", "Low Diversity Protein Group", 
                                      "High Diversity Protein Group", "Highest Diversity Protein Group")
levels(sampled_data$protein_pool) <- c("Super Low", "Lowest Diversity Protein Group", "Low Diversity Protein Group", 
                                       "High Diversity Protein Group", "Highest Diversity Protein Group", "ambiguous")

levels(actual_data$protein_pool) <- seq(1:5)
levels(sampled_data$protein_pool) <- c(seq(1:5), NA)

# sampled_data <- sampled_data %>% filter(trial == 2)
# actual_data <- actual_data %>% filter(trial == 2)
# observed_data <- observed_data %>% filter(trial == 2)


# # Check that the ambigious number of sequences is different acro --------

sampled_data %>% 
  filter(organism != 'ambiguous') %>%
  # filter(protein_pool == 'Super Low') %>% 
  group_by(trial, protein_pool, organism) %>% 
  summarize(number_per_unique = n()) %>% 
  ggplot(aes(y = number_per_unique, x = protein_pool)) + 
  geom_boxplot() +
  # scale_y_log10() +
  geom_jitter()

sampled_data %>% 
  filter(organism != 'ambiguous') %>%
  # filter(protein_pool == 'Super Low') %>% 
  group_by(trial, protein_pool) %>% 
  summarize(number_per_unique = n()) %>% 
  ggplot(aes(y = number_per_unique, x = protein_pool)) + 
  geom_boxplot() +
  # scale_y_log10() +
  geom_jitter()


# # summarizing dataset ---------------------------------------------------

trial_s_1 <- sampled_data %>% #filter(trial == 1) %>% 
  group_by(protein_pool, organism, trial) %>% 
  summarize(mean_ion_int = mean(ion_intensity),
            sum_ion_int = sum(ion_intensity),
            median_ion_int = median(ion_intensity),
            number_of_peps = n())

actual_s_1 <- actual_data %>% #filter(trial == 1) %>% 
  group_by(protein_pool, organism, trial) %>% 
  summarize(mean_ion_int = mean(ion_intensity),
            sum_ion_int = sum(ion_intensity),
            median_ion_int = median(ion_intensity))

actual_s_1 <- actual_s_1 %>% dplyr::rename(mean_ion_int_actual = mean_ion_int,
                                    sum_ion_int_actual = sum_ion_int,
                                    median_ion_int_actual = median_ion_int)

# combining the simulated true and simulated sampled
tester <- inner_join(trial_s_1, actual_s_1, 
                     by = c('protein_pool', 'organism', 'trial'))

tester$protein_pool <- as.numeric(tester$protein_pool)

sum_cor <- tester %>% 
  filter(organism != 'ambiguous') %>%
  ggplot(aes(y = sum_ion_int, 
             x = sum_ion_int_actual)) + 
  geom_point(alpha = 0.5, 
             aes(colour = as.numeric(protein_pool),
                 size = number_of_peps)) +
  # geom_bin2d() +
  # facet_grid(~protein_pool) + 
  scale_x_log10() + scale_y_log10() +
  theme_bw() + 
  scale_size_continuous('Number of Unique Peptides') +
  # geom_smooth(method = 'lm', aes(group = protein_pool, colour = protein_pool)) +
  geom_abline(intercept = 0) +
  scale_colour_continuous('Protein Pool Diversity') +
  # stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~"))) +
  ylab('Observed Summed Ion Intensity \nper Protein Group') + 
  xlab('Actual Summed Ion Intensity per Protein Group');sum_cor


# plotting root mean sq error ---------------------------------------------
# 
# rmse_sum_dist <- tester %>%
#   filter(organism != 'ambiguous',
#          protein_pool != 'Low Diversity Protein Group', 
#          protein_pool != 'High Diversity Protein Group',
#          protein_pool != 'ambiguous') %>%
#   ggplot(aes(sqrt((sum_ion_int - sum_ion_int_actual)^2),
#              fill = protein_pool)) +
#   geom_histogram(alpha = 0.8) +
#   # geom_bin2d() +
#   # facet_grid(~protein_pool) + 
#   scale_x_log10() +
#   scale_y_log10() +
#   theme_bw() +
#   scale_fill_discrete('Protein Pool') +
#   facet_wrap(~protein_pool, nrow = 4) +
#   ylab('Peptide Count') + 
#   theme(legend.position = "none",
#         strip.background = element_rect(fill = 'white'),
#         strip.text = element_text(size = 12),
#         axis.title = element_text(size = 12)) +  xlab('RMSE of Summed Peptide Intensities per Group');rmse_sum_dist
# 
# rmse_mean_dist <- tester %>%
#   filter(organism != 'ambiguous',
#          protein_pool != 'Low Diversity Protein Group', protein_pool != 'High Diversity Protein Group') %>%
#   ggplot(aes(sqrt((mean_ion_int - mean_ion_int_actual)^2),
#              fill = protein_pool)) + 
#   geom_histogram() +
#   # geom_bin2d() +
#   # facet_grid(~protein_pool) + 
#   scale_x_log10() +
#   scale_y_log10() +
#   theme_bw() +
#   scale_fill_discrete('Protein Pool') +
#   facet_wrap(~protein_pool, nrow = 4) +
#   ylab('Peptide Count') + 
#   theme(legend.position = "none",
#         strip.background = element_rect(fill = 'white'),
#         strip.text = element_text(size = 12),
#         axis.title = element_text(size = 12)) +
#   xlab('RMSE of Mean Peptide Intensities per Group');rmse_mean_dist
# 
# 
# grid.arrange(rmse_mean_dist, rmse_sum_dist, nrow = 1)
# 
# 

# plotting correlation results --------------------------------------------

# mean_cor <- tester %>% 
#   filter(organism != 'ambiguous') %>%
#   ggplot(aes(y = mean_ion_int, 
#              x = mean_ion_int_actual)) + 
#   geom_point(alpha = 1, aes(colour = protein_pool)) +
#   # facet_grid(~protein_pool) + 
#   scale_x_log10() + 
#   scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   scale_colour_discrete('Protein Pool') +
#   # geom_smooth(method = 'lm', aes(group = protein_pool, colour = protein_pool)) +
#   # stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~", group = protein_pool))) +
#   ylab('Observed Mean Ion Intensity \nper Protein Group') + 
#   xlab('Actual Mean Ion Intensity per Protein Group') + 
#   ggtitle('Protein Groups ranging from low (1) to high (3) sequence diversity\nMean Peptide Intensity');mean_cor



# ggsave(sum_cor,
#        filename = 'simulated_peptide_diversity_counts_observability.png',
#        width = 8.76, height = 5.7)

# grid.arrange(mean_cor, sum_cor)



# implementing a correction factor ----------------------------------------
# 
# # what is the average ion intensity, or what is the total sum
# mean_mean_int <- trial_s_1 %>% 
#   group_by(trial) %>% 
#   summarize(mean_int = mean(mean_ion_int))
# 
# sum_mean_int <- trial_s_1 %>% 
#   group_by(trial) %>% 
#   summarize(sum_int = sum(mean_ion_int))
# 
# # correction factor based on organism abundance
# cor_factors_organism <- trial_s_1 %>% 
#   filter(organism != 'ambiguous') %>% 
#   group_by(organism, trial) %>% 
#   summarize(cor_factor_sum_org = sum(mean_ion_int), # this is the sum of peptides attributable to a taxa
#             cor_factor_mean_org = mean(mean_ion_int)) %>% # this is a weighting factor centered around 1
#   ungroup() %>% 
#   inner_join(mean_mean_int, by = 'trial') %>% 
#   inner_join(sum_mean_int, by = 'trial') %>% 
#   mutate(cor_factor_sum_org_prop = cor_factor_sum_org/sum_int,
#          cor_factor_mean_org_prop = cor_factor_mean_org/mean_int)
# 
# # correction factor based on protein pool grouping
# cor_factors_pool <- trial_s_1 %>% 
#   filter(protein_pool != 'ambiguous') %>% 
#   group_by(protein_pool, trial) %>% 
#   dplyr::summarize(cor_factor_sum_pool = sum(mean_ion_int), # this is the sum of peptides attributable to a taxa
#                    cor_factor_mean_pool = mean(mean_ion_int)) %>% # this is a weighting factor centered around 1
#   ungroup() %>% 
#   inner_join(mean_mean_int, by = 'trial') %>% 
#   inner_join(sum_mean_int, by = 'trial') %>% 
#   mutate(cor_factor_sum_pool_prop = cor_factor_sum_pool/sum_int,
#          cor_factor_mean_pool_prop = cor_factor_mean_pool/mean_int)
# 
# # aggregating the data
# model_out_cor <- inner_join(tester, cor_factors_organism, by = c('organism', 'trial'))
# model_out_cor <- inner_join(model_out_cor, cor_factors_pool, by = c('protein_pool', 'trial'))
# 
# # calculating the number of total peptides per protein pool per trial
# per_pool_per_trial <- model_out_cor %>% 
#   group_by(trial, protein_pool) %>% 
#   summarize(n_per_pool = sum(number_of_peps))
# 
# model_out_cor <- inner_join(model_out_cor, per_pool_per_trial, by = c('trial', 'protein_pool'))
# 
# 
# 
# # plotting correctino factors ---------------------------------------------
# 
# 
# mean_sum_cor_factor <- model_out_cor %>% 
#   filter(organism != 'ambiguous') %>%
#   mutate(mean_ion_int_cor = mean_ion_int*cor_factor_sum_org_prop) %>% 
#   ggplot(aes(y = mean_ion_int_cor, 
#              x = mean_ion_int_actual)) + 
#   geom_point(alpha = 0.6) +
#   facet_grid(~protein_pool) + 
#   scale_x_log10() + scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~"))) +
#   ylab('Corrected Mean Ion Intensity \nper Protein Group') + 
#   xlab('Actual Mean Ion Intensity \nper Protein Group') + 
#   ggtitle('\nTax Sum Correction Applied');mean_sum_cor_factor
# 
# mean_mean_cor_factor <- model_out_cor %>% 
#   filter(organism != 'ambiguous') %>%
#   mutate(mean_ion_int_cor = mean_ion_int*cor_factor_mean_org_prop) %>%
#   ggplot(aes(y = mean_ion_int_cor, 
#              x = mean_ion_int_actual)) + 
#   geom_point(alpha = 0.6) +
#   facet_grid(~protein_pool) + 
#   scale_x_log10() + scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~"))) +
#   ylab('Corrected Mean Ion Intensity \nper Protein Group') + 
#   xlab('Actual Mean Ion Intensity \nper Protein Group') + 
#   ggtitle('\nTax Mean Correction Applied');mean_mean_cor_factor
# 
# no_cor_factor <- model_out_cor %>% 
#   filter(organism != 'ambiguous') %>%
#   # mutate(mean_ion_int_cor = mean_ion_int*cor_factor_mean_org) %>% 
#   ggplot(aes(y = mean_ion_int, 
#              x = mean_ion_int_actual, 
#              colour = protein_pool)) + 
#   geom_point(alpha = 0.3) +
#   facet_grid(~protein_pool) + 
#   scale_x_log10() + scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~")), colour = 'black') +
#   ylab('Corrected Mean Ion Intensity \nper Protein Group') + 
#   theme(strip.background = element_rect(fill = 'white'),
#         legend.position = 'none') +
#   xlab('Actual Mean Ion Intensity per Protein Group') + 
#   ggtitle('\nNo Correction Factor');no_cor_factor
# 
# no_cor_factor_sum <- model_out_cor %>% 
#   filter(organism != 'ambiguous') %>%
#   # mutate(mean_ion_int_cor = mean_ion_int*cor_factor_mean_org) %>% 
#   ggplot(aes(y = sum_ion_int, 
#              x = sum_ion_int_actual, 
#              colour = protein_pool)) + 
#   geom_point(alpha = 0.9) +
#   # facet_grid(~protein_pool) + 
#   scale_x_log10() + scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   # stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~")), colour = 'black') +
#   ylab('Observed Summed Ion Intensity per Protein Group') + 
#   theme(strip.background = element_rect(fill = 'white'),
#         legend.position = 'none') +
#   xlab('Actual Summed Ion Intensity per Protein Group') + 
#   ggtitle('\nNo Correction Factor');no_cor_factor_sum
# 
# grid.arrange(no_cor_factor, 
#              mean_mean_cor_factor) 
#              # mean_sum_cor_factor)
# 
# 
# ### looking at protein pool correction factor
# 
# mean_sum_cor_factor_pool <- model_out_cor %>% 
#   filter(organism != 'ambiguous') %>%
#   mutate(mean_ion_int_cor = mean_ion_int*cor_factor_sum_pool_prop) %>% 
#   ggplot(aes(y = mean_ion_int_cor, 
#              x = mean_ion_int_actual)) + 
#   geom_point(alpha = 0.6) +
#   facet_grid(~protein_pool) + 
#   scale_x_log10() + scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~"))) +
#   ylab('Corrected Mean Ion Intensity \nper Protein Group') + 
#   xlab('Actual Mean Ion Intensity \nper Protein Group') + 
#   ggtitle('\nPool Sum Correction Applied');mean_sum_cor_factor_pool
# 
# mean_mean_cor_factor_pool <- model_out_cor %>% 
#   filter(organism != 'ambiguous') %>%
#   mutate(mean_ion_int_cor = mean_ion_int*cor_factor_mean_pool_prop) %>%
#   ggplot(aes(y = mean_ion_int_cor, 
#              x = mean_ion_int_actual)) + 
#   geom_point(alpha = 0.6) +
#   facet_grid(~protein_pool) + 
#   scale_x_log10() + scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~"))) +
#   ylab('Corrected Mean Ion Intensity \nper Protein Group') + 
#   xlab('Actual Mean Ion Intensity per Protein Group') + 
#   ggtitle('\nPool Mean Correction Applied');mean_mean_cor_factor_pool
# 
# 
# mean_sum_cor_factor_pool <- model_out_cor %>% 
#   filter(organism != 'ambiguous') %>%
#   mutate(mean_ion_int_cor = sum_ion_int*cor_factor_sum_pool_prop) %>% 
#   ggplot(aes(y = mean_ion_int_cor, 
#              x = mean_ion_int_actual)) + 
#   geom_point(alpha = 0.6) +
#   facet_grid(~protein_pool) + 
#   scale_x_log10() + scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~"))) +
#   ylab('Corrected Mean Ion Intensity \nper Protein Group') + 
#   xlab('Actual Mean Ion Intensity \nper Protein Group') + 
#   ggtitle('\nPool Sum Correction Applied');mean_sum_cor_factor_pool
# 
# mean_mean_cor_factor_pool <- model_out_cor %>% 
#   filter(organism != 'ambiguous') %>%
#   mutate(mean_ion_int_cor = sum_ion_int*cor_factor_mean_pool_prop) %>%
#   ggplot(aes(y = mean_ion_int_cor, 
#              x = mean_ion_int_actual)) + 
#   geom_point(alpha = 0.6) +
#   facet_grid(~protein_pool) + 
#   scale_x_log10() + scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~"))) +
#   ylab('Corrected Mean Ion Intensity \nper Protein Group') + 
#   xlab('Actual Mean Ion Intensity per Protein Group') + 
#   ggtitle('\nPool Mean Correction Applied');mean_mean_cor_factor_pool
# 
# grid.arrange(no_cor_factor, 
#              mean_mean_cor_factor_pool)
# 
# 
# # looking at both
# # correction factors calculated from both the sum of pool and sum of org applied to mean ion intensity
# mean_sum_cor_factor_pool_org <- model_out_cor %>% 
#   filter(organism != 'ambiguous') %>%
#   mutate(mean_ion_int_cor = mean_ion_int*cor_factor_sum_pool_prop*cor_factor_sum_org_prop) %>% 
#   ggplot(aes(y = mean_ion_int_cor, 
#              x = mean_ion_int_actual)) + 
#   geom_point(alpha = 1) +
#   facet_grid(~protein_pool) + 
#   scale_x_log10() + scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~"))) +
#   ylab('Corrected Mean Ion Intensity \nper Protein Group') + 
#   xlab('Actual Mean Ion Intensity \nper Protein Group') + 
#   ggtitle('\nPool/Tax Sum Correction Applied');mean_sum_cor_factor_pool_org
# 
# # correction factors calculated from both the mean of pool and mean of org applied to mean ion intensity
# mean_mean_cor_factor_pool_org <- model_out_cor %>% 
#   filter(organism != 'ambiguous') %>%
#   mutate(mean_ion_int_cor = mean_ion_int*cor_factor_mean_pool_prop*cor_factor_mean_org_prop) %>%
#   ggplot(aes(y = mean_ion_int_cor, 
#              x = mean_ion_int_actual,
#              colour = protein_pool)) + 
#   facet_grid(~protein_pool) + 
#   scale_x_log10() + scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~")), colour = 'black') +
#   ylab('Corrected Mean Ion Intensity \nper Protein Group') + 
#   xlab('Actual Mean Ion Intensity per Protein Group') + 
#   geom_point(alpha = 0.3) +
#   theme(strip.background = element_rect(fill = 'white'),
#         legend.position = 'none') +
#   ggtitle('\nCorrection Factor Applied');mean_mean_cor_factor_pool_org
# 
# 
# sum_cor_factor_pool_org <- model_out_cor %>% 
#   filter(organism != 'ambiguous') %>%
#   mutate(sum_ion_int_cor = sum_ion_int*cor_factor_sum_pool_prop*cor_factor_sum_org_prop) %>% 
#   ggplot(aes(y = sum_ion_int_cor, 
#              x = sum_ion_int_actual, colour = protein_pool)) + 
#   geom_point(alpha = 0.9) +
#   # facet_grid(~protein_pool) + 
#   scale_x_log10() + scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   # stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~")), colour = 'black') +
#   ylab('Corrected Summed Ion Intensity per Protein Group') + 
#   xlab('Actual Summed Ion Intensity per Protein Group') + 
#   theme(strip.background = element_rect(fill = 'white'),
#         legend.position = 'none') +
#   ggtitle('\nPool/Tax Sum Correction Applied');sum_cor_factor_pool_org
# 
# grid.arrange(no_cor_factor, 
#              mean_mean_cor_factor_pool_org)
# 
# cor_factor_comparison_with_sum <- ggarrange(no_cor_factor_sum, 
#              sum_cor_factor_pool_org)
# 
# ggsave(cor_factor_comparison_with_sum, filename = 'figures/correction_factor_comparison_with_sum.png')
# 
# sum_histogram_cor <- model_out_cor %>% 
#   filter(organism != 'ambiguous') %>%
#   mutate(sum_ion_int_cor = sum_ion_int*cor_factor_sum_pool_prop*cor_factor_sum_org_prop,
#          rmse = sqrt((sum_ion_int_cor - sum_ion_int_actual)^2)) %>% 
#   ggplot(aes(rmse)) +
#   geom_histogram() +
#   scale_x_log10() +
#   facet_wrap(~protein_pool, nrow = 5)
# 
# sum_histogram_no_cor <- model_out_cor %>% 
#   filter(organism != 'ambiguous') %>%
#   mutate(sum_ion_int_cor = sum_ion_int*cor_factor_sum_pool_prop*cor_factor_sum_org_prop,
#          rmse = sqrt((sum_ion_int - sum_ion_int_actual)^2)) %>% 
#   ggplot(aes(rmse)) +
#   geom_histogram() +
#   scale_x_log10() +
#   facet_wrap(~protein_pool, nrow = 5)
# 
# 
# grid.arrange(sum_histogram_no_cor, sum_histogram_cor, ncol = 2)
# 
# 
# ## implementing weighted correction
# mean_cor_factors_org <- cor_factors_organism %>%
#   group_by(trial, organism) %>% 
#   summarize(min_cor_factor_mean_org = mean(cor_factor_mean_org_prop))
# 
# mean_cor_factors_pool <- cor_factors_pool %>%
#   group_by(trial, protein_pool) %>% 
#   summarize(min_cor_factor_mean_pool = mean(cor_factor_mean_pool_prop))
# 
# 
# # weighted_cor <- tester %>%
# #   group_by(organism) %>%
# #   dplyr::summarize(weighted_up_mean_sum = (1/2)*mean_cor_factor_sum*(1/log(n())) + (1/2)*sum(mean_ion_int)/total_mean_int*(log(n())/log(n()) + 1),
# #                    weighted_up_mean_mean = (1/2)*mean_cor_factor_mean*(1/log(n())) + (1/2)*mean(mean_ion_int)/mean_mean_int*(log(n())/log(n()) + 1),
# #                    weighted_up_min_sum = (1/2)*min_cor_factor_sum*(1/log(n())) + (1/2)*sum(mean_ion_int)/total_mean_int*(log(n())/log(n()) + 1),
# #                    weighted_up_min_mean = (1/2)*min_cor_factor_mean*(1/log(n())) + (1/2)*mean(mean_ion_int)/mean_mean_int*(log(n())/log(n()) + 1))
# 
# weighted_cor_pool <- model_out_cor %>% 
#   group_by(organism, protein_pool, trial) %>% 
#   inner_join(mean_cor_factors_org, by = c("trial", "organism")) %>% 
#   inner_join(mean_cor_factors_pool, by = c("trial", "protein_pool")) %>%
#   dplyr::summarize(pool_correction = (1/2)*min_cor_factor_mean_pool*(1/number_of_peps) + 
#                      (1/2)*(cor_factor_mean_pool_prop)*(number_of_peps/(number_of_peps + 1)),
#                    org_correction = (1/2)*min_cor_factor_mean_org*(1/number_of_peps) + 
#                      (1/2)*(cor_factor_mean_pool_prop)*(number_of_peps/(number_of_peps + 1)),
#                    pool_org_weighted = pool_correction*org_correction)
# 
# weighted_cor_pool_data <- inner_join(model_out_cor, 
#                                      weighted_cor_pool, 
#                                      by = c("trial", "organism", "protein_pool"))
# 
# # cor_factor_mean_pool_prop*cor_factor_mean_org_prop
# 
# exp_cor_pool <- model_out_cor %>% 
#   mutate(cor_pool_above_or_below = ifelse(cor_factor_mean_pool_prop > 1, 1, -1),
#          cor_org_above_or_below = ifelse(cor_factor_mean_org_prop > 1, 1, -1),
#          exp_mean_cor_factor_pool = 1 + cor_pool_above_or_below*abs(1 - cor_factor_mean_pool_prop)/sqrt(number_of_peps),
#          exp_mean_cor_factor_org = 1 + cor_org_above_or_below*abs(1 - cor_factor_mean_org_prop)/sqrt(number_of_peps))
#          
#          # exp_mean_cor_factor_pool = cor_factor_mean_pool_prop*exp(0.05*cor_pool_above_or_below*number_of_peps),
#          # exp_mean_cor_factor_org = cor_factor_mean_org_prop*exp(0.05*cor_org_above_or_below*number_of_peps))
#   
# exp_cor_pool_data <- inner_join(model_out_cor, 
#                                 exp_cor_pool, 
#                                      by = c("trial", "organism", "protein_pool"))
# 
# cor_mean_mean <- exp_cor_pool_data %>% 
#   filter(organism != 'ambiguous') %>%
#   mutate(mean_ion_int_cor = mean_ion_int.x*exp_mean_cor_factor_org*exp_mean_cor_factor_pool) %>%  
#   ggplot(aes(y = mean_ion_int_cor, 
#              x = mean_ion_int_actual.x)) + 
#   geom_point(alpha = 1) +
#   facet_grid(~protein_pool) + 
#   scale_x_log10() + scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~"))) +
#   ylab('Corrected Mean Ion Intensity \nper Protein Group') + 
#   xlab('Actual Mean Ion Intensity \nper Protein Group') + 
#   ggtitle('\n Mean, Mean, Weighted Correction Applied');cor_mean_mean
# 
# 
# # cor_mean_sum <- weighted_cor_pool %>% 
# #   filter(organism != 'ambiguous') %>%
# #   mutate(mean_ion_int_cor = mean_ion_int*weighted_up_mean_sum) %>% 
# #   ggplot(aes(y = mean_ion_int_cor, 
# #              x = mean_ion_int_actual)) + 
# #   geom_point(alpha = 1) +
# #   facet_grid(~protein_pool) + 
# #   scale_x_log10() + scale_y_log10() +
# #   theme_bw() + 
# #   stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~"))) +
# #   geom_abline(intercept = 0) +
# #   ylab('Corrected Mean Ion Intensity \nper Protein Group') + 
# #   xlab('Actual Mean Ion Intensity \nper Protein Group') + 
# #   ggtitle('\n Mean, Sum, Weighted Correction Applied');cor_mean_sum
# 
# cor_mean_mean <- weighted_cor_pool_data %>% 
#   filter(organism != 'ambiguous') %>%
#   mutate(mean_ion_int_cor = mean_ion_int*pool_org_weighted) %>%  
#   ggplot(aes(y = mean_ion_int_cor, 
#              x = mean_ion_int_actual)) + 
#   geom_point(alpha = 1) +
#   facet_grid(~protein_pool) + 
#   scale_x_log10() + scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~"))) +
#   ylab('Corrected Mean Ion Intensity \nper Protein Group') + 
#   xlab('Actual Mean Ion Intensity \nper Protein Group') + 
#   ggtitle('\n Mean, Mean, Weighted Correction Applied');cor_mean_mean
# 
# cor_min_sum <- model_out_weighted %>% 
#   filter(organism != 'ambiguous') %>%
#   mutate(mean_ion_int_cor = mean_ion_int*weighted_up_min_sum) %>% 
#   ggplot(aes(y = mean_ion_int_cor, 
#              x = mean_ion_int_actual)) + 
#   geom_point(alpha = 1) +
#   facet_grid(~protein_pool) + 
#   scale_x_log10() + scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~"))) +
#   ylab('Corrected Mean Ion Intensity \nper Protein Group') + 
#   xlab('Actual Mean Ion Intensity \nper Protein Group') + 
#   ggtitle('\n Min, Sum, Weighted Correction Applied');cor_min_sum
# cor_min_mean <- model_out_weighted %>% 
#   filter(organism != 'ambiguous') %>%
#   mutate(mean_ion_int_cor = mean_ion_int*weighted_up_min_mean) %>% 
#   ggplot(aes(y = mean_ion_int_cor, 
#              x = mean_ion_int_actual)) + 
#   geom_point(alpha = 1) +
#   facet_grid(~protein_pool) + 
#   scale_x_log10() + scale_y_log10() +
#   theme_bw() + 
#   geom_abline(intercept = 0) +
#   # stat_cor() +
#   stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"))) +
#   ylab('Corrected Mean Ion Intensity \nper Protein Group') + 
#   xlab('Actual Mean Ion Intensity \nper Protein Group') + 
#   ggtitle('\n Min, Mean, Weighted Correction Applied');cor_min_mean
# 
# grid.arrange(cor_mean_sum, cor_mean_mean, 
#              cor_min_sum, cor_min_mean)
# 
# 
# model_out_weighted %>% 
#   ggplot(aes(mean_ion_int*weighted_up_min_mean, 
#              mean_ion_int*weighted_up_mean_mean)) +
#   geom_point()
# 
# 
# 
# ggsave(sum_cor, file = 'sum_cor.png')
# ggsave(mean_cor, file = 'mean_cor.png')
# # ggsave(rmse_sum_int, file = 'rmse_sum_int.png')
# # ggsave(rmse_mean_int, file = 'rmse_mean_int.png')
# ggsave(rmse_sum_dist, file = 'rmse_sum.png')
# ggsave(rmse_mean_dist, file = 'rmse_mean.png')
