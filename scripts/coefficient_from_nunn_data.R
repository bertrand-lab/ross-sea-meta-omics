# estimating the coefficient of variation across fixed conditions, which means that
# all of the proteome should be 'fixed'

library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)

nunn_data <- read.csv("../../mn-fe-allocation/data/culture_proteomes/Nunn2013_Table_S1_annotated_formatted.csv")

# getting the sum of spectral counts per coarse grained protein group for + Fe
plus_fe_cv <- nunn_data %>%
  filter(Primary.vs.Secondary.IDs.homologous.protein.sequences. == 'primary') %>% 
  rowwise() %>% 
  mutate(tp_fe_mean = mean(c(tp_fe1, tp_fe2, tp_fe3, tp_fe4), na.rm = TRUE),
         tp_fe_sd = sd(c(tp_fe1, tp_fe2, tp_fe3, tp_fe4), na.rm = TRUE),
         tp_fe_cv = tp_fe_sd/tp_fe_mean)

minus_fe_cv <- nunn_data %>%
  filter(Primary.vs.Secondary.IDs.homologous.protein.sequences. == 'primary') %>% 
  rowwise() %>% 
  mutate(tp_nofe_mean = mean(c(tp_nofe1, tp_nofe2, tp_nofe3, tp_nofe4), na.rm = TRUE),
         tp_nofe_sd = sd(c(tp_nofe1, tp_nofe2, tp_nofe3, tp_nofe4), na.rm = TRUE),
         tp_nofe_cv = tp_nofe_sd/tp_nofe_mean)

# histogram of all cv for both minus and plus cv
histogram_vals_cv <- data.frame(cv_proteins = c(plus_fe_cv$tp_fe_cv,
                                                minus_fe_cv$tp_nofe_cv),
                                mean_proteins = c(plus_fe_cv$tp_fe_mean,
                                                  minus_fe_cv$tp_nofe_mean),
                                cv_condition = c(rep('Fe', nrow(plus_fe_cv)),
                                                 rep('noFe', nrow(minus_fe_cv))))

nunn_protein_histogram_cv <- histogram_vals_cv %>% 
  filter(cv_proteins > 0) %>% 
  ggplot(aes(x = cv_proteins)) +
  geom_histogram() +
  theme_bw() +
  ylab('Count') +
  xlab('Protein-specific Coefficient of Variation\nacross constant conditions');nunn_protein_histogram_cv
  # facet_wrap(~cv_condition, nrow = 2)

mean_cv_cor_nunn <- histogram_vals_cv %>% 
  filter(cv_proteins > 0) %>% 
  ggplot(aes(x = mean_proteins, y = cv_proteins)) +
  geom_point(alpha = 0.8) +
  theme_bw() +
  ylab('Protein-specific Coefficient of Variation') +
  xlab('Mean Protein Abundance Value');mean_cv_cor_nunn


nunn_supp_cv_plot <- ggarrange(nunn_protein_histogram_cv, 
                               mean_cv_cor_nunn, nrow = 1,
                               labels = c('a', 'b'))

ggsave(nunn_supp_cv_plot, filename = '../figures/nunn_supp_cv_plot.png',
       width = 9.02, height = 5.75)


histogram_vals_cv$cv_proteins[histogram_vals_cv$cv_proteins > 0] %>% quantile(c(0.4, 0.5, 0.6))


