
## analysis of peptide discoverability between ribosomes and photosynthetic machinery

library(ggplot2)
library(magrittr)
library(dplyr)


##### must first load the R script 'getting_coarse_grained_all_taxa_....R'

# get all taxon specific, protein pools, protein ID's

all_greedy_discoverability <- all_greedy
all_greedy_discoverability$coarse_grains <- dplyr::recode(all_greedy_discoverability$coarse_grains,
                                                         P = "Photosynthetic Proteins",
                                                         R = "Ribosomal Proteins")
number_peps_discoverability <- all_greedy_discoverability %>% 
  dplyr::filter(tax_assign %in% selected_taxa_big,
                coarse_grains %in% c('Photosynthetic Proteins', 
                                     'Ribosomal Proteins')) %>% 
  group_by(coarse_grains, tax_assign, day) %>% 
  summarize(number_unique_peptides = n()) %>% 
  ggplot(aes(x = tax_assign, 
             y = number_unique_peptides)) +
  geom_point(size = 4) +
  facet_grid(~coarse_grains) +
  ylab('Number of Unique Peptides') +
  xlab('Taxonomic Assignment') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1));number_peps_discoverability

ggsave('../figures/number_peptides_coarse_grain_discoverability.png',
       width = 6.8, height = 3.85)

### must run the plotting_metapep_sim.R script first

arranged_discover <- ggarrange(sum_cor, number_peps_discoverability, labels = c('a', 'b'),
          nrow = 2, ncol = 1)

ggsave(arranged_discover, 
       filename = '../figures/discoverability_sim_observed.png',
       width = 8.76, height = 8.45)
