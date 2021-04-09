# poster plotting normalization

library(ggplot2)
library(dplyr)
library(testthat)
library(magrittr)

# bf_norm <- read.csv("../data/bf_normalization_factors.csv")
# db_id_norm <- read.csv("../data/db-search-id-only-normalization.csv")
# db_norm <- read.csv("../data/db-dependent-normalization-factors.csv")
pp_norm <- read.csv("../data/pp_normalization_factors.csv")

tfg_t0_norm <- read.csv("../data/normalization-data/norm-tfg-t0.csv")
tfg_all_norm <- read.csv("../data/normalization-data/norm-tfg-all-database.csv")
one_sample_norm <- read.csv("../data/normalization-data/norm-one-sample.csv")
pooled_norm <- read.csv("../data/normalization-data/norm-pooled-database.csv")
specific_norm <- read.csv("../data/normalization-data/norm-specific-database.csv")

length(pp_norm$file_name) == length(tfg_t0_norm$rep_file_list)

filter_size <- rep(c(0.1, 0.1, 0.8, 0.8, 0.8, 3.0, 3.0, 3.0), 4)
time_point <- c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8))

# format_norm_out <- function(norm_df){
#   expect_is(nrow(norm_df), 32)
#   norm_df$filter_size <- filter_size
#   norm_df$time_point <- time_point
#   norm_df$db_used
# }

aggregate_norm <- data.frame(TIC = pp_norm$tic,
                             `TFG All` = tfg_all_norm$norm_factors,
                             `TFG T=0` = tfg_t0_norm$norm_factors,
                             `One Sample` = one_sample_norm$norm_factors,
                             `Pooled by Size Fraction` = pooled_norm$norm_factors,
                             `Sample Specific` = specific_norm$norm_factors,
                             filter_size = filter_size,
                             time_point = time_point)

library(GGally)

norm01 <- aggregate_norm %>% 
  dplyr::filter(filter_size == 0.1) %>% 
  select(-c(filter_size, time_point)) %>% 
  ggpairs() +
  theme_bw() +
  ggtitle('Filter Size: 0.1um')
norm08 <- aggregate_norm %>% 
  dplyr::filter(filter_size == 0.8) %>% 
  select(-c(filter_size, time_point)) %>% 
  ggpairs() +
  theme_bw() +
  ggtitle('Filter Size: 0.8um')
norm30 <- aggregate_norm %>% 
  dplyr::filter(filter_size == 3.0) %>% 
  select(-c(filter_size, time_point)) %>% 
  ggpairs() +
  theme_bw() +
  ggtitle('Filter Size: 3.0um')

ggsave(norm01, filename = '../figures/normalization-factor-comp-01.png',
       width = 12.1, height = 9.61)
ggsave(norm08, filename = '../figures/normalization-factor-comp-08.png',
       width = 12.1, height = 9.61)
ggsave(norm30, filename = '../figures/normalization-factor-comp-30.png',
       width = 12.1, height = 9.61)
