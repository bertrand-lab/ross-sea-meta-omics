## formatting cofragmentation results that are peptide specific

library(dplyr)
library(ggplot2)
library(stringr)

filter_01_week1 <- read.csv('../data/digested-databases/GOS_927_0_1_nonredun_single_tryptic_peptide_rt_oligo_lc-retention-times_mi-0.00833333_ipw-1.44_para-7_co-sim.csv')
filter_08_week1 <- read.csv('../data/digested-databases/GOS_927_0_8_nonredun_single_tryptic_peptide_rt_oligo_lc-retention-times_mi-0.00833333_ipw-1.44_para-7_co-sim.csv')
filter_30_week1 <- read.csv('../data/digested-databases/GOS_927_3_0_nonredun_single_tryptic_peptide_rt_oligo_lc-retention-times_mi-0.00833333_ipw-1.44_para-7_co-sim.csv')

filter_01_week2 <- read.csv('../data/digested-databases/GOS_930_0_1_nonredun_single_tryptic_peptide_rt_oligo_lc-retention-times_mi-0.00833333_ipw-1.44_para-7_co-sim.csv')
filter_08_week2 <- read.csv('../data/digested-databases/GOS_930_0_8_nonredun_single_tryptic_peptide_rt_oligo_lc-retention-times_mi-0.00833333_ipw-1.44_para-7_co-sim.csv')
filter_30_week2 <- read.csv('../data/digested-databases/GOS_930_3_0_nonredun_single_tryptic_peptide_rt_oligo_lc-retention-times_mi-0.00833333_ipw-1.44_para-7_co-sim.csv')


cofrag_files <- dir(path = '../data/digested-databases/')[grepl(pattern = 'co-sim', 
                                                                x = dir(path = '../data/digested-databases/'))]
extract_filter_date <- function(co_sim_file_name){
  filter_date <- str_sub(string = co_sim_file_name, start = 1, end = 7)
  if(filter_date == 'GOS_927'){filter_week <- 1}
  if(filter_date == 'GOS_930'){filter_week <- 2}
  if(filter_date == 'GOS_933'){filter_week <- 3}
  if(filter_date == 'GOS_935'){filter_week <- 4}
  return(filter_week)
}
extract_filter_size <- function(co_sim_file_name){
  filter_size <- str_sub(string = co_sim_file_name, start = 9, end = 11)
  filter_size_period <- sub(pattern = '_', replacement = '.', x = filter_size, fixed = TRUE) %>% as.numeric()
  return(filter_size_period)
}

extract_filter_size(cofrag_files[1])
extract_filter_date(cofrag_files[10])


read_in_cofrag_output <- function(cofrag_out_location = '../data/digested-databases/',
                                  file_names = cofrag_files){
  
  empty_df <- data.frame(X = numeric(),
                         mean_cofrag_score = numeric(),
                         median_cofrag_score = numeric(),
                         sd_cofrag_score = numeric(),
                         pep_seq = numeric(),
                         day = numeric(),
                         filter_size = numeric())
  for(i in 1:length(cofrag_files)){
    temp_df <- read.csv(file = paste0(cofrag_out_location, cofrag_files[i]))
    week_to_add <- rep(extract_filter_date(cofrag_files[i]), nrow(temp_df))
    filter_to_add <- rep(extract_filter_size(cofrag_files[i]), nrow(temp_df))
    temp_df$day <- week_to_add
    temp_df$filter_size <- filter_to_add
    empty_df <- rbind(empty_df, temp_df)
  }
  return(empty_df)
}

cofrag_analysis_output <- read_in_cofrag_output()

cofrag_analysis_output %>%
  group_by(day, pep_seq) %>% 
  summarize(average_across_filters = mean(mean_cofrag_score, na.rm = TRUE)) %>% 
  ggplot(aes(x = pep_seq, y = average_across_filters)) +
  geom_point(aes(shape = day %>% as.factor()), 
             size = 4, alpha = 0.6) +
  theme_bw()

### this command is to find the taxon mappings, which needs data from 'getting_coarse_grained_all_taxa...R'
# pep_abundance_greedy %>% 
#   filter(grepl("plastocyanin", kegg_desc)) %>% 
#   group_by(peptide, tax_assign) %>% summarize(mean_blah = mean(db_abund_avg))
peptide_to_taxa <- data.frame(pep_seq = c('GGPHNVVFVEDAIPK',
                                          'GGPHNVVFDEDTLPAGVSQEK',
                                          'MGSDSGQLVFVPDEVK', 
                                          'VVQMGTDSGGLK'),
                              tax_assign = rep(c('Haptophytes', 'Diatoms'), 2))

cofrag_analysis_output_summarized <- cofrag_analysis_output %>%
  group_by(day, pep_seq) %>% 
  summarize(average_across_filters = mean(mean_cofrag_score, na.rm = TRUE)) %>% 
  inner_join(peptide_to_taxa, by = 'pep_seq') %>% 
  ungroup() %>% 
  group_by(day, tax_assign) %>% 
  summarize(mean_cofrag_score = mean(average_across_filters, na.rm = TRUE))

write.csv(x = cofrag_analysis_output_summarized,
          file = '../data/cobia-analysis/cobia_analysis_summarized.csv',
          row.names = FALSE)
