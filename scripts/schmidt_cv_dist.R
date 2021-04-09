
library(ggplot2)
library(magrittr)
library(readxl)
library(ggridges)
library(dplyr)
library(tidyr)
library(reshape2)

cv_schmidt <- read_xlsx('../data/culture-proteomics-data/41587_2016_BFnbt3418_MOESM18_ESM.xlsx', sheet = 6, skip = 2)

cv_schmidt_sub <- cv_schmidt[ , c(1, 51:72)]

quantile(cv_schmidt_sub$`Glycerol + AA...53`/100, 0.75)
quantile(cv_schmidt_sub$LB...52/100, 0.75)
quantile(cv_schmidt_sub$`Glycerol + AA...53`/100, 0.75)
quantile(cv_schmidt_sub$`Glycerol + AA...53`/100, 0.75)


test <- melt(cv_schmidt_sub)

tquant <- test %>% 
  group_by(variable) %>% 
  mutate(value2 = value/100) %>% 
  dplyr::summarize(quant50 = quantile(value2, probs = c(0.5)),
                   quant75 = quantile(value2, probs = c(0.75)),
                   quant90 = quantile(value2, probs = c(0.90)),
                   quant95 = quantile(value2, probs = c(0.95)))

test %>% 
  ggplot(aes(y = variable, x = value/100)) +
  ggridges::geom_density_ridges()

cv_dist_schmidt <- test %>% 
  ggplot(aes(x = value/100)) +
  geom_histogram(binwidth = 0.01) +
  ylab('Count') +
  xlab('Protein-level coefficient of variation across constant conditions') +
  theme_bw() +
  geom_vline(xintercept = tquant$quant75 %>% mean())


tquant$quant75 %>% mean()

cv_schmidt %>% 
  ggplot(aes(x = `Glucose...7`,
             y= `Glucose...51`)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

cv_schmidt %>% 
  ggplot(aes(x = `Xylose...24`,
             y= `Xylose...68`)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

cv_schmidt_mean <- cv_schmidt[,c(1, 7:28)]

schmidt_mean <- melt(cv_schmidt_mean) %>% extract(variable, "condition")
schmidt_cv <- test %>% extract(variable, "condition")
names(schmidt_cv)[3] <- c("cv")

mean_cv_df <- data.frame(uniprot_acc = schmidt_mean$`Uniprot Accession`,
                         condition = schmidt_mean$condition, 
                         mean_val = schmidt_mean$value,
                         cv_val = schmidt_cv$cv)

mean_cv_relationship <- mean_cv_df %>% 
  ggplot(aes(x = cv_val, y = mean_val)) +
  geom_point(alpha = 0.2) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  ylab('Mean Protein Abundance') +
  xlab('Protein Abundance Coefficient of Variation')

combined_plot <- ggarrange(cv_dist_schmidt, mean_cv_relationship, nrow = 1, ncol = 2, labels = c('a', 'b'))

ggsave(combined_plot, 
       filename = '../figures/mean_cv_relationship_schmidt.png', width = 10.4, height = 6.57)

mean_cv_df %>% 
  filter(cv_val < 5 & mean_val < 1) %>% 
  ggplot(aes(x = cv_val, 
             y = mean_val)) +
  geom_point(alpha = 0.9) +
  scale_y_log10() +
  scale_x_log10() +
  geom_text(aes(label = uniprot_acc))

mean_cv_df %>% 
  filter(cv_val < 5 & mean_val < 1) %>% 
  ggplot(aes(x = cv_val, 
             y = mean_val)) +
  geom_point(alpha = 0.9) +
  scale_y_log10() +
  scale_x_log10() +
  geom_text(aes(label = uniprot_acc))


weirdos <- mean_cv_df %>% 
  filter(cv_val < 5 & mean_val < 1)

accession_numbers <- weirdos$uniprot_acc %>% unique()

cv_schmidt %>% filter(`Uniprot Accession` %in% accession_numbers) %>% as.data.frame() %>% select(Description)
