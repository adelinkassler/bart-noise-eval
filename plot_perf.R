library(tidyverse)

perf <- read_csv('eval_compare.csv') %>%
  filter(!(method %in% c('t1_rBARTimpT', 'ofpsbart', 't1_H-TMLE'))) %>%
  filter(variable=='Overall' & is.na(year))

perf <- rbind(perf, rep(NA, ncol(perf)))
perf$method[nrow(perf)] <- 'Oracle'
perf$rmse[nrow(perf)] <- 0
perf$width[nrow(perf)] <- 0
perf$cover[nrow(perf)] <- .9

perf <- perf %>%
  mutate(Method = ifelse(method == 'curia', 'Curia', 
                         ifelse(method == 'Oracle', 'Oracle', 'Other')))

method_colors <- c(Curia = "#F8766D", Oracle = "#00BFC4", Other = "darkgrey")

meta <- readRDS('real_submissions.RDS')

curia_rmse <- perf %>%
  filter(Method=='Curia') %>%
  pull(rmse)

perf %>%
  summarise(percent_better = mean(rmse < curia_rmse, na.rm=T))

perf %>%
  ggplot(aes(x=rmse, y=cover, color=Method)) + 
  geom_hline(yintercept = 0.9, color='darkgrey') +
  geom_point() + 
  xlab('RMSE, overall SATT') +
  ylab('Coverage of overall SATT 90% interval') +
  scale_color_manual(values=method_colors)
ggsave('rmse_cover.png', height=3.5, width=4.5)

perf %>%
  ggplot(aes(x=width, y=cover, color=Method)) + 
  geom_hline(yintercept = 0.9, color='darkgrey')+
  geom_point() + 
  xlab('Width of overall SATT 90% interval') +
  ylab('Coverage of overall SATT 90% interval') +
  scale_color_manual(values=method_colors)
ggsave('width_cover.png', height=3.5, width=4.5)
