#' pre-process the ACIC winners


# Setup -------------------------------------------------------------------

library(tidyverse)

# Load/combine ------------------------------------------------------------

acic_winners <- readRDS("Data/eval_simple_tosend.RDS")

acic_winners.df <- names(acic_winners) %>% 
  #keep(~ grepl("^t2", .x)) %>% 
  map(~ acic_winners[[.x]]$eval %>% mutate(method = .x)) %>%
  list_rbind() %>% 
  filter(!(method %in% c('t1_rBARTimpT', 'ofpsbart', 't1_H-TMLE')))
# select(-DGP_name, -id.practice) %>% 
# mutate(noise = FALSE)

saveRDS(acic_winners.df, "Save/acic_winners_performance.rds")

# DGPs --------------------------------------------------------------------

estimand_truths <- read_csv("Data/ACIC_estimand_truths.csv")

dgp_levels <- full_join(
  estimand_truths %>% 
    select(dataset.num:`Idiosyncrasy of Impacts`) %>% 
    distinct(),
  acic_winners.df %>% 
    distinct(DGP_name, dataset.num)
)
  
dgp_levels %>% 
  #filter(dataset.num %in% dataset_nums) %>% 
  distinct(pick(-dataset.num))
