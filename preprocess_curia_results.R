#' pre-process curia estimates we were given into a form that's useful

# Setup -------------------------------------------------------------------

library(tidyverse)


# Background data ---------------------------------------------------------

dataset_nums <- read_csv("Data/curia_data/dataset_nums.csv") %>% 
  pull(dataset_num)


# Load and combine --------------------------------------------------------

curia_estimates_1x <- dir("Data/curia_results/1x/", full.names = TRUE) %>%
  # grep("acic_[a-z]+_([0-9]{4})\\.csv", ., value=T) %>%
  map(~ read_csv(.x)) %>%
  list_rbind() %>%
  mutate(dataset.num = sprintf("%04d", dataset.num)) %>% 
  mutate(noise = FALSE)

curia_estimates_10x <- dir("Data/curia_results/satt_results_10x/", full.names = TRUE) %>%
  grep("acic_[a-z]+_([0-9]{4}b)\\.csv", ., value=T) %>% 
  map(~ read_csv(.x)) %>%
  list_rbind() %>% 
  mutate(dataset.num = sprintf("%04d", dataset.num)) %>% 
  mutate(noise = TRUE)

curia_estimates <- bind_rows(curia_estimates_1x, curia_estimates_10x) %>% 
  rename(satt = model_satt, se = model_se)

saveRDS(curia_estimates, "Save/all_curia_estimates_1x_10x.rds")

