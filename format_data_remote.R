#  A script to use on remote nodes to generate noise and format data
# Adelin Kassler - Nov 2024


# Setup -------------------------------------------------------------------

if (interactive()) {
  dataset_num <- '0001'
} else {
  dataset_num <- commandArgs(TRUE)[1]
}

path_p <- "Data/track2_20220404/practice/acic_practice_%s.csv"
path_py <- "Data/track2_20220404/practice_year/acic_practice_year_%s.csv"
path_fmt <- "Data/formatted/practice_%s.csv"
path_noise <- "Data/formatted/practice_%s_noise_%s.csv"

library(tidyverse)
library(BART)
library(caret)


# Load data ---------------------------------------------------------------

data_py.raw <- read.csv(sprintf(path_py, dataset_num), stringsAsFactors = TRUE)
data_p.raw <- read.csv(sprintf(path_p, dataset_num), stringsAsFactors = TRUE)


# Format data -------------------------------------------------------------

data_py.wide <- data_py.raw %>%
  pivot_longer(!c(id.practice, year, Z)) %>%
  group_by(id.practice, name) %>%
  arrange(id.practice, name, year) %>% 
  mutate(diff1 = value - lag(value),
         diff2 = value - lag(lag(value)),
         diff3 = value - lag(lag(lag(value)))) %>% 
  ungroup() %>%
  pivot_wider(id_cols = c(id.practice, Z),
              names_from = c(name, year),
              values_from = c(value, diff1, diff2, diff3),
              names_glue = "{name}{ifelse(.value=='value','','_diff')}.y{year}{ifelse(.value=='value','',paste0('_y',year-parse_number(.value,na='value')))}") %>% 
  select(!starts_with("post") & 
           !ends_with(c('0', '-1', '-2')) & 
           !starts_with(c("Y_diff.y3", "Y_diff.y4"))) %>% 
  # TODO: outcome variable should be Y_avg.post - Y_avg.pre
  mutate(Y_avg.post = (Y.y3 * n.patients.y3 + Y.y4 * n.patients.y4) /
           (n.patients.y3 + n.patients.y4)) %>%
  select(-Y.y4, -Y.y3)

data.main <- inner_join(data_py.wide, data_p.raw, by = "id.practice")
class(data.main) <- "data.frame"

rm(data_py.raw, data_p.raw, data_py.wide)


# Add noise variables -----------------------------------------------------

# todo: make sure noise generation is at the same level (practice, p-y) as curia is doing

# Adding variables is done within a python script provided by Curia

write_csv(data.main, file = sprintf(path_fmt, dataset_num))
# Call python code
# TODO: Incorporate script into process, get random seed.
# Output files as path_noise, with %s dataset num and %s noise type


# Load python-generated noise ---------------------------------------------

# Fit BART/DART to data+noise ---------------------------------------------

# Call fit functions on all data sets, sparse and not
# Calculate partial dependence functions
# Calculate pdep averages for all subgroups? (could also do that locally)
# Bring down to local and aggregate

# Performance evals in another script/notebook, combined with curia evals &c.