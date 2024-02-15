#' This script is used for pre-processing data for analysis with BART and DART.
#' Author: Adelin Kassler, 2024/01/20


# Setup -------------------------------------------------------------------

.arguments <- commandArgs(TRUE)

.version.params <- list(
  name = "default",
  runtime = Sys.time(),
  notes = NULL,
  sdy = TRUE,
  pscores = TRUE,
  weight_cols = c('n.patients.y4', 'n.patients.y3'),
  bart = list(
    nskip = 200,
    keepevery = 1
  ),
  dart = list(
    nskip = 2000,
    keepevery = 5
  )
)

suppressPackageStartupMessages({
  # QoL
  library(yaml)
  library(testthat)
  library(progress)
  library(spsUtils)
  
  # Computational tools
  library(BART)
  library(caret)
  
  # General use
  library(magrittr)
  library(glue)
  library(tidyverse)
})

options(readr.show_col_types = FALSE)

# Convenience functions ---------------------------------------------------

my_timer <- function (expr, gcFirst = TRUE) {
  ppt <- function(y) {
    if (!is.na(y[4L])) 
      y[1L] <- y[1L] + y[4L]
    if (!is.na(y[5L])) 
      y[2L] <- y[2L] + y[5L]
    paste(formatC(y[1L:3L]), collapse = " ")
  }
  if (gcFirst) 
    gc(FALSE)
  time <- Sys.time()
  on.exit(message("Timing stopped at: ", ppt(Sys.time() - 
                                               time)))
  expr
  new.time <- Sys.time()
  on.exit()
  structure(new.time - time, class = "difftime")
}

parse_dataset_num <- function(path) {
  str_match(basename(path), '[0-9]{4}')[1,1]
}

# Processing curia data ---------------------------------------------------

dedup_cols <- function(x) {
  dup_cols <- duplicated(t(x))
  x[!dup_cols]
}

get_practice_level_vars <- function(x) {
  # expect_in(c('id.practice', 'Z', 'year'), names(x))
  
  x %>% 
    group_by(id.practice) %>% 
    summarise(
      n = n(),
      across(!n, n_distinct)
    ) %>% 
    select(-id.practice, -Z) %>% 
    {as.data.frame(t(rbind(
      mean = summarise(., across(everything(), mean)) %>% round(1),
      sd = summarise(., across(everything(), sd)) %>% round(3)
    )))} %>% 
    filter(mean == 1, sd == 0) %>% 
    rownames()
}

select_practice_level <- function(x) {
  expect_in(c('id.practice', 'Z'), names(x))
  practice_level_vars <- get_practice_level_vars(x)
  select(x, id.practice, Z, all_of(practice_level_vars))
}

select_practice_year_level <- function(x) {
  expect_in(c('id.practice', 'Z'), names(x))
  practice_level_vars <- get_practice_level_vars(x)
  select(x, id.practice, Z, !all_of(practice_level_vars))
}



# Format data -------------------------------------------------------------

strings_to_factors <- function(x) {
  mutate(x, across(where(is.character), as.factor))
}

add_sd_y_practice <- function(x, dsnum) {
  sdy <- readRDS(glue("~/Documents/Consulting/BCBS/Data/formatted/track1_sdy/acic_practice_{dsnum}_sdy.rds"))
  left_join(x, sdy, by = 'id.practice')
}

format_practice_year_wide <- function(practice_year) {
  expect_in(c('id.practice', 'Z', 'year', 'Y', 'n.patients'), 
            names(practice_year))
  
  practice_year %>%
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
             (n.patients.y3 + n.patients.y4),
           Y_avg.pre = (Y.y1 * n.patients.y1 + Y.y2 * n.patients.y2) /
             (n.patients.y1 + n.patients.y2),
           Y_avg_diff.post_pre = Y_avg.post - Y_avg.pre) %>%
    select(-Y.y4, -Y.y3)
}

add_pscores <- function(x, method = 'BART') {
  data <- as_rhs(x)
  rhs_names <- setdiff(names(data), 'Z')
  # lbart and pbart don't accept weights
  # weights <- rowSums(x[.version.params$weight_cols])
  # weights %<>% {. / sum(.)}
  
  # Extensible to other methods
  if (method == 'BART') {
    fit <- lbart(data[rhs_names], data$Z, sparse = F, nskip = .version.params$bart$nskip,
                 keepevery = .version.params$bart$keepevery)
  } else if (method == 'DART') {
    fit <- lbart(data[rhs_names], data$Z, sparse = T, nskip = .version.params$dart$nskip, 
                 keepevery = .version.params$dart$keepevery)
  } else if (method == 'linear') {
    .fmla = as.formula(glue("Z ~ ", glue_collapse(rhs_names, sep = ' + ')))
    fit <- lm(data = data, formula = .fmla)
  }
  
  x %>% 
    mutate(pscore = colMeans(fit$yhat.train))
}

join_practice_level <- function(practice, practice_year) {
  # x <- inner_join(practice, practice_year, by = c('id.practice', 'Z'))
  x <- inner_join(practice, practice_year)
  class(x) <- 'data.frame'
  #TODO: figure out what this join is producing 4x repeated entries
  distinct(x)
}

get_rhs_names <- function(x) {
  setdiff(names(x), c("Y_outcome", 'id.practice', 'dataset.num', 'Y_avg.post', 
                      'Y_avg.pre', 'Y_avg_diff.post_pre'))
}

as_rhs <- function(x) {
  if (!(class(x)[1] %in% c('data.frame', 'matrix'))) 
    class(x) <- 'data.frame'
  
  x[get_rhs_names(x)] %>% 
    strings_to_factors()
}

# Data loading and pre-processing -----------------------------------------

load_clean_and_split_curia_dataset <- function(path) {
  merged_data <- read_csv(path) %>% 
    dedup_cols() %>% 
    strings_to_factors() %>% 
    rename(Y = label, n.patients = n_patients_in_practice) %>% 
    select(-mo)
  
  list(
    practice = select_practice_level(merged_data),
    practice_year = select_practice_year_level(merged_data)
  )
}

preproc_data <-
  function(acic_practice_path = NULL,
           acic_practice_year_path = NULL,
           curia_merged_path = NULL,
           sdy = .version.params$sdy,
           pscores = .version.params$pscores
  ) {
    if (is.null(acic_practice_path) & is.null(acic_practice_year_path) &
        !is.null(curia_merged_path)) {
      dsname <- "curia_10x"
      dsnum <- parse_dataset_num(curia_merged_path)
      split <- load_clean_and_split_curia_dataset(curia_merged_path)
      
      practice <- split$practice
      practice_year <- split$practice_year
      rm(split)
    } else if (!is.null(acic_practice_path) & !is.null(acic_practice_year_path) &
               is.null(curia_merged_path)) {
      dsname <- "acic_base"
      dsnum <- parse_dataset_num(acic_practice_path)
      expect_equal(dsnum, parse_dataset_num(acic_practice_year_path))
      
      practice <- read_csv(acic_practice_path) %>% 
        strings_to_factors()
      practice_year <- read_csv(acic_practice_year_path) %>% 
        strings_to_factors()
    } else {
      stop(glue("Expecting either acic practice and practice-year paths, OR curia merged paths.",
                " Got the following arguments: {as.list(match.call()[-1])}"))
    }
    
    data <- join_practice_level(
      practice_year = practice_year %>% format_practice_year_wide(),
      practice = practice
    )
    if (sdy) data <- data %>% add_sd_y_practice(dsnum)
    if (pscores) data <- data %>% add_pscores()
    attr(data, "dataset.num") <- dsnum
    
    saveRDS(data, glue("{dir_path}/patient_wide_preproc_{dsnum}_{dsname}.rds"))
    # write_csv(data, glue("{dir_path}/patient_wide_preproc_{dsnum}_{dsname}.csv"))
    return(data)
}


# Run main ----------------------------------------------------------------

dir_path <- glue("~/Documents/Consulting/BCBS/Data/model_ready/{.version.params$name}")
if (!dir.exists(dir_path)) dir.create(dir_path)
write_yaml(.version.params, file.path(dir_path, "PARAMETERS.yml"))

focus_datasets <- read_csv("Data/curia_data/dataset_nums.csv", show_col_types = FALSE)

acic_dsnums <- focus_datasets$dataset_num[focus_datasets$p == '1x']
pb <- progress_bar$new(
  format = ":current/:total in :elapsed [:bar] eta: :eta\n",
  total = length(acic_dsnums)
)
pb$tick(0)
iwalk(acic_dsnums, function(dsnum, i) {
  quiet(preproc_data(
    acic_practice_path = glue("~/Documents/Consulting/BCBS/Data/track2_20220404/practice/acic_practice_{dsnum}.csv"),
    acic_practice_year_path = glue("~/Documents/Consulting/BCBS/Data/track2_20220404/practice_year/acic_practice_year_{dsnum}.csv")
  ))
  pb$tick()
  # message(glue("({i}/{length(acic_dsnums)}) Processed ACIC {dsnum} in {format(elapsed)}."))
})

curia_10x_dsnums <- focus_datasets$dataset_num[focus_datasets$p == '10x']
pb <- progress_bar$new(
  format = ":current/:total in :elapsed [:bar] eta: :eta\n",
  total = length(curia_10x_dsnums)
)
pb$tick(0)
iwalk(curia_10x_dsnums, function(dsnum, i) {
  quiet(preproc_data(
    curia_merged_path = glue("~/Documents/Consulting/BCBS/Data/curia_data/df_merged_raw_practice_level_10x/merged_dataset_practice_{dsnum}b.csv")
  ))
  pb$tick()
  # message(glue("({i}/{length(curia_10x_dsnums)}) Processed Curia 10x {dsnum} in {format(elapsed)}."))
})
