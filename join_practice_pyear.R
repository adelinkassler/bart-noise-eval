# Join practice, practice-year datasets

library(glue)


# Join practice/practice-year ---------------------------------------------

for (i in 1:3400) {
  dsnum <- str_pad(i,4,'left','0')
  
  acic_practice_path = glue("~/Documents/Consulting/BCBS/Data/track2_20220404/practice/acic_practice_{dsnum}.csv")
  acic_practice_year_path = glue("~/Documents/Consulting/BCBS/Data/track2_20220404/practice_year/acic_practice_year_{dsnum}.csv")

  practice <- read_csv(acic_practice_path) %>% 
    strings_to_factors()
  practice_year <- read_csv(acic_practice_year_path) %>% 
    strings_to_factors()
  
  t2_joined <- full_join(practice, practice_year)
  saveRDS(t2_joined, glue("Data/t2_joined/acic_ppy_{dsnum}.rds"))
}


# Add born/died -----------------------------------------------------------

for (i in 1:3400) {
  dsnum <- str_pad(i,4,'left','0')
  
  t2_joined <- readRDS(glue("Data/t2_joined/acic_ppy_{dsnum}.rds"))
  t2_joined_jd <- t2_joined %>% 
    mutate(joined = NA, died = NA)
  saveRDS(t2_joined_jd, glue("Data/t2_joined/acic_ppy_{dsnum}_jd.rds"))
}
