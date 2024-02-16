library(tidyverse)
library(glue)
source('acic_bcf_from_d.r')

impacts1 <- run_bcf(dataset=1)
impacts2 <- run_bcf(dataset=2)
impacts3 <- run_bcf(dataset=3)
