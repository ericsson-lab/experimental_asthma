library(tidyverse)
library(dplyr)

metadata <- read_delim("data/metadata.txt")
metadata %>% count(Sample)
metadata.filt <- metadata %>% filter(Time_Point == 0 |
                      Time_Point == 6 |
                      Time_Point == 36) %>% 
  mutate(Time_Point = case_when(Time_Point == 0 ~ "Health", 
                                Time_Point == 6 ~ "Acute",
                                Time_Point == 36 ~ "Chronic"),
         Sample = case_when(Sample == "BAL" ~ "BALF",
                            Sample == "OP" ~ "OP",
                            Sample == "Fecal" ~ "Rectal"))
write_tsv(metadata.filt, "data/metadata_filt.tsv")
