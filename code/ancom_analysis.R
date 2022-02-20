library(readr)
library(tidyverse)
library(tidyr)
library(compositions)


## Read in metadata and table
metadata <- read.delim("data/metadata.txt")
metadata$Time_Point <- as.factor(metadata$Time_Point)

table <- read_tsv("data/feature_table.1488.txt", col_names = TRUE) %>% 
  rownames_to_column(var = "ASV") %>% 
  as.data.frame()
rownames(table) <- paste("ASV", table$ASV, sep = "")

# Generate Taxonomy file
taxonomy <- table %>% dplyr::select(featureid) %>% 
  separate(col = "featureid", 
           into = c("Kingdom", "Phylum", "Class", 
                    "Order", "Family", "Genus"), 
           sep = ";",
           fill = "right") %>% 
  rownames_to_column(var = "ASV")

## Table with Taxa
table.tax <- table %>% dplyr::select(-featureid, -ASV) %>% 
  rownames_to_column(var = "ASV") %>% 
  pivot_longer(-ASV) %>% 
  inner_join(., taxonomy, by = "ASV") %>% 
  pivot_longer(cols = c("Kingdom", "Phylum", "Class", 
                        "Order", "Family", "Genus"),
               names_to = "level",
               values_to = "taxon") %>% 
  as.data.frame()

# Family level abundance
family.table <- table.tax %>% 
  filter(level == "Family") %>% 
  as.data.frame() %>% 
  group_by(name, taxon) %>% 
  summarise(family_count = sum(value)) %>% 
  ungroup() %>% 
  mutate(taxon = gsub("D_4__", "", taxon))

#######################
# Healthy vs Chronic
## filter down metadata to chronic and healthy
metadata.filt.chronic <- metadata %>% 
  dplyr::select(sampleid, Time_Point, Sample) %>% 
  filter(Sample == "BAL") %>% 
  filter(Time_Point == 0 | 
           Time_Point == 36) %>% 
  as.data.frame() %>% 
  mutate(Time_Point = case_when(Time_Point == 0 ~ "Healthy",
                                Time_Point == 36 ~ "Chronic"))


## filter down table and rename rows to family ID
table.family.filt.chronic <- family.table %>% 
  filter(name %in% metadata.filt.chronic$sampleid) %>% 
  arrange(desc(name)) %>% 
  pivot_wider(names_from = name,
              values_from = family_count)
table.family.filt.chronic

# Removes any NA row names and sets to "Undefined Family"
table.family.filt.chronic[is.na(table.family.filt.chronic)]<-"Undefined Family"
# Move rowname to col 1
table.family.filt.chronic <-column_to_rownames(table.family.filt.chronic, var = "taxon")
## psuedocount 1 to every cell in table
table.family.filt.chronic.pseudo <- table.family.filt.chronic + 1


## Folowing ANCOM Standard Analysis
# https://github.com/FrederickHuangLin/ANCOM#standard-analysis 
otu_data.chronic = table.family.filt.chronic.pseudo

meta_data.chronic = metadata.filt.chronic
meta_data.chronic = meta_data.chronic %>% rename(Sample.ID = sampleid)

source("code/ANCOM-master/scripts/ancom_v2.1.R")

# Step 1: Data preprocessing

feature_table = otu_data.chronic; sample_var = "Sample.ID"; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data.chronic, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "Time_Point"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL

res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)

res$out <- res$out %>% rename(family = taxa_id)

write_csv(res$out, "data/ANCOM/h_vs_c.csv")

### Combine ANCOM and ALDEx2 data
aldex2_output <- read_csv("data/ALDEx2/h_vs_c_ALDEx2_output.csv")
aldex2_ancom_h_v_c <- left_join(res$out, aldex2_output, by = "family")

# Write combined output
write_csv(aldex2_ancom_h_v_c, "data/aldex2_ancom_h_v_c_output.csv")
