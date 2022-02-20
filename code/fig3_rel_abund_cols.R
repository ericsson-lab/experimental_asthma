library(ggtext)
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(tidyr)
library(stringr)

## get metadata
metadata <- read_delim("data/metadata.txt")

## get feature table and assign taxonomy as rownames
orig.table <- read_delim("data/feature_table.1488.txt", col_types = NULL) %>% 
  as.data.frame()
head(orig.table)

## Filter down metadata to points of interest
metadata.filtered <- metadata %>% 
  filter(Time_Point == 0 |
           Time_Point == 6 |
           Time_Point == 36) %>% 
  select(-Count) %>% 
  filter(Sample == "BAL") %>% 
  # Add "stage" column for facet wrap in ggplot
  mutate(stage = case_when(Time_Point == 0 ~ "Healthy",
                          Time_Point == 6 ~ "Acute",
                          Time_Point == 36 ~ "Chronic"))

# Add ASV column to parse taxonomy out
table <- orig.table %>% 
  rownames_to_column("ASV")

# Generate taxonomy and separate identifiers
taxonomy <- table %>% select(ASV, featureid) %>% 
  separate(featureid,
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
           sep = ";",
           fill = "right")

# Make feature table longer
table.counts <- table %>% 
  select(ASV, metadata.filtered$sampleid) %>% 
  pivot_longer(-c(ASV),
               names_to = "sampleid",
               values_to = "Count")

# Join metadata to feature table and taxonomy
table.rel.abund <- inner_join(metadata.filtered, 
                              table.counts, 
                              by = "sampleid") %>% 
  inner_join(., taxonomy, by = "ASV") %>% 
  group_by(sampleid) %>% 
  # calculate % relative abundance
  mutate(rel_abund = 100*Count / sum(Count)) %>% 
  ungroup() %>% 
  # Break up taxonomy
  pivot_longer(cols = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
               names_to = "level",
               values_to = "taxon") %>% 
  # Select only family level
  filter(level == "Family",) %>%
  # Remove string header from silva/gg(?)
  mutate(taxon = gsub("D_4__", "", taxon),
         taxon = str_replace(taxon,
                             "^(\\S*)$","__*\\1*__")) %>% 
  ungroup() %>% 
  group_by(sampleid, taxon) %>% 
  # sum up mean rel abund for each taxa within a sample
  summarise(rel_abund = sum(rel_abund), .groups = "drop") %>% 
  inner_join(., metadata.filtered, by = "sampleid")
  
# Generate condition for Other pool.
# Logical condition evaluates if the max mean relabund in any sample
# in any sample is < 1. If yes, taxa is grouped into other.
taxon.pool <- table.rel.abund %>% 
  group_by(taxon) %>% 
  summarise(pool = max(mean(rel_abund)) < 1, .groups = "drop")

# Count # of max(mean(rel_abund)) taxa > 1% to set color pallete in plot
taxon.pool %>% count(pool)
  
# Join pool table to rel.abund and rename family to other is mean(max(rel_abund)) < 1%
bal.table <- inner_join(table.rel.abund, taxon.pool, by = "taxon") %>% 
    mutate(taxon = if_else(pool, "Other", taxon)) 

# Clean up data for plotting
time = bal.table$stage
value = bal.table$rel_abund
taxa = bal.table$taxon
sampleid = bal.table$`Cat `

# Plot mean relabund col
data.frame(sampleid, time, value, taxa) %>% 
    group_by(sampleid) %>% 
    ggplot(., aes(x=sampleid, y = value, fill = taxa)) +
    geom_col() +
    facet_wrap(factor(time, levels = c("Healthy", "Acute", "Chronic")), 
               strip.position = "bottom") +
    scale_y_continuous(expand = c(0,0)) +
  
    scale_fill_manual(name = NULL,
                      values = c(brewer.pal(9, "Set3"),
                                 "#808080")) +
    
    labs(x = NULL,
         y = "Mean Relative Abundance (%)") +
    scale_x_continuous(breaks = 1:3,
                       labels = paste0(c("Healthy",
                                         "Acute",
                                         "Chronic"))) +
    theme_cowplot() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          
          axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "bold"),
          
          legend.title = element_blank(),
          legend.text = element_markdown(face = "bold"),
          
          strip.background = NULL,
          strip.placement = "outside",
          strip.text = element_text(face = "bold")) +
    guides(position = guide_legend(override.aes = list(size = 4.5,
                                                    shape = c(15,16,17))))

# And save...
ggsave("plots/bal_rel_abund_col_family.png",
       dpi = 600,
       width = 7,
       height = 4,
       units = c("in"),
       bg = "white")

#


########################################

# Join metadata to feature table and taxonomy
table.rel.abund <- inner_join(metadata.filtered, 
                              table.counts, 
                              by = "sampleid") %>% 
  inner_join(., taxonomy, by = "ASV") %>% 
  group_by(sampleid) %>% 
  # calculate % relative abundance
  mutate(rel_abund = 100*Count / sum(Count)) %>% 
  ungroup() %>% 
  # Break up taxonomy
  pivot_longer(cols = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
               names_to = "level",
               values_to = "taxon") %>% 
  # Select only family level
  filter(level == "Phylum",) %>%
  # Remove string header from silva/gg(?)
  mutate(taxon = gsub("D_1__", "", taxon),
         taxon = str_replace(taxon,
                             "^(\\S*)$","__*\\1*__")) %>% 
  ungroup() %>% 
  group_by(sampleid, taxon) %>% 
  # sum up mean rel abund for each taxa within a sample
  summarise(rel_abund = sum(rel_abund), .groups = "drop") %>% 
  inner_join(., metadata.filtered, by = "sampleid")

# Generate condition for Other pool.
# Logical condition evaluates if the max mean relabund in any sample
# in any sample is < 1. If yes, taxa is grouped into other.
taxon.pool <- table.rel.abund %>% 
  group_by(taxon) %>% 
  summarise(pool = max(mean(rel_abund)) < 1, .groups = "drop")

# Count # of max(mean(rel_abund)) taxa > 1% to set color pallete in plot
taxon.pool %>% count(pool)

# Join pool table to rel.abund and rename family to other is mean(max(rel_abund)) < 1%
bal.table <- inner_join(table.rel.abund, taxon.pool, by = "taxon") %>% 
  mutate(taxon = if_else(pool, "Other", taxon)) 

# Clean up data for plotting
time = bal.table$stage
value = bal.table$rel_abund
taxa = bal.table$taxon
sampleid = bal.table$`Cat `

# Plot mean relabund col
data.frame(sampleid, time, value, taxa) %>% 
  group_by(sampleid) %>% 
  ggplot(., aes(x=sampleid, y = value, fill = taxa)) +
  geom_col() +
  facet_wrap(factor(time, levels = c("Healthy", "Acute", "Chronic")), 
             strip.position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  
  scale_fill_manual(name = NULL,
                    values = c(brewer.pal(9, "Set3"),
                               "#808080")) +
  
  labs(x = NULL,
       y = "Mean Relative Abundance (%)") +
  scale_x_continuous(breaks = 1:3,
                     labels = paste0(c("Healthy",
                                       "Acute",
                                       "Chronic"))) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        
        legend.title = element_blank(),
        legend.text = element_markdown(face = "bold"),
        
        strip.background = NULL,
        strip.placement = "outside",
        strip.text = element_text(face = "bold")) +
  guides(position = guide_legend(override.aes = list(size = 4.5,
                                                     shape = c(15,16,17))))

#And save...
ggsave("plots/bal_rel_abund_col_phylum.png",
       dpi = 600,
       width = 7,
       height = 4,
       units = c("in"),
       bg = "white")

#



