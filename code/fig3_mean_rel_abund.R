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
  select(-Count)

table <- orig.table %>% 
  rownames_to_column("ASV")

taxonomy <- table %>% select(ASV, featureid) %>% 
  separate(featureid,
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
           sep = ";",
           fill = "right")
taxonomy

table.counts <- table %>% 
  select(ASV, metadata.filtered$sampleid) %>% 
  pivot_longer(-c(ASV),
               names_to = "sampleid",
               values_to = "Count")

table.rel.abund <- inner_join(metadata.filtered, table.counts, by = "sampleid") %>% 
  inner_join(., taxonomy, by = "ASV") %>% 
  group_by(sampleid) %>% 
  mutate(rel_abund = Count / sum(Count)) %>% 
  ungroup() %>% 
  pivot_longer(cols = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
               names_to = "level",
               values_to = "taxon")

rel_abund_plot_family <- function(sample, other_level) {
  bal.table.mean <- table.rel.abund %>% filter(Sample == sample) %>% 
    filter(level == "Family",) %>%
   mutate(taxon = gsub("D_4__", "", taxon),
           taxon = str_replace(taxon,
                              "^(\\S*)$","__*\\1*__")) %>% 
    group_by(Time_Point, sampleid, taxon) %>% 
    summarise(rel_abund = sum(rel_abund), .groups = "drop") %>% 
    group_by(Time_Point, taxon) %>% 
    summarise(mean_rel_abund = 100*mean(rel_abund), .groups = "drop")

  taxon.pool <- bal.table.mean %>% 
    group_by(taxon) %>% 
    summarise(pool = max(mean_rel_abund) < other_level, .groups = "drop")

  bal.table <- inner_join(bal.table.mean, taxon.pool, by = "taxon") %>% 
    mutate(taxon = if_else(pool, "Other", taxon)) %>% 
    group_by(Time_Point, taxon) %>% 
    summarise(mean_rel_abund = sum(mean_rel_abund)) %>% 
    mutate(order = case_when(Time_Point == 0 ~ 1,
                             Time_Point == 6 ~ 2,
                             Time_Point == 36 ~ 3))

  time = bal.table$order
  value = bal.table$mean_rel_abund
  taxa = bal.table$taxon

  data.frame(time, value, taxa) %>% 
    ggplot(., aes(x=time, y = value, fill = taxa)) +
    geom_area() +

   scale_y_continuous(expand = c(0,0))+
   scale_fill_manual(name = NULL,
                      values = c(brewer.pal(11, "Set3"),
                                 "#808080")) +
  
    labs(x = NULL,
         y = "Mean Relative Abundance (%)") +
    scale_x_continuous(breaks = 1:3,
                      labels = paste0(c("Healthy",
                                        "Acute",
                                        "Chronic"))) +
    theme_cowplot() +
    theme(axis.text = element_text(face = "bold"),
         axis.title = element_text(face = "bold"),
        
         legend.title = element_blank(),
         legend.text = element_markdown(face = "bold")) 
}

rel_abund_plot_family("BAL", 1)


# ggsave("plots/bal_rel_abund_family.png",
#        dpi = 600,
#        width = 6,
#        height = 4,
#        units = c("in"),
#        bg = "white")

################## Phylum

bal.table.mean <- table.rel.abund %>% filter(Sample == "BAL") %>% 
  filter(level == "Phylum",) %>%
  mutate(taxon = gsub("D_1__", "", taxon),
         taxon = str_replace(taxon,
                             "^(\\S*)$","__*\\1*__")) %>% 
  group_by(Time_Point, sampleid, taxon) %>% 
  summarise(rel_abund = sum(rel_abund), .groups = "drop") %>% 
  group_by(Time_Point, taxon) %>% 
  summarise(mean_rel_abund = 100*mean(rel_abund), .groups = "drop")

taxon.pool <- bal.table.mean %>% 
  group_by(taxon) %>% 
  summarise(pool = max(mean_rel_abund) < 1, .groups = "drop")

bal.table <- inner_join(bal.table.mean, taxon.pool, by = "taxon") %>% 
  mutate(taxon = if_else(pool, "Other", taxon)) %>% 
  group_by(Time_Point, taxon) %>% 
  summarise(mean_rel_abund = sum(mean_rel_abund)) %>% 
  mutate(order = case_when(Time_Point == 0 ~ 1,
                           Time_Point == 6 ~ 2,
                           Time_Point == 36 ~ 3))

time = bal.table$order
value = bal.table$mean_rel_abund
taxa = bal.table$taxon

data.frame(time, value, taxa) %>% 
  ggplot(., aes(x=time, y = value, fill = taxa)) +
  geom_area() +
  
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(name = NULL,
                    values = c(brewer.pal(11, "Set3"),
                               "#808080")) +
  
  labs(x = NULL,
       y = "Mean Relative Abundance (%)") +
  scale_x_continuous(breaks = 1:3,
                     labels = paste0(c("Healthy",
                                       "Acute",
                                       "Chronic"))) +
  theme_cowplot() +
  theme(axis.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        
        legend.title = element_blank(),
        legend.text = element_markdown(face = "bold")) 

ggsave("plots/bal_rel_abund_phylum.png",
       dpi = 600,
       width = 6,
       height = 4,
       units = c("in"),
       bg = "white")
