library(ggtext)
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(tidyr)
library(stringr)
library(forcats)
library(ggthemes)
library(scales)




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


## Plotting Group by phylum


#rel_abund_plot_family <- function(sample, other_level) {

  bal.table.mean <- table.rel.abund %>% filter(Sample == "BAL") %>% 
    filter(level == "Family",) %>%
   mutate(taxon = gsub("D_4__", "", taxon)) %>% 
           # taxon = str_replace(taxon,
           #                    "^(\\S*)$","__*\\1*__")) %>% 
    group_by(Time_Point, sampleid, taxon) %>% 
    summarise(rel_abund = sum(rel_abund), .groups = "drop") %>% 
    group_by(Time_Point, taxon) %>% 
    summarise(mean_rel_abund = 100*mean(rel_abund), .groups = "drop")
unique(bal.table.mean$taxon)
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
  taxon = bal.table$taxon

df <- data.frame(time, value, taxon)



# Gather all of the family to be plotted and group by phylum
#https://stackoverflow.com/questions/62627480/how-to-creat-a-bar-graph-of-microbiota-data-with-one-color-for-higher-taxonomic
phylums <- c('Actinobacteria','Bacteroidetes','Firmicutes', "Proteobacteria")

proteobacteria <- c("Beijerinckiaceae", "Pseudomonadaceae",
                    "Xanthobacteraceae", "Burkholderiaceae",
                    "Moraxellaceae", "Xanthomonadaceae",
                    "Sphingomonadaceae")
bacteroidetes <- c("Sphingobacteriaceae","Chitinophagaceae")  
actinobacteria <- c("Propionibacteriaceae")
firmicutes <- c("Lachnospiraceae")

df <- df %>% 
  mutate(family = gsub("__\\*",  "", taxon),
         family = gsub("\\*__", "", family),
         phylum = case_when(family %in% proteobacteria ~ "Proteobacteria",
                            family %in% bacteroidetes ~ "Bacteroidetes",
                            family %in% actinobacteria ~ "Actinobacteria",
                            family %in% firmicutes ~ "Firmicutes", 
                            family == "Other" ~ "Other")) %>% 
  mutate(phylum=factor(phylum, levels=c(phylums, "Other")),
         family=fct_reorder(family, 10*as.integer(phylum) + grepl("Other", family)))


palette <- c("#3366cc","#dc3912" ,"#0099c6", 
         "#109618", "#990099", "#ff9900" ,
         "#dd4477", "#66aa00" ,"#b82e2e" ,
         "#316395", "#329262", "#808080")



  

df %>% 
    ggplot(., aes(x=time, y = ,value, fill = family)) +
    geom_area(aes(alpha = 0.5)) +
  scale_alpha(guide = "none")+

   scale_y_continuous(expand = c(0,0))+
   scale_fill_manual(name = NULL,
                    values = alpha(palette, 0.7)) +
  
    labs(x = NULL,
         y = "Mean Relative Abundance (%)") +
    scale_x_continuous(breaks = 1:3,
                      labels = paste0(c("Health",
                                        "Acute",
                                        "Chronic"))) +
    theme_cowplot() +
    theme(axis.text = element_text(face = "bold"),
         axis.title = element_text(face = "bold"),
        
         legend.title = element_blank(),
         legend.text = element_markdown(face = "bold"))



ggsave("plots/bal_rel_abund_family.png",
       dpi = 600,
       width = 6,
       height = 4,
       units = c("in"),
       bg = "white")

################## Phylum

bal.table.mean <- table.rel.abund %>% filter(Sample == "BAL") %>% 
  filter(level == "Phylum",) %>%
  mutate(taxon = gsub("D_1__", "", taxon)) %>% 
         # taxon = str_replace(taxon,
         #                     "^(\\S*)$","__*\\1*__")) %>% 
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
taxon = bal.table$taxon


plot_color <- c("#ff9900","#3366cc","#dc3912", "#66aa00", "#808080")

data.frame(time, value, taxon) %>% 
  ggplot(., aes(x=time, y = value, fill = factor(taxon, levels = c("Actinobacteria",
                                                                   "Bacteroidetes",
                                                                   "Firmicutes",
                                                                   "Proteobacteria",
                                                                   "Other")))) +
  geom_area(aes(alpha = 0.5)) +
  scale_alpha(guide = "none")+
  
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(name = NULL,
                    values = alpha(plot_color, 0.7)) +
  
  labs(x = NULL,
       y = "Mean Relative Abundance (%)") +
  scale_x_continuous(breaks = 1:3,
                     labels = paste0(c("Health",
                                       "Acute",
                                       "Chronic"))) +
  theme_cowplot() +
  theme(axis.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        
        legend.title = element_blank(),
        legend.text = element_markdown(face = "bold")) 

ggsave("plots/bal_rel_abund_phylum.png",
       dpi = 600,
       width = 5.75,
       height = 4,
       units = c("in"),
       bg = "white")
