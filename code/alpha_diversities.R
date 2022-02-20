library(readr)
library(tidyverse)
library(tidyr)
library(vegan)
library(cowplot)
library(ggpubr)
library(dplyr)

## Read in metadata and table
metadata <- read.delim("data/metadata.txt") %>% 
  mutate(Time_Point = case_when(Time_Point == 0 ~ "Healthy",
                                Time_Point == 6 ~ "Acute",
                                Time_Point == 36 ~ "Chronic")) %>% 
  mutate(Sample = case_when(Sample == "Fecal" ~ "Feces",
                            Sample == "OP" ~ "OP",
                            Sample == "BAL" ~ "BAL")) %>% 
  drop_na()

table <- read_tsv("data/feature_table.1488.txt", col_names = TRUE) %>% 
  rownames_to_column(var = "ASV") %>% 
  as.data.frame()
rownames(table) <- paste("ASV", table$ASV, sep = "")
metadata

alpha_metrics <- table %>% 
  dplyr::select(-c(ASV,featureid)) %>% 
  rownames_to_column(var = "ASV") %>% 
  pivot_longer(-ASV,
               names_to = "sampleid",
               values_to = "counts") %>% 
  group_by(sampleid) %>% 
  mutate(obs_richness = specnumber(counts)) %>%

  mutate(shannon_e = diversity(counts, 
                               index = "shannon",
                               base = exp(1))) %>% 
  
  mutate(simpson = diversity(counts, 
                             index = "simpson")) %>% 
  dplyr::select(-c(ASV, counts)) %>% 
  pivot_longer(-sampleid, names_to = "alpha_metric", values_to = "value") %>% 
  left_join(metadata, .)


alpha.sum.stats <- alpha_metrics %>% 
  group_by(Time_Point, alpha_metric, Sample) %>% 
  summarize(min = min(value),
            q1 = quantile(value, 0.25),
            median = median(value),
            mean = mean(value),
            q3 = quantile(value, 0.75),
            max = max(value),
            sd = sd(value)) %>% 
  arrange(alpha_metric, Sample)

write.csv(alpha.sum.stats, "data/alpha_summary_stats.csv",
          row.names = F)
  
## Richness Plots
alpha_metrics %>% 
  group_by(Sample, Time_Point) %>% 
  ggplot(aes(x = factor(Time_Point, levels = c("Healthy", "Acute", "Chronic")),
             y = value, 
             fill = Time_Point)) +
  geom_boxplot(aes(alpha = 0.5),
               color = "black") +

  facet_wrap(factor(Sample, levels = c("BAL", "OP", "Feces")) ~alpha_metric, 
             scales = "free",
             strip.position = "top",

  ) +

  scale_fill_manual(values = c("yellow", "red", "green")) +
  
  theme_bw() +
  
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black",
                                 face = "bold"),
        
        legend.position = "none",
        
        panel.background = element_blank(),
        panel.grid = element_blank(),

        strip.background = element_blank(),
        strip.text = element_blank()
        )

# ggsave("plots/alpha_diversity_labs.png",
#        dpi = 600,
#        width = 6,
#        height = 6,
#        units = c("in"))




