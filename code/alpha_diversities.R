library(readr)
library(tidyverse)
library(tidyr)
library(vegan)
library(cowplot)
library(ggpubr)
library(dplyr)
library(ggthemes)
library(rstatix)
library(readxl)
library(facetscales)


## Read in metadata and table
metadata <- read.delim("data/metadata.txt") %>% 
  mutate(Time_Point = case_when(Time_Point == 0 ~ "Health",
                                Time_Point == 6 ~ "Acute",
                                Time_Point == 36 ~ "Chronic")) %>% 
  mutate(Sample = case_when(Sample == "Fecal" ~ "Feces",
                            Sample == "OP" ~ "OP",
                            Sample == "BAL" ~ "BALF")) %>% 
  drop_na()
metadata
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

#write.csv(alpha.sum.stats, "data/alpha_summary_stats.csv",
 #         row.names = F)

obs_richness <- alpha_metrics %>% 
  filter(alpha_metric == "obs_richness") %>% 
  select(Cat, Time_Point, Sample, value) %>% 
  unique() %>% 
  pivot_wider(names_from = Sample, values_from = value)
simpson <- alpha_metrics %>% 
  filter(alpha_metric == "simpson") %>% 
  select(Cat, Time_Point, Sample, value) %>% 
  unique() %>% 
  pivot_wider(names_from = Sample, values_from = value)
shannon <- alpha_metrics %>% 
  filter(alpha_metric == "shannon_e") %>% 
  select(Cat, Time_Point, Sample, value) %>% 
  unique() %>% 
  pivot_wider(names_from = Sample, values_from = value)
write.csv(obs_richness, "data/alpha_metrics/obs_richness.csv", row.names = F)
write.csv(simpson, "data/alpha_metrics/simpson.csv", row.names = F)
write.csv(shannon, "data/alpha_metrics/shannon.csv", row.names = F)


scales_y <- list(
  `Observed ASVs` = scale_y_continuous(limits = c(0, 125)),
  `Simpson Index` = scale_y_continuous(limits = c(0.1, 1.1)),
  `Shannon Index` = scale_y_continuous(limits = c(0,4))
)

alpha_metrics
## Richness Plots
alpha_metrics %>%
  mutate(alpha_metric = case_when(alpha_metric == "obs_richness" ~ "Observed ASVs",
                                  alpha_metric == "shannon_e" ~ "Shannon Index",
                                  alpha_metric == "simpson" ~ "Simpson Index"),
         Sample = case_when(Sample == "OP" ~ "OP",
                            Sample == "BALF" ~ "BALF",
                            Sample == "Feces" ~ "Rectal")) %>% 
  #filter(alpha_metric == "obs_richness") %>% 
  group_by(Sample, Time_Point) %>% 
  ggplot(aes(x = factor(Time_Point, levels = c("Health", "Acute", "Chronic")),
             y = value, 
             fill = Time_Point)) +
  geom_boxplot(aes(alpha = 0.5),
               color = "black") +
  
  geom_point(aes(color = Time_Point)) +

  # facet_wrap(alpha_metric ~factor(Sample, levels = c("OP", "BALF", "Feces") ), 
  #            scales = "free",
  #            strip.position = "top"

  # ) +

  scale_fill_manual(values = gdocs_pal()(3)) +
  scale_color_manual(values = gdocs_pal()(3)) +
  #ylim(c(0,115)) +
  facet_grid_sc(rows = vars(alpha_metric),
                cols = vars(factor(Sample, levels = c("OP", "BALF", "Rectal"))),
                scales = list(y = scales_y), 
                switch = 'y') +

  
  theme_bw() +
  
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black",
                                 face = "bold", 
                                 size = 10),
        
        legend.position = "none",
        
        panel.background = element_blank(),
        panel.grid = element_blank(),

        strip.background = element_blank(),
        strip.text = element_text(face = "bold",
                                  color = "black",
                                  size = 13),
        strip.placement = "outside"
        #strip.text = element_blank()
        )


# ggsave("plots/richness.png",
#               dpi = 600,
#               width = 6,
#               height = 4,
#               units = c("in"))

# ggsave("plots/alpha_diversity_labs.png",
#        dpi = 600,
#        width = 6,
#        height = 6,
#        units = c("in"))




