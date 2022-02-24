# install.packages("devtools")
# devtools::install_github("ggloor/ALDEx_bioc")

library(ALDEx2)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggtext)
library(cowplot)
library(ggrepel)
library(EnhancedVolcano)
library(ggthemes)
library(scales)


## Read in metadata and table
metadata <- read.delim("data/metadata.txt")
metadata$Time_Point <- as.factor(metadata$Time_Point)

table <- read_tsv("data/feature_table.1488.txt", col_names = TRUE) %>% 
  rownames_to_column(var = "ASV") %>% 
  as.data.frame()
rownames(table) <- paste("ASV", table$ASV, sep = "")
table
# Generate Taxonomy file
taxonomy <- table %>% dplyr::select(featureid) %>% 
  separate(col = "featureid", 
           into = c("Kingdom", "Phylum", "Class", 
                    "Order", "Family", "Genus"), 
           sep = ";",
           fill = "right") %>% 
  rownames_to_column(var = "ASV")


metadata
taxonomy

table.tax <- table %>% dplyr::select(-featureid, -ASV) %>% 
  rownames_to_column(var = "ASV") %>% 
  pivot_longer(-ASV) %>% 
  inner_join(., taxonomy, by = "ASV") %>% 
  pivot_longer(cols = c("Kingdom", "Phylum", "Class", 
                        "Order", "Family", "Genus"),
               names_to = "level",
               values_to = "taxon")

family.table <- table.tax %>% 
  filter(level == "Family") %>% 
  group_by(name, taxon) %>% 
  summarise(family_count = sum(value)) %>% 
  ungroup() %>% 
  mutate(taxon = gsub("D_4__", "", taxon))


family.table


#######################################
# Healthy (0) vs Acute (6)

## filter down metadata to chronic and healthy
metadata.filt.acute <- metadata %>% 
  dplyr::select(sampleid, Time_Point, Sample) %>% 
  filter(Sample == "BAL") %>% 
  filter(Time_Point == 0 | 
           Time_Point == 6) %>% 
  as.data.frame()

## filter down table and rename rows to family ID
table.family.filt.acute <- family.table %>% 
  filter(name %in% metadata.filt.acute$sampleid) %>% 
  arrange(desc(name)) %>% 
  pivot_wider(names_from = name,
              values_from = family_count)

# Removes any NA row names and sets to "Undefined Family"
table.family.filt.acute[is.na(table.family.filt.acute)]<-"Undefined Family"
# Move rowname to col 1
table.family.filt.acute <-column_to_rownames(table.family.filt.acute, var = "taxon")
## psuedocount 1 to every cell in table
table.family.filt.acute.pseudo <- table.family.filt.acute + 1

## Mutation to calculate log2 fold change
log2FC.acute <- table.family.filt.acute.pseudo %>%
  # move family col back to rowname
  rownames_to_column(var = "family") %>% 
  # lengthen table to perform mutations
  pivot_longer(-family, names_to = "sampleid") %>% 
  left_join(., metadata.filt.acute, by = "sampleid") %>% 
  # Group by family and time point to find mean within each family at
  # both time points
  group_by(family, Time_Point) %>% 
  summarise(mean_fam = mean(value), .groups = "drop") %>% 
  ## widen table to split up time points
  pivot_wider(names_from = Time_Point, values_from = mean_fam)  %>% 
  # group by family for future plotting and perform log2 FC calculation
  # + FC = up in Chronic (right); - FC = up in Healthy (left)
  group_by(family) %>% 
  summarise(log2FC = log2(`6`/`0`)) %>% 
  arrange(desc(log2FC))

## DESEq2 Analysis, following 
# https://microbiome.github.io/OMA/differential-abundance.html#aldex2
# Build condition table for DESEq2
conds.acute <- c(rep("Healthy",8), rep("Acute",8))

x.acute <- aldex.clr(
  reads = table.family.filt.acute.pseudo,
  conds = conds.acute, 
  # 128 recommened for ttest, 1000 for rigorous effect size calculation
  mc.samples = 128, 
  denom = "all",
  verbose = FALSE
)

x_tt.acute <- aldex.ttest(
  x.acute, 
  paired.test = FALSE, 
  verbose = FALSE)

aldex_out.acute <- data.frame(x_tt.acute) %>% 
  rownames_to_column(var = "family") %>% 
  left_join(., log2FC.acute)


## Plots generated using Benjami-Hochberg corrected Welch's T test
## and 20 fold change in mean family counts between groups
h_vs_a.plot <- EnhancedVolcano(aldex_out.acute,
                               lab = aldex_out.acute$family,
                               x = 'log2FC',
                               y = 'we.eBH',
                               ylim = c(0,6),
                               pCutoff = 5e-02,
                               FCcutoff = log2(20),
                               drawConnectors = TRUE,
                               title = NULL,
                               subtitle = NULL,
                               caption = NULL,
                               axisLabSize = 12,
                               legendPosition = "right",
                               legendLabels = c("NS", 
                                                "Log2 FC", 
                                                "p < 0.05",
                                                "p < 0.05 and Log2 FC"),
                               legendLabSize = 10,
                               labSize = 4,
                               labFace = "bold",
                               colAlpha = 0.9) +
  theme_set(
    theme(axis.text = element_text(color = "black",
                                   face = "bold"),
          axis.title = element_text(color = "black",
                                    face = "bold"),
          legend.text = element_text(color = "black",
                                     face = "bold"),
          legend.text.align	= 0)
  )
h_vs_a.plot

ggsave("plots/h_vs_a_aldex2_we.BH_volcano.png",
       dpi = 600,
       width = 6,
       height = 4,
       units = c("in"),
       bg = "white")

#######################################
# Healthy (0) vs Chronic (36)

## filter down metadata to chronic and healthy
metadata.filt.chronic <- metadata %>% 
  dplyr::select(sampleid, Time_Point, Sample) %>% 
  filter(Sample == "BAL") %>% 
  filter(Time_Point == 0 | 
           Time_Point == 36) %>% 
  as.data.frame()

## filter down table and rename rows to family ID
table.family.filt.chronic <- family.table %>% 
  filter(name %in% metadata.filt.chronic$sampleid) %>% 
  arrange(desc(name)) %>% 
  pivot_wider(names_from = name,
              values_from = family_count)

# Removes any NA row names and sets to "Undefined Family"
table.family.filt.chronic[is.na(table.family.filt.chronic)]<-"Undefined Family"
# Move rowname to col 1
table.family.filt.chronic <-column_to_rownames(table.family.filt.chronic, var = "taxon")
## psuedocount 1 to every cell in table
table.family.filt.chronic.pseudo <- table.family.filt.chronic + 1

## Mutation to calculate log2 fold change
log2FC.chronic <- table.family.filt.chronic.pseudo %>%
  # move family col back to rowname
  rownames_to_column(var = "family") %>% 
  # lengthen table to perform mutations
  pivot_longer(-family, names_to = "sampleid") %>% 
  left_join(., metadata.filt.chronic, by = "sampleid") %>% 
  # Group by family and time point to find mean within each family at
  # both time points
  group_by(family, Time_Point) %>% 
  summarise(mean_fam = mean(value), .groups = "drop") %>% 
  ## widen table to split up time points
  pivot_wider(names_from = Time_Point, values_from = mean_fam)  %>% 
  # group by family for future plotting and perform log2 FC calculation
  # + FC = up in Chronic (right); - FC = up in Healthy (left)
  group_by(family) %>% 
  summarise(log2FC = log2(`36`/`0`)) %>% 
  arrange(desc(log2FC))

## DESEq2 Analysis, following 
# https://microbiome.github.io/OMA/differential-abundance.html#aldex2
# Build condition table for DESEq2
conds.chronic <- c(rep("Healthy",8), rep("Chronic",8))

x.chronic <- aldex.clr(
  reads = table.family.filt.chronic.pseudo,
  conds = conds.chronic, 
  # 128 recommened for ttest, 1000 for rigorous effect size calculation
  mc.samples = 128, 
  denom = "all",
  verbose = FALSE
)

x_tt.chronic <- aldex.ttest(
  x.chronic, 
  paired.test = FALSE, 
  verbose = FALSE)

aldex_out.chronic <- data.frame(x_tt.chronic) %>% 
  rownames_to_column(var = "family") %>% 
  left_join(., log2FC.chronic)

#write_csv(aldex_out.chronic, file = "data/ALDEx2/h_vs_c_ALDEx2_output.csv")



## Plots generated using Benjami-Hochberg corrected Welch's T test
## and 20 fold change in mean family counts between groups

h_vs_c.plot <- EnhancedVolcano(aldex_out.chronic,
                lab = aldex_out.chronic$family,
                x = 'log2FC',
                y = 'we.eBH',
                ylim = c(0,6),
                pCutoff = 5e-02,
                FCcutoff = log2(20),
                drawConnectors = TRUE,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                axisLabSize = 12,
                legendPosition = "bottom",
                legendLabSize = 10,
                legendLabels = c("NS", 
                                 "Log2 FC", 
                                 "p < 0.05",
                                 "p < 0.05 and Log2 FC"),
                legendDropLevels = TRUE,
                labSize = 4,
                labFace = "bold",
                colAlpha = 0.9) +
  scale_fill_brewer(breaks = c("NS", 
                               "p < 0.05",
                               "p < 0.05 and Log2 FC")) +
  theme_set(
    theme(axis.text = element_text(color = "black",
                                   face = "bold"),
          axis.title = element_text(color = "black",
                                    face = "bold"),
          legend.text = element_text(color = "black",
                                     face = "bold"),
          legend.text.align	= 0
  ))
h_vs_c.plot

ggsave("plots/h_vs_c_aldex2_we.BH_volcano.png",
       dpi = 600,
       width = 6,
       height = 4,
       units = c("in"),
       bg = "white")


######################3
# Boxplots of >20 fold change
sigfc <- log2FC.chronic %>% 
  mutate(value = abs(log2FC) > log2(20))
sigfc


show_col(gdocs_pal()(3))
gdocs_pal()(3)
table.family.filt.chronic
sig_dif_boxplots <- table.family.filt.chronic %>% 
  rownames_to_column(var = "family") %>% 
  left_join(., sigfc, by = "family") %>% 
  filter(value == TRUE) %>% 
  select(-c(value, log2FC)) %>% 
  pivot_longer(-family, 
               names_to = "sampleid", 
               values_to = "count") %>% 
  left_join(., metadata.filt.chronic, by = "sampleid") %>% 
  select(-Sample) %>% 
  mutate(Time_Point = case_when(Time_Point == 0 ~ "Health",
                                Time_Point == 36 ~ "Chronic")) %>% 
  ggplot(aes(x = factor(Time_Point, level = c("Health", "Chronic")),
             y = count)) +
  geom_boxplot(aes(fill = Time_Point, alpha = 0.5),
               color = "black") +
  geom_dotplot(aes(fill = Time_Point),
               binaxis = 'y',
               stackdir = "center",
               color = NA) +
  facet_wrap(~factor(family, levels=c('Pseudomonadaceae',
                                     'Sphingomonadaceae',
                                     'Sphingobacteriaceae',
                                     "Moraxellaceae",
                                     "Xanthomonadaceae",
                                     "Xanthobacteraceae")),
         scales = "free",
         as.table = T) +
  
  ylab("Family Feature Count") +
  
  scale_fill_manual(values = c("#dc3912", "#ff9900")) +
  theme_cowplot() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text = element_text(face = "bold",
                                 size = 9),
        
        legend.position = "none",
        
        strip.background = NULL,
        strip.text = element_text(face = "bold",
                                  color = "black",
                                  size = 9))
sig_dif_boxplots

show_col(gdocs_pal()(2))
ggsave("plots/h_vs_c_sig.fam_boxplots.png",
       dpi = 600,
       width = 6,
       height = 4,
       unit = c("in"),
       bg = "white")  



plot_grid(h_vs_c.plot,
          sig_dif_boxplots, 
          nrow = 2,
          ncol = 1,
          rel_heights = c(0.5,0.5))

ggsave("plots/h_vs_c_volcano_barplots.png",
       height = 8,
       width = 6,
       dpi = 600,
       bg = "white")


