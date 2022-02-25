library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(dplyr)
library(plyr)
library(vegan) #2.5-7
library(ape)
library(ggthemes)
library(scales)



set.seed(1851)
permutations = 9999

## get metadata
metadata <- read_delim("data/metadata.txt") 

## get feature table and assign taxonomy as rownames
orig.table <- read_delim("data/feature_table.1488.txt") %>% 
  as.data.frame(row.names = featureid)
rownames(orig.table) <- orig.table$featureid

## Filter down metadata to points of interest
metadata.filtered <- metadata %>% 
  filter(Time_Point == 0 |
         Time_Point == 6 |
         Time_Point == 36) %>% 
  mutate(Sample = case_when(Sample == "BAL" ~ "BALF",
                            Sample == "OP" ~ "OP",
                            Sample == "Fecal" ~ "Rectal"),
         Time_Point = case_when(Time_Point == 0 ~ "Health",
                                Time_Point == 6 ~ "Acute",
                                Time_Point == 36 ~ "Chronic"))

metadata.filtered

## Filter feature table to only samples in metadata.filtered
## duplicate rownames and featureid in col 1, but select will solve that
table <- orig.table %>% 
  select(metadata.filtered$sampleid)
# Drop 0 rows
table <- table[rowSums(table[])>0,]
# Transpose table
table <- table %>% t()


## PCoA plots for all sample types

#quarter root transformation of table
table.transformed <- table^1/4

# Calculate jaccard & bray-curtis distances
bc_dist <- vegdist(table.transformed, method= "bray")

# Calculate PCoA with BC/J distances
bc_pcoa <- pcoa(bc_dist, correction = "cailliez") %>% 
  mutate(percent_var = (.$values$Eigenvalues/.$trace)*100)

## Pull PCo1/PCo2 vectors for plotting

# Rename values
bc_pcoa_vectors <- bc_pcoa$vectors %>% as_tibble(rownames = "sampleid") %>% 
  select(sampleid, Axis.1,Axis.2)
colnames(bc_pcoa_vectors) <- c("sampleid", "bc_PCo1", "bc_PCo2")

bc_hull <- as.data.frame(bc_pcoa_vectors) %>% 
  left_join(metadata.filtered) %>% 
  group_by(Sample) %>% 
  slice(chull(bc_PCo1, bc_PCo2))

bc_hull

## Pull variance represented for axis lables
bc_variance_rep <- round(bc_pcoa$percent_var[1:2],2)



#Join with metadata for plotting
pcoa.vectors.metadata <- inner_join(metadata.filtered, bc_pcoa_vectors, by = "sampleid")

## Colors used in following PCoA plot
display.brewer.pal(3, "Set2")
pcoa.vectors.metadata$Sample <- factor(pcoa.vectors.metadata$Sample, levels = c("OP", "BALF", "Rectal"))
# plot bc pcoa for all samples

show_col(gdocs_pal()(3))
pcoa.vectors.metadata %>% 
  ggplot(aes(x = bc_PCo1, y = bc_PCo2,
             color = Sample)) +
  
  geom_polygon(data = bc_hull,
               aes(fill = Sample,
                   color = Sample),
               alpha = 0.1,
               show.legend = FALSE) +

   # stat_ellipse(show.legend = F) +
  geom_point(aes(shape = as.factor(Time_Point)),
             alpha = 0.9,
             size = 3) +
  
  labs(x = paste("PCo1 - ",bc_variance_rep[1], "%", sep = ""),
       y = paste("PCo2 - ",bc_variance_rep[2], "%", sep = "")) +
  
  scale_shape_manual(values=c(0,1,2), labels = c("Health",
                                                 "Acute",
                                                 "Chronic"))+
  scale_color_manual(values=c("#dc3912","#ff9900", "#3366cc"), 
                     breaks = c("Rectal", "OP", "BALF"),
                     labels = c("Rectal", "OP", "BALF")) +
  scale_fill_manual(values=c("#dc3912","#ff9900", "#3366cc"), 
                    breaks = c("Rectal", "OP", "BALF"),
                    labels = c("Rectal", "OP", "BALF")) +
  
  theme_cowplot() +
  
  theme(legend.title = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text (face = "bold"),
        legend.text = element_text(face = "bold")) +
  guides(alpha = "none",
         shape = guide_legend(override.aes = list(size = 4,
                                                  shape = c(0,1,2),
                                                  stroke = 1.05)),
         color = guide_legend(override.aes = list(size = 4.5,
                                                  shape = 15)))

ggsave("plots/bc_pcoa_all_sample_types_convex_hulls.png",
       dpi = 600,
       bg = "white",
       width = 6,
       height = 4,
       units = c("in"))
# ggsave("plots/bc_pcoa_all_sample_types_ellipses.png",
#        dpi = 600,
#        bg = "white",
#        width = 6,
#        height = 4,
#        units = c("in"))
metadata.filtered
adonis.all.sample <- adonis(bc_dist ~ Time_Point * Sample, 
       permutations = permutations, data = metadata.filtered)


# Call:
# adonis(formula = bc_dist ~ Time_Point * Sample, data = metadata.filtered,      permutations = permutations) 
# 
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
#                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Time_Point          1     2.0893  2.0893  17.082 0.08581  1e-04 ***
#   Sample            2     2.1458  6.0729  49.653 0.49882  1e-04 ***
#   Time_Point:Sample 2     2.8976  1.4488  11.846 0.11901  1e-04 ***
#   Residuals         59    7.2161  0.1223         0.29636           
#   Total             64   24.3488                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


