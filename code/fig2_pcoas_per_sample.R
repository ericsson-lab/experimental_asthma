library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(dplyr)
library(plyr)
library(vegan) #2.5-7
library(ape)
library(scales)
library(ggthemes)
library(EcolUtils)

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

# Runs BC PCoA Analysis for plotting
run_pcoa <- function(sample_type) {
  metadata.filtered <- metadata.filtered %>% filter(Sample == sample_type)

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
  
  # Calculate PCoA with BC distances
  bc_pcoa <- pcoa(bc_dist, correction = "cailliez") %>%
    mutate(percent_var = (.$values$Eigenvalues/.$trace)*100)

  ## Pull PCo1/PCo2 vectors for plotting
  # Rename values
  bc_pcoa_vectors <- bc_pcoa$vectors %>% as_tibble(rownames = "sampleid") %>%
    select(sampleid, Axis.1,Axis.2)
  colnames(bc_pcoa_vectors) <- c("sampleid", "bc_PCo1", "bc_PCo2")

  ## Pull variance represented for axis lables
  bc_variance_rep <- round(bc_pcoa$percent_var[1:2],2)

  #Join with metadata for plotting
  pcoa.vectors.metadata <- inner_join(metadata.filtered, bc_pcoa_vectors, by = "sampleid")

  ## Generate hulls
  hull <- as.data.frame(pcoa.vectors.metadata) %>%
    group_by(as.factor(Time_Point)) %>%
    slice(chull(bc_PCo1, bc_PCo2))

  output <- list(pcoa.vectors.metadata, hull, bc_variance_rep, bc_dist, metadata.filtered)

}
 
# Create pcoa plots for each sample location
bal.pcoa <- run_pcoa("BALF")
fecal.pcoa <- run_pcoa("Rectal")
op.pcoa <- run_pcoa("OP")

# Colors for plotting    
# bal
# brewer.pal(3,"Set2")[1]
# feces
# brewer.pal(3,"Set2")[2]
# op
# brewer.pal(3,"Set2")[3]

## Colors used in following PCoA plot
show_col(gdocs_pal()(3))
bal.pcoa

# Plot BC pcoa results
plot_pcoa <- function(pcoa_object, color) {
  # plot bc pcoa for all samples
  ## Take vectors for each sample
  pcoa_object[[1]] %>% 
  ggplot(aes(x = bc_PCo1, y = bc_PCo2,
             color = as.factor(Time_Point))) +
  
  # Plots hulls for samples grouped by time
  # geom_polygon(data = pcoa_object[[2]],
  #              aes(fill = as.factor(Time_Point),
  #                  color = as.factor(Time_Point)),
  #              alpha = 0.1,
  #              show.legend = FALSE) +
   stat_ellipse(show.legend = F) +
  
  # Plot points for each sample above convex hull
  geom_point(aes(shape = as.factor(Time_Point)),alpha = 0.9,
             size = 3) +
    
  # Add labels from pcoa calculation
  labs(x = paste("PCo1 - ",pcoa_object[[3]][1], "%", sep = ""),
       y = paste("PCo2 - ",pcoa_object[[3]][2], "%", sep = "")) +
  
  # Changes point shapes to open symbols
  scale_shape_manual(values = c(0,1,2),
                     labels = c("Health",
                                "Acute",
                                "Chronic"))+
  # Change color to match fig 1
   scale_color_manual(values=c(rep(gdocs_pal()(3)[color],3)),
                      labels = c("Health",
                                 "Acute",
                                 "Chronic")) +
   scale_fill_manual(values=c(rep(gdocs_pal()(3)[color],3))) +
   
  theme_cowplot() +
  
  theme(legend.title = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text (face = "bold"),
        legend.text = element_text(face = "bold"))+
   guides(alpha = "none",
          shape = guide_legend(override.aes = list(size = 4.5,
                                                   shape = c(15,16,17))))
 
}

## Plot pcoa results
bal.plot <- plot_pcoa(bal.pcoa, 1)
fecal.plot <- plot_pcoa(fecal.pcoa, 2)
op.plot <- plot_pcoa(op.pcoa, 3)

# Get legend from one of the plots
leg <- get_legend(bal.plot + guides(shape = guide_legend(override.aes = list(color = "black",
                                                                             stroke = 1.05)) ))
## Plot pcoas together without legends.
pcoa.plots <- plot_grid(op.plot + theme(legend.position = "none"),
          bal.plot + theme(legend.position = "none"),
          fecal.plot + theme(legend.position = "none"),
          nrow = 1,
          labels = "AUTO")
## Plot pcoa plots and legend
plot_grid(pcoa.plots, 
          leg, 
          rel_widths = c(0.88, 0.12))

# Save plot
# ggsave("plots/bc_pcoa_split_ellipses.png",
#        dpi = 600,
#        bg = "white",
#        width = 10,
#        height = 4,
#        units = c("in"))
# ggsave("plots/bc_pcoa_split_convex_hulls.png",
#        dpi = 600,
#        bg = "white",
#        width = 10,
#        height = 4,
#        units = c("in"))


### PERMANOVA

bal.permanova <- adonis(bal.pcoa[4][[1]] ~ Time_Point, data = bal.pcoa[5][[1]], 
       permutations = permutations)
bal.permanova
# 
# Call:
#   adonis(formula = bal.pcoa[4][[1]] ~ Time_Point, data = bal.pcoa[5][[1]],      permutations = permutations) 
# 
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
#                 Df  SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#   Time_Point    1      3.6471  3.6471  45.039 0.67183  1e-04 ***
#    Residuals    22     1.7815  0.0810         0.32817           
#     Total      23      5.4285                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# BAL Pairwise
adonis.pair(bal.pcoa[4][[1]], as.factor(bal.pcoa[5][[1]]$Time_Point), 
            nper = permutations, corr.method = "BH")

# combination SumsOfSqs  MeanSqs   F.Model        R2 P.value P.value.corrected
# 1     0 <-> 6  0.862562 0.862562  10.42165 0.4267382   4e-04             4e-04
# 2    0 <-> 36  3.845041 3.845041 374.79133 0.9639910   4e-04             4e-04
# 3    6 <-> 36  1.528647 1.528647  17.26329 0.5521904   3e-04             4e-04

OP.permanova <- adonis(op.pcoa[4][[1]] ~ Time_Point, data = op.pcoa[5][[1]], 
                        permutations = permutations)
OP.permanova
# Call:
#   adonis(formula = op.pcoa[4][[1]] ~ Time_Point, data = op.pcoa[5][[1]],      permutations = permutations) 
# 
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
#               Df SumsOfSqs MeanSqs F.Model    R2 Pr(>F)  
#   Time_Point  1   0.35788 0.35788  2.9348 0.155 0.0106 *
#   Residuals  16   1.95108 0.12194         0.845         
#   Total      17   2.30895                 1.000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis.pair(op.pcoa[4][[1]], as.factor(op.pcoa[5][[1]]$Time_Point), 
            nper = permutations, corr.method = "BH")

# combination  SumsOfSqs    MeanSqs   F.Model         R2 P.value P.value.corrected
# 1     0 <-> 6 0.07274222 0.07274222 0.6457784 0.06694933  0.7024            0.7024
# 2    0 <-> 36 0.33431143 0.33431143 2.6219132 0.17931396  0.0087            0.0261
# 3    6 <-> 36 0.19559915 0.19559915 1.4358455 0.13758784  0.1712            0.2568


feces.permanova <- adonis(fecal.pcoa[4][[1]] ~ Time_Point, data = fecal.pcoa[5][[1]], 
                        permutations = permutations)
feces.permanova
# Call:
#   adonis(formula = fecal.pcoa[4][[1]] ~ Time_Point, data = fecal.pcoa[5][[1]],      permutations = permutations) 
# 
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#   Time_Point  1    0.9760 0.97603  5.8838 0.21886  3e-04 ***
#   Residuals  21    3.4836 0.16589         0.78114           
#   Total      22    4.4596                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis.pair(fecal.pcoa[4][[1]], as.factor(fecal.pcoa[5][[1]]$Time_Point), 
            nper = permutations, corr.method = "BH")
# combination SumsOfSqs   MeanSqs  F.Model         R2 P.value P.value.corrected
# 1     0 <-> 6 0.2246406 0.2246406 1.321253 0.08623664  0.2121           0.21210
# 2    0 <-> 36 0.9221754 0.9221754 6.208020 0.32319936  0.0013           0.00390
# 3    6 <-> 36 0.6203754 0.6203754 3.532139 0.21365286  0.0103           0.01545
