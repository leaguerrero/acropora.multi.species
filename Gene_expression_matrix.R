# Multi Species TagSeq
# DESeq2 to make gene expression matricies
# For the heated samples, control samples, and the difference

# Load Libraries
library(DESeq2)
library(tidyverse)

# File Set-up

my.dir <- "/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/HTseq_out"
my.files <- grep(".txt", list.files(my.dir), value=TRUE)
my.metadata <- read.csv("/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/tagseq_meta_data.csv", header = TRUE)

# Create sample table
sampleNames <- my.metadata$seq.sample.id
my.sampleTable <- data.frame(sampleName = sampleNames, 
                             fileName = my.files, 
                             condition = my.metadata)
# Convert variables to factors
factorVars <- c("condition.species", "condition.colony.id", "condition.treatment")
my.sampleTable[, factorVars] <- lapply(my.sampleTable[, factorVars], factor)
colnames(my.sampleTable)

# Create DESeqDataSet Object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = my.sampleTable,
  directory = my.dir,
    design = ~ condition.species +
               condition.treatment +
               condition.species:condition.treatment)
# Pre-filter
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]


# Differential Expression Analysis
dds <- DESeq(ddsHTSeq) 
resultsNames(dds)

# Report: Number of transcripts
num_transcripts <- nrow(dds)
num_transcripts #22090

# Variance stabilizing transformation
vtd <- vst(dds)

# PCA: heat stress and species
deseq_PCA_heat_stress_species <- plotPCA(vtd,intgroup=c( "condition.treatment"), returnData = TRUE) 
deseq_PCA_heat_stress_species$species <- my.metadata$species
#deseq_PCA_heat_stress_species$species[deseq_PCA_heat_stress_species$species == "Acropora.unknown"] <- "Acropora.nasuta"
official.species.colors <- c("Acropora.tutuilensis" = "#F7C11E",
                             "Acropora.lutkeni" = "#FF9A17",
                             "Acropora.pulchra" = "#fb5607",
                             "Acropora.retusa" = "#ff006e",
                             "Acropora.nasuta" = "#c11cad",
                             "Acropora.robusta" = "#8338ec",
                             "Acropora.hyacinthus" = "#3E6FCB",
                             "Acropora.abrotanoides" = "#3a86ff") 

species.order <- c("Acropora.tutuilensis", "Acropora.lutkeni",
                   "Acropora.pulchra", "Acropora.retusa","Acropora.nasuta",
                   "Acropora.robusta","Acropora.hyacinthus","Acropora.abrotanoides")
deseq_PCA_heat_stress_species$species <- factor(deseq_PCA_heat_stress_species$species, levels = species.order)


# Plot PCA Heat Stress
legend.names <- c("Acropora tutuilensis", "Acropora lutkeni",
                                   "Acropora pulchra", "Acropora retusa","Acropora nasuta","Acropora robusta","Acropora hyacinthus","Acropora abrotanoides") 
shape.names <- c('29 °C','36 °C') 
PCA_plot_hs_species <- ggplot(deseq_PCA_heat_stress_species,aes(x=PC1,y=PC2, color=species, shape=condition.treatment)) +
  theme_classic(base_size = 20) +
  labs(title="PCA of gene expression colored by species") +
  geom_point(size = 5) +
  scale_color_manual(values = official.species.colors, labels = legend.names) +  
  #scale_shape_manual(labels = shape.names) +
  labs(color = "Species", shape = "Treatment") 
  
 

PCA_plot_hs_species + distance.v.ind.ed.50
ggsave(filename = "/Users/tillandsia/Lab Notebook/Exit seminar photos/multispecies_PCA.jpeg", 
       plot = PCA_plot_hs_species, 
       width = 10, 
       height = 8, 
       dpi = 800,
       device = "jpeg"
)

ggsave(filename = "/Users/tillandsia/Lab Notebook/Exit seminar photos/multispecies_PCA_transparent.jpeg", 
       plot = transparent_pca, 
       width = 10, 
       height = 8, 
       dpi = 800,
       device = "jpeg"
)

ggsave(filename = "/Users/tillandsia/Lab Notebook/Exit seminar photos/multispecies_PCA_euc_distance.jpeg", 
       plot = PCA_plot_hs_species +distance.v.ind.ed.50 , 
       width = 18, 
       height = 8, 
       dpi = 800,
       device = "jpeg"
)

#rm(list = ls(all.names = TRUE))


# Gene Expression: Heat stressed samples data frame
dds.hs <- dds[,dds@colData@listData$condition.treatment == "heat"]
hs.gene.expression.matrix <- counts(dds.hs, normalized=TRUE)

# Gene Expression: Control samples data frame
dds.cont <- dds[,dds@colData@listData$condition.treatment == "control"]
cont.gene.expression.matrix <- counts(dds.cont, normalized=TRUE)

## Heated - Control Matrix
diff.expression.matrix <- hs.gene.expression.matrix - cont.gene.expression.matrix

# Rename columns of each matrix
# Heated: Adjust column names
hs.metadata <- my.metadata[my.metadata$treatment=="heat",]
hs.species.colony.names <- paste(hs.metadata$species, hs.metadata$colony.id,"H", sep = ".")
head(hs.gene.expression.matrix)
colnames(hs.gene.expression.matrix) <- hs.species.colony.names
head(hs.gene.expression.matrix)
write.table(hs.gene.expression.matrix, 
            file = "/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/heat.exp.counts.txt",
            sep = '\t',
            quote = FALSE)



# Control: Adjust column names
cont.metadata <- my.metadata[my.metadata$treatment=="control",]
cont.species.colony.names <- paste(cont.metadata$species, cont.metadata$colony.id,"C", sep = ".")
head(cont.gene.expression.matrix)
colnames(cont.gene.expression.matrix) <- cont.species.colony.names
head(cont.gene.expression.matrix)
write.table(cont.gene.expression.matrix, 
            file = "/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/control.exp.counts.txt",
            sep = '\t',
            quote = FALSE)

# Diff: Adjust column names
diff.species.colony.names <- paste(cont.metadata$species, cont.metadata$colony.id,"diff", sep = ".")
colnames(diff.expression.matrix) <- diff.species.colony.names
head(diff.expression.matrix)
write.table(diff.expression.matrix, 
            file = "/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/difference.exp.counts.txt",
            sep = '\t',
            quote = FALSE)

# Final Check
head(hs.gene.expression.matrix)
head(cont.gene.expression.matrix)
head(diff.expression.matrix)


# Calculating the euclidean distance between control and heat stress samples for each colony in the PCA
deseq_PCA_heat_stress_species$colony <- my.metadata$colony.id[match((deseq_PCA_heat_stress_species$name), my.metadata$seq.sample.id)]
# Function to calculate Euclidean distance
euclidean_distance <- function(x1, y1, x2, y2) {
  return(sqrt((x2 - x1)^2 + (y2 - y1)^2))
}

# Initialize an empty dataframe to store distances
distance_df <- data.frame(colony = character(),
                          distance = numeric(),
                          stringsAsFactors = FALSE)

# Loop through unique colonies
for (colony in unique(deseq_PCA_heat_stress_species$colony)) {
  # Subset data for the current colony
  colony_data <- deseq_PCA_heat_stress_species[deseq_PCA_heat_stress_species$colony == colony, ]
  
  # Extract control and heat points
  control_points <- subset(colony_data, group == "control")[, c("PC1", "PC2")]
  heat_points <- subset(colony_data, group == "heat")[, c("PC1", "PC2")]
  
  # Calculate distances for each pair of control and heat points
  for (i in 1:nrow(control_points)) {
    for (j in 1:nrow(heat_points)) {
      # Calculate Euclidean distance
      dist <- euclidean_distance(control_points[i, 1], control_points[i, 2],
                                 heat_points[j, 1], heat_points[j, 2])
      
      # Append the distance to the distance dataframe
      distance_df <- rbind(distance_df, data.frame(colony = colony, distance = dist))
    }
  }
}

# Print the resulting dataframe
distance_df
distance_df$species <- my.metadata$species[match((distance_df$colony), my.metadata$colony.id)]
species.ed.50 <- read.csv("~/Lab Notebook/Chapter3/Thermal_stress_analysis/Phylogeny resolved ED50 curves/individual.specied.ED50.csv")
distance_df$species.ed.50 <-species.ed.50$mean.therm.tol[match((distance_df$colony), species.ed.50$colony.id)]
distance_df$ind.ed.50 <- species.ed.50$ind.rel.therm.tol[match((distance_df$colony), species.ed.50$colony.id)]
#distance_df$species[distance_df$species =="Acropora.unknown"] <- "Acropora nasuta"
distance_df$species <- factor(distance_df$species, levels = species.order)

library(ggpubr)
# distance.v.ind.ed.50 <- ggplot(distance_df, aes(x = ind.ed.50, y = distance)) +
#   theme_classic(base_size = 15) +
#   labs(title = "PC distance between\nheat-stressed and control samples",
#        x = "Estimated colony thermal tolerance",
#        y = "Euclidean distance") +
#   geom_point(aes(color = species), size = 5) +
#   #geom_smooth(method = "lm", se = FALSE, aes(group = 1)) +
#   geom_smooth(method = "lm", se = FALSE) +
#   scale_color_manual(values = official.species.colors, labels = legend.names) +  
#   labs(color = "Species") 

# remove the legend:
distance.v.ind.ed.50 <- ggplot(distance_df, aes(x = ind.ed.50, y = distance)) +
  theme_classic(base_size = 15) +
  labs(title = "PC distance between\nheat-stressed and control samples",
       x = "Estimated colony thermal tolerance",
       y = "Euclidean distance") +
  geom_point(aes(color = species), size = 5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = official.species.colors, labels = legend.names) +  
  theme(legend.position = "none")

# Calculate the linear regression model
lm_model <- lm(distance ~ ind.ed.50, data = distance_df)

# Extract the p-value for the slope coefficient
p_value <- summary(lm_model)$coefficients[2, 4]

# Add the p-value to the plot
distance.v.ind.ed.50 <- distance.v.ind.ed.50 + annotate("text", x = max(distance_df$ind.ed.50), 
                                                        y = max(distance_df$distance), 
                                                        label = paste("p-value =", round(p_value, 4)), 
                                                        hjust = 1, vjust = 1, size = 4) 

distance.v.ind.ed.50

# Calculate the Spearman correlation
spearman_result <- cor.test(distance_df$ind.ed.50, distance_df$distance, method = "spearman")

# Extract the correlation coefficient and p-value
spearman_cor <- round(spearman_result$estimate, 2)
spearman_p <- format.pval(spearman_result$p.value, digits = 2)

# Create the plot and add the correlation result
distance.v.ind.ed.50 <- ggplot(distance_df, aes(x = ind.ed.50, y = distance)) +
  theme_classic(base_size = 20) +
  labs(title = "PC distance between\nheat-stressed and control samples",
       x = "Estimated colony thermal tolerance",
       y = "Euclidean distance") +
  geom_point(aes(color = species), size = 5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = official.species.colors) +  
  theme(legend.position = "none") +
  labs(color = "Species") +
  annotate("text", x = Inf, y = Inf, label = paste("Spearman's ρ =", spearman_cor, "\np-value =", spearman_p), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black", fontface = "italic")

# Print the plot
print(distance.v.ind.ed.50)

library(patchwork)
PCA_plot_hs_species +
  inset_element(distance.v.ind.ed.50, 0.6, 0.6, 1, 1) +
  theme_classic()
PCA_plot_hs_species + distance.v.ind.ed.50 + plot_annotation(tag_levels = list(c('A.', 'B.')))

## Distance no abro
# Exclude abro samples
abro_samples <- c("RBAMS01", "RBAMS02", "RBAMS03", "RBAMS04", "RBAMS05", "RBAMS06",
                  "RBAMS07", "RBAMS08", "RBAMS09", "RBAMS10")

no_abro <- !colnames(vtd) %in% abro_samples

vtd_no_abro <- vtd[, no_abro]

deseq_PCA_heat_stress_species_no_abro <- plotPCA(vtd_no_abro,intgroup=c( "condition.treatment"), returnData = TRUE) 
deseq_PCA_heat_stress_species_no_abro$species <- my.metadata$species[-(1:10)]
official.species.colors.no.abro <- c("Acropora.tutuilensis" = "#F7C11E",
                             "Acropora.lutkeni" = "#FF9A17",
                             "Acropora.pulchra" = "#fb5607",
                             "Acropora.retusa" = "#ff006e",
                             "Acropora.nasuta" = "#c11cad",
                             "Acropora.robusta" = "#8338ec",
                             "Acropora.hyacinthus" = "#3E6FCB")

species.order.no.abro <- c("Acropora.tutuilensis", "Acropora.lutkeni",
                   "Acropora.pulchra", "Acropora.retusa","Acropora.nasuta",
                   "Acropora.robusta","Acropora.hyacinthus")

deseq_PCA_heat_stress_species_no_abro$species <- factor(deseq_PCA_heat_stress_species_no_abro$species, levels = species.order)


# Plot PCA Heat Stress
legend.names.no.abro <- c("Acropora tutuilensis", "Acropora lutkeni",
                  "Acropora pulchra", "Acropora retusa","Acropora nasuta","Acropora robusta","Acropora hyacinthus") 
shape.names <- c('29 °C','36 °C') 
PCA_plot_hs_species_no_abro <- ggplot(deseq_PCA_heat_stress_species_no_abro,aes(x=PC1,y=PC2, color=species, shape=condition.treatment)) +
  theme_classic(base_size = 20) +
  labs(title="PCA of gene expression colored by species") +
  geom_point(size = 5) +
  scale_color_manual(values = official.species.colors.no.abro, labels = legend.names.no.abro) +  
  #scale_shape_manual(labels = shape.names) +
  labs(color = "Species", shape = "Treatment") 



PCA_plot_hs_species_no_abro
ggsave(filename = "/Users/tillandsia/Lab Notebook/Chapter2/Data_analysis/TagSeq/Figures/Exploratory_pcas/PCA_hs_day.jpg", 
       plot = PCA_plot_hs_species_no_abro, 
       width = 8, 
       height = 6, 
       dpi = 300,
       device = "jpeg"
)

library(patchwork)


# Calculate distance:
head(deseq_PCA_heat_stress_species_no_abro)
my.metadata.no.abro <- my.metadata[-(1:10),]
deseq_PCA_heat_stress_species_no_abro$colony <- my.metadata.no.abro$colony.id[match((deseq_PCA_heat_stress_species_no_abro$name), my.metadata.no.abro$seq.sample.id)]
#Subset the dataframe into control and heat data
control_data <- subset(deseq_PCA_heat_stress_species_no_abro, condition.treatment == "control")
heat_data <- subset(deseq_PCA_heat_stress_species_no_abro, condition.treatment == "heat")

merged_data <- merge(control_data, heat_data, by = "colony", suffixes = c(".control", ".heat"))

# Calculate the distance in PC1 between control and heat for each colony
merged_data$PC1_distance <- abs(merged_data$PC1.control - merged_data$PC1.heat)

# Step 3: Create a new dataframe with the desired results
distance_dataframe <- merged_data[, c("colony", "PC1_distance")]
distance_dataframe$species <- my.metadata.no.abro$species[match((distance_dataframe$colony), my.metadata.no.abro$colony.id)]
distance_dataframe$ind.ed.50 <- species.ed.50$ind.rel.therm.tol[match((distance_dataframe$colony), species.ed.50$colony.id)]
distance_dataframe$species <- factor(distance_dataframe$species, levels = species.order.no.abro)


spearman_result <- cor.test(distance_dataframe$ind.ed.50, distance_dataframe$PC1_distance, method = "spearman")

# Extract the correlation coefficient and p-value
spearman_cor <- round(spearman_result$estimate, 2)
spearman_p <- format.pval(spearman_result$p.value, digits = 2)

distance.v.ind.ed.50.no.abro <- ggplot(distance_dataframe, aes(x = ind.ed.50, y = PC1_distance)) +
  theme_classic(base_size = 20) +
  labs(title = "PC1 distance between\nheat-stressed and control samples",
       x = "Estimated colony thermal tolerance",
       y = "PC1 distance") +
  geom_point(aes(color = species), size = 5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = official.species.colors.no.abro, labels = legend.names.no.abro) +  
  labs(color = "Species") +
  annotate("text", x = Inf, y = Inf, label = paste("Spearman's ρ =", spearman_cor, "\np-value =", spearman_p), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black", fontface = "italic")

# Print the plot
print(distance.v.ind.ed.50.no.abro)



# Switching the axes:
# Calculate the Spearman correlation
spearman_result <- cor.test(distance_df$distance,distance_df$ind.ed.50, method = "spearman")

# Extract the correlation coefficient and p-value
spearman_cor <- round(spearman_result$estimate, 2)
spearman_p <- format.pval(spearman_result$p.value, digits = 2)

# Create the plot and add the correlation result
distance.v.ind.ed.50 <- ggplot(distance_df, aes(x = distance ,y = ind.ed.50,)) +
  theme_classic(base_size = 20) +
  labs(title = "PC distance between\nheat-stressed and control samples",
       y = "Estimated colony thermal tolerance",
       x = "Euclidean distance") +
  geom_point(aes(color = species), size = 5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = official.species.colors) +  
  theme(legend.position = "none") +
  labs(color = "Species") +
  annotate("text", x = Inf, y = Inf, label = paste("Spearman's ρ =", spearman_cor, "\np-value =", spearman_p), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black", fontface = "italic")

# Print the plot
print(distance.v.ind.ed.50)

distance.v.ind.ed.50

deseq_PCA_heat_stress_species$alpha <- ifelse(deseq_PCA_heat_stress_species$colony == "c017" &
                                                deseq_PCA_heat_stress_species$condition.treatment %in% c("heat", "control"), 
                                              1,  # Full opacity for specific points
                                              0.2)  # Transparency for all other points



transparent_pca <- ggplot(deseq_PCA_heat_stress_species,aes(x=PC1,y=PC2, color=species, shape=condition.treatment)) +
  theme_classic(base_size = 20) +
  labs(title="PCA of gene expression colored by species") +
  geom_point(aes(alpha = alpha), size = 5) +
  scale_alpha_identity() + 
  scale_color_manual(values = official.species.colors, labels = legend.names) +  
  #scale_shape_manual(labels = shape.names) +
  labs(color = "Species", shape = "Treatment") 



