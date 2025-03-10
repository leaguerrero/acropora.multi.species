# Variance partitioning eigengenes?
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("variancePartition")
library("variancePartition")

###### Gene Expression

## NOTE 04/03: Add 'site' as a variable?

###### Eigengenes
# Load expression data:
# eigengene.exp = as.data.frame(t(MEs[, -c(1)]))
# write.csv(eigengene.exp, file = "~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/eigengene.exp.varpar.csv")
eigengene.exp <- read.csv("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/eigengene.exp.varpar.csv", header = TRUE)
rownames(eigengene.exp) <- eigengene.exp$X
eigengene.exp$X <- NULL  # Remove the 'X' column
# Remove grey module
eigengene.exp <- eigengene.exp[!(rownames(eigengene.exp) == "MEgrey"), ]
# load metadata:
my.metadata <- read.csv("/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/tagseq_meta_data.csv", header = TRUE)
ind.ed50 <- read.csv("~/Lab Notebook/Chapter3/Thermal_stress_analysis/Phylogeny resolved ED50 curves/individual.ED50.csv", header = TRUE)
my.metadata$thermal.treatment <-ifelse(my.metadata$treatment == "control", 29, 
                                       ifelse(my.metadata$treatment == "heat", 36, NA)) 

my.metadata$species.class <-as.numeric(match(my.metadata$species, 
                                             unique(my.metadata$species)))

my.metadata$ind.ed50 <- ind.ed50$ind.rel.therm.tol[match(my.metadata$colony.id, ind.ed50$colony.id)]
my.metadata$batch <- as.factor(my.metadata$batch)
my.metadata$collection.site <- as.factor(my.metadata$collection.site)
head(my.metadata)
row.names(my.metadata) <- my.metadata$seq.sample.id

# first run
# form <- ~ (1| treatment) + (1 | species) + (1 | colony.id) + ind.ed50 + (1 | batch) + (1 | species:treatment)
# varPart <- fitExtractVarPartModel(eigengene.exp, form, my.metadata)
# vp <- sortCols(varPart)
# plotPercentBars(vp[1:27, ])
# plotVarPart(vp)

# second run:
form <- ~ (1| treatment) + (1 | species) + (1 | colony.id) + ind.ed50 + (1 | batch) + (1 | collection.site)+ (1 | species:treatment)
varPart <- fitExtractVarPartModel(eigengene.exp, form, my.metadata)
vp <- sortCols(varPart)
plotPercentBars(vp[1:27, ])
plotVarPart(vp)

# Sort by variance explained by treatment
varPart.txt <- varPart[order(varPart$'species:treatment', decreasing = TRUE), ]
vp.txt <- sortCols(varPart.txt)
plotPercentBars(vp.txt[1:27, ])

# Sort by variance explained by residual
varPart.res <- varPart[order(varPart$Residuals, decreasing = FALSE), ]
vp.res <- sortCols(varPart.res)
plotPercentBars(vp.res[1:27, ])

#### High variance explained by residuals means the model can be improved?
### Does the residual persist when looking at all genes individually?
gene.counts.matrix <- read.csv("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/acro.counts.matrix.csv", header = TRUE)
row.names(gene.counts.matrix) <- gene.counts.matrix$gene
gene.counts.matrix <- gene.counts.matrix[,-1]
form <- ~ (1| treatment) + (1 | species) + (1 | colony.id) + ind.ed50 + (1 | batch) +(1 | species:treatment)
varPart.gene <- fitExtractVarPartModel(gene.counts.matrix, form, my.metadata)
vp.gene <- sortCols(varPart.gene)
plotPercentBars(vp.gene[1:50, ])
plotVarPart(vp.gene)

### Modelling species and treatment as fixed effects and batch and colony.id as
### random effects

form <- ~ treatment + species + ind.ed50 + batch + collection.site + species:treatment
varPart <- fitExtractVarPartModel(eigengene.exp, form, my.metadata)
vp <- sortCols(varPart)
plotPercentBars(vp[1:27, ])
plotVarPart(vp)


# Sort by variance explained by treatment
varPart.txt <- varPart[order(varPart$treatment, decreasing = TRUE), ]
vp.txt <- sortCols(varPart.txt)
plotPercentBars(vp.txt[1:27, ])

# Sort by variance explained by residual
varPart.res <- varPart[order(varPart$Residuals, decreasing = FALSE), ]
vp.res <- sortCols(varPart.res)
plotPercentBars(vp.res[1:27, ])

# Sort by variance explained by treatment:species
varPart.txt.spe <- varPart[order(varPart$'treatment:species', decreasing = TRUE), ]
vp.txt.spe <- sortCols(varPart.txt.spe)
plotPercentBars(vp.txt.spe[1:27, ])

#### Accounting for batch effects
# Analysis of residuals
library("limma")
# subtract out effect of Batch
fit <- lmFit(eigengene.exp, model.matrix(~batch, my.metadata))
res <- residuals(fit, eigengene.exp)

# fit model on residuals
form <- ~ treatment + species + ind.ed50 + collection.site + species:treatment
varPartResid <- fitExtractVarPartModel(res, form, my.metadata)
vp.res <- sortCols(varPartResid)
plotPercentBars(vp.res[1:27, ])
plotVarPart(vp.res)

# Sort by variance explained by treatment:species
varPart.res.txt.spe <- varPartResid[order(varPartResid$'treatment:species', decreasing = TRUE), ]
vp.res.txt.spe <- sortCols(varPart.res.txt.spe)
plotPercentBars(vp.res.txt.spe[1:27, ])

# Sort by variance explained by residuals
varPart.res <- varPartResid[order(varPartResid$Residuals, decreasing = FALSE), ]
vp.res <- sortCols(varPart.res)
plotPercentBars(vp.res[1:27, ])

## Removing colony ID
form <- ~ treatment + species + collection.site + species:treatment
varPartResid <- fitExtractVarPartModel(res, form, my.metadata)
vp.res <- sortCols(varPartResid)
plotPercentBars(vp.res[1:27, ])
plotVarPart(vp.res)

varPart.res.txt.spe <- varPartResid[order(varPartResid$'treatment:species', decreasing = TRUE), ]
vp.res.txt.spe <- sortCols(varPart.res.txt.spe)
plotPercentBars(vp.res.txt.spe[1:27, ])

# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
C <- canCorPairs(form, my.metadata)

# Plot correlation matrix
# between all pairs of variables
plotCorrMatrix(C)

# Adding a line that represents the average variance explained by treatment:species
library(tidyr)
library(dplyr)
vp.res.txt.spe$module <- row.names(vp.res.txt.spe)
long_df <- pivot_longer(vp.res.txt.spe, 
                        cols = c('treatment:species',treatment, species, collection.site, Residuals), 
                        names_to = "Variable", 
                        values_to = "Value")
long_df_sorted <- long_df %>%
  arrange(ifelse(Variable == "treatment:species", -Value, Value))

cutoff.var.exp.txt.spec <- mean(varPart.res.txt.spe$`treatment:species`)+sd(varPart.res.txt.spe$`treatment:species`)
mean.var.exp.txt.spec <- mean(varPart.res.txt.spe$`treatment:species`)

# violin plot
ggplot(long_df_sorted, aes(x=Variable, y=Value, fill = Variable))+
  #geom_violin( width=2) +
  #geom_boxplot(width=0.1, fill = "#e9ecef" ) +
  geom_boxplot() +
  scale_fill_manual(values = c("#e7cfbc","#f2b880",  "#006792","#c98686", "#966b9d"))+
  labs(y = "Variance Explained (%)") +
  theme_classic() +
  theme(legend.position="none", axis.title.x = element_blank()) +
  scale_x_discrete(limits=c("treatment:species", "species", "treatment", "collection.site", "Residuals"))

# stacked bar plot
ggplot(long_df_sorted, aes( x = reorder(module, ifelse(Variable == "treatment:species", Value, -Value)), y = Value, fill = factor(Variable, levels=c("Residuals", "treatment", "species", "collection.site", "treatment:species" )))) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#e7cfbc","#f2b880", "#006792", "#c98686", "#966b9d"))+
  geom_hline(aes(yintercept = mean.var.exp.txt.spec ), linetype = "dashed", color = "white") +
  labs(y = "Variance Explained (%)", x = "Module") +
  guides(fill = guide_legend(title = NULL)) +
  coord_flip() +
  theme_classic()
  

## Adding symbiont type as variable
eigengene.exp <- read.csv("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/eigengene.exp.varpar.csv", header = TRUE)
rownames(eigengene.exp) <- eigengene.exp$X
eigengene.exp$X <- NULL  # Remove the 'X' column
# Remove grey module
eigengene.exp <- eigengene.exp[!(rownames(eigengene.exp) == "MEgrey"), ]
# load metadata:
my.metadata <- read.csv("/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/tagseq_meta_data.csv", header = TRUE)
ind.ed50 <- read.csv("~/Lab Notebook/Chapter3/Thermal_stress_analysis/Phylogeny resolved ED50 curves/ind.ed50.symbionts.csv", header = TRUE)
my.metadata$thermal.treatment <-ifelse(my.metadata$treatment == "control", 29, 
                                       ifelse(my.metadata$treatment == "heat", 36, NA)) 

my.metadata$species.class <-as.numeric(match(my.metadata$species, 
                                             unique(my.metadata$species)))

my.metadata$ind.ed50 <- ind.ed50$ind.rel.therm.tol[match(my.metadata$colony.id, ind.ed50$colony.id)]
my.metadata$batch <- as.factor(my.metadata$batch)
my.metadata$collection.site <- as.factor(my.metadata$collection.site)
my.metadata$maj.sym <- ind.ed50$maj.sym[match(my.metadata$colony.id, ind.ed50$colony.id)]
head(my.metadata)
row.names(my.metadata) <- my.metadata$seq.sample.id

# first run:
form <- ~ (1| treatment) + (1 | species) + (1 | colony.id) + ind.ed50 + (1|maj.sym) +(1 | batch) + (1 | collection.site)+ (1 | species:treatment)
varPart <- fitExtractVarPartModel(eigengene.exp, form, my.metadata)
vp <- sortCols(varPart)
plotPercentBars(vp[1:27, ])
plotVarPart(vp)

# Sort by variance explained by treatment
varPart.txt <- varPart[order(varPart$'species:treatment', decreasing = TRUE), ]
vp.txt <- sortCols(varPart.txt)
plotPercentBars(vp.txt[1:27, ])

# Sort by variance explained by residual
varPart.res <- varPart[order(varPart$Residuals, decreasing = FALSE), ]
vp.res <- sortCols(varPart.res)
plotPercentBars(vp.res[1:27, ])

#### High variance explained by residuals means the model can be improved?
### Does the residual persist when looking at all genes individually?
gene.counts.matrix <- read.csv("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/acro.counts.matrix.csv", header = TRUE)
row.names(gene.counts.matrix) <- gene.counts.matrix$gene
gene.counts.matrix <- gene.counts.matrix[,-1]
form <- ~ (1| treatment) + (1 | species) + (1 | colony.id) + ind.ed50 + (1 | batch) +(1 | species:treatment)
varPart.gene <- fitExtractVarPartModel(gene.counts.matrix, form, my.metadata)
vp.gene <- sortCols(varPart.gene)
plotPercentBars(vp.gene[1:50, ])
plotVarPart(vp.gene)

### Modelling species and treatment as fixed effects and batch and colony.id as
### random effects

form <- ~ treatment + species + ind.ed50 + batch + collection.site + species:treatment
varPart <- fitExtractVarPartModel(eigengene.exp, form, my.metadata)
vp <- sortCols(varPart)
plotPercentBars(vp[1:27, ])
plotVarPart(vp)


# Sort by variance explained by treatment
varPart.txt <- varPart[order(varPart$treatment, decreasing = TRUE), ]
vp.txt <- sortCols(varPart.txt)
plotPercentBars(vp.txt[1:27, ])

# Sort by variance explained by residual
varPart.res <- varPart[order(varPart$Residuals, decreasing = FALSE), ]
vp.res <- sortCols(varPart.res)
plotPercentBars(vp.res[1:27, ])

# Sort by variance explained by treatment:species
varPart.txt.spe <- varPart[order(varPart$'treatment:species', decreasing = TRUE), ]
vp.txt.spe <- sortCols(varPart.txt.spe)
plotPercentBars(vp.txt.spe[1:27, ])

#### Accounting for batch effects
# Analysis of residuals
library("limma")
# subtract out effect of Batch
fit <- lmFit(eigengene.exp, model.matrix(~batch, my.metadata))
res <- residuals(fit, eigengene.exp)

# fit model on residuals
form <- ~ treatment + species + ind.ed50 + collection.site + species:treatment
varPartResid <- fitExtractVarPartModel(res, form, my.metadata)
vp.res <- sortCols(varPartResid)
plotPercentBars(vp.res[1:27, ])
plotVarPart(vp.res)

# Sort by variance explained by treatment:species
varPart.res.txt.spe <- varPartResid[order(varPartResid$'treatment:species', decreasing = TRUE), ]
vp.res.txt.spe <- sortCols(varPart.res.txt.spe)
plotPercentBars(vp.res.txt.spe[1:27, ])

# Sort by variance explained by residuals
varPart.res <- varPartResid[order(varPartResid$Residuals, decreasing = FALSE), ]
vp.res <- sortCols(varPart.res)
plotPercentBars(vp.res[1:27, ])

## Removing colony ID
form <- ~ treatment + species + collection.site + species:treatment
varPartResid <- fitExtractVarPartModel(res, form, my.metadata)
vp.res <- sortCols(varPartResid)
plotPercentBars(vp.res[1:27, ])
plotVarPart(vp.res)

varPart.res.txt.spe <- varPartResid[order(varPartResid$'treatment:species', decreasing = TRUE), ]
vp.res.txt.spe <- sortCols(varPart.res.txt.spe)
plotPercentBars(vp.res.txt.spe[1:27, ])

# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
C <- canCorPairs(form, my.metadata)

# Plot correlation matrix
# between all pairs of variables
plotCorrMatrix(C)

# Adding a line that represents the average variance explained by treatment:species
library(tidyr)
library(dplyr)
vp.res.txt.spe$module <- row.names(vp.res.txt.spe)
long_df <- pivot_longer(vp.res.txt.spe, 
                        cols = c('treatment:species',treatment, species, collection.site, Residuals), 
                        names_to = "Variable", 
                        values_to = "Value")
long_df_sorted <- long_df %>%
  arrange(ifelse(Variable == "treatment:species", -Value, Value))

cutoff.var.exp.txt.spec <- mean(varPart.res.txt.spe$`treatment:species`)+sd(varPart.res.txt.spe$`treatment:species`)
mean.var.exp.txt.spec <- mean(varPart.res.txt.spe$`treatment:species`)

# violin plot
ggplot(long_df_sorted, aes(x=Variable, y=Value, fill = Variable))+
  #geom_violin( width=2) +
  #geom_boxplot(width=0.1, fill = "#e9ecef" ) +
  geom_boxplot() +
  scale_fill_manual(values = c("#e7cfbc","#f2b880",  "#006792","#c98686", "#966b9d"))+
  labs(y = "Variance Explained (%)") +
  theme_classic() +
  theme(legend.position="none", axis.title.x = element_blank()) +
  scale_x_discrete(limits=c("treatment:species", "species", "treatment", "collection.site", "Residuals"))

# stacked bar plot
variance_partitioning <- ggplot(long_df_sorted, aes( x = reorder(module, ifelse(Variable == "treatment:species", Value, -Value)), y = Value, fill = factor(Variable, levels=c("Residuals", "treatment", "species", "collection.site", "treatment:species" )))) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#e7cfbc","#f2b880", "#006792", "#c98686", "#966b9d"))+
  geom_hline(aes(yintercept = mean.var.exp.txt.spec ), linetype = "dashed", color = "white") +
  labs(y = "Variance Explained (%)", x = "Module") +
  guides(fill = guide_legend(title = NULL)) +
  coord_flip() +
  theme_classic(base_size = 15)

ggsave("~/Lab Notebook/Chapter3/Figures/Supplements/Var_par.png", plot = variance_partitioning, dpi = 600)
