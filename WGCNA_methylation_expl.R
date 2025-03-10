# Chapter 3: WGCNA on methylation data
# WGCNA is designed to identify modules or groups of genes that
# are co-methylated.
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.18")

# BiocManager::install("preprocessCore")
# BiocManager::install("WGCNA", force = TRUE)

library(WGCNA)
library(methylKit)
library(tidyr)
library(tidyverse)
library(broom)
library(purrr)
library(rtracklayer)
library(GenomicRanges)
library(ggpubr)
options(stringsAsFactors = FALSE)

#Enable multithread
enableWGCNAThreads(nThreads = 6)

# format the methylation data
my.dir <- "/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/"
meta.file.name <- "../meth_meta.csv"
meta.path <- file.path(my.dir, meta.file.name)
meta <- read.csv(meta.path, header = TRUE)
meta$ID <- as.factor(meta$ID)
meta$species <- as.factor(meta$species)
meta$species.class <- as.numeric(meta$species.class)

# Offical_species_colors
official.species.colors <- c("Acropora tutuilensis" = "#F7C11E",
                             "Acropora lutkeni" = "#FF9A17",
                             "Acropora pulchra" = "#fb5607",
                             "Acropora retusa" = "#ff006e",
                             "Acropora unknown" = "#c11cad",
                             "Acropora robusta" = "#8338ec",
                             "Acropora hyacinthus" = "#3a86ff",
                             "Acropora abrotanoides" = "#3E6FCB")

# Files for methylKit sorted in order to group the species
file.list <- list("C033_S22.rmdup_CpG.methylKit",
                  "C075_S29.rmdup_CpG.methylKit",
                  "C076_S28.rmdup_CpG.methylKit",
                  "C077_S30.rmdup_CpG.methylKit",
                  "C078_S26.rmdup_CpG.methylKit",
                  "C011_S25.rmdup_CpG.methylKit",
                  "C017_S23.rmdup_CpG.methylKit",
                  "C018_S27.rmdup_CpG.methylKit",
                  "C038_S31.rmdup_CpG.methylKit",
                  "C051_S24.rmdup_CpG.methylKit",
                  "C047_S36.rmdup_CpG.methylKit",
                  "C196_S35.rmdup_CpG.methylKit",
                  "C210_S41.rmdup_CpG.methylKit",
                  "C004_S34.rmdup_CpG.methylKit",
                  "C007_S39.rmdup_CpG.methylKit",
                  "C013_S40.rmdup_CpG.methylKit",
                  "C016_S38.rmdup_CpG.methylKit",
                  "C030_S37.rmdup_CpG.methylKit",
                  "C032_S33.rmdup_CpG.methylKit",
                  "C040_S32.rmdup_CpG.methylKit",
                  "C094_S42.rmdup_CpG.methylKit",
                  "C090_S44.rmdup_CpG.methylKit",
                  "C123_S45.rmdup_CpG.methylKit",
                  "C167_S46.rmdup_CpG.methylKit",
                  "C211_S47.rmdup_CpG.methylKit",
                  "C014_S48.rmdup_CpG.methylKit",
                  "C024_S49.rmdup_CpG.methylKit",
                  "C066_S50.rmdup_CpG.methylKit",
                  "C071_S51.rmdup_CpG.methylKit",
                  "C122_S52.rmdup_CpG.methylKit",
                  "C137_S53.rmdup_CpG.methylKit",
                  "C164_S54.rmdup_CpG.methylKit",
                  "C037_S55.rmdup_CpG.methylKit",
                  "C012_S56.rmdup_CpG.methylKit",
                  "C120_S58.rmdup_CpG.methylKit",
                  "C009_S59.rmdup_CpG.methylKit")


# Samples
sample.id=list("c033",
               "c075",
               "c076",
               "c077",
               "c078",
               "c011",
               "c017",
               "c018",
               "c038",
               "c051",
               "c047",
               "c196",
               "c210",
               "c004",
               "c007",
               "c013",
               "c016",
               "c030",
               "c032",
               "c040",
               "c094",
               "c090",
               "c123",
               "c167",
               "c211",
               "c014",
               "c024",
               "c066",
               "c071",
               "c122",
               "c137",
               "c164",
               "c037",
               "c012",
               "c120",
               "c009")


# Methylation object creation
setwd("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
myobj=methRead(file.list,                   
               sample.id=sample.id,       
               assembly="Ahyacinthus.chrsV1",
               treatment=c(0,0,0,0,0,1,1,1,1,1,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,6,6,6,7,7,7,3),
               context="CpG",
               dbtype = "tabix",
               dbdir = "methylDB")
filtered.myobj=filterByCoverage(myobj,
                                lo.count=10,
                                lo.perc=NULL,
                                hi.count=NULL,
                                hi.perc=99.9)



# Create GR object to overlay with methylation data
gr.obj <- import("~/Lab Notebook/Chapter3/WGBS_analysis/Ahyacinthuns.genes.gff") 
# Select ranges from chr1:chr14
chr.regions <-gr.obj[seqnames(gr.obj) %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14")]
gr.df <- data.frame(chr.regions)
names(gr.df)
# Make a column with chr.start.end
gr.df$concat <- paste(gr.df$seqnames,gr.df$start,gr.df$end,sep=".")

# Region counts
genes <- methylKit::regionCounts(filtered.myobj,chr.regions,save.db=FALSE) # Gene counts for individuals

# plot number of genes/sample and number of genes/species
num_observations <- c()
sample_ids <- c()

# Iterate over each element of the list
for (i in seq_along(genes)) {
  # Extract number of observations and sample ID
  num_obs <- nrow(genes[[i]])
  sample_id <- attr(genes[[i]], "sample.id")
  
  # Append to the vectors
  num_observations <- c(num_observations, num_obs)
  sample_ids <- c(sample_ids, sample_id)
}

# Create a data frame with the results
result_df <- data.frame(sample_id = sample_ids, num_observations = num_observations)
result_df$prop_obs <- result_df$num_observations/24511
result_df$species <- meta$species[match((result_df$sample_id),meta$ID)]
 
ggplot(result_df, aes(x = sample_id, y = prop_obs)) +
  geom_bar(stat = "identity", fill = "#f4a261", color = "#ffbf69") +
  theme_minimal() +
  ylim(0,1)+
  labs(x = "Sample ID", y = "Proportion of Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggplot(result_df, aes(x = species, y = num_observations)) +
  geom_bar(stat = "identity", fill = "#c8b6ff", color = "#ffd6ff") +
  theme_minimal() +
  labs(x = "Species", y = "Number of Observations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
 
genes.unite <- methylKit::unite(genes,
                                destrand = FALSE, #Combine
                                min.per.group = 2L, # in at least 2 samples per species
                                # 2L is a good balance between gene and sample retention when
                                # assessing missing data.
                                save.db = F)

gene.pm.df <- as.data.frame(percMethylation(genes.unite, rowids = TRUE))
gene.pm.df$concat <- rownames(gene.pm.df)

gene.names <- gr.df$ID[match((gene.pm.df$concat),gr.df$concat)]

rownames(gene.pm.df) <- gene.names

gene.pm.df <- gene.pm.df[,!names(gene.pm.df) %in% c("concat")] 

dim(gene.pm.df) # 4376 genes, 36 samples
## Rename columns to reflect species and colony ID
spec.and.colony.colnames <- c("c033-abro", "c075-abro", "c076-abro", "c077-abro", "c078-abro", "c011-hya", "c017-hya", "c018-hya", "c038-hya", "c051-hya", "c047-lut", "c196-lut",
                              "c210-lut", "c004-pul", "c007-pul", "c013-pul", "c016-pul", "c030-ret", "c032-ret", "c040-ret", "c094-ret", "c090-rob", "c123-rob", "c167-rob",
                              "c211-rob", "c014-tut", "c024-tut", "c066-tut", "c071-tut", "c122-tut", "c137-tut", "c164-tut", "c037-unk", "c012-unk", "c120-unk", "c009-pul")

colnames(gene.pm.df) <- spec.and.colony.colnames

# Looking at the pre-imputed distribution of the percent methylation across genes
gene.pm.long <- gather(gene.pm.df, sample_id, percent_methylation)
gene.pm.long$species <- meta$species[match((gene.pm.long$sample_id),meta$col.spec.id)]

# Ridgeline plot- nonimputed data
library(ggridges)
gene.pm.long$imputation_status <- "non imputed"
gene.pm.long$species <- factor(gene.pm.long$species, levels =c("Acropora abrotanoides",
                                                                       "Acropora hyacinthus",
                                                                       "Acropora robusta",
                                                                       "Acropora unknown",
                                                                       "Acropora retusa" ,
                                                                       "Acropora pulchra",
                                                                       "Acropora lutkeni",
                                                                       "Acropora tutuilensis") )
non_imp_ridgeline <- ggplot(gene.pm.long, aes(x = percent_methylation, y =species , fill = species)) + 
  stat_density_ridges(scale = 1, quantile_lines = TRUE, na.rm = TRUE, quantiles = 2) +
  scale_fill_manual(values = official.species.colors) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Non-imputed % methylation of genes") +
  xlab("% methylation")
non_imp_ridgeline

# Reshape data for WGCNA
perc.meth.df <- tibble::rownames_to_column(gene.pm.df, "Genes")
datMeth0 = as.data.frame(t(perc.meth.df [, -c(1)]))
names(datMeth0) = perc.meth.df$Genes;
rownames(datMeth0) = names(perc.meth.df)[-c(1)]
sample.names <- names(perc.meth.df)[-c(1)]

gsg = goodSamplesGenes(datMeth0, verbose = 3)
gsg$allOK

# 1 samples and 4 genes have too many NAs
# Removing bad genes and samples
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datMeth0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datMeth0)[!gsg$goodSamples], collapse = ", ")))
  datMeth0 = datMeth0[gsg$goodSamples, gsg$goodGenes]
}

# Samples removed: c033-abro

# Retry
gsg = goodSamplesGenes(datMeth0, verbose = 3)
gsg$allOK
rm(gsg)

# Check Sample Outliters
sampleTree.1 = hclust(dist(datMeth0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "~/Lab Notebook/Chapter3/WGBS_analysis/Figures/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree.1, main = "Sample clustering to detect outliers", 
     sub="",
     xlab="", 
     cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 250, col = "red")
dev.off()
# No outliers
datMeth <- datMeth0 
rm(datMeth0)

# Impute the data here:
# Summarize NAs by column (sample)
na_by_gene <- colSums(is.na(datMeth))

# Summarize NAs by row (gene)
na_by_sample <- rowSums(is.na(datMeth))

summary(na_by_gene)
summary(na_by_sample)
# Histogram of NAs by sample
hist(na_by_sample, xlab = "Number of NAs",ylab = "Number of samples", main = "Distribution of NAs by Sample")

# Bar plot of NAs by gene
barplot(na_by_gene, xlab = "Gene", ylab = "Number of NAs", axisnames = FALSE, main = "NAs by Gene")

# missMDA to impute missing values
# https://martha-labbook.netlify.app/posts/example-of-data-impuation-with-missmda/

library(missMDA)
library(FactoMineR)
nb = estim_ncpPCA(datMeth)

res.comp = imputePCA(datMeth,ncp=nb$ncp)
impute.na.by.gene <- colSums(is.na(res.comp$fittedX))
summary(impute.na.by.gene)

# Comparing percent methylation density plots

# long data frame
gene.pm.imp.df <- res.comp$fittedX
gene.pm.imp.df= as.data.frame(t(gene.pm.imp.df))
gene.pm.imp.long <- gather(gene.pm.imp.df, sample_id, percent_methylation)
gene.pm.imp.long$species <- meta$species[match((gene.pm.imp.long$sample_id),meta$col.spec.id)]
gene.pm.imp.long$imputation_status <- "imputed"

# visualize
gene.pm.imp.long$species <- factor(gene.pm.imp.long$species, levels =c("Acropora abrotanoides",
                                                                       "Acropora hyacinthus",
                                                                       "Acropora robusta",
                                                                       "Acropora unknown",
                                                                       "Acropora retusa" ,
                                                                       "Acropora pulchra",
                                                                       "Acropora lutkeni",
                                                                       "Acropora tutuilensis") )

imputed_ridgeline <- ggplot(gene.pm.imp.long, aes(x = percent_methylation, y =species , fill = species)) +
  stat_density_ridges(scale = 1, quantile_lines = TRUE, quantiles = 2) +
  scale_fill_manual(values = official.species.colors) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Imputed % methylation of genes") +
  xlab("% methylation")

library(gridExtra)
grid.arrange(non_imp_ridgeline, imputed_ridgeline, ncol=2)

# Join imputed and non-imputed data for better visualization
combined.gene.pm.df <- rbind(gene.pm.imp.long, gene.pm.long)
ggplot(combined.gene.pm.df, aes(x = percent_methylation, y =imputation_status , fill = species)) +
  stat_density_ridges(scale = 0.5, quantile_lines = TRUE, quantiles = 2) +
  geom_boxplot(aes(x = percent_methylation, y = imputation_status, fill = species), 
               width = 0.18, outlier.shape = NA, position = position_nudge(y = -0.2)) +
  scale_fill_manual(values = official.species.colors) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("% methylation of genes") +
  xlab("% methylation") +
  ylab("Imputation status") +
  facet_wrap(~species)

# Assess difference in means
## Data is not normally distributed
wilcox.res<- wilcox.test(gene.pm.long$percent_methylation, gene.pm.imp.long$percent_methylation)
wilcox.res


# Add significance stars to the plot
ggboxplot(combined.gene.pm.df, x = "imputation_status", y = "percent_methylation", 
          fill = "imputation_status", palette = c("#FFA585", "#FFEDA0"), outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  labs(title = "Boxplot of Percent Methylation",
       x = "Imputation Status",
       y = "Percent Methylation") +
  theme_minimal() + 
  theme(legend.position = "none") # + facet_wrap(~species)
  

# there is a significant difference in the means of percent of methylation
# between the imputed and non imputed data

# Assess SD
sd.imputed <- sd(gene.pm.imp.long$percent_methylation, na.rm = TRUE)
sd.nonimputed <- sd(gene_pm_long$percent_methylation, na.rm = TRUE)
sd.imputed
sd.nonimputed

# SD is slightly smaller in the imputed data than the non-imputed data


# Eigengene analysis
traitData = meta[-c(1),] # Removed c033 in GoodSamplesGenes
traitData$col.spec <- spec.and.colony.colnames[-1]
dim(traitData)
names(traitData)

# Form a data frame analogous to methylation data that will hold the species names.
meth.Samples = rownames(datMeth)
traitRows = match(meth.Samples, traitData$col.spec)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]
collectGarbage()

# Picking beta for expression data
# Choose a set of soft-thresholding powers
powers.e = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function

# Associating Modules and Species
#Relating modules to characteristics and identifying important genes
#Defining the number of genes and samples
nGenes = ncol(datMeth.imp) #4372
nSamples = nrow(datMeth.imp) #35

traitData = meta
dim(traitData)
names(traitData)
traitData.wide <- traitData %>%
  pivot_wider(
    id_cols = c(ID, col.spec.id),
    names_from = species,
    values_from = species.class,
    values_fill = 0
  )
# Order columns based on species relative thermal tolerances (rtt)
# Lowest rtt to highest:
# A. tut
# A. lut
# A. pul
# A. ret
# A. unk
# A. rob
# A. hya
# A. abr

ordered.col.names <- c("ID",
                       "col.spec.id",
                       "Acropora tutuilensis",
                       "Acropora lutkeni",
                       "Acropora pulchra",
                       "Acropora retusa", 
                       "Acropora unknown",
                       "Acropora robusta",
                       "Acropora hyacinthus",
                       "Acropora abrotanoides")
traitData.wide <- traitData.wide[,ordered.col.names]

# Form a data frame analogous to expression data that will hold the clinical traits.
datMeth.imp <- res.comp$fittedX
wgbs.Samples = rownames(datMeth.imp)
traitRows = match(wgbs.Samples, traitData.wide$col.spec.id)
datTraits = as.data.frame(traitData.wide[traitRows, -1])
traitDataSubset <- as.vector(traitData.wide[traitRows, 1])
rownames(datTraits) <- traitDataSubset$col.spec.id
collectGarbage()

powers.e = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft.e = pickSoftThreshold(datMeth.imp, powerVector = powers.e, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft.e$fitIndices[,1], -sign(sft.e$fitIndices[,3])*sft.e$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft.e$fitIndices[,1], -sign(sft.e$fitIndices[,3])*sft.e$fitIndices[,2],
     labels=powers.e,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.82,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft.e$fitIndices[,1], sft.e$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft.e$fitIndices[,1], sft.e$fitIndices[,5], labels=powers.e, cex=cex1,col="red")
abline(h=0.90,col="red")

# power = 3 and minModuleSize = 100
me.network= blockwiseModules(datMeth.imp, power = 3,
                             TOMType = "signed", minModuleSize = 100,
                             reassignThreshold = 0, mergeCutHeight = 0.25,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             saveTOMs = FALSE,
                             verbose = 3)
table(me.network$colors)     

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(me.network$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(me.network$dendrograms[[1]], mergedColors[me.network$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)                        

moduleLabels.e = me.network$colors
moduleColors.e = labels2colors(me.network$colors)
ME.es = me.network$MEs;
geneTree.e = me.network$dendrograms[[1]]        

softPower = 3 #Chosen in the graphs before
adjacency = adjacency(datMeth.imp, power = softPower, type = "signed") 
#Transforming the adjacency matrix in a topological overlap
TOM = TOMsimilarity(adjacency) #Calculating the topological overlap matrix
dissTOM = 1-TOM 

geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 100;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datMeth.imp, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.1
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datMeth.imp, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleColors = mergedColors

#Building numeric labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

# Associating Modules and Phenotypes
#Relating modules to characteristics and identifying important genes

#Recalculating MEs with label colors
MEs0 = moduleEigengenes(datMeth.imp, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs,datTraits[,-c(1)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)




#####
# Variance partitioning of the eigengenes
library("variancePartition")

# 
# ###### Eigengenes
# # Load methyltaion data:
eigengene.meth = as.data.frame(t(MEs[, -c(1)]))
# write.csv(eigengene.meth, file = "~/Lab Notebook/Chapter3/WGBS_analysis/eigengene.meth.csv")
# Remove grey module: No gray module?
# load metadata:


# remove c033 (didn't pass WGCNA test for missing data):

meta <- meta[-c(1),]
row.names(meta) <- meta$col.spec.id
form <- ~ (1| species) 
varPart <- fitExtractVarPartModel(eigengene.meth, form, meta)
vp <- sortCols(varPart)
plotPercentBars(vp)
plotVarPart(vp)

# Sort by variance explained by treatment
varPart.spec <- varPart[order(varPart$'species', decreasing = TRUE), ]
vp.spec <- sortCols(varPart.spec)
plotPercentBars(vp.spec)

# What genes are in the methylation eigengenes (modules?)
head(row.names(datMeth))

# Extract gene names in the turquoise, brown, blue, green and red methylation modules
# MEturquoise.meth.genes <- colnames(datMeth.imp)[moduleColors=="turquoise"] # 845
# MEbrown.meth.genes  <-colnames(datMeth.imp)[moduleColors=="brown"] #313
# MEblue.meth.genes  <-colnames(datMeth.imp)[moduleColors=="blue"] #701
# MEgreen.meth.genes <-colnames(datMeth.imp)[moduleColors=="green"] #298
# MEred.meth.genes <-colnames(datMeth.imp)[moduleColors=="red"] #196

# Visualize the gene level methylation across species:
# log of percent methylation of each species for these modules
# Turquoise module



official.species.colors <- c("Acropora tutuilensis" = "#F7C11E",
                             "Acropora lutkeni" = "#FF9A17",
                             "Acropora pulchra" = "#fb5607",
                             "Acropora retusa" = "#ff006e",
                             "Acropora nasuta" = "#c11cad",
                             "Acropora robusta" = "#8338ec",
                             "Acropora hyacinthus" = "#3a86ff",
                             "Acropora abrotanoides" = "#3E6FCB")

species.order <- c("Acropora tutuilensis", "Acropora lutkeni",
                   "Acropora pulchra", "Acropora retusa","Acropora nasuta",
                   "Acropora robusta","Acropora hyacinthus","Acropora abrotanoides")
all.species.pm.means.turquoise$species <- factor(all.species.pm.means.turquoise$species, levels = species.order)

# Plot
# rain_height <- .1
# pm.mean.plot.all.samples.turq <- ggplot(all.species.pm.means.turquoise, aes(x = species, y = pm.mean, fill = species )) +
#   # # clouds
#   # introdataviz::geom_flat_violin(trim=FALSE,
#   #                                position = position_nudge(x = rain_height+.05),
#   #                                width = 5) +
#   # rain
#   geom_point(aes(color = species), size = 1, alpha = 0.4,  show.legend = FALSE, 
#              position = position_jitter(width = rain_height, height = 0)) +
#   # boxplots
#   geom_boxplot(width = rain_height, show.legend = FALSE, 
#                outlier.shape = NA,
#                position = position_nudge(x = -rain_height*2)) +
#   #coord_flip() +
#   scale_fill_manual(values = official.species.colors) +
#   scale_color_manual(values = official.species.colors)  +
#   
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
# 
#   labs(x = NULL) +
#   ggtitle("Percent methylation\nof genes in turquoise module")
# pm.mean.plot.all.samples.turq


# Test if the eigenvalue of each module is significantly associated with ED50?

# Isolate matrix of average eigenvector values:
eigenvector.meth.df <- as.data.frame(MEs)
eigenvector.meth.df$seq.sample.id <- rownames(eigenvector.meth.df)
eigenvector.meth.df$species<- meta[match(eigenvector.meth.df$seq.sample.id, meta$col.spec.id),3]
eigenvector.meth.df$colony.id<- meta[match(eigenvector.meth.df$seq.sample.id, meta$col.spec.id),1]

# Add species level ED50 average:
# Read in ED50 data
ED50.data <- read.csv("~/Lab Notebook/Chapter3/Thermal_stress_analysis/Phylogeny resolved ED50 curves/individual.specied.ED50.csv", 
                      row.names = NULL)
eigenvector.meth.df$mean.ed50 <- ED50.data[match(eigenvector.meth.df$colony.id, ED50.data$colony.id),2]
eigenvector.meth.df$indv.ed50 <- ED50.data[match(eigenvector.meth.df$colony.id, ED50.data$colony.id),3]

yellow.plot <- ggplot(eigenvector.meth.df, aes(x = indv.ed50, y = MEyellow, color = species)) +
  geom_point(size = 3) +  # Adjust the size of the points
  labs(x = "indv.ed50", y = "MEyellow") + 
  scale_color_manual(values = official.species.colors) +# Label axes and legend
  theme_classic()   # Use a minimal theme
yellow.plot

blue.plot <- ggplot(eigenvector.meth.df, aes(x = indv.ed50, y = MEblue, color = species)) +
  geom_point(size = 3) +  # Adjust the size of the points
  labs(x = "indv.ed50", y = "MEblue") + 
  scale_color_manual(values = official.species.colors) +# Label axes and legend
  theme_classic()   # Use a minimal theme
blue.plot

brown.plot <- ggplot(eigenvector.meth.df, aes(x = indv.ed50, y = MEbrown)) +
  geom_point(aes(color = species), size = 3) +  # Adjust the size of the points and set color within aes
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "indv.ed50", y = "MEbrown") + 
  scale_color_manual(values = official.species.colors) + # Label axes and legend
  theme_classic()   # Use a minimal theme
brown.plot

# Perform linear regression
model <- lm(MEbrown ~ indv.ed50, data = eigenvector.meth.df)

# Get the summary of the model
model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- model_summary$coefficients[2, 4]
r_squared <- model_summary$r.squared

# Print the results
cat("P-value:", p_value, "\n")
cat("R-squared:", r_squared, "\n")

# Create the plot
brown.plot <- ggplot(eigenvector.meth.df, aes(x = indv.ed50, y = MEbrown)) +
  geom_point(aes(color = species), size = 3) +  # Adjust the size of the points and set color within aes
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "indv.ed50", y = "MEbrown", 
       title = paste("Linear Regression\nP-value:", round(p_value, 4), "R-squared:", round(r_squared, 4))) + 
  scale_color_manual(values = official.species.colors) + # Label axes and legend
  theme_classic()   # Use a minimal theme
brown.plot


## Yellow meth module
# Perform linear regression
model <- lm(MEyellow ~ indv.ed50, data = eigenvector.meth.df)

# Get the summary of the model
model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- model_summary$coefficients[2, 4]
r_squared <- model_summary$r.squared

# Print the results
cat("P-value:", p_value, "\n")
cat("R-squared:", r_squared, "\n")

# Create the plot
yellow.plot <- ggplot(eigenvector.meth.df, aes(x = indv.ed50, y = MEyellow)) +
  geom_point(aes(color = species), size = 3) +  # Adjust the size of the points and set color within aes
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "indv.ed50", y = "MEyellow", 
       title = paste("Linear Regression\nP-value:", round(p_value, 4), "R-squared:", round(r_squared, 4))) + 
  scale_color_manual(values = official.species.colors) + # Label axes and legend
  theme_classic()   # Use a minimal theme
yellow.plot


## blue meth module
# Perform linear regression
model <- lm(MEblue ~ indv.ed50, data = eigenvector.meth.df)

# Get the summary of the model
model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- model_summary$coefficients[2, 4]
r_squared <- model_summary$r.squared

# Print the results
cat("P-value:", p_value, "\n")
cat("R-squared:", r_squared, "\n")

# Create the plot
blue.plot <- ggplot(eigenvector.meth.df, aes(x = indv.ed50, y = MEblue)) +
  geom_point(aes(color = species), size = 3) +  # Adjust the size of the points and set color within aes
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "indv.ed50", y = "MEblue", 
       title = paste("Linear Regression\nP-value:", round(p_value, 4), "R-squared:", round(r_squared, 4))) + 
  scale_color_manual(values = official.species.colors) + # Label axes and legend
  theme_classic()   # Use a minimal theme
blue.plot


exp.and.perc.meth.abro.brown <- exp.and.perc.meth.abro %>%
  filter(gene %in% MEbrown.meth.genes)
exp.and.perc.meth.abro.brown$species <- "Acropora abrotanoides"
exp.and.perc.meth.abro.brown$pm.mean <- exp.and.perc.meth.abro.brown$pm.mean.abro
exp.and.perc.meth.abro.brown <- exp.and.perc.meth.abro.brown[,c(-9)]

exp.and.perc.meth.hya.brown <- exp.and.perc.meth.hya %>%
  filter(gene %in% MEbrown.meth.genes)
exp.and.perc.meth.hya.brown$species <- "Acropora hyacinthus"
exp.and.perc.meth.hya.brown$pm.mean <- exp.and.perc.meth.hya.brown$pm.mean.hya
exp.and.perc.meth.hya.brown <- exp.and.perc.meth.hya.brown[,c(-9)]

exp.and.perc.meth.lut.brown <- exp.and.perc.meth.lut %>%
  filter(gene %in% MEbrown.meth.genes)
exp.and.perc.meth.lut.brown$species <- "Acropora lutkeni"
exp.and.perc.meth.lut.brown$pm.mean <- exp.and.perc.meth.lut.brown$pm.mean.lut
exp.and.perc.meth.lut.brown <- exp.and.perc.meth.lut.brown[,c(-9)]

exp.and.perc.meth.nas.brown <- exp.and.perc.meth.nas %>%
  filter(gene %in% MEbrown.meth.genes)
exp.and.perc.meth.nas.brown$species <- "Acropora nasuta"
exp.and.perc.meth.nas.brown$pm.mean <- exp.and.perc.meth.nas.brown$pm.mean.nas
exp.and.perc.meth.nas.brown <- exp.and.perc.meth.nas.brown[,c(-9)]

exp.and.perc.meth.pul.brown <- exp.and.perc.meth.pul %>%
  filter(gene %in% MEbrown.meth.genes)
exp.and.perc.meth.pul.brown$species <- "Acropora pulchra"
exp.and.perc.meth.pul.brown$pm.mean <- exp.and.perc.meth.pul.brown$pm.mean.pul
exp.and.perc.meth.pul.brown <- exp.and.perc.meth.pul.brown[,c(-9)]

exp.and.perc.meth.ret.brown <- exp.and.perc.meth.ret %>%
  filter(gene %in% MEbrown.meth.genes)
exp.and.perc.meth.ret.brown$species <- "Acropora retusa"
exp.and.perc.meth.ret.brown$pm.mean <- exp.and.perc.meth.ret.brown$pm.mean.ret
exp.and.perc.meth.ret.brown <- exp.and.perc.meth.ret.brown[,c(-9)]

exp.and.perc.meth.rob.brown <- exp.and.perc.meth.rob %>%
  filter(gene %in% MEbrown.meth.genes)
exp.and.perc.meth.rob.brown$species <- "Acropora robusta"
exp.and.perc.meth.rob.brown$pm.mean <- exp.and.perc.meth.rob.brown$pm.mean.rob
exp.and.perc.meth.rob.brown <- exp.and.perc.meth.rob.brown[,c(-9)]

exp.and.perc.meth.tut.brown <- exp.and.perc.meth.tut %>%
  filter(gene %in% MEbrown.meth.genes)
exp.and.perc.meth.tut.brown$species <- "Acropora tutuilensis"
exp.and.perc.meth.tut.brown$pm.mean <- exp.and.perc.meth.tut.brown$pm.mean.tut
exp.and.perc.meth.tut.brown <- exp.and.perc.meth.tut.brown[,c(-9)]


all.species.pm.means.brown<- rbind(exp.and.perc.meth.tut.brown,
                                       exp.and.perc.meth.rob.brown,
                                       exp.and.perc.meth.ret.brown,
                                       exp.and.perc.meth.pul.brown,
                                       exp.and.perc.meth.nas.brown, 
                                       exp.and.perc.meth.lut.brown, 
                                       exp.and.perc.meth.hya.brown, 
                                       exp.and.perc.meth.abro.brown)


# Test if brown module gene percent methylation is significantly associated with basemean exp of the brown module genes?
# Perform linear regression
model <- lm(log10_baseMean ~ pm.mean, data = all.species.pm.means.brown)

# Get the summary of the model
model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- model_summary$coefficients[2, 4]
r_squared <- model_summary$r.squared

# Print the results
cat("P-value:", p_value, "\n")
cat("R-squared:", r_squared, "\n")

# Create the plot
basemean.pm.meth.brown <- ggplot(all.species.pm.means.brown, aes(x = pm.mean, y = log10_baseMean)) +
  geom_point(aes(color = species), size = 3, alpha = 0.7) +  # Adjust the size of the points and set color within aes
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "pm.mean", y = "log(basemean)", 
       title = paste("Linear Regression\nP-value:", round(p_value, 4), "R-squared:", round(r_squared, 4))) + 
  scale_color_manual(values = official.species.colors) + # Label axes and legend
  theme_classic()   # Use a minimal theme
basemean.pm.meth.brown # not much different than looking at all the genes


# Test if brown module gene percent methylation is significantly associated with the log2FC of HRGs?
# Perform linear regression
model <- lm(abs(log2FoldChange) ~ pm.mean, data = all.species.pm.means.brown)

# Get the summary of the model
model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- model_summary$coefficients[2, 4]
r_squared <- model_summary$r.squared

# Print the results
cat("P-value:", p_value, "\n")
cat("R-squared:", r_squared, "\n")

# Create the plot
l2fc.pm.meth.brown <- ggplot(all.species.pm.means.brown, aes(x = pm.mean, y = abs(log2FoldChange))) +
  geom_point(aes(color = species), size = 3, alpha = 0.7) +  # Adjust the size of the points and set color within aes
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "pm.mean", y = "log2FoldChange", 
       title = paste("Linear Regression\nP-value:", round(p_value, 4), "R-squared:", round(r_squared, 4))) + 
  scale_color_manual(values = official.species.colors) + # Label axes and legend
  theme_classic()   # Use a minimal theme
l2fc.pm.meth.brown # not much different than looking at all the genes

all.species.pm.means.brown$species <- factor(all.species.pm.means.brown$species, levels = species.order)

# Plot
rain_height <- .1

l2fc.plot.all.samples.brown<- ggplot(all.species.pm.means.brown, aes(x = species, y = log2FoldChange, fill = species )) +
  # # clouds
  # introdataviz::geom_flat_violin(trim=FALSE,
  #                                position = position_nudge(x = rain_height+.05),
  #                                width = 1) +
  # # rain
  geom_point(aes(color = species), size = 1, alpha = 0.4,  show.legend = FALSE, 
             position = position_jitter(width = rain_height, height = 0)) +
  # boxplots
  geom_boxplot(width = rain_height, show.legend = FALSE, 
               outlier.shape = NA,
               position = position_nudge(x = -rain_height*2)) +
  #coord_flip() +
  scale_fill_manual(values = official.species.colors) +
  scale_color_manual(values = official.species.colors)  +
  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  
  # stat_compare_means(method = "wilcox.test",
  #                    comparisons = list(c("TRUE", "FALSE")),
  #                    label = "p.format",
  #                    step.increase = 0.05) +
  theme(legend.position = "none") +
  labs(x = NULL) +
  ggtitle("Heat-stress Log2FC\nof genes in brown meth Module")
l2fc.plot.all.samples.brown

# overlap between salmon and yellow expression modules and the brown methylation module?
combined.genes <- union(salmon.genes$x, yellow.genes$x)
all.species.pm.means.brown.salm.yello.exp.mod <- all.species.pm.means.brown[all.species.pm.means.brown$gene %in% combined.genes, ]

# Perform linear regression
model <- lm(log10_baseMean ~ pm.mean, data = all.species.pm.means.brown.salm.yello.exp.mod )

# Get the summary of the model
model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- model_summary$coefficients[2, 4]
r_squared <- model_summary$r.squared

# Print the results
cat("P-value:", p_value, "\n")
cat("R-squared:", r_squared, "\n")

# Create the plot
basemean.pm.meth.brown <- ggplot(all.species.pm.means.brown.salm.yello.exp.mod, aes(x = pm.mean, y = log10_baseMean)) +
  geom_point(aes(color = species), size = 3, alpha = 0.7) +  # Adjust the size of the points and set color within aes
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "pm.mean", y = "log(basemean)", 
       title = paste("Linear Regression\nP-value:", round(p_value, 4), "R-squared:", round(r_squared, 4))) + 
  scale_color_manual(values = official.species.colors) + # Label axes and legend
  theme_classic()   # Use a minimal theme
basemean.pm.meth.brown

# Perform linear regression
model <- lm(abs(log2FoldChange) ~ pm.mean, data = all.species.pm.means.brown.salm.yello.exp.mod)

# Get the summary of the model
model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- model_summary$coefficients[2, 4]
r_squared <- model_summary$r.squared

# Print the results
cat("P-value:", p_value, "\n")
cat("R-squared:", r_squared, "\n")

# Create the plot
l2fc.pm.meth.brown <- ggplot(all.species.pm.means.brown.salm.yello.exp.mod, aes(x = pm.mean, y = abs(log2FoldChange))) +
  geom_point(aes(color = species), size = 3, alpha = 0.7) +  # Adjust the size of the points and set color within aes
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "pm.mean", y = "log2FoldChange", 
       title = paste("Linear Regression\nP-value:", round(p_value, 4), "R-squared:", round(r_squared, 4))) + 
  scale_color_manual(values = official.species.colors) + # Label axes and legend
  theme_classic()   # Use a minimal theme
l2fc.pm.meth.brown # not much different than looking at all the genes

l2fc.plot.all.samples.brown<- ggplot(all.species.pm.means.brown.salm.yello.exp.mod, aes(x = species, y = log2FoldChange, fill = species )) +
  # # clouds
  # introdataviz::geom_flat_violin(trim=FALSE,
  #                                position = position_nudge(x = rain_height+.05),
  #                                width = 1) +
  # # rain
  geom_point(aes(color = species), size = 1, alpha = 0.4,  show.legend = FALSE, 
             position = position_jitter(width = rain_height, height = 0)) +
  # boxplots
  geom_boxplot(width = rain_height, show.legend = FALSE, 
               outlier.shape = NA,
               position = position_nudge(x = -rain_height*2)) +
  #coord_flip() +
  scale_fill_manual(values = official.species.colors) +
  scale_color_manual(values = official.species.colors)  +
  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  
  # stat_compare_means(method = "wilcox.test",
  #                    comparisons = list(c("TRUE", "FALSE")),
  #                    label = "p.format",
  #                    step.increase = 0.05) +
  theme(legend.position = "none") +
  labs(x = NULL) +
  ggtitle("Heat-stress Log2FC\nof genes in brown meth Module")
l2fc.plot.all.samples.brown


# percent methylation of genes in yellow and salmon modules

  






# topGO of the modules




