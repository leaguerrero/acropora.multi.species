# Chapter 3: WGCNA
# WGCNA is designed to identify modules or groups of genes that
# are co-expressed, meaning their expression levels tend to
# change in a coordinated manner across samples.
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.18")

# BiocManager::install("preprocessCore")
# BiocManager::install("WGCNA", force = TRUE)

library(WGCNA)
library(DESeq2)
library(tidyr)
library(tidyverse)
library(broom)
library(purrr)
options(stringsAsFactors = FALSE)

#Enable multithread
enableWGCNAThreads()

# Load the expression data
# File Set-up
my.dir <- "/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/HTseq_out"
my.files <- grep(".txt", list.files(my.dir), value=TRUE)
my.metadata <- read.csv("/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/tagseq_meta_data.csv", header = TRUE)

my.metadata$thermal.treatment <-ifelse(my.metadata$treatment == "control", 29, 
                                ifelse(my.metadata$treatment == "heat", 36, NA)) 

my.metadata$species.class <-as.numeric(match(my.metadata$species, 
                                             unique(my.metadata$species)))


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

acro.counts.matrix <- ddsHTSeq@assays@data@listData$counts

# Transform matrix
acro.counts.mat <- varianceStabilizingTransformation(acro.counts.matrix , blind = TRUE)
# write.csv(acro.counts.mat, file = "~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/acro.counts.matrix.csv")
# Transpose the expression data
counts.df <- as.data.frame(acro.counts.mat)
counts.df <- tibble::rownames_to_column(counts.df, "Genes")
datExpr.0 = as.data.frame(t(counts.df[, -c(1)]))
names(datExpr.0) = counts.df$Genes;
rownames(datExpr.0) = names(counts.df)[-c(1)]
sample.names <- names(counts.df)[-c(1)]

gsg.1 = goodSamplesGenes(datExpr.0, verbose = 3)
gsg.1$allOK
rm(gsg.1)
# Check Sample Outliters
sampleTree.1 = hclust(dist(datExpr.0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
#pdf(file = "~/Lab Notebook/Chapter3/Tag_seq_analysis/Figures/sampleClustering.pdf", width = 12, height = 9);
#par(cex = 0.6);
#par(mar = c(0,4,2,0))
plot(sampleTree.1, main = "Sample clustering to detect outliers",
     sub="",
     xlab="",
     cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 250, col = "red")
# dev.off()
# # No outliers
datExpr <- datExpr.0
# 
# 
# ######################################
# traitData = my.metadata
# dim(traitData)
# names(traitData)
# 
# 
# # Form a data frame analogous to expression data that will hold the clinical traits.
# rna.Samples = rownames(datExpr)
# traitRows = match(rna.Samples, traitData$seq.sample.id)
# datTraits = traitData[traitRows, -1]
# rownames(datTraits) = traitData[traitRows, 1]
# collectGarbage()
# 
# sampleTree.2 = hclust(dist(datExpr), method = "average")
# # Convert traits to a color representation: white means low, red means high, grey means missing entry
# traitColors = numbers2colors(datTraits[,c(5,6)], signed = FALSE);
# # Plot the sample dendrogram and the colors underneath.
# pdf(file = "~/Lab Notebook/Chapter3/Tag_seq_analysis/Figures/preliminary.acro.dendro.pdf", wi = 9, he = 6)
# plotDendroAndColors(sampleTree.2, traitColors,
#                     groupLabels = c("Thermal treatment", "Species ID"),
#                     main = "Sample dendrogram and trait heatmap")
# dev.off()
# 
# Picking beta for expression data
# Choose a set of soft-thresholding powers
powers.e = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft.e = pickSoftThreshold(datExpr, powerVector = powers.e, verbose = 5)
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

# # based on plot, I will pick a threshold of 5
# 
# # Constructing the gene network
# ge.network= blockwiseModules(datExpr, power = 5,
#                              TOMType = "unsigned", minModuleSize = 30,
#                              reassignThreshold = 0, mergeCutHeight = 0.25,
#                              numericLabels = TRUE, pamRespectsDendro = FALSE,
#                              saveTOMs = FALSE,
#                              verbose = 3)
# 
# # Restart R Session
# # library(WGCNA)
# # search(): ".GlobalEnv", "package:WGCNA" 
# table(ge.network$colors)                         
# 
# 
# sizeGrWindow(12, 9)
# # Convert labels to colors for plotting
# mergedColors = labels2colors(ge.network$colors)
# # Plot the dendrogram and the module colors underneath
# plotDendroAndColors(ge.network$dendrograms[[1]], mergedColors[ge.network$blockGenes[[1]]],
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)                        
# 
# moduleLabels.e = ge.network$colors
# moduleColors.e = labels2colors(ge.network$colors)
# ME.es = ge.network$MEs;
# geneTree.e = ge.network$dendrograms[[1]]                        
# 
# softPower = 5 #Chosen in the graphs before
# adjacency = adjacency(datExpr, power = softPower, type = "unsigned") 
# #Transforming the adjacency matrix in a topological overlap
# TOM = TOMsimilarity(adjacency) #Calculating the topological overlap matrix
# # rm(adjacency)
# 
# dissTOM = 1-TOM 
# 
# # rm(TOM)
# geneTree = hclust(as.dist(dissTOM), method = "average");
# # Plot the resulting clustering tree (dendrogram)
# sizeGrWindow(12,9)
# plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
#      labels = FALSE, hang = 0.04);
# 
# # We like large modules, so we set the minimum module size relatively high:
# minModuleSize = 200;
# # Module identification using dynamic tree cut:
# dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
#                             deepSplit = 2, pamRespectsDendro = FALSE,
#                             minClusterSize = minModuleSize);
# table(dynamicMods)
# 
# 
# # Convert numeric labels into colors
# dynamicColors = labels2colors(dynamicMods)
# table(dynamicColors)
# # Plot the dendrogram and colors underneath
# sizeGrWindow(8,6)
# plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05,
#                     main = "Gene dendrogram and module colors")
# ## 2.b.5 Merging of modules whose expression profiles are very similar
# 
# # Calculate eigengenes
# MEList = moduleEigengenes(datExpr, colors = dynamicColors)
# MEs = MEList$eigengenes
# # Calculate dissimilarity of module eigengenes
# MEDiss = 1-cor(MEs);
# # Cluster module eigengenes
# METree = hclust(as.dist(MEDiss), method = "average");
# # Plot the result
# sizeGrWindow(7, 6)
# plot(METree, main = "Clustering of module eigengenes",
#      xlab = "", sub = "")
# MEDissThres = 0.1
# # Plot the cut line into the dendrogram
# abline(h=MEDissThres, col = "red")
# # Call an automatic merging function
# merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# # The merged module colors
# mergedColors = merge$colors;
# # Eigengenes of the new merged modules:
# mergedMEs = merge$newMEs;
# 
# 
# sizeGrWindow(12, 9)
# #pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
# plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
#                     c("Dynamic Tree Cut", "Merged dynamic"),
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# #dev.off()
# 
# moduleColors = mergedColors
# 
# #Building numeric labels corresponding to the colors
# colorOrder = c("grey", standardColors(50))
# moduleLabels = match(moduleColors, colorOrder)-1
# MEs = mergedMEs
# 
# 
# 
# # Associating Modules and Phenotypes
# #Relating modules to characteristics and identifying important genes
# #Defining the number of genes and samples
# nGenes = ncol(datExpr)
# nSamples = nrow(datExpr)
# 
# #Recalculating MEs with label colors
# MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
# MEs = orderMEs(MEs0)
# moduleTraitCor = cor(MEs,datTraits[,c(5,6)], use = "p")
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# 
# 
# #sizeGrWindow(8,4)
# 
# #Displaying correlations and its p-values
# textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
#                     signif(moduleTraitPvalue, 1), ")", sep = "")
# dim(textMatrix) = dim(moduleTraitCor)
# par(mar = c(6, 8, 4, 6))
# 
# #Displaying the correlation values in a heatmap plot
# pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
# labeledHeatmap(Matrix = moduleTraitCor,
#                xLabels = c("thermal treatment", "species"),
#                yLabels = names(MEs),
#                ySymbols = names(MEs),
#                colorLabels = FALSE,
#                colors = blueWhiteOrange(50),
#                textMatrix = textMatrix,
#                setStdMargins = FALSE,
#                cex.text = 0.5,
#                zlim = c(-1,1),
#                main = paste("Module-trait relationships"))

#############################################################################
# Species and heatstress?
# So, have one column be a.abro.hs and a.abro.control, etc for all species?
# Module-sample relationship 

traitData = my.metadata
dim(traitData)
names(traitData)
traitData.wide <- traitData %>%
  pivot_wider(
    id_cols = c(seq.sample.id, colony.id,file.prefix),
    names_from = c(species, treatment),
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
ordered.col.names <- c("seq.sample.id",
                       "colony.id",
                       "file.prefix", 
                       "Acropora.tutuilensis_control",
                       "Acropora.tutuilensis_heat",
                       "Acropora.lutkeni_control",
                       "Acropora.lutkeni_heat",
                       "Acropora.pulchra_control",
                       "Acropora.pulchra_heat",
                       "Acropora.retusa_control",
                       "Acropora.retusa_heat", 
                       "Acropora.nasuta_control",
                       "Acropora.nasuta_heat",
                       "Acropora.robusta_control",
                       "Acropora.robusta_heat",
                       "Acropora.hyacinthus_control",
                       "Acropora.hyacinthus_heat",
                       "Acropora.abrotanoides_control",
                       "Acropora.abrotanoides_heat")
traitData.wide <- traitData.wide[,ordered.col.names]

# Form a data frame analogous to expression data that will hold the clinical traits.
rna.Samples = rownames(datExpr)
traitRows = match(rna.Samples, traitData.wide$seq.sample.id)
datTraits = as.data.frame(traitData.wide[traitRows, -1])
traitDataSubset <- as.vector(traitData.wide[traitRows, 1])
rownames(datTraits) <- traitDataSubset$seq.sample.id
collectGarbage()

sampleTree.2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits[,c(5,6)], signed = TRUE);
# Plot the sample dendrogram and the colors underneath.
#pdf(file = "~/Lab Notebook/Chapter3/Tag_seq_analysis/Figures/preliminary.acro.dendro.pdf", wi = 9, he = 6)
#plotDendroAndColors(sampleTree.2, traitColors,
 #                   groupLabels = c("Thermal treatment", "Species ID"),
  #                  main = "Sample dendrogram and trait heatmap")
#dev.off()

# Picking beta for expression data
# Choose a set of soft-thresholding powers
powers.e = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft.e = pickSoftThreshold(datExpr, powerVector = powers.e, verbose = 5)
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

# based on plot, I will pick a threshold of 5

# Constructing the gene network
ge.network= blockwiseModules(datExpr, power = 5,
                             TOMType = "unsigned", minModuleSize = 30, # check signed vs unsigned and minModuleSize
                             reassignThreshold = 0, mergeCutHeight = 0.25,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             saveTOMs = FALSE,
                             verbose = 3)
table(ge.network$colors)                         


sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(ge.network$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(ge.network$dendrograms[[1]], mergedColors[ge.network$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)                        

moduleLabels.e = ge.network$colors
moduleColors.e = labels2colors(ge.network$colors)
ME.es = ge.network$MEs;
geneTree.e = ge.network$dendrograms[[1]]                        

softPower = 5 #Chosen in the graphs before
adjacency = adjacency(datExpr, power = softPower, type = "unsigned") # unsigned matrix
#Transforming the adjacency matrix in a topological overlap
TOM = TOMsimilarity(adjacency) #Calculating the topological overlap matrix
dissTOM = 1-TOM 

# Clean up vector space
rm(acro.counts.mat, acro.counts.matrix, counts.df, ddsHTSeq,adjacency)

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
## 2.b.5 Merging of modules whose expression profiles are very similar

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
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
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

moduleColors = mergedColors

#Building numeric labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs



# Associating Modules and Phenotypes
#Relating modules to characteristics and identifying important genes
#Defining the number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculating MEs with label colors
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs,datTraits[,-c(1,2)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


#sizeGrWindow(8,4)

#Displaying correlations and its p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 8, 6, 8))
par(oma = c(8, 0, 4, 0))  # Adjust the outer margin

#Displaying the correlation values in a heatmap plot
#pdf(file = "~/Lab Notebook/Chapter3/Tag_seq_analysis/Figures/module_sample_relationship.pdf", wi = 15, he = 8.5)
#par(mar=c(9,8,4,1)+.1)
#my_palette <- colorRampPalette(c('#FF0000',"white", '#0000FF'))(n = 50)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits[,-c(1,2)]),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = my_palette,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-sample relationships"))
dev.off()

# Looking at the distribution of the correlation values between the genes 
# of the module and the character

#Defining the variable "Acropora.tutuilensis_control" containing the column "Acropora.tutuilensis_control" of datTrait
Acropora.tutuilensis_control = as.data.frame(datTraits$Acropora.tutuilensis_control)
names(Acropora.tutuilensis_control) = "Acropora.tutuilensis_control"

#names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr, Acropora.tutuilensis_control, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(Acropora.tutuilensis_control), sep="")
names(GSPvalue) = paste("p.GS.", names(Acropora.tutuilensis_control), sep="")

module = "white" #########################putting the color below the plot
column = match(module, modNames)
moduleGenes = moduleColors==module

#sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for \n Acropora.tutuilensis_control",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#Display the gene names inside the module
# colnames(datExpr)[moduleColors=="white"] 
length(colnames(datExpr)[moduleColors=="white"]) #147
length(colnames(datExpr)[moduleColors=="royalblue"]) #237
length(colnames(datExpr)[moduleColors=="skyblue"]) #142
#Identifying most important genes for one determined characteristic inside of the cluster
geneInfo0 = data.frame(EST = colnames(datExpr),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)


modOrder = order(-abs(cor(MEs, Acropora.tutuilensis_control, use = "p")))
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Acropora.tutuilensis_control))
geneInfo = geneInfo0[geneOrder, ]
#if you want to write the information in a csv file, just uncomment line below
#write.csv(geneInfo, file = "~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/A.tut.control.GeneInfo.csv")


## Tutorial: https://rstudio-pubs-static.s3.amazonaws.com/687551_ed469310d8ea4652991a2e850b0018de.html

## Testing correlations between modules and species level ED50:

# Isolate matrix of average eigenvector values:
eigenvector.df <- as.data.frame(MEs)
eigenvector.df$seq.sample.id <- rownames(eigenvector.df)
eigenvector.df<- merge(eigenvector.df, my.metadata[, c("seq.sample.id", "species","colony.id", "treatment")], all.x = TRUE)

# Add species level ED50 average:
# Read in ED50 data
ED50.data <- read.csv("~/Lab Notebook/Chapter3/Thermal_stress_analysis/Phylogeny resolved ED50 curves/individual.specied.ED50.csv", 
                      row.names = NULL)
eigenvector.df$mean.ed50 <- ED50.data[match(eigenvector.df$colony.id, ED50.data$colony.id),2]
# Need to create a difference matrix: control eigenvector value - heat eigenvector value
## Control and heat eigenvector DF
control.eigenvector.df <-subset(eigenvector.df, treatment == "control")
heat.eigenvector.df <- subset(eigenvector.df, treatment == "heat")

# Difference in the eigenvector values
cols.to.subtract <- colnames(MEs)

# Subtract corresponding elements
difference.ev.df <- control.eigenvector.df[,-33]
difference.ev.df[, cols.to.subtract] <- control.eigenvector.df[, cols.to.subtract] - heat.eigenvector.df[, cols.to.subtract]


# Subtract the absolute values
abs.control.eigenvector.df<-abs(control.eigenvector.df[,cols.to.subtract])
abs.heat.eigenvector.df<-abs(heat.eigenvector.df[,cols.to.subtract])
difference.av.ev.df <- control.eigenvector.df[,-33]
difference.av.ev.df[, cols.to.subtract] <- abs.control.eigenvector.df[, cols.to.subtract] - abs.heat.eigenvector.df[, cols.to.subtract]

# Test the relationship between species level changes in gene expression and ed50
# lm(MEcyan.diff ~ ED50 + Species)
# 
# ^^ by individual
# 
# cor.test(MEcyan.speciesavg.diff ~ ED50.speciesavg)

# create data frame with avg species difference of each module
difference.av.ev.df$mean.ed50 <- as.factor(difference.av.ev.df$mean.ed50)
ME.cols <- difference.av.ev.df[,2:31]
species.avg.ev.df <- aggregate(ME.cols, by = list(species = ME.cols$species), FUN = mean)
species.avg.ev.df <- species.avg.ev.df[,-31]
species.avg.ev.df$mean.ed50 <- as.numeric(as.character(difference.av.ev.df[match(species.avg.ev.df$species, difference.av.ev.df$species),33]))
# write.csv(species.avg.ev.df,
#           file = "~/Lab Notebook/Chapter3/Thermal_stress_analysis/species.avg.ev.df.csv",
#           quote = FALSE,
#           row.names = FALSE)

# Loop through correlations:
columns.to.correlate <- cols.to.subtract
# Initialize an empty data frame to store results
correlation_results_df <- data.frame(Column_Name = character(),
                                     Correlation_Coefficient = numeric(),
                                     P_Value = numeric(),
                                     stringsAsFactors = FALSE)

# Loop through each column combination
for (col in columns.to.correlate) {
  correlation_result <- cor.test(species.avg.ev.df[[col]], species.avg.ev.df$mean.ed50, method = "pearson")
  
  # Extract correlation coefficient and p-value
  corr_coef <- correlation_result$estimate
  p_value <- correlation_result$p.value
  
  # Append the results to the data frame
  correlation_results_df <- rbind(correlation_results_df, 
                                  data.frame(Column_Name = col,
                                             Correlation_Coefficient = corr_coef,
                                             P_Value = p_value,
                                             stringsAsFactors = FALSE))
}

bubble_plot <- ggplot(correlation_results_df, aes(x = Column_Name, y = Correlation_Coefficient, size = P_Value, alpha = ifelse(P_Value < 0.05, 0.7, 0.3))) +
  geom_point(color = "#219ebc") +
  scale_size_continuous(range = c(2,10)) +  # Adjust the range of bubble sizes
  labs(x = "Module", y = "Correlation Coefficient between \n eigenvalue difference and species level ED50", size = "P-Value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility
bubble_plot

# one plot for each sig module that plots eigen difference v. Ed50 colored by species
sig.eigen.diff.ed50.species <- species.avg.ev.df[,c("species","MEroyalblue", "MEskyblue", "MEwhite", "mean.ed50")]
royal.blue.plot <- ggplot(sig.eigen.diff.ed50.species, aes(x = mean.ed50, y = MEroyalblue, color = species)) +
  geom_point(size = 3) +  # Adjust the size of the points
  labs(x = "mean.ed50", y = "MEroyalblue", color = "Species") +  # Label axes and legend
  theme_classic()   # Use a minimal theme
  royal.blue.plot

 sky.blue.plot <- ggplot(sig.eigen.diff.ed50.species, aes(x = mean.ed50, y = MEskyblue)) +
    geom_point(size = 3, aes(color = species)) +
   geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence interval
    labs(x = "mean.ed50", y = "MEroyalblue", color = "Species") +  # Label axes and legend
    theme_classic()   # Use a minimal theme
  sky.blue.plot  

  white_plot <- ggplot(sig.eigen.diff.ed50.species, aes(x = mean.ed50, y = MEwhite)) +
    geom_point(size = 4, aes(color = species)) +  # Adjust the size of the points
    geom_smooth(method = "lm", se = TRUE) +  # Add a linear regression line without confidence intervals
    labs(x = "Species Average ED50", y = "MEwhite", color = "Species") +  # Label axes and legend
    scale_color_discrete(labels = c("A. abrotanoides",
                                    "A. hyacinthus",
                                    "A. lutkeni",
                                    "A. pulchra",
                                    "A. retusa",
                                    "A. robusta",
                                    "A. tutuilensis",
                                    "A. spp")) + 
    theme_classic()+   # Use a classic theme
    theme(axis.text = element_text(size = 14),  
          axis.title = element_text(size = 16))
  
  white_plot
  # How many genes in the white module:
  length(colnames(datExpr)[moduleColors=="white"])
  # Gene names:
  genes_in_white_module <- colnames(datExpr)[moduleColors=="white"]
  # Calculate the linear regression model
  lm_model <- lm(MEwhite ~ mean.ed50, data = sig.eigen.diff.ed50.species)
  
  # Extract the p-value for the slope coefficient
  p_value <- summary(lm_model)$coefficients[2, 4]
  
  # Add the p-value to the plot
  white_plot <- white_plot + annotate("text", x = max(sig.eigen.diff.ed50.species$mean.ed50), 
                                      y = max(sig.eigen.diff.ed50.species$MEwhite), 
                                      label = paste("p-value =", round(p_value, 4)), 
                                      hjust = 1, vjust = 1, size = 6) 
  
  white_plot
  ggsave(filename = "whitemodule_vs_speciesED50_v2.jpg", plot = white_plot, width = 10, height = 6, dpi = 300)
  
  
  
skyblue_plot <- ggplot(sig.eigen.diff.ed50.species, aes(x = mean.ed50, y = MEskyblue)) +
    geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
    geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
    labs(x = "mean.ed50", y = "MEskyblue", color = "Species") +  # Label axes and legend
    theme_classic()   # Use a classic theme
  
  # Calculate the linear regression model
  lm_model <- lm(MEskyblue ~ mean.ed50, data = sig.eigen.diff.ed50.species)
  
  # Extract the p-value for the slope coefficient
  p_value <- summary(lm_model)$coefficients[2, 4]
  
  # Add the p-value to the plot
  skyblue_plot <- skyblue_plot + annotate("text", x = max(sig.eigen.diff.ed50.species$mean.ed50), 
                                      y = max(sig.eigen.diff.ed50.species$MEskyblue), 
                                      label = paste("p-value =", round(p_value, 4)), 
                                      hjust = 1, vjust = 1, size = 4) 
  
  skyblue_plot
  
  
  royalblue_plot <- ggplot(sig.eigen.diff.ed50.species, aes(x = mean.ed50, y = MEroyalblue)) +
    geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
    geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
    labs(x = "mean.ed50", y = "MEroyalblue", color = "Species") +  # Label axes and legend
    theme_classic()   # Use a classic theme
  
  # Calculate the linear regression model
  lm_model <- lm(MEroyalblue ~ mean.ed50, data = sig.eigen.diff.ed50.species)
  
  # Extract the p-value for the slope coefficient
  p_value <- summary(lm_model)$coefficients[2, 4]
  
  # Add the p-value to the plot
  royalblue_plot <- royalblue_plot + annotate("text", x = max(sig.eigen.diff.ed50.species$mean.ed50), 
                                          y = max(sig.eigen.diff.ed50.species$MEroyalblue), 
                                          label = paste("p-value =", round(p_value, 4)), 
                                          hjust = 1, vjust = 1, size = 4) 
  
  royalblue_plot
###### Tangent but doesn't touch on 
# the hypothesis of interest
# # Visualize correlation matrix
# # install.packages("corrplot")
# library(corrplot)
# 
# # compute correlation matrix:
# # Make dataframe for computing correlation matrix
# species.avg.ev.df.for.correlations <- species.avg.ev.df
# row.names(species.avg.ev.df.for.correlations) <- species.avg.ev.df.for.correlations$mean.ed50
# ed50.eigen.dif.res <- cor(species.avg.ev.df.for.correlations[,-1])
# round(ed50.eigen.dif.res, 3)
#####

# Species level
cor.test(species.avg.ev.df$MEcyan, species.avg.ev.df$mean.ed50,  method = "pearson")

# By individual:
library(performance)
library(qqplotr)

# lm(MEcyan.diff ~ ED50 + Species)
difference.av.ev.df$ind.ed50 <- ED50.data[match(difference.av.ev.df$colony.id, ED50.data$colony.id),3]
# write.csv(difference.av.ev.df,
#           file = "~/Lab Notebook/Chapter3/Thermal_stress_analysis/Diff_eigenvalues.csv",
#           quote = FALSE,
#           row.names = FALSE)


cyan.model <- lm(difference.av.ev.df$MEcyan ~ difference.av.ev.df$ind.ed50+ difference.av.ev.df$species)
summary(cyan.model)
coef(summary(cyan.model))
plot(check_model(cyan.model))
anova(cyan.model)
summary(cyan.model)$coefficients[2, 1]
summary(cyan.model)$coefficients[2, 4]

# Loop to test the relationship between the modules and individual ed.50
# Initialize a list to store the models

###### Code Needs assessment
coef_list <- list()
coef_p_value_list <- list()


# Loop through each response variable
for (response_var in columns.to.correlate) {
  # Perform linear regression
  model <- lm(paste(response_var, " ~ ind.ed50 + species"), data = difference.av.ev.df)
  
  # Extract coefficients and their p-values
  coef_summary <- summary(model)
  coef_list[[response_var]] <- coef_summary$coefficients[2, 1]  # Coefficient
  coef_p_value_list[[response_var]] <- coef_summary$coefficients[2, 4]  # P-value for coefficient
  
  # Extract p-value for the model
  model_summary <- anova(model)
  model_p_value_list[[response_var]] <- coef_summary$fstatistic["Pr(>F)"]  # P-value for the model
}


# Create a data frame from the results
results_df <- data.frame(Response_Variable = columns.to.correlate,
                         Coefficient = unlist(coef_list),
                         Coefficient_P_Value = unlist(coef_p_value_list))

results_df

# Plot results
ind.ed50.module.relationship <- ggplot(data = results_df, aes(x = Response_Variable, y = Coefficient, size = Coefficient_P_Value, alpha = ifelse(Coefficient_P_Value < 0.05, 0.7, 0.2))) +
  geom_point(shape = 21, fill = "#7826e3") +  # Use shape 21 for filled circles
  scale_size_continuous(range = c(3, 10)) +  # Adjust the range of bubble sizes
  scale_alpha_continuous(range = c(0.3, 0.9)) +  # Adjust the range of alpha values
  labs(x = "Module", y = "Coefficient", size = "P-Value", alpha = "P-Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility

ind.ed50.module.relationship

## Looking at each module in each individual colored by species
difference.av.ev.df
sig.eigen.diff.ed50.ind <- difference.av.ev.df[,c("seq.sample.id","MEbrown", "MEmidnightblue", "MEwhite","MEsalmon","MEturquoise","MEyellow" ,"species", "ind.ed50")]

brown_plot <- ggplot(sig.eigen.diff.ed50.ind, aes(x = ind.ed50, y = MEbrown)) +
  geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "ind.ed50", y = "MEbrown", color = "Species") +  # Label axes and legend
  theme_classic()   # Use a classic theme

# Calculate the linear regression model
lm_model <- lm(MEbrown ~ ind.ed50, data = sig.eigen.diff.ed50.ind)

# Extract the p-value for the slope coefficient
p_value <- summary(lm_model)$coefficients[2, 4]

# Add the p-value to the plot
brown_plot <- brown_plot + annotate("text", x = max(sig.eigen.diff.ed50.ind$ind.ed50), 
                                        y = max(sig.eigen.diff.ed50.ind$MEbrown), 
                                        label = paste("p-value =", round(p_value, 4)), 
                                        hjust = 1, vjust = 1, size = 4) 

brown_plot


midnightblue_plot <- ggplot(sig.eigen.diff.ed50.ind, aes(x = ind.ed50, y = MEmidnightblue)) +
  geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "ind.ed50", y = "MEmidnightblue", color = "Species") +  # Label axes and legend
  theme_classic()   # Use a classic theme

# Calculate the linear regression model
lm_model <- lm(MEmidnightblue ~ ind.ed50, data = sig.eigen.diff.ed50.ind)

# Extract the p-value for the slope coefficient
p_value <- summary(lm_model)$coefficients[2, 4]

# Add the p-value to the plot
midnightblue_plot <- midnightblue_plot + annotate("text", x = max(sig.eigen.diff.ed50.ind$ind.ed50), 
                                    y = max(sig.eigen.diff.ed50.ind$MEmidnightblue), 
                                    label = paste("p-value =", round(p_value, 4)), 
                                    hjust = 1, vjust = 1, size = 4) 

midnightblue_plot


salmon_plot <- ggplot(sig.eigen.diff.ed50.ind, aes(x = ind.ed50, y = MEsalmon)) +
  geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "ind.ed50", y = "MEsalmon", color = "Species") +  # Label axes and legend
  theme_classic()   # Use a classic theme

# Calculate the linear regression model
lm_model <- lm(MEsalmon ~ ind.ed50, data = sig.eigen.diff.ed50.ind)

# Extract the p-value for the slope coefficient
p_value <- summary(lm_model)$coefficients[2, 4]

# Add the p-value to the plot
salmon_plot <- salmon_plot + annotate("text", x = max(sig.eigen.diff.ed50.ind$ind.ed50), 
                                                  y = max(sig.eigen.diff.ed50.ind$MEsalmon), 
                                                  label = paste("p-value =", round(p_value, 4)), 
                                                  hjust = 1, vjust = 1, size = 4) 

salmon_plot


white_plot <- ggplot(sig.eigen.diff.ed50.ind, aes(x = ind.ed50, y = MEwhite)) +
  geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "ind.ed50", y = "MEwhite", color = "Species") +  # Label axes and legend
  theme_classic()   # Use a classic theme

# Calculate the linear regression model
lm_model <- lm(MEwhite ~ ind.ed50, data = sig.eigen.diff.ed50.ind)

# Extract the p-value for the slope coefficient
p_value <- summary(lm_model)$coefficients[2, 4]

# Add the p-value to the plot
white_plot <- white_plot + annotate("text", x = max(sig.eigen.diff.ed50.ind$ind.ed50), 
                                                  y = max(sig.eigen.diff.ed50.ind$MEwhite), 
                                                  label = paste("p-value =", round(p_value, 4)), 
                                                  hjust = 1, vjust = 1, size = 4) 

white_plot


yellow_plot <- ggplot(sig.eigen.diff.ed50.ind, aes(x = ind.ed50, y = MEyellow)) +
  geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "ind.ed50", y = "MEyellow", color = "Species") +  # Label axes and legend
  theme_classic()   # Use a classic theme

# Calculate the linear regression model
lm_model <- lm(MEyellow ~ ind.ed50, data = sig.eigen.diff.ed50.ind)

# Extract the p-value for the slope coefficient
p_value <- summary(lm_model)$coefficients[2, 4]

# Add the p-value to the plot
yellow_plot <- yellow_plot + annotate("text", x = max(sig.eigen.diff.ed50.ind$ind.ed50), 
                                                  y = max(sig.eigen.diff.ed50.ind$MEyellow), 
                                                  label = paste("p-value =", round(p_value, 4)), 
                                                  hjust = 1, vjust = 1, size = 4) 

yellow_plot

turquoise_plot <- ggplot(sig.eigen.diff.ed50.ind, aes(x = ind.ed50, y = MEturquoise)) +
  geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "ind.ed50", y = "MEturquoise", color = "Species") +  # Label axes and legend
  theme_classic()   # Use a classic theme

# Calculate the linear regression model
lm_model <- lm(MEturquoise ~ ind.ed50, data = sig.eigen.diff.ed50.ind)

# Extract the p-value for the slope coefficient
p_value <- summary(lm_model)$coefficients[2, 4]

# Add the p-value to the plot
turquoise_plot <- turquoise_plot + annotate("text", x = max(sig.eigen.diff.ed50.ind$ind.ed50), 
                                                  y = max(sig.eigen.diff.ed50.ind$MEturquoise), 
                                                  label = paste("p-value =", round(p_value, 4)), 
                                                  hjust = 1, vjust = 1, size = 4) 

turquoise_plot


######
# Test eignengene diff with red channel intensity
#####
# load in channel intensity data
color.data.ed50 <- read.csv("~/Lab Notebook/Chapter3/Thermal_stress_analysis/color.data.ed50.csv")
species.avg.ev.df <- read.csv("~/Lab Notebook/Chapter3/Thermal_stress_analysis/species.avg.ev.df.csv")

# Calculate the species level slope of Change in R channel intensity
# Greater slope indicates a greater response to heat stress.


# Calculate the species level slopes for the R Channel Intensity
# Acropora abrotanoides
a.abr.model <- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora abrotanoides"))
a.abro.slope <- coef(a.abr.model)[2]


# Acropora hycinthus
a.hya.model <- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora hyacinthus"))
a.hya.slope <- coef(a.hya.model)[2]

# Acropora lutkeni
a.lut.model <- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora lutkeni"))
a.lut.slope <- coef(a.lut.model)[2]

# Acropora pulchra
a.pul.model <- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora pulchra"))
a.pul.slope <- coef(a.pul.model)[2]

# Acropora retusa
a.ret.model <- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora retusa"))
a.ret.slope <- coef(a.ret.model)[2]

# Acropora robusta
a.rob.model<- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora robusta"))
a.rob.slope <- coef(a.rob.model)[2]

# Acropora nasuta
a.nas.model<- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora nasuta"))
a.nas.slope <- coef(a.nas.model)[2]

# Acropora tutuilensis
a.tut.model<- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora tutuilensis"))
a.tut.slope <- coef(a.tut.model)[2]

# Make a data frame for all the slopes
slopes.df <- data.frame(
  species = c("Acropora.abrotanoides", "Acropora.hyacinthus", "Acropora.lutkeni",
              "Acropora.pulchra", "Acropora.retusa", "Acropora.robusta",
              "Acropora.nasuta", "Acropora.tutuilensis"),
  slope = c(a.abro.slope, a.hya.slope, a.lut.slope,
            a.pul.slope, a.ret.slope, a.rob.slope,
            a.nas.slope, a.tut.slope)
)


# Merge slopes with ED.50 data ### changed species spelling and now don't match
# color.data.ed50 <- merge(color.data.ed50, slopes.df , by = "species")
# unique(color.data.ed50$species)

# write.csv(color.data.ed50,
#           file = "~/Lab Notebook/Chapter3/Thermal_stress_analysis/color.data.ed50.rchan.slopes.csv",
#           quote = FALSE,
#           row.names = FALSE)

# Merge the eignevector matrix R channel slope
species.avg.ev.df<- species.avg.ev.df %>%
  mutate(slope = slopes.df$slope[match(species, slopes.df$species)])


## Since large slope indicates heat sensitivity, lets take the inverse of the 
## slopes to test the relationship between eigengene value and R channel
species.avg.ev.df$inv.slope <- 1/species.avg.ev.df$slope

# Test the relationship between species level changes in gene expression and 
# red channel intensity slope. 

# Loop through correlations:
columns.to.correlate <- colnames(MEs)
# Initialize an empty data frame to store results
correlation_results_df <- data.frame(Column_Name = character(),
                                     Correlation_Coefficient = numeric(),
                                     P_Value = numeric(),
                                     stringsAsFactors = FALSE)

# Loop through each column combination
for (col in columns.to.correlate) {
  correlation_result <- cor.test(species.avg.ev.df[[col]], species.avg.ev.df$inv.slope, method = "pearson")
  
  # Extract correlation coefficient and p-value
  corr_coef <- correlation_result$estimate
  p_value <- correlation_result$p.value
  
  # Append the results to the data frame
  correlation_results_df <- rbind(correlation_results_df, 
                                  data.frame(Column_Name = col,
                                             Correlation_Coefficient = corr_coef,
                                             P_Value = p_value,
                                             stringsAsFactors = FALSE))
}

bubble_plot <- ggplot(correlation_results_df, aes(x = Column_Name, y = Correlation_Coefficient, size = P_Value, alpha = ifelse(P_Value < 0.05, 0.7, 0.3))) +
  geom_point(color = "#ff595e") +
  scale_size_continuous(range = c(2,10)) +  # Adjust the range of bubble sizes
  labs(x = "Module", y = "Correlation Coefficient between \n eigenvalue difference and species level R channel slope", size = "P-Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility
bubble_plot


## Looking at the individual modules:
# MEblue
MEblue_plot <- ggplot(species.avg.ev.df, aes(x = inv.slope, y = MEblue)) +
  geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "Inverse of R channel slope", y = "MEblue", color = "Species") +  # Label axes and legend
  theme_classic()   # Use a classic theme

# Calculate the linear regression model
lm_model <- lm(MEblue ~ inv.slope, data = species.avg.ev.df)

# Extract the p-value for the slope coefficient
p_value <- summary(lm_model)$coefficients[2, 4]

# Add the p-value to the plot
MEblue_plot <- MEblue_plot + annotate("text", x = max(species.avg.ev.df$inv.slope), 
                                            y = max(species.avg.ev.df$MEblue), 
                                            label = paste("p-value =", round(p_value, 4)), 
                                            hjust = 1, vjust = 1, size = 4) 

MEblue_plot

# MEwhite
MEwhite_plot <- ggplot(species.avg.ev.df, aes(x = inv.slope, y = MEwhite)) +
  geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "Inverse of R channel slope", y = "MEwhite", color = "Species") +  # Label axes and legend
  theme_classic()   # Use a classic theme

# Calculate the linear regression model
lm_model <- lm(MEwhite ~ slope, data = species.avg.ev.df)

# Extract the p-value for the slope coefficient
p_value <- summary(lm_model)$coefficients[2, 4]

# Add the p-value to the plot
MEwhite_plot <- MEwhite_plot + annotate("text", x = max(species.avg.ev.df$inv.slope), 
                                      y = max(species.avg.ev.df$MEwhite), 
                                      label = paste("p-value =", round(p_value, 4)), 
                                      hjust = 1, vjust = 1, size = 4) 

MEwhite_plot


# MEdarkorange
MEdarkorange_plot <- ggplot(species.avg.ev.df, aes(x = slope, y = MEdarkorange)) +
  geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "R channel slope", y = "MEdarkorange", color = "Species") +  # Label axes and legend
  theme_classic()   # Use a classic theme

# Calculate the linear regression model
lm_model <- lm(MEdarkorange ~ slope, data = species.avg.ev.df)

# Extract the p-value for the slope coefficient
p_value <- summary(lm_model)$coefficients[2, 4]

# Add the p-value to the plot
MEdarkorange_plot <- MEdarkorange_plot + annotate("text", x = max(species.avg.ev.df$slope), 
                                        y = max(species.avg.ev.df$MEdarkorange), 
                                        label = paste("p-value =", round(p_value, 4)), 
                                        hjust = 1, vjust = 1, size = 4) 

MEdarkorange_plot


## Save genes in MEwhite, MEroyalblue, and MEskyblue
skyblue.genes <- colnames(datExpr)[moduleColors=="skyblue"]
skyblue.id <- rep("MEskyblue", length(skyblue.genes))
skyblue.genes.df <- as.data.frame(cbind(skyblue.genes,skyblue.id))
colnames(skyblue.genes.df) <- c("gene", "module")

royalblue.genes <- colnames(datExpr)[moduleColors=="royalblue"]
royalblue.id <- rep("MEroyalblue", length(royalblue.genes))
royalblue.genes.df <- as.data.frame(cbind(royalblue.genes,royalblue.id))
colnames(royalblue.genes.df) <- c("gene", "module")

white.genes <- colnames(datExpr)[moduleColors=="white"]
white.id <- rep("MEwhite", length(white.genes))
white.genes.df <- as.data.frame(cbind(white.genes,white.id))
colnames(white.genes.df) <- c("gene", "module")


yellow.genes <- colnames(datExpr)[moduleColors=="yellow"]
yellow.id <- rep("MEyellow", length(yellow.genes))
yellow.genes.df <- as.data.frame(cbind(yellow.genes,yellow.id))
colnames(yellow.genes.df) <- c("gene", "module")

salmon.genes <- colnames(datExpr)[moduleColors=="salmon"]
salmon.id <- rep("MEsalmon", length(salmon.genes))
salmon.genes.df <- as.data.frame(cbind(salmon.genes,salmon.id))
colnames(salmon.genes.df) <- c("gene", "module")

# merge genes from 3 modules:
genes.of.interest <- bind_rows(skyblue.genes.df,royalblue.genes.df)
genes.of.interest <- bind_rows(genes.of.interest,white.genes.df)

# write.csv(genes.of.interest,
#           file = "~/Lab Notebook/Chapter3/Thermal_stress_analysis/genes.of.interest.eigenvector.csv",
#           quote = FALSE,
#           row.names = FALSE)

# Heat stress eigengenes: What eigengenes have significantly different eigenvalues in heat vs control across all species?
control.eigenvector.df
heat.eigenvector.df

control.eigenvector.values <- control.eigenvector.df[,-c(1,31:34)]
heat.eigenvetor.values <- heat.eigenvector.df[,-c(1,31:34)]

# t-test for difference in each module between control and heat stress
t.test.results <- sapply(colnames(control.eigenvector.values), function(module) {
  t_test <- t.test(control.eigenvector.values[[module]], heat.eigenvetor.values[[module]])
  return(list(module = module, p_value = t_test$p.value))
  })

t_test_results_df <- as.data.frame(do.call(rbind, t.test.results), stringsAsFactors = FALSE)

modules <- t_test_results_df$V1[seq(1, nrow(t_test_results_df), by = 2)]
p_values <- as.numeric(t_test_results_df$V1[seq(2, nrow(t_test_results_df), by = 2)])

t_test_results_df <- data.frame(module = modules, p_value = p_values, stringsAsFactors = FALSE)


p_value_plot <- ggplot(t_test_results_df, aes(x = p_value, y = module)) +
  geom_point(aes(alpha = p_value < 0.05), size = 4, color = "#76c893") +
  scale_alpha_manual(values = c(0.5, 1), guide = FALSE) +  # Adjust alpha
  labs(title = "Module t.test p-values comparing \nheat stress and control samples",
       x = "P-value",
       y = "Module") +
  theme_minimal()

# Display the plot
print(p_value_plot)

# Modules that are significant are yellow, tan, salmon, royal blue, red, 
# midnight blue, lightyellow, grey60, green yellow, darkred, dark orange, black.
# This seems like a lot. Maybe this isn't the best test.


# Using lm and anova to test difference in eigengene value
eigenvector.df

# Fit linear models for each module
lm_results <- lapply(names(eigenvector.df)[2:30], function(module) {
  lm_model <- lm(paste(module, "~ treatment"), data = eigenvector.df)
  lm_summary <- summary(lm_model)
  lm_tidy <- tidy(lm_model)
  lm_tidy$Module <- module
  list(tidy = lm_tidy)
})

# Combine results into a single data frame
lm_df <- map_dfr(lm_results, ~{
  module <- .x$tidy$Module
  slope <- .x$tidy$estimate[2]  # Assuming the slope is the second estimate
  p_value <- .x$tidy$p.value[2]  # Assuming the p-value of the slope is the second p-value
  tibble(Module = module, Slope = slope, P_Value = p_value)
})
lm_df <- lm_df %>% distinct()

anova_results <- lapply(names(eigenvector.df)[2:30], function(module) {
  lm_model <- lm(paste(module, "~ treatment"), data = eigenvector.df)
  anova_result <- anova(lm_model)
  return(anova_result)
})

# Combine reults into a dataframe to plot the slope and p-pvalue of the
# eigengene value difference between control and heat

# Multiple regression adding the species effect 
# NOTE: Comparing all species to a. abrotanoides
lm_results_species_interaction <- lapply(names(eigenvector.df)[2:30], function(module) {
  lm_model_species_interaction <- lm(paste(module, "~ treatment * species"), data = eigenvector.df)
  lm_summary_species_interaction <- summary(lm_model_species_interaction)
  lm_tidy_species_interaction <- tidy(lm_model_species_interaction)
  lm_tidy_species_interaction$Module <- module
  list(summary = lm_summary_species_interaction, tidy = lm_tidy_species_interaction)
})



anova_results_species <- lapply(names(eigenvector.df)[2:30], function(module) {
  lm_model_species <- lm(paste(module, "~ treatment * species"), data = eigenvector.df)
  anova_result_species <- anova(lm_model_species)
  return(anova_result_species)
})


# Combine results into a single data frame
lm_spec_txt_df <- map_dfr(lm_results_species_interaction, ~{
  term <- .x$tidy$term
  module <- .x$tidy$Module
  slope <- .x$tidy$estimate 
  p_value <- .x$tidy$p.value  
  tibble(Term = term, Module = module, Slope = slope, P_Value = p_value)
})
lm_spec_txt_df<-lm_spec_txt_df %>% distinct()
lm_spec_txt_df <- lm_spec_txt_df %>%
  mutate(Color = gsub("^ME", "", Module))
View(lm_spec_txt_df)

# Remove intercept rows:
lm_spec_txt_df <- filter(lm_spec_txt_df, Term != "(Intercept)")

# Plot lm results:

ggplot(data = lm_spec_txt_df, aes(x = Term, y = Slope, color = Color, alpha = ifelse(P_Value < 0.05, 0.7, 0.3))) +
  geom_point(size = 3, color = "black") +  # Adding black outlines
  scale_color_identity() +
  scale_alpha_continuous(guide = "none") +  # Remove the legend for alpha
  theme_minimal() +
  facet_wrap(~ Module, scales = "free_y", ncol = 5) +
  labs(x = "Term", y = "Slope", title = "Slope vs Term") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# filter for treatmentheat p-value < 0.05
modules_to_keep <- lm_spec_txt_df %>%
  filter(Term == 'treatmentheat' & P_Value < 0.05) %>%
  pull(Module)

lm_spec_txt_df_subset <- lm_spec_txt_df %>%
  filter(Module %in% modules_to_keep)

ggplot(data = lm_spec_txt_df_subset, aes(x = Term, y = Slope, color = Color, alpha = ifelse(P_Value < 0.05, 0.7, 0.3))) +
  geom_point(size = 3, color = "black") +  # Adding black outlines
  scale_color_identity() +
  scale_alpha_continuous(guide = "none") +  # Remove the legend for alpha
  theme_minimal() +
  facet_wrap(~ Module, scales = "free_y", ncol = 5) +
  labs(x = "Term", y = "Slope", title = "Slope vs Term") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Anova results:
anova_species_interaction <- map_dfr(anova_results_species, ~{
  term <- .x$tidy$term
  module <- .x$tidy$Module
  slope <- .x$tidy$estimate 
  p_value <- .x$tidy$p.value  
  tibble(Term = term, Module = module, Slope = slope, P_Value = p_value)
})

ggplot(lm_df, aes(x = Slope, y = Module, alpha = ifelse(P_Value < 0.05, 0.7, 0.3), size = P_Value)) +
  geom_point(shape = 21, fill = "blue") +
  scale_size_continuous(range = c(2, 10)) +
  labs(x = "Slope", y = "Module", size = "P_Value", fill = "P_Value", alpha = "Alpha") +
  theme_minimal()


## Looking at the just one significant module: MEroyalblue:
control.royalblue <- data.frame(eigenvalue = control.eigenvector.df$MEroyalblue,
                                treatment = rep("control", nrow(control.eigenvector.df)),
                                species = control.eigenvector.df$species)
heat.royalblue <- data.frame(eigenvalue = heat.eigenvector.df$MEroyalblue,
                             treatment = rep("heat", nrow(heat.eigenvector.df)),
                             species = heat.eigenvector.df$species)
royalblue.control.v.heat <- rbind(control.royalblue,heat.royalblue )

library(ggpubr)
ggplot(royalblue.control.v.heat, aes(x = treatment, y = eigenvalue)) +
  geom_jitter(aes(color = species), width = 0.1, height = 0) +  # Add jittered points
 # geom_boxplot(width = 0.4, position = position_dodge(width = 0.5), outlier.shape = NA) +  # Add half of a boxplot
  labs(x = "Treatment", y = "MEroyalblue Eigenvalue", color = "Species", fill = "Treatment") +
  theme_minimal() +
  stat_compare_means(method = "t.test", comparisons = list(c("control", "heat")))


# Let's see what the gene expression is like for the genes in the 
# RoyalBlue module

# Gene names in the royalblue module:
royalblue.genes 

# Select columns (genes) in the vector royalblue.genes
datExpr.0[1:5,1:5]

royal.blue.counts <-datExpr.0[, royalblue.genes]
royal.blue.counts$species <- eigenvector.df[match(eigenvector.df$seq.sample.id,
                                                  row.names(royal.blue.counts)),
                                            31]
royal.blue.counts$treatment <- eigenvector.df[match(eigenvector.df$seq.sample.id,
                                                    row.names(royal.blue.counts)),
                                              33]

head(royal.blue.counts[,c(1:5, 238:239)])

# Convert row names to a column
royal.blue.counts$Sample <- rownames(royal.blue.counts)

# Reshape the data frame from wide to long format
royal.blue.counts_long <- tidyr::gather(royal.blue.counts, Gene, Count, -Sample, -species, -treatment)

# Filter the rows containing Ahyacinthus genes
royal.blue.counts_ahyacinthus <- royal.blue.counts_long[grep("Ahyacinthus", royal.blue.counts_long$Gene), ]

# Plot
ggplot(royal.blue.counts_ahyacinthus, aes(x = treatment, y = Count)) +
  geom_boxplot() +
  #facet_wrap(~Gene, scales = "free_y") +
  labs(x = "Treatment", y = "Counts", color = "Species") +
  theme_minimal() +
  stat_compare_means(method = "t.test", comparisons = list(c("control", "heat")))

ggplot(royal.blue.counts_ahyacinthus, aes(x = treatment, y = Count, color = species)) +
  geom_violin() +
  geom_boxplot(width = 0.09) +
  stat_compare_means(method = "t.test", comparisons = list(c("control", "heat")), label.y = 10) + # Adjust label.y as needed
  facet_wrap(~species, scales = "free_y") +
  labs(x = "Treatment", y = "Counts", color = "Species") +
  theme_minimal() +
  theme(legend.position = "none")

# ##### There is a species effect with how eigengenes respond to heat stress.
# #### To run a single model on all the eigengenes, I need to reshape the data:
eigenvector_long <- pivot_longer(eigenvector.df,
                                 cols = MEcyan:MEgrey,
                                 names_to = "module",
                                 values_to = "eigengene")
# 
# # Run the model
lm_model_treatment_species <- lm(eigengene ~ treatment * species, data = eigenvector_long)
summary(lm_model_treatment_species)
anova(lm_model_treatment_species)
# 
## Visualize heat vs control eigenvalues for each species
eigenvector_long<- eigenvector_long %>%
  mutate(color = gsub("^ME", "", module))

ggplot(eigenvector_long, aes(x = treatment, y = eigengene, color = species)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +  # Add points for each eigengene
  geom_line(aes(group = species), position = position_dodge(width = 0.5)) +  # Add lines connecting points
  labs(x = "Treatment", y = "Eigengene Value", color = "Species") +
  facet_wrap(~ module, scales = "free_y", ncol = 5) +  # Facet by module
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Species level average:
species_txt_eigenvector_mean <- eigenvector_long %>%
  group_by(species, treatment, module) %>%
  summarise_at(vars(eigengene),
               list(mean_eigengene = mean))

ggplot(species_txt_eigenvector_mean, aes(x = treatment, y = mean_eigengene, color = species)) +
  geom_point(size = 2) +  # Add points for each eigengene
  geom_line(aes(group = species)) +  # Add lines connecting points
  labs(x = "Treatment", y = "Eigengene Value", color = "Species") +
  facet_wrap(~ module, scales = "free_y", ncol = 5) +  # Facet by module
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

## To Find heat response eigengenes:
# 1) treatment:species significant and 2) ED50 correlated with response.
# Sig treatment:species: MESalmon, MElightgreen, MEyellow, MEdarkred,
# MEgreen

#2) ED50 correlated with response
## Plot the species averate ED50 and the difference of the eigengene value:
sig.eigen.diff.ed50.ind <- difference.av.ev.df[,c("seq.sample.id","MEsalmon", "MElightgreen", "MEyellow","MEdarkred", "MEgreen", "species", "mean.ed50","ind.ed50")]

## Salmon plot
official.species.colors <- c("Acropora.tutuilensis" = "#F7C11E",
                             "Acropora.lutkeni" = "#FF9A17",
                             "Acropora.pulchra" = "#fb5607",
                             "Acropora.retusa" = "#ff006e",
                             "Acropora.nasuta" = "#c11cad",
                             "Acropora.robusta" = "#8338ec",
                             "Acropora.hyacinthus" = "#3E6FCB",
                             "Acropora.abrotanoides" = "#3a86ff")

salmon_plot <- ggplot(sig.eigen.diff.ed50.ind, aes(x = ind.ed50, y = MEsalmon)) +
  geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
  scale_color_manual(values = official.species.colors)+
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "Colony thermal tolerance", y = "MEsalmon", color = "Species") +  # Label axes and legend
  theme_classic() +  # Use a classic theme
  theme(legend.position = "none")
# Calculate the linear regression model
lm_model <- lm(MEsalmon ~ ind.ed50, data = sig.eigen.diff.ed50.ind)

# Extract the p-value for the slope coefficient
p_value <- summary(lm_model)$coefficients[2, 4]

# Add the p-value to the plot

salmon_plot <- salmon_plot + annotate("text", x = max(sig.eigen.diff.ed50.ind$ind.ed50), 
                                    y = max(sig.eigen.diff.ed50.ind$MEsalmon), 
                                    label = paste("p-value =", round(p_value, 4)), 
                                    hjust = 1, vjust = 1, size = 4) 

salmon_plot

## lightgreen plot
lightgreen_plot <- ggplot(sig.eigen.diff.ed50.ind, aes(x = ind.ed50, y = MElightgreen)) +
  geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "Colony thermal tolerance", y = "MElightgreen", color = "Species") +  # Label axes and legend
  theme_classic()   # Use a classic theme

# Calculate the linear regression model
lm_model <- lm(MElightgreen ~ ind.ed50, data = sig.eigen.diff.ed50.ind)

# Extract the p-value for the slope coefficient
p_value <- summary(lm_model)$coefficients[2, 4]

# Add the p-value to the plot
lightgreen_plot <- lightgreen_plot + annotate("text", x = max(sig.eigen.diff.ed50.ind$ind.ed50), 
                                      y = max(sig.eigen.diff.ed50.ind$MElightgreen), 
                                      label = paste("p-value =", round(p_value, 4)), 
                                      hjust = 1, vjust = 1, size = 4) 

lightgreen_plot

## yellow plot
yellow_plot <- ggplot(sig.eigen.diff.ed50.ind, aes(x = ind.ed50, y = MEyellow)) +
  geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
  scale_color_manual(values = official.species.colors)+
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "Colony thermal tolerance", y = "MEyellow", color = "Species") +  # Label axes and legend
  theme_classic()   # Use a classic theme

# Calculate the linear regression model
lm_model <- lm(MEyellow ~ ind.ed50, data = sig.eigen.diff.ed50.ind)

# Extract the p-value for the slope coefficient
p_value <- summary(lm_model)$coefficients[2, 4]

# Add the p-value to the plot
yellow_plot <- yellow_plot + annotate("text", x = max(sig.eigen.diff.ed50.ind$ind.ed50), 
                                              y = max(sig.eigen.diff.ed50.ind$MEyellow), 
                                              label = paste("p-value =", round(p_value, 4)), 
                                              hjust = 1, vjust = 1, size = 4) 

yellow_plot

library(patchwork)
salmon_plot + yellow_plot
## darkred plot
darkred_plot <- ggplot(sig.eigen.diff.ed50.ind, aes(x = ind.ed50, y = MEdarkred)) +
  geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "ind.ed50", y = "MEdarkred", color = "Species") +  # Label axes and legend
  theme_classic()   # Use a classic theme

# Calculate the linear regression model
lm_model <- lm(MEdarkred ~ ind.ed50, data = sig.eigen.diff.ed50.ind)

# Extract the p-value for the slope coefficient
p_value <- summary(lm_model)$coefficients[2, 4]

# Add the p-value to the plot
darkred_plot <- darkred_plot + annotate("text", x = max(sig.eigen.diff.ed50.ind$ind.ed50), 
                                      y = max(sig.eigen.diff.ed50.ind$MEdarkred), 
                                      label = paste("p-value =", round(p_value, 4)), 
                                      hjust = 1, vjust = 1, size = 4) 

darkred_plot

## green plot
green_plot <- ggplot(sig.eigen.diff.ed50.ind, aes(x = ind.ed50, y = MEgreen)) +
  geom_point(size = 3, aes(color = species)) +  # Adjust the size of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "ind.ed50", y = "MEgreen", color = "Species") +  # Label axes and legend
  theme_classic()   # Use a classic theme

# Calculate the linear regression model
lm_model <- lm(MEgreen ~ ind.ed50, data = sig.eigen.diff.ed50.ind)

# Extract the p-value for the slope coefficient
p_value <- summary(lm_model)$coefficients[2, 4]

# Add the p-value to the plot
green_plot <- green_plot + annotate("text", x = max(sig.eigen.diff.ed50.ind$ind.ed50), 
                                        y = max(sig.eigen.diff.ed50.ind$MEgreen), 
                                        label = paste("p-value =", round(p_value, 4)), 
                                        hjust = 1, vjust = 1, size = 4) 

green_plot


### YELLOW AND SALMON ARE OUR HEAT RESPONSE EIGENGENES ###

#### Signed module analysis
# Originall version produces no gray module meaning no unassigned genes. 
# ge.network= blockwiseModules(datExpr, power = 5,
#                              TOMType = "signed", minModuleSize = 30, # check signed vs unsigned and minModuleSize
#                              reassignThreshold = 0, mergeCutHeight = 0.25,
#                              numericLabels = TRUE, pamRespectsDendro = FALSE,
#                              saveTOMs = FALSE,
#                              verbose = 3)

ge.network = blockwiseModules(datExpr, power = 5,
                              TOMType = "signed", minModuleSize = 30, 
                              reassignThreshold = 0.05, mergeCutHeight = 0.25,
                              numericLabels = TRUE, pamRespectsDendro = FALSE,
                              saveTOMs = FALSE,
                              verbose = 3)

table(ge.network$colors)                         


sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(ge.network$colors) 
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(ge.network$dendrograms[[1]], mergedColors[ge.network$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)                        

moduleLabels.e = ge.network$colors
moduleColors.e = labels2colors(ge.network$colors)
ME.es = ge.network$MEs;
geneTree.e = ge.network$dendrograms[[1]]                        

softPower = 5 #Chosen in the graphs before
adjacency = adjacency(datExpr, power = softPower, type = "signed") # signed matrix
#Transforming the adjacency matrix in a topological overlap
TOM = TOMsimilarity(adjacency) #Calculating the topological overlap matrix
dissTOM = 1-TOM 

# Clean up vector space
rm(acro.counts.mat, acro.counts.matrix, counts.df, ddsHTSeq,adjacency)

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
## 2.b.5 Merging of modules whose expression profiles are very similar

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
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
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

moduleColors = mergedColors

#Building numeric labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs



# Associating Modules and Phenotypes
#Relating modules to characteristics and identifying important genes
#Defining the number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculating MEs with label colors
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs,datTraits[,-c(1,2)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


#sizeGrWindow(8,4)

#Displaying correlations and its p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 8, 6, 8))
par(oma = c(8, 0, 4, 0))  # Adjust the outer margin

#Displaying the correlation values in a heatmap plot
pdf(file = "~/Lab Notebook/Chapter3/Tag_seq_analysis/Figures/module_sample_relationship.pdf", wi = 15, he = 8.5)
par(mar=c(9,8,4,1)+.1)
my_palette <- colorRampPalette(c('#FF0000',"white", '#0000FF'))(n = 50)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits[,-c(1,2)]),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = my_palette,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-sample relationships"))
dev.off()

# Isolate matrix of average eigenvector values:
eigenvector.df <- as.data.frame(MEs)
eigenvector.df$seq.sample.id <- rownames(eigenvector.df)
eigenvector.df<- merge(eigenvector.df, my.metadata[, c("seq.sample.id", "species","colony.id", "treatment")], all.x = TRUE)

# Add species level ED50 average:
# Read in ED50 data
ED50.data <- read.csv("~/Lab Notebook/Chapter3/Thermal_stress_analysis/Phylogeny resolved ED50 curves/individual.specied.ED50.csv", 
                      row.names = NULL)
eigenvector.df$mean.ed50 <- ED50.data[match(eigenvector.df$colony.id, ED50.data$colony.id),2]
# Need to create a difference matrix: control eigenvector value - heat eigenvector value
## Control and heat eigenvector DF
control.eigenvector.df <-subset(eigenvector.df, treatment == "control")
heat.eigenvector.df <- subset(eigenvector.df, treatment == "heat")

# Difference in the eigenvector values
cols.to.subtract <- colnames(MEs)

# Subtract corresponding elements
difference.ev.df <- control.eigenvector.df[,-33]
difference.ev.df[, cols.to.subtract] <- control.eigenvector.df[, cols.to.subtract] - heat.eigenvector.df[, cols.to.subtract]


# Subtract the absolute values
abs.control.eigenvector.df<-abs(control.eigenvector.df[,cols.to.subtract])
abs.heat.eigenvector.df<-abs(heat.eigenvector.df[,cols.to.subtract])
difference.av.ev.df <- control.eigenvector.df[,-33]
difference.av.ev.df[, cols.to.subtract] <- abs.control.eigenvector.df[, cols.to.subtract] - abs.heat.eigenvector.df[, cols.to.subtract]

eigengene.exp = as.data.frame(t(MEs[, -c(1)]))


