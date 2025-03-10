# DESeq2 of genes in eigen genes with significant relationship with ED 50
# Load Libraries
library(DESeq2)
library(tidyr)
library(tidyverse)

# Load in data
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
  design = ~ condition.treatment)
# Pre-filter
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]

# Estimate normalization factors
ddsHTSeq <- estimateSizeFactors(ddsHTSeq)
ddsHTSeq <- estimateDispersions(ddsHTSeq)

# Differential expression
dds <- DESeq(ddsHTSeq) 
resultsNames(dds)
num_transcripts <- nrow(dds)
num_transcripts # 22090

# Heat stress response
names(dds@colData) # condition.treatment
resultsNames(dds)
hs.results <- lfcShrink(dds,
                        contrast = c("condition.treatment", "heat", "control"),
                        type = "ashr")

heat.response.results <- hs.results[complete.cases(hs.results),]
heat.response.results.df <- as.data.frame(heat.response.results)
# Isolate genes in the yellow and salmon modules
# MEsalmon.genes$gene <- colnames(datExpr)[moduleColors=="salmon"]
# MEyellow.genes$gene <-colnames(datExpr)[moduleColors=="yellow"]

# save the salmon and yellow genes 
# salmon.genes <- write.csv(MEsalmon.genes$gene,
#                           file = "~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/MEsalmon.genes$gene.csv",
#                           row.names = FALSE)
# yellow.genes <- write.csv(MEyellow.genes$gene,
#                           file = "~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/MEyellow.genes$gene.csv",
#                           row.names = FALSE)



MEsalmon.genes <- read.csv("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/MEsalmon.genes.csv")
MEyellow.genes <- read.csv("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/MEyellow.genes.csv")

me.salmon <- rownames(heat.response.results.df) %in% MEsalmon.genes$gene
me.yellow <- rownames(heat.response.results.df) %in% MEyellow.genes$gene

sum(me.salmon)
sum(me.yellow)

# Add the me.salmon column to the dataframe
heat.response.results.df$me.salmon <- me.salmon
heat.response.results.df$me.yellow <- me.yellow


## Create dds object for each species where design is treatment only
### A. abrotanoides
abro.files.list <- paste(filter(my.metadata, species == "Acropora.abrotanoides")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
abro.metadata <- filter(my.metadata, species == "Acropora.abrotanoides")
abro.sampleNames <- abro.metadata$seq.sample.id

# Create sample table
abro.sampleTable <- data.frame(sampleName = abro.sampleNames, 
                             fileName = abro.files.list, 
                             condition = abro.metadata)
# Convert variables to factors
factorVars <- c("condition.species", "condition.colony.id", "condition.treatment")
abro.sampleTable[, factorVars] <- lapply(abro.sampleTable[, factorVars], factor)
colnames(abro.sampleTable)

# Create DESeqDataSet Object
abro.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = abro.sampleTable,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep <- rowSums(counts(abro.ddsHTSeq)) >= 10
abro.ddsHTSeq <- abro.ddsHTSeq[keep,]

# Estimate normalization factors
abro.ddsHTSeq <- estimateSizeFactors(abro.ddsHTSeq)
abro.ddsHTSeq <- estimateDispersions(abro.ddsHTSeq)

# Differential expression
abro.dds <- DESeq(abro.ddsHTSeq) 
resultsNames(abro.dds)
nrow(abro.dds)
# 15182

# Heat stress response
names(abro.dds@colData) # condition.treatment
resultsNames(abro.dds)
abro.hs.results <- lfcShrink(abro.dds,
                        contrast = c("condition.treatment", "heat", "control"),
                        type = "ashr")

abro.heat.response.results <- abro.hs.results[complete.cases(abro.hs.results),]
abro.heat.response.results.df <- as.data.frame(abro.heat.response.results)
# Isolate genes in the yellow and salmon modules
# MEsalmon.genes$gene <- colnames(datExpr)[moduleColors=="salmon"]
# MEyellow.genes$gene <-colnames(datExpr)[moduleColors=="yellow"]

abro.me.salmon <- rownames(abro.heat.response.results.df) %in% MEsalmon.genes$gene
abro.me.yellow <- rownames(abro.heat.response.results.df) %in% MEyellow.genes$gene

sum(abro.me.salmon) #157
sum(abro.me.yellow) #580
# Add the me.salmon column to the dataframe
abro.heat.response.results.df$me.salmon <- abro.me.salmon
abro.heat.response.results.df$me.yellow <- abro.me.yellow
abro.heat.response.results.df$species <- "Acropora abrotanoides"

### A. hyacinthus
hya.files.list <- paste(filter(my.metadata, species == "Acropora.hyacinthus")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
hya.metadata <- filter(my.metadata, species == "Acropora.hyacinthus")
hya.sampleNames <- hya.metadata$seq.sample.id

# Create sample table
hya.sampleTable <- data.frame(sampleName = hya.sampleNames, 
                               fileName = hya.files.list, 
                               condition = hya.metadata)
# Convert variables to factors
factorVars <- c("condition.species", "condition.colony.id", "condition.treatment")
hya.sampleTable[, factorVars] <- lapply(hya.sampleTable[, factorVars], factor)
colnames(hya.sampleTable)

# Create DESeqDataSet Object
hya.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = hya.sampleTable,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep <- rowSums(counts(hya.ddsHTSeq)) >= 10
hya.ddsHTSeq <- hya.ddsHTSeq[keep,]

# Estimate normalization factors
hya.ddsHTSeq <- estimateSizeFactors(hya.ddsHTSeq)
hya.ddsHTSeq <- estimateDispersions(hya.ddsHTSeq)

# Differential expression
hya.dds <- DESeq(hya.ddsHTSeq) 
resultsNames(hya.dds)
nrow(hya.dds)
# 17248

# Heat stress response
names(hya.dds@colData) # condition.treatment
resultsNames(hya.dds)
hya.hs.results <- lfcShrink(hya.dds,
                             contrast = c("condition.treatment", "heat", "control"),
                             type = "ashr")

hya.heat.response.results <- hya.hs.results[complete.cases(hya.hs.results),]
hya.heat.response.results.df <- as.data.frame(hya.heat.response.results)
# Isolate genes in the yellow and salmon modules
# MEsalmon.genes$gene <- colnames(datExpr)[moduleColors=="salmon"]
# MEyellow.genes$gene <-colnames(datExpr)[moduleColors=="yellow"]

hya.me.salmon <- rownames(hya.heat.response.results.df) %in% MEsalmon.genes$gene
hya.me.yellow <- rownames(hya.heat.response.results.df) %in% MEyellow.genes$gene

sum(hya.me.salmon) #203
sum(hya.me.yellow) #931
# Add the me.salmon column to the dataframe
hya.heat.response.results.df$me.salmon <- hya.me.salmon
hya.heat.response.results.df$me.yellow <- hya.me.yellow
hya.heat.response.results.df$species <- "Acropora hyacinthus"


### A. pulchra
pul.files.list <- paste(filter(my.metadata, species == "Acropora.pulchra")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
pul.metadata <- filter(my.metadata, species == "Acropora.pulchra")
pul.sampleNames <- pul.metadata$seq.sample.id

# Create sample table
pul.sampleTable <- data.frame(sampleName = pul.sampleNames, 
                              fileName = pul.files.list, 
                              condition = pul.metadata)
# Convert variables to factors
factorVars <- c("condition.species", "condition.colony.id", "condition.treatment")
pul.sampleTable[, factorVars] <- lapply(pul.sampleTable[, factorVars], factor)
colnames(pul.sampleTable)

# Create DESeqDataSet Object
pul.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = pul.sampleTable,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep <- rowSums(counts(pul.ddsHTSeq)) >= 10
pul.ddsHTSeq <- pul.ddsHTSeq[keep,]

# Estimate normalization factors
pul.ddsHTSeq <- estimateSizeFactors(pul.ddsHTSeq)
pul.ddsHTSeq <- estimateDispersions(pul.ddsHTSeq)

# Differential expression
pul.dds <- DESeq(pul.ddsHTSeq) 
resultsNames(pul.dds)
nrow(pul.dds)
# 16967

# Heat stress response
names(pul.dds@colData) # condition.treatment
resultsNames(pul.dds)
pul.hs.results <- lfcShrink(pul.dds,
                            contrast = c("condition.treatment", "heat", "control"),
                            type = "ashr")

pul.heat.response.results <- pul.hs.results[complete.cases(pul.hs.results),]
pul.heat.response.results.df <- as.data.frame(pul.heat.response.results)
# Isolate genes in the yellow and salmon modules
# MEsalmon.genes$gene <- colnames(datExpr)[moduleColors=="salmon"]
# MEyellow.genes$gene <-colnames(datExpr)[moduleColors=="yellow"]

pul.me.salmon <- rownames(pul.heat.response.results.df) %in% MEsalmon.genes$gene
pul.me.yellow <- rownames(pul.heat.response.results.df) %in% MEyellow.genes$gene

sum(pul.me.salmon) #81
sum(pul.me.yellow) #269
# Add the me.salmon column to the dataframe
pul.heat.response.results.df$me.salmon <- pul.me.salmon
pul.heat.response.results.df$me.yellow <- pul.me.yellow
pul.heat.response.results.df$species <- "Acropora pulchra"

### A. lutkeni
lut.files.list <- paste(filter(my.metadata, species == "Acropora.lutkeni")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
lut.metadata <- filter(my.metadata, species == "Acropora.lutkeni")
lut.sampleNames <- lut.metadata$seq.sample.id

# Create sample table
lut.sampleTable <- data.frame(sampleName = lut.sampleNames, 
                              fileName = lut.files.list, 
                              condition = lut.metadata)
# Convert variables to factors
factorVars <- c("condition.species", "condition.colony.id", "condition.treatment")
lut.sampleTable[, factorVars] <- lapply(lut.sampleTable[, factorVars], factor)
colnames(lut.sampleTable)

# Create DESeqDataSet Object
lut.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = lut.sampleTable,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep <- rowSums(counts(lut.ddsHTSeq)) >= 10
lut.ddsHTSeq <- lut.ddsHTSeq[keep,]

# Estimate normalization factors
lut.ddsHTSeq <- estimateSizeFactors(lut.ddsHTSeq)
lut.ddsHTSeq <- estimateDispersions(lut.ddsHTSeq)

# Differential expression
lut.dds <- DESeq(lut.ddsHTSeq) 
resultsNames(lut.dds)
nrow(lut.dds)
# 15251

# Heat stress response
names(lut.dds@colData) # condition.treatment
resultsNames(lut.dds)
lut.hs.results <- lfcShrink(lut.dds,
                            contrast = c("condition.treatment", "heat", "control"),
                            type = "ashr")

lut.heat.response.results <- lut.hs.results[complete.cases(lut.hs.results),]
lut.heat.response.results.df <- as.data.frame(lut.heat.response.results)


lut.me.salmon <- rownames(lut.heat.response.results.df) %in% MEsalmon.genes$gene
lut.me.yellow <- rownames(lut.heat.response.results.df) %in% MEyellow.genes$gene

sum(lut.me.salmon) #230
sum(lut.me.yellow) #1175
# Add the me.salmon column to the dataframe
lut.heat.response.results.df$me.salmon <- lut.me.salmon
lut.heat.response.results.df$me.yellow <- lut.me.yellow
lut.heat.response.results.df$species <- "Acropora lutkeni"

### Acropora retusa
ret.files.list <- paste(filter(my.metadata, species == "Acropora.retusa")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
ret.metadata <- filter(my.metadata, species == "Acropora.retusa")
ret.sampleNames <- ret.metadata$seq.sample.id

# Create sample table
ret.sampleTable <- data.frame(sampleName = ret.sampleNames, 
                              fileName = ret.files.list, 
                              condition = ret.metadata)
# Convert variables to factors
factorVars <- c("condition.species", "condition.colony.id", "condition.treatment")
ret.sampleTable[, factorVars] <- lapply(ret.sampleTable[, factorVars], factor)
colnames(ret.sampleTable)

# Create DESeqDataSet Object
ret.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = ret.sampleTable,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep <- rowSums(counts(ret.ddsHTSeq)) >= 10
ret.ddsHTSeq <- ret.ddsHTSeq[keep,]

# Estimate normalization factors
ret.ddsHTSeq <- estimateSizeFactors(ret.ddsHTSeq)
ret.ddsHTSeq <- estimateDispersions(ret.ddsHTSeq)

# Differential expression
ret.dds <- DESeq(ret.ddsHTSeq) 
resultsNames(ret.dds)
nrow(ret.dds)
# 16212

# Heat stress response
names(ret.dds@colData) # condition.treatment
resultsNames(ret.dds)
ret.hs.results <- lfcShrink(ret.dds,
                            contrast = c("condition.treatment", "heat", "control"),
                            type = "ashr")

ret.heat.response.results <- ret.hs.results[complete.cases(ret.hs.results),]
ret.heat.response.results.df <- as.data.frame(ret.heat.response.results)


ret.me.salmon <- rownames(ret.heat.response.results.df) %in% MEsalmon.genes$gene
ret.me.yellow <- rownames(ret.heat.response.results.df) %in% MEyellow.genes$gene

sum(ret.me.salmon) #238
sum(ret.me.yellow) #1181
# Add the me.salmon column to the dataframe
ret.heat.response.results.df$me.salmon <- ret.me.salmon
ret.heat.response.results.df$me.yellow <- ret.me.yellow
ret.heat.response.results.df$species <- "Acropora retusa"

### Acropora robusta
rob.files.list <- paste(filter(my.metadata, species == "Acropora.robusta")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
rob.metadata <- filter(my.metadata, species == "Acropora.robusta")
rob.sampleNames <- rob.metadata$seq.sample.id

# Create sample table
rob.sampleTable <- data.frame(sampleName = rob.sampleNames, 
                              fileName = rob.files.list, 
                              condition = rob.metadata)
# Convert variables to factors
factorVars <- c("condition.species", "condition.colony.id", "condition.treatment")
rob.sampleTable[, factorVars] <- lapply(rob.sampleTable[, factorVars], factor)
colnames(rob.sampleTable)

# Create DESeqDataSet Object
rob.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = rob.sampleTable,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep <- rowSums(counts(rob.ddsHTSeq)) >= 10
rob.ddsHTSeq <- rob.ddsHTSeq[keep,]

# Estimate normalization factors
rob.ddsHTSeq <- estimateSizeFactors(rob.ddsHTSeq)
rob.ddsHTSeq <- estimateDispersions(rob.ddsHTSeq)

# Differential expression
rob.dds <- DESeq(rob.ddsHTSeq) 
resultsNames(rob.dds)
nrow(rob.dds)
# 16308

# Heat stress response
names(rob.dds@colData) # condition.treatment
resultsNames(rob.dds)
rob.hs.results <- lfcShrink(rob.dds,
                            contrast = c("condition.treatment", "heat", "control"),
                            type = "ashr")

rob.heat.response.results <- rob.hs.results[complete.cases(rob.hs.results),]
rob.heat.response.results.df <- as.data.frame(rob.heat.response.results)


rob.me.salmon <- rownames(rob.heat.response.results.df) %in% MEsalmon.genes$gene
rob.me.yellow <- rownames(rob.heat.response.results.df) %in% MEyellow.genes$gene

sum(rob.me.salmon) #224
sum(rob.me.yellow) #1069
# Add the me.salmon column to the dataframe
rob.heat.response.results.df$me.salmon <- rob.me.salmon
rob.heat.response.results.df$me.yellow <- rob.me.yellow
rob.heat.response.results.df$species <- "Acropora robusta"

### Acropora tutuilensis
tutu.files.list <- paste(filter(my.metadata, species == "Acropora.tutuilensis")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
tutu.metadata <- filter(my.metadata, species == "Acropora.tutuilensis")
tutu.sampleNames <- tutu.metadata$seq.sample.id

# Create sample table
tutu.sampleTable <- data.frame(sampleName = tutu.sampleNames, 
                              fileName = tutu.files.list, 
                              condition = tutu.metadata)
# Convert variables to factors
factorVars <- c("condition.species", "condition.colony.id", "condition.treatment")
tutu.sampleTable[, factorVars] <- lapply(tutu.sampleTable[, factorVars], factor)
colnames(tutu.sampleTable)

# Create DESeqDataSet Object
tutu.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = tutu.sampleTable,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep <- rowSums(counts(tutu.ddsHTSeq)) >= 10
tutu.ddsHTSeq <- tutu.ddsHTSeq[keep,]

# Estimate normalization factors
tutu.ddsHTSeq <- estimateSizeFactors(tutu.ddsHTSeq)
tutu.ddsHTSeq <- estimateDispersions(tutu.ddsHTSeq)

# Differential expression
tutu.dds <- DESeq(tutu.ddsHTSeq) 
resultsNames(tutu.dds)
nrow(tutu.dds)
# 18121

# Heat stress response
names(tutu.dds@colData) # condition.treatment
resultsNames(tutu.dds)
tutu.hs.results <- lfcShrink(tutu.dds,
                            contrast = c("condition.treatment", "heat", "control"),
                            type = "ashr")

tutu.heat.response.results <- tutu.hs.results[complete.cases(tutu.hs.results),]
tutu.heat.response.results.df <- as.data.frame(tutu.heat.response.results)


tutu.me.salmon <- rownames(tutu.heat.response.results.df) %in% MEsalmon.genes$gene
tutu.me.yellow <- rownames(tutu.heat.response.results.df) %in% MEyellow.genes$gene

sum(tutu.me.salmon) #248
sum(tutu.me.yellow) #1327
# Add the me.salmon column to the dataframe
tutu.heat.response.results.df$me.salmon <- tutu.me.salmon
tutu.heat.response.results.df$me.yellow <- tutu.me.yellow
tutu.heat.response.results.df$species <- "Acropora tutuilensis"

### Acrpora nasuta
nas.files.list <- paste(filter(my.metadata, species == "Acropora.nasuta")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
nas.metadata <- filter(my.metadata, species == "Acropora.nasuta")
nas.sampleNames <- nas.metadata$seq.sample.id

# Create sample table
nas.sampleTable <- data.frame(sampleName = nas.sampleNames, 
                               fileName = nas.files.list, 
                               condition = nas.metadata)
# Convert variables to factors
factorVars <- c("condition.species", "condition.colony.id", "condition.treatment")
nas.sampleTable[, factorVars] <- lapply(nas.sampleTable[, factorVars], factor)
colnames(nas.sampleTable)

# Create DESeqDataSet Object
nas.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = nas.sampleTable,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep <- rowSums(counts(nas.ddsHTSeq)) >= 10
nas.ddsHTSeq <- nas.ddsHTSeq[keep,]

# Estimate normalization factors
nas.ddsHTSeq <- estimateSizeFactors(nas.ddsHTSeq)
nas.ddsHTSeq <- estimateDispersions(nas.ddsHTSeq)

# Differential expression
nas.dds <- DESeq(nas.ddsHTSeq) 
resultsNames(nas.dds)
nrow(nas.dds)
# 14427

# Heat stress response
names(nas.dds@colData) # condition.treatment
resultsNames(nas.dds)
nas.hs.results <- lfcShrink(nas.dds,
                             contrast = c("condition.treatment", "heat", "control"),
                             type = "ashr")

nas.heat.response.results <- nas.hs.results[complete.cases(nas.hs.results),]
nas.heat.response.results.df <- as.data.frame(nas.heat.response.results)


nas.me.salmon <- rownames(nas.heat.response.results.df) %in% MEsalmon.genes$gene
nas.me.yellow <- rownames(nas.heat.response.results.df) %in% MEyellow.genes$gene

sum(nas.me.salmon) #152
sum(nas.me.yellow) #614
# Add the me.salmon column to the dataframe
nas.heat.response.results.df$me.salmon <- nas.me.salmon
nas.heat.response.results.df$me.yellow <- nas.me.yellow
nas.heat.response.results.df$species <- "Acropora nasuta"


# Combine data frames for plotting salmon gene expression
salmon.abro <- abro.heat.response.results.df[abro.heat.response.results.df$me.salmon,]
salmon.hya <- hya.heat.response.results.df[hya.heat.response.results.df$me.salmon,]
salmon.pul <- pul.heat.response.results.df[pul.heat.response.results.df$me.salmon,]
salmon.lut <- lut.heat.response.results.df[lut.heat.response.results.df$me.salmon,]
salmon.ret <- ret.heat.response.results.df[ret.heat.response.results.df$me.salmon,]
salmon.rob <- rob.heat.response.results.df[rob.heat.response.results.df$me.salmon,]
salmon.tut <- tutu.heat.response.results.df[tutu.heat.response.results.df$me.salmon,]
salmon.nas <- nas.heat.response.results.df[nas.heat.response.results.df$me.salmon, ]

fil.salmon.abro <- salmon.abro[salmon.abro$padj < 0.05, ]
fil.salmon.hya <- salmon.hya[salmon.hya$padj < 0.05, ]
fil.salmon.pul <- salmon.pul[salmon.pul$padj < 0.05, ]
fil.salmon.lut <- salmon.lut[salmon.lut$padj < 0.05, ]
fil.salmon.ret <- salmon.ret[salmon.ret$padj < 0.05, ]
fil.salmon.rob <- salmon.rob[salmon.rob$padj < 0.05, ]
fil.salmon.tut <- salmon.tut[salmon.tut$padj < 0.05, ]
fil.salmon.nas <- salmon.nas[salmon.nas$padj < 0.05, ]


# Combine data frames for plotting yellow gene expression
yellow.abro <- abro.heat.response.results.df[abro.heat.response.results.df$me.yellow,]
yellow.hya <- hya.heat.response.results.df[hya.heat.response.results.df$me.yellow,]
yellow.pul <- pul.heat.response.results.df[pul.heat.response.results.df$me.yellow,]
yellow.lut <- lut.heat.response.results.df[lut.heat.response.results.df$me.yellow,]
yellow.ret <- ret.heat.response.results.df[ret.heat.response.results.df$me.yellow,]
yellow.rob <- rob.heat.response.results.df[rob.heat.response.results.df$me.yellow,]
yellow.tut <- tutu.heat.response.results.df[tutu.heat.response.results.df$me.yellow,]
yellow.nas <- nas.heat.response.results.df[nas.heat.response.results.df$me.yellow, ]

fil.yellow.abro <- yellow.abro[yellow.abro$padj < 0.05, ]
fil.yellow.hya <- yellow.hya[yellow.hya$padj < 0.05, ]
fil.yellow.pul <- yellow.pul[yellow.pul$padj < 0.05, ]
fil.yellow.lut <- yellow.lut[yellow.lut$padj < 0.05, ]
fil.yellow.ret <- yellow.ret[yellow.ret$padj < 0.05, ]
fil.yellow.rob <- yellow.rob[yellow.rob$padj < 0.05, ]
fil.yellow.tut <- yellow.tut[yellow.tut$padj < 0.05, ]
fil.yellow.nas <- yellow.nas[yellow.nas$padj < 0.05, ]

yellow.df.list <- list(yellow.abro, yellow.hya, yellow.pul, yellow.lut,yellow.ret, yellow.rob, yellow.tut, yellow.nas)
all.species.yellow.df <- do.call(rbind, yellow.df.list)
official.species.colors <- c("Acropora tutuilensis" = "#F7C11E",
                             "Acropora lutkeni" = "#FF9A17",
                             "Acropora pulchra" = "#fb5607",
                             "Acropora retusa" = "#ff006e",
                             "Acropora nasuta" = "#c11cad",
                             "Acropora robusta" = "#8338ec",
                             "Acropora hyacinthus" = "#3E6FCB",
                             "Acropora abrotanoides" = "#3a86ff")

species.order <- c("Acropora tutuilensis", "Acropora lutkeni",
                   "Acropora pulchra", "Acropora retusa","Acropora nasuta",
                  "Acropora robusta","Acropora hyacinthus","Acropora abrotanoides")
all.species.yellow.df$species <- factor(all.species.yellow.df$species, levels = species.order)


# Plot
rain_height <- .1
l2fc.plot.all.samples.yel <- ggplot(all.species.yellow.df, aes(x = species, y = log2FoldChange, fill = species )) +
  # clouds
  # introdataviz::geom_flat_violin(trim=FALSE,
  #                                position = position_nudge(x = rain_height+.05),
  #                                width = 5) +
  # rain
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
  # theme(legend.position = "none") +
   labs(x = NULL) +
   ggtitle("Heat-stress Log2FC\nof genes in Yellow Module")
  l2fc.plot.all.samples.yel

## salmon module
  salmon.df.list <- list(salmon.abro, salmon.hya, salmon.pul, salmon.lut,salmon.ret, salmon.rob, salmon.tut, salmon.nas)
  all.species.salmon.df <- do.call(rbind, salmon.df.list)
  official.species.colors <- c("Acropora tutuilensis" = "#F7C11E",
                               "Acropora lutkeni" = "#FF9A17",
                               "Acropora pulchra" = "#fb5607",
                               "Acropora retusa" = "#ff006e",
                               "Acropora nasuta" = "#c11cad",
                               "Acropora robusta" = "#8338ec",
                               "Acropora hyacinthus" = "#3E6FCB",
                               "Acropora abrotanoides" = "#3a86ff")
  
  species.order <- c("Acropora tutuilensis", "Acropora lutkeni",
                     "Acropora pulchra", "Acropora retusa","Acropora nasuta",
                     "Acropora robusta","Acropora hyacinthus","Acropora abrotanoides")
  all.species.salmon.df$species <- factor(all.species.salmon.df$species, levels = species.order)
  
  
  # Plot
  rain_height <- .1
  
  l2fc.plot.all.samples.salmon <- ggplot(all.species.salmon.df, aes(x = species, y = log2FoldChange, fill = species )) +
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
    ggtitle("Heat-stress Log2FC\nof genes in Salmon Module")
  l2fc.plot.all.samples.salmon
  
  
# plot the mean abs log2FC ~ Species ED 50
ed.50.data <- read.csv("~/Lab Notebook/Chapter3/Thermal_stress_analysis/Phylogeny resolved ED50 curves/individual.specied.ED50.csv")
ed.50.data <- ed.50.data %>%
  mutate(species = ifelse(species == "Acropora sp", "Acropora nasuta", species))

all.species.yellow.df$abs.l2fc <- abs(all.species.yellow.df$log2FoldChange)
species_mean_abs_l2fc <- all.species.yellow.df %>%
  dplyr::group_by(species) %>%
  dplyr::summarize(mean_abs_l2fc = mean(abs.l2fc, na.rm = TRUE))

species_mean_abs_l2fc <- species_mean_abs_l2fc %>%
  left_join(ed.50.data %>% dplyr::select(species, mean.therm.tol), by = "species")
species_mean_abs_l2fc_yellow_mod_genes <- species_mean_abs_l2fc %>%
  distinct()

ind.therm.toler.species.l2fc.yellow <- ed.50.data %>%
  left_join(species_mean_abs_l2fc_yellow_mod_genes[,-3], by = "species")


# Spearmans cor
spearman_result_yellow <- cor.test(species_mean_abs_l2fc_yellow_mod_genes$mean.therm.tol, 
                            species_mean_abs_l2fc_yellow_mod_genes$mean_abs_l2fc, 
                            method = "spearman")


# Pearson's
pearson_result_yellow <- cor.test(ind.therm.toler.species.l2fc.yellow$ind.rel.therm.tol,
                           ind.therm.toler.species.l2fc.yellow$mean_abs_l2fc,
                           method = "pearson")

# Extract the correlation coefficient and p-value
pearson_cor_yellow <- round(pearson_result_yellow$estimate, 3)
pearson_p_yellow <- format.pval(pearson_result_yellow$p.value, digits = 3)

# # Extract the correlation coefficient and p-value
# spearman_cor_yellow <- round(spearman_result_yellow$estimate, 3)
# spearman_p_yellow <- format.pval(spearman_result_yellow$p.value, digits = 3)
# 
# yellow.module.genes.ed50 <- ggplot(species_mean_abs_l2fc_yellow_mod_genes, aes(x = mean.therm.tol, y = mean_abs_l2fc)) +
#   geom_point(size = 5, aes(color = species)) +  # Adjust the size of the points
#   scale_color_manual(values = official.species.colors) +
#   geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
#   labs(x = "Species thermal tolerance", y = "Mean abs(l2FC)", color = "Species") +  # Label axes and legend
#   theme_classic(base_size = 15) +  # Use a classic theme
#   ggtitle("Yellow module") +
#   theme(legend.position = "none") +
#   xlim(34.6, NA) +
#   annotate("text", x = Inf, y = Inf, label = paste("Spearman's ρ =", spearman_cor_yellow, "\np-value =", spearman_p_yellow),
#            hjust = 1.1, vjust = 1.1, size = 5, color = "black", fontface = "italic")
# 
# # Print the plot
# yellow.module.genes.ed50
 

yellow.module.genes.ed50 <- ggplot(ind.therm.toler.species.l2fc.yellow, aes(x = mean_abs_l2fc, y = ind.rel.therm.tol)) +
  geom_point(size = 5, aes(color = species)) +  # Adjust the size of the points
  scale_color_manual(values = official.species.colors) +
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "Mean abs(l2FC)", y ="Colony thermal tolerance", color = "Species") +  # Label axes and legend
  theme_classic() +  # Use a classic theme
  ggtitle("Yellow module") +
  theme(legend.position = "none") +
 annotate("text", x = Inf, y = Inf, label = paste("Pearsons's cor =", pearson_cor_yellow, "\np-value =", pearson_p_yellow),
          hjust = 1.1, vjust = 1.1, size = 5, color = "black", fontface = "italic")

# Print the plot
yellow.module.genes.ed50



 
# test the cor without pulchra
spearman_result_yellow_no_pul <- cor.test(species_mean_abs_l2fc_yellow_mod_genes$mean.therm.tol[-3], 
                                   species_mean_abs_l2fc_yellow_mod_genes$mean_abs_l2fc[-3], 
                                   method = "spearman")



all.species.salmon.df$abs.l2fc <- abs(all.species.salmon.df$log2FoldChange)
species_mean_abs_l2fc <- all.species.salmon.df %>%
  dplyr::group_by(species) %>%
  dplyr::summarize(mean_abs_l2fc = mean(abs.l2fc, na.rm = TRUE))

species_mean_abs_l2fc <- species_mean_abs_l2fc %>%
  left_join(ed.50.data %>% dplyr::select(species, mean.therm.tol), by = "species")
species_mean_abs_l2fc_salmon_mod_genes <- species_mean_abs_l2fc %>%
  distinct()



ind.therm.toler.species.l2fc.salmon <- ed.50.data %>%
  left_join(species_mean_abs_l2fc_salmon_mod_genes[,-3], by = "species")


# Spearmans cor

# 
# spearman.res.2.salmon <- cor.test(ind.therm.toler.species.l2fc.salmon$ind.rel.therm.tol,
#                            ind.therm.toler.species.l2fc.salmon$mean_abs_l2fc,
#                            method = "spearman")
# 
# # Extract the correlation coefficient and p-value
# spearman_cor_salmon <- round(spearman.res.2.salmon$estimate, 3)
# spearman_p_salmon <- format.pval(spearman.res.2.salmon$p.value, digits = 3)


# Spearmans cor
spearman_result_salmon <- cor.test(species_mean_abs_l2fc_salmon_mod_genes$mean.therm.tol,
                                   species_mean_abs_l2fc_salmon_mod_genes$mean_abs_l2fc,
                                   method = "spearman")

# Extract the correlation coefficient and p-value
spearman_cor_salmon <- round(spearman_result_salmon$estimate, 2)
spearman_p_salmon <- format.pval(spearman_result_salmon$p.value, digits = 2)

# Pearsons cor
pearson_result_salmon <- cor.test(ind.therm.toler.species.l2fc.salmon$ind.rel.therm.tol,
                                  ind.therm.toler.species.l2fc.salmon$mean_abs_l2fc,
                                   method = "pearson")

# Extract the correlation coefficient and p-value
pearson_cor_salmon <- round(pearson_result_salmon$estimate, 3)
pearson_p_salmon <- format.pval(pearson_result_salmon$p.value, digits = 3)


salmon.module.genes.ed50 <- ggplot(ind.therm.toler.species.l2fc.salmon, aes(x = mean_abs_l2fc, y = ind.rel.therm.tol)) +
  geom_point(size = 5, aes(color = species)) +  # Adjust the size of the points
  scale_color_manual(values = official.species.colors) +
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line without confidence intervals
  labs(x = "Mean abs(l2FC)", y = "Colony thermal tolerance", color = "Species") +  # Label axes and legend
  theme_classic() +  # Use a classic theme
  ggtitle("Salmon module") +
  theme(legend.position = "none") +
  annotate("text", x = Inf, y = Inf, label = paste("Pearsons's cor =", pearson_cor_salmon, "\np-value =", pearson_p_salmon), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black", fontface = "italic")

# Print the plot
salmon.module.genes.ed50 

yellow.module.genes.ed50 + salmon.module.genes.ed50 
  
library(patchwork)
l2fc.plot.all.samples.yel + l2fc.plot.all.samples.salmon + yellow.module.genes.ed50 + salmon.module.genes.ed50 + plot_annotation(tag_levels = list(c('A.', 'B.', 'C.', 'D.')))
  
## Species base mean expression in yellow and salmon modules
bmean.plot.all.samples.salmon <- ggplot(all.species.salmon.df, aes(x = species, y = log(baseMean), fill = species )) +
  # clouds
  introdataviz::geom_flat_violin(trim=FALSE,
                                 position = position_nudge(x = rain_height+.05),
                                 width = 1) +
  # # rain
  geom_point(aes(color = species), size = 1, alpha = 0.4,  show.legend = FALSE, 
             position = position_jitter(width = rain_height, height = 0)) +
  # boxplots
  geom_boxplot(width = rain_height, show.legend = FALSE, 
               outlier.shape = NA,
               position = position_nudge(x = -rain_height*2)) +
 # coord_flip() +
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
  ggtitle("Log of Base Mean\nExpression Salmon Module")
bmean.plot.all.samples.salmon

bmean.plot.all.samples.yellow <- ggplot(all.species.yellow.df, aes(x = species, y = log(baseMean), fill = species )) +
  # clouds
  introdataviz::geom_flat_violin(trim=FALSE,
                                 position = position_nudge(x = rain_height+.05),
                                 width = 1) +
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
  ggtitle("Log of Base Mean\nExpression Yellow Module")
bmean.plot.all.samples.yellow

(l2fc.plot.all.samples.yel | l2fc.plot.all.samples.salmon ) /
(bmean.plot.all.samples.yellow | bmean.plot.all.samples.salmon) 

ggsave( filename = "~/Lab Notebook/Chapter3/Figures/figure6.png", plot = ((l2fc.plot.all.samples.yel | l2fc.plot.all.samples.salmon ) /
(yellow.module.genes.ed50 | salmon.module.genes.ed50) /
  (bmean.plot.all.samples.yellow | bmean.plot.all.samples.salmon) + plot_annotation(tag_levels = list(c('A.', 'B.', 'C.', 'D.', 'E.', 'F.')))), 
width = 12, height = 12, dpi = 1500)
   


anova_result <- aov(baseMean ~ species, data = all.species.yellow.df)
summary(anova_result)
tukey_hsd_result <- TukeyHSD(anova_result)
tukey_hsd_result
plot(tukey_hsd_result)   
library(multcompView)
tukey_df <- as.data.frame(tukey_hsd_result$species)

# Add a column for comparisons
tukey_df$Comparison <- rownames(tukey_df)

# Add a column for alpha value based on CI overlap with zero
tukey_df$alpha <- ifelse(tukey_df$lwr <= 0 & tukey_df$upr >= 0, 0.2, 1)

# Order the data frame for better plotting
tukey_df <- tukey_df %>%
  arrange(desc(diff)) %>%
  mutate(Comparison = factor(Comparison, levels = Comparison))

# Plot the Tukey HSD results using ggplot2
species.comp.basemean.yellow.module <- ggplot(tukey_df, aes(x = Comparison, y = diff, alpha = alpha)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  geom_point(size = 3) +
  coord_flip() +
  theme_minimal(base_size = 15) +
  labs(title = "Diffference in mean expression: Yellow Module",
       x = "Species Comparison",
       y = "Difference in Means",
       caption = "Error bars represent 95% confidence intervals") +
  scale_alpha_identity() + # Use the alpha column for transparency
  geom_hline(yintercept = 0, color = "red") 
species.comp.basemean.yellow.module


anova_result <- aov(baseMean ~ species, data = all.species.salmon.df)
summary(anova_result)   
# not significant
tukey_hsd_result <- TukeyHSD(anova_result)
tukey_hsd_result
plot(tukey_hsd_result)   
library(multcompView)
tukey_df <- as.data.frame(tukey_hsd_result$species)

# Add a column for comparisons
tukey_df$Comparison <- rownames(tukey_df)

# Add a column for alpha value based on CI overlap with zero
tukey_df$alpha <- ifelse(tukey_df$lwr <= 0 & tukey_df$upr >= 0, 0.2, 1)

# Order the data frame for better plotting
tukey_df <- tukey_df %>%
  arrange(desc(diff)) %>%
  mutate(Comparison = factor(Comparison, levels = Comparison))

# Plot the Tukey HSD results using ggplot2
species.comp.basemean.salmon.module <- ggplot(tukey_df, aes(x = Comparison, y = diff, alpha = alpha)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  geom_point(size = 3) +
  coord_flip() +
  theme_minimal(base_size = 15) +
  labs(title = "Diffference in mean expression: Salmon Module",
       x = "Species Comparison",
       y = "Difference in Means",
       caption = "Error bars represent 95% confidence intervals") +
  scale_alpha_identity() +
  geom_hline(yintercept = 0, color = "red") 
species.comp.basemean.salmon.module

species.comp.basemean.yellow.module/ species.comp.basemean.salmon.module
ggsave("~/Lab Notebook/Chapter3/Figures/Tukey.basemean.species.comp.jpeg", plot = (species.comp.basemean.yellow.module/ species.comp.basemean.salmon.module), width = 14, height = 16, units = "in", dpi = 600)

# all.species.yellow.df$species.ed.50 <- NA
# all.species.salmon.df$species.ed.50 <- NA
# # Match species between the two data frames
# matching_indices <- match(all.species.yellow.df$species, ed.50.data$species)
# matching_indices <- match(all.species.salmon.df$species, ed.50.data$species)
# 
# # Replace NA values in 'species.ed.50' column with corresponding values from ed.50.data
# all.species.yellow.df$species.ed.50 <- ifelse(is.na(matching_indices), NA, ed.50.data$mean.therm.tol[matching_indices])
# all.species.salmon.df$species.ed.50 <- ifelse(is.na(matching_indices), NA, ed.50.data$mean.therm.tol[matching_indices])
# 
# 
# # Test if the more thermally tolerant species has higher baseMean expression:
# anova(lm(baseMean~species.ed.50, data = all.species.yellow.df))
# anova(lm(baseMean~species.ed.50, data = all.species.salmon.df))
# summary(lm(baseMean~species, data = all.species.yellow.df))

anova_result <- aov(abs.l2fc ~ species, data = all.species.yellow.df)
summary(anova_result)
tukey_hsd_result <- TukeyHSD(anova_result)
tukey_hsd_result
plot(tukey_hsd_result)   
library(multcompView)
tukey_df <- as.data.frame(tukey_hsd_result$species)

# Add a column for comparisons
tukey_df$Comparison <- rownames(tukey_df)

# Add a column for alpha value based on CI overlap with zero
tukey_df$alpha <- ifelse(tukey_df$lwr <= 0 & tukey_df$upr >= 0, 0.2, 1)

# Order the data frame for better plotting
tukey_df <- tukey_df %>%
  arrange(desc(diff)) %>%
  mutate(Comparison = factor(Comparison, levels = Comparison))

# Plot the Tukey HSD results using ggplot2
species.comp.abs.l2fc.yellow.module <- ggplot(tukey_df, aes(x = Comparison, y = diff, alpha = alpha)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  geom_point(size = 3) +
  coord_flip() +
  theme_minimal(base_size = 15) +
  labs(title = "Diffference in GE Response Magnitude: Yellow Module",
       x = "Species Comparison",
       y = "Difference in Means",
       caption = "Error bars represent 95% confidence intervals") +
  scale_alpha_identity() + # Use the alpha column for transparency
  geom_hline(yintercept = 0, color = "red") 
species.comp.abs.l2fc.yellow.module


anova_result <- aov(abs.l2fc ~ species, data = all.species.salmon.df)
summary(anova_result)   
# not significant
tukey_hsd_result <- TukeyHSD(anova_result)
tukey_hsd_result
plot(tukey_hsd_result)   
library(multcompView)
tukey_df <- as.data.frame(tukey_hsd_result$species)

# Add a column for comparisons
tukey_df$Comparison <- rownames(tukey_df)

# Add a column for alpha value based on CI overlap with zero
tukey_df$alpha <- ifelse(tukey_df$lwr <= 0 & tukey_df$upr >= 0, 0.2, 1)

# Order the data frame for better plotting
tukey_df <- tukey_df %>%
  arrange(desc(diff)) %>%
  mutate(Comparison = factor(Comparison, levels = Comparison))

# Plot the Tukey HSD results using ggplot2
species.comp.absl2fc.salmon.module <- ggplot(tukey_df, aes(x = Comparison, y = diff, alpha = alpha)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  geom_point(size = 3) +
  coord_flip() +
  theme_minimal(base_size = 15) +
  labs(title = "Diffference in GE response magnitude: Salmon Module",
       x = "Species Comparison",
       y = "Difference in Means",
       caption = "Error bars represent 95% confidence intervals") +
  scale_alpha_identity() +
  geom_hline(yintercept = 0, color = "red") 
species.comp.absl2fc.salmon.module

species.comp.abs.l2fc.yellow.module/ species.comp.absl2fc.salmon.module  

#ggsave("~/Lab Notebook/Chapter3/Figures/Tukey.absl2FC.species.comp.jpeg", plot = (species.comp.abs.l2fc.yellow.module/ species.comp.absl2fc.salmon.module), width = 14, height = 16, units = "in", dpi = 600)

# Testing Front loading with control samples only
## Create dds object for each species with no design
### A. abrotanoides
abro.files.list <- paste(filter(my.metadata, species == "Acropora.abrotanoides" & treatment == "control")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
abro.metadata <- filter(my.metadata, species == "Acropora.abrotanoides" & treatment == "control")
abro.sampleNames <- abro.metadata$seq.sample.id

# Create sample table
abro.sampleTable <- data.frame(sampleName = abro.sampleNames, 
                               fileName = abro.files.list, 
                               condition = abro.metadata)
# Convert variables to factors
colnames(abro.sampleTable)

# Create DESeqDataSet Object
abro.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = abro.sampleTable,
  directory = my.dir,
  design = ~ 1)
# Pre-filter
keep <- rowSums(counts(abro.ddsHTSeq)) >= 10
abro.ddsHTSeq <- abro.ddsHTSeq[keep,]

# Estimate normalization factors
abro.ddsHTSeq <- estimateSizeFactors(abro.ddsHTSeq)
abro.ddsHTSeq <- estimateDispersions(abro.ddsHTSeq)

# Differential expression
abro.dds <- DESeq(abro.ddsHTSeq) 
resultsNames(abro.dds)
nrow(abro.dds)
# 13658

# Heat stress response
names(abro.dds@colData) # condition.treatment
resultsNames(abro.dds)
abro.basemean.results <- results(abro.dds)

abro.basemean.results.df <-abro.basemean.results [complete.cases(abro.basemean.results ),]
abro.basemean.results.df <- as.data.frame(abro.basemean.results)
# Isolate genes in the yellow and salmon modules
# MEsalmon.genes$gene <- colnames(datExpr)[moduleColors=="salmon"]
# MEyellow.genes$gene <-colnames(datExpr)[moduleColors=="yellow"]

abro.me.salmon <- rownames(abro.basemean.results.df) %in% MEsalmon.genes$gene
abro.me.yellow <- rownames(abro.basemean.results.df) %in% MEyellow.genes$gene

sum(abro.me.salmon) #295
sum(abro.me.yellow) #1482
# Add the me.salmon column to the dataframe
abro.basemean.results.df$me.salmon <- abro.me.salmon
abro.basemean.results.df$me.yellow <- abro.me.yellow
abro.basemean.results.df$species <- "Acropora abrotanoides"

### A. hyacinthus
hya.files.list <- paste(filter(my.metadata, species == "Acropora.hyacinthus" & treatment == "control")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
hya.metadata <- filter(my.metadata, species == "Acropora.hyacinthus" & treatment == "control")
hya.sampleNames <- hya.metadata$seq.sample.id

# Create sample table
hya.sampleTable <- data.frame(sampleName = hya.sampleNames, 
                               fileName = hya.files.list, 
                               condition = hya.metadata)
# Convert variables to factors
colnames(hya.sampleTable)

# Create DESeqDataSet Object
hya.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = hya.sampleTable,
  directory = my.dir,
  design = ~ 1)
# Pre-filter
keep <- rowSums(counts(hya.ddsHTSeq)) >= 10
hya.ddsHTSeq <- hya.ddsHTSeq[keep,]

# Estimate normalization factors
hya.ddsHTSeq <- estimateSizeFactors(hya.ddsHTSeq)
hya.ddsHTSeq <- estimateDispersions(hya.ddsHTSeq)

# Differential expression
hya.dds <- DESeq(hya.ddsHTSeq) 
resultsNames(hya.dds)
nrow(hya.dds)
# 15615

# Heat stress response
names(hya.dds@colData) # condition.treatment
resultsNames(hya.dds)
hya.basemean.results <- results(hya.dds)

hya.basemean.results.df <-hya.basemean.results [complete.cases(hya.basemean.results ),]
hya.basemean.results.df <- as.data.frame(hya.basemean.results)
# Isolate genes in the yellow and salmon modules
# MEsalmon.genes$gene <- colnames(datExpr)[moduleColors=="salmon"]
# MEyellow.genes$gene <-colnames(datExpr)[moduleColors=="yellow"]

hya.me.salmon <- rownames(hya.basemean.results.df) %in% MEsalmon.genes$gene
hya.me.yellow <- rownames(hya.basemean.results.df) %in% MEyellow.genes$gene

sum(hya.me.salmon) #313
sum(hya.me.yellow) #1639
# Add the me.salmon column to the dataframe
hya.basemean.results.df$me.salmon <- hya.me.salmon
hya.basemean.results.df$me.yellow <- hya.me.yellow
hya.basemean.results.df$species <- "Acropora hyacinthus"


### A. pulchra
pul.files.list <- paste(filter(my.metadata, species == "Acropora.pulchra" & treatment == "control")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
pul.metadata <- filter(my.metadata, species == "Acropora.pulchra" & treatment == "control")
pul.sampleNames <- pul.metadata$seq.sample.id

# Create sample table
pul.sampleTable <- data.frame(sampleName = pul.sampleNames, 
                               fileName = pul.files.list, 
                               condition = pul.metadata)
# Convert variables to factors
colnames(pul.sampleTable)

# Create DESeqDataSet Object
pul.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = pul.sampleTable,
  directory = my.dir,
  design = ~ 1)
# Pre-filter
keep <- rowSums(counts(pul.ddsHTSeq)) >= 10
pul.ddsHTSeq <- pul.ddsHTSeq[keep,]

# Estimate normalization factors
pul.ddsHTSeq <- estimateSizeFactors(pul.ddsHTSeq)
pul.ddsHTSeq <- estimateDispersions(pul.ddsHTSeq)

# Differential expression
pul.dds <- DESeq(pul.ddsHTSeq) 
resultsNames(pul.dds)
nrow(pul.dds)
# 14772

# Heat stress response
names(pul.dds@colData) # condition.treatment
resultsNames(pul.dds)
pul.basemean.results <- results(pul.dds)

pul.basemean.results.df <-pul.basemean.results [complete.cases(pul.basemean.results ),]
pul.basemean.results.df <- as.data.frame(pul.basemean.results)
# Isolate genes in the yellow and salmon modules
# MEsalmon.genes$gene <- colnames(datExpr)[moduleColors=="salmon"]
# MEyellow.genes$gene <-colnames(datExpr)[moduleColors=="yellow"]

pul.me.salmon <- rownames(pul.basemean.results.df) %in% MEsalmon.genes$gene
pul.me.yellow <- rownames(pul.basemean.results.df) %in% MEyellow.genes$gene

sum(pul.me.salmon) #297
sum(pul.me.yellow) #1603
# Add the me.salmon column to the dataframe
pul.basemean.results.df$me.salmon <- pul.me.salmon
pul.basemean.results.df$me.yellow <- pul.me.yellow
pul.basemean.results.df$species <- "Acropora pulchra"


### A. lutkeni
lut.files.list <- paste(filter(my.metadata, species == "Acropora.lutkeni" & treatment == "control")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
lut.metadata <- filter(my.metadata, species == "Acropora.lutkeni" & treatment == "control")
lut.sampleNames <- lut.metadata$seq.sample.id

# Create sample table
lut.sampleTable <- data.frame(sampleName = lut.sampleNames, 
                               fileName = lut.files.list, 
                               condition = lut.metadata)
# Convert variables to factors
colnames(lut.sampleTable)

# Create DESeqDataSet Object
lut.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = lut.sampleTable,
  directory = my.dir,
  design = ~ 1)
# Pre-filter
keep <- rowSums(counts(lut.ddsHTSeq)) >= 10
lut.ddsHTSeq <- lut.ddsHTSeq[keep,]

# Estimate normalization factors
lut.ddsHTSeq <- estimateSizeFactors(lut.ddsHTSeq)
lut.ddsHTSeq <- estimateDispersions(lut.ddsHTSeq)

# Differential expression
lut.dds <- DESeq(lut.ddsHTSeq) 
resultsNames(lut.dds)
nrow(lut.dds)
# 13568

# Heat stress response
names(lut.dds@colData) # condition.treatment
resultsNames(lut.dds)
lut.basemean.results <- results(lut.dds)

lut.basemean.results.df <-lut.basemean.results [complete.cases(lut.basemean.results ),]
lut.basemean.results.df <- as.data.frame(lut.basemean.results)
# Isolate genes in the yellow and salmon modules
# MEsalmon.genes$gene <- colnames(datExpr)[moduleColors=="salmon"]
# MEyellow.genes$gene <-colnames(datExpr)[moduleColors=="yellow"]

lut.me.salmon <- rownames(lut.basemean.results.df) %in% MEsalmon.genes$gene
lut.me.yellow <- rownames(lut.basemean.results.df) %in% MEyellow.genes$gene

sum(lut.me.salmon) #282
sum(lut.me.yellow) #1417
# Add the me.salmon column to the dataframe
lut.basemean.results.df$me.salmon <- lut.me.salmon
lut.basemean.results.df$me.yellow <- lut.me.yellow
lut.basemean.results.df$species <- "Acropora lutkeni"

### Acropora retusa
ret.files.list <- paste(filter(my.metadata, species == "Acropora.retusa" & treatment == "control")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
ret.metadata <- filter(my.metadata, species == "Acropora.retusa" & treatment == "control")
ret.sampleNames <- ret.metadata$seq.sample.id

# Create sample table
ret.sampleTable <- data.frame(sampleName = ret.sampleNames, 
                               fileName = ret.files.list, 
                               condition = ret.metadata)
# Convert variables to factors
colnames(ret.sampleTable)

# Create DESeqDataSet Object
ret.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = ret.sampleTable,
  directory = my.dir,
  design = ~ 1)
# Pre-filter
keep <- rowSums(counts(ret.ddsHTSeq)) >= 10
ret.ddsHTSeq <- ret.ddsHTSeq[keep,]

# Estimate normalization factors
ret.ddsHTSeq <- estimateSizeFactors(ret.ddsHTSeq)
ret.ddsHTSeq <- estimateDispersions(ret.ddsHTSeq)

# Differential expression
ret.dds <- DESeq(ret.ddsHTSeq) 
resultsNames(ret.dds)
nrow(ret.dds)
# 14320

# Heat stress response
names(ret.dds@colData) # condition.treatment
resultsNames(ret.dds)
ret.basemean.results <- results(ret.dds)

ret.basemean.results.df <-ret.basemean.results [complete.cases(ret.basemean.results ),]
ret.basemean.results.df <- as.data.frame(ret.basemean.results)
# Isolate genes in the yellow and salmon modules
# MEsalmon.genes$gene <- colnames(datExpr)[moduleColors=="salmon"]
# MEyellow.genes$gene <-colnames(datExpr)[moduleColors=="yellow"]

ret.me.salmon <- rownames(ret.basemean.results.df) %in% MEsalmon.genes$gene
ret.me.yellow <- rownames(ret.basemean.results.df) %in% MEyellow.genes$gene

sum(ret.me.salmon) #291
sum(ret.me.yellow) #1515
# Add the me.salmon column to the dataframe
ret.basemean.results.df$me.salmon <- ret.me.salmon
ret.basemean.results.df$me.yellow <- ret.me.yellow
ret.basemean.results.df$species <- "Acropora retusa"

### Acropora robusta
rob.files.list <- paste(filter(my.metadata, species == "Acropora.robusta" & treatment == "control")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
rob.metadata <- filter(my.metadata, species == "Acropora.robusta" & treatment == "control")
rob.sampleNames <- rob.metadata$seq.sample.id

# Create sample table
rob.sampleTable <- data.frame(sampleName = rob.sampleNames, 
                               fileName = rob.files.list, 
                               condition = rob.metadata)
# Convert variables to factors
colnames(rob.sampleTable)

# Create DESeqDataSet Object
rob.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = rob.sampleTable,
  directory = my.dir,
  design = ~ 1)
# Pre-filter
keep <- rowSums(counts(rob.ddsHTSeq)) >= 10
rob.ddsHTSeq <- rob.ddsHTSeq[keep,]

# Estimate normalization factors
rob.ddsHTSeq <- estimateSizeFactors(rob.ddsHTSeq)
rob.ddsHTSeq <- estimateDispersions(rob.ddsHTSeq)

# Differential expression
rob.dds <- DESeq(rob.ddsHTSeq) 
resultsNames(rob.dds)
nrow(rob.dds)
# 14437

# Heat stress response
names(rob.dds@colData) # condition.treatment
resultsNames(rob.dds)
rob.basemean.results <- results(rob.dds)

rob.basemean.results.df <-rob.basemean.results [complete.cases(rob.basemean.results ),]
rob.basemean.results.df <- as.data.frame(rob.basemean.results)
# Isolate genes in the yellow and salmon modules
# MEsalmon.genes$gene <- colnames(datExpr)[moduleColors=="salmon"]
# MEyellow.genes$gene <-colnames(datExpr)[moduleColors=="yellow"]

rob.me.salmon <- rownames(rob.basemean.results.df) %in% MEsalmon.genes$gene
rob.me.yellow <- rownames(rob.basemean.results.df) %in% MEyellow.genes$gene

sum(rob.me.salmon) #300
sum(rob.me.yellow) #1532
# Add the me.salmon column to the dataframe
rob.basemean.results.df$me.salmon <- rob.me.salmon
rob.basemean.results.df$me.yellow <- rob.me.yellow
rob.basemean.results.df$species <- "Acropora robusta"

### Acropora tutuilensis
tut.files.list <- paste(filter(my.metadata, species == "Acropora.tutuilensis" & treatment == "control")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
tut.metadata <- filter(my.metadata, species == "Acropora.tutuilensis" & treatment == "control")
tut.sampleNames <- tut.metadata$seq.sample.id

# Create sample table
tut.sampleTable <- data.frame(sampleName = tut.sampleNames, 
                               fileName = tut.files.list, 
                               condition = tut.metadata)
# Convert variables to factors
colnames(tut.sampleTable)

# Create DESeqDataSet Object
tut.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = tut.sampleTable,
  directory = my.dir,
  design = ~ 1)
# Pre-filter
keep <- rowSums(counts(tut.ddsHTSeq)) >= 10
tut.ddsHTSeq <- tut.ddsHTSeq[keep,]

# Estimate normalization factors
tut.ddsHTSeq <- estimateSizeFactors(tut.ddsHTSeq)
tut.ddsHTSeq <- estimateDispersions(tut.ddsHTSeq)

# Differential expression
tut.dds <- DESeq(tut.ddsHTSeq) 
resultsNames(tut.dds)
nrow(tut.dds)
# 16630

# Heat stress response
names(tut.dds@colData) # condition.treatment
resultsNames(tut.dds)
tut.basemean.results <- results(tut.dds)

tut.basemean.results.df <-tut.basemean.results [complete.cases(tut.basemean.results ),]
tut.basemean.results.df <- as.data.frame(tut.basemean.results)
# Isolate genes in the yellow and salmon modules
# MEsalmon.genes$gene <- colnames(datExpr)[moduleColors=="salmon"]
# MEyellow.genes$gene <-colnames(datExpr)[moduleColors=="yellow"]

tut.me.salmon <- rownames(tut.basemean.results.df) %in% MEsalmon.genes$gene
tut.me.yellow <- rownames(tut.basemean.results.df) %in% MEyellow.genes$gene

sum(tut.me.salmon) #318
sum(tut.me.yellow) #1746
# Add the me.salmon column to the dataframe
tut.basemean.results.df$me.salmon <- tut.me.salmon
tut.basemean.results.df$me.yellow <- tut.me.yellow
tut.basemean.results.df$species <- "Acropora tutuilensis"


### Acrpora nasuta
nas.files.list <- paste(filter(my.metadata, species == "Acropora.nasuta" & treatment == "control")[,5],"_gene_count.txt", sep = "")

# Filter metadata and create sample table
nas.metadata <- filter(my.metadata, species == "Acropora.nasuta" & treatment == "control")
nas.sampleNames <- nas.metadata$seq.sample.id

# Create sample table
nas.sampleTable <- data.frame(sampleName = nas.sampleNames, 
                               fileName = nas.files.list, 
                               condition = nas.metadata)
# Convert variables to factors
colnames(nas.sampleTable)

# Create DESeqDataSet Object
nas.ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = nas.sampleTable,
  directory = my.dir,
  design = ~ 1)
# Pre-filter
keep <- rowSums(counts(nas.ddsHTSeq)) >= 10
nas.ddsHTSeq <- nas.ddsHTSeq[keep,]

# Estimate normalization factors
nas.ddsHTSeq <- estimateSizeFactors(nas.ddsHTSeq)
nas.ddsHTSeq <- estimateDispersions(nas.ddsHTSeq)

# Differential expression
nas.dds <- DESeq(nas.ddsHTSeq) 
resultsNames(nas.dds)
nrow(nas.dds)
# 12835

# Heat stress response
names(nas.dds@colData) # condition.treatment
resultsNames(nas.dds)
nas.basemean.results <- results(nas.dds)

nas.basemean.results.df <-nas.basemean.results [complete.cases(nas.basemean.results ),]
nas.basemean.results.df <- as.data.frame(nas.basemean.results)
# Isolate genes in the yellow and salmon modules
# MEsalmon.genes$gene <- colnames(datExpr)[moduleColors=="salmon"]
# MEyellow.genes$gene <-colnames(datExpr)[moduleColors=="yellow"]

nas.me.salmon <- rownames(nas.basemean.results.df) %in% MEsalmon.genes$gene
nas.me.yellow <- rownames(nas.basemean.results.df) %in% MEyellow.genes$gene

sum(nas.me.salmon) #265
sum(nas.me.yellow) #1374
# Add the me.salmon column to the dataframe
nas.basemean.results.df$me.salmon <- nas.me.salmon
nas.basemean.results.df$me.yellow <- nas.me.yellow
nas.basemean.results.df$species <- "Acropora nasuta"


# Combine data frames for plotting salmon gene expression
salmon.abro <- abro.basemean.results.df[abro.basemean.results.df$me.salmon,]
salmon.hya <- hya.basemean.results.df[hya.basemean.results.df$me.salmon,]
salmon.pul <- pul.basemean.results.df[pul.basemean.results.df$me.salmon,]
salmon.lut <- lut.basemean.results.df[lut.basemean.results.df$me.salmon,]
salmon.ret <- ret.basemean.results.df[ret.basemean.results.df$me.salmon,]
salmon.rob <- rob.basemean.results.df[rob.basemean.results.df$me.salmon,]
salmon.tut <- tut.basemean.results.df[tut.basemean.results.df$me.salmon,]
salmon.nas <- nas.basemean.results.df[nas.basemean.results.df$me.salmon, ]


# Combine data frames for plotting yellow gene expression
yellow.abro <- abro.basemean.results.df[abro.basemean.results.df$me.yellow,]
yellow.hya <- hya.basemean.results.df[hya.basemean.results.df$me.yellow,]
yellow.pul <- pul.basemean.results.df[pul.basemean.results.df$me.yellow,]
yellow.lut <- lut.basemean.results.df[lut.basemean.results.df$me.yellow,]
yellow.ret <- ret.basemean.results.df[ret.basemean.results.df$me.yellow,]
yellow.rob <- rob.basemean.results.df[rob.basemean.results.df$me.yellow,]
yellow.tut <- tut.basemean.results.df[tut.basemean.results.df$me.yellow,]
yellow.nas <- nas.basemean.results.df[nas.basemean.results.df$me.yellow, ]



yellow.df.list <- list(yellow.abro, yellow.hya, yellow.pul, yellow.lut,yellow.ret, yellow.rob, yellow.tut, yellow.nas)
all.species.yellow.df <- do.call(rbind, yellow.df.list)
official.species.colors <- c("Acropora tutuilensis" = "#F7C11E",
                             "Acropora lutkeni" = "#FF9A17",
                             "Acropora pulchra" = "#fb5607",
                             "Acropora retusa" = "#ff006e",
                             "Acropora nasuta" = "#c11cad",
                             "Acropora robusta" = "#8338ec",
                             "Acropora hyacinthus" = "#3E6FCB",
                             "Acropora abrotanoides" = "#3a86ff")

species.order <- c("Acropora tutuilensis", "Acropora lutkeni",
                   "Acropora pulchra", "Acropora retusa","Acropora nasuta",
                   "Acropora robusta","Acropora hyacinthus","Acropora abrotanoides")
all.species.yellow.df$species <- factor(all.species.yellow.df$species, levels = species.order)

salmon.df.list <- list(salmon.abro, salmon.hya, salmon.pul, salmon.lut,salmon.ret, salmon.rob, salmon.tut, salmon.nas)
all.species.salmon.df <- do.call(rbind, salmon.df.list)
official.species.colors <- c("Acropora tutuilensis" = "#F7C11E",
                             "Acropora lutkeni" = "#FF9A17",
                             "Acropora pulchra" = "#fb5607",
                             "Acropora retusa" = "#ff006e",
                             "Acropora nasuta" = "#c11cad",
                             "Acropora robusta" = "#8338ec",
                             "Acropora hyacinthus" = "#3E6FCB",
                             "Acropora abrotanoides" = "#3a86ff")

species.order <- c("Acropora tutuilensis", "Acropora lutkeni",
                   "Acropora pulchra", "Acropora retusa","Acropora nasuta",
                   "Acropora robusta","Acropora hyacinthus","Acropora abrotanoides")
all.species.salmon.df$species <- factor(all.species.salmon.df$species, levels = species.order)


## Species base mean expression in yellow and salmon modules
bmean.plot.all.samples.yellow <- ggplot(all.species.yellow.df, aes(x = species, y = log(baseMean), fill = species )) +
  # clouds
  introdataviz::geom_flat_violin(trim=FALSE,
                                 position = position_nudge(x = rain_height+.05),
                                 width = 1) +
  # # rain
  geom_point(aes(color = species), size = 1, alpha = 0.4,  show.legend = FALSE, 
             position = position_jitter(width = rain_height, height = 0)) +
  # boxplots
  geom_boxplot(width = rain_height, show.legend = FALSE, 
               outlier.shape = NA,
               position = position_nudge(x = -rain_height*2)) +
 # coord_flip() +
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
  ggtitle("Log of Base Mean\nExpression Yellow Module")
bmean.plot.all.samples.yellow

bmean.plot.all.samples.salmon <- ggplot(all.species.salmon.df, aes(x = species, y = log(baseMean), fill = species )) +
  # clouds
  introdataviz::geom_flat_violin(trim=FALSE,
                                 position = position_nudge(x = rain_height+.05),
                                 width = 1) +
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
  ggtitle("Log of Base Mean\nExpression Salmon Module")
bmean.plot.all.samples.salmon




anova_result <- aov(baseMean ~ species, data = all.species.yellow.df)
summary(anova_result)
tukey_hsd_result <- TukeyHSD(anova_result)
tukey_hsd_result
plot(tukey_hsd_result)   
library(multcompView)
tukey_df <- as.data.frame(tukey_hsd_result$species)

# Add a column for comparisons
tukey_df$Comparison <- rownames(tukey_df)

# Add a column for alpha value based on CI overlap with zero
tukey_df$alpha <- ifelse(tukey_df$lwr <= 0 & tukey_df$upr >= 0, 0.2, 1)

# Order the data frame for better plotting
tukey_df <- tukey_df %>%
  arrange(desc(diff)) %>%
  mutate(Comparison = factor(Comparison, levels = Comparison))

# Plot the Tukey HSD results using ggplot2
species.comp.basemean.yellow.module <- ggplot(tukey_df, aes(x = Comparison, y = diff, alpha = alpha)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  geom_point(size = 3) +
  coord_flip() +
  theme_minimal(base_size = 15) +
  labs(title = "Diffference in mean expression: Yellow Module",
       x = "Species Comparison",
       y = "Difference in Means",
       caption = "Error bars represent 95% confidence intervals") +
  scale_alpha_identity() + # Use the alpha column for transparency
  geom_hline(yintercept = 0, color = "red") 
species.comp.basemean.yellow.module


anova_result <- aov(baseMean ~ species, data = all.species.salmon.df)
summary(anova_result)   
# not significant
tukey_hsd_result <- TukeyHSD(anova_result)
tukey_hsd_result
plot(tukey_hsd_result)   
library(multcompView)
tukey_df <- as.data.frame(tukey_hsd_result$species)

# Add a column for comparisons
tukey_df$Comparison <- rownames(tukey_df)

# Add a column for alpha value based on CI overlap with zero
tukey_df$alpha <- ifelse(tukey_df$lwr <= 0 & tukey_df$upr >= 0, 0.2, 1)

# Order the data frame for better plotting
tukey_df <- tukey_df %>%
  arrange(desc(diff)) %>%
  mutate(Comparison = factor(Comparison, levels = Comparison))

# Plot the Tukey HSD results using ggplot2
species.comp.basemean.salmon.module <- ggplot(tukey_df, aes(x = Comparison, y = diff, alpha = alpha)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  geom_point(size = 3) +
  coord_flip() +
  theme_minimal(base_size = 15) +
  labs(title = "Diffference in mean expression: Salmon Module",
       x = "Species Comparison",
       y = "Difference in Means",
       caption = "Error bars represent 95% confidence intervals") +
  scale_alpha_identity() +
  geom_hline(yintercept = 0, color = "red") 
species.comp.basemean.salmon.module

species.comp.basemean.yellow.module/ species.comp.basemean.salmon.module
ggsave("~/Lab Notebook/Chapter3/Figures/Tukey.basemean.species.comp.jpeg", plot = (species.comp.basemean.yellow.module/ species.comp.basemean.salmon.module), width = 14, height = 16, units = "in", dpi = 600)
