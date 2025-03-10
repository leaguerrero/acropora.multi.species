# Chapter 3 GO terms of salmon and yellow module genes
library(GenomicRanges)
library(topGO)
library(tidyverse)
library(readxl)
library(DESeq2)

# Gene Ontology (GO) terms for yellow and salmon module genes

# Set-up
# Read in annotation table and reformat GO column
annos <- read_excel("~/Lab Notebook/Chapter3/Bioinformatics/Eggnog_outs/out.emapper.annotations.xlsx")
# annos$GOlist <- str_replace_all(annos$gene_ontology_blast,"(GO:[0-9]*)\\^.*?\\`","\\1,")
# annos$GOlist <- str_replace_all(annos$GOlist,"\\^.*?$","")
# GOmap <- annos[,c("transcript_id","GOlist")]
# # write.table(GOmap,"~/Lab Notebook/Chapter2/Data_analysis/Amil_GOmap.txt",quote=F,row.names=F,col.names=F,sep="\t")
# geneID2GO <- readMappings(file="~/Lab Notebook/Chapter2/Data_analysis/Amil_GOmap.txt",sep="\t",IDsep=",")

annos$transcript_id <- annos$query
annos$transcript_id <- gsub("-RA$", "", annos$transcript_id)

GOmap <- annos[,c("transcript_id","GOs")]
GOmap <- GOmap %>%
  mutate(GOs = na_if(GOs, "-"))
# write.table(GOmap,"~/Lab Notebook/Chapter3/Tag_seq_analysis/Ahya_GOmap.txt",quote=F,row.names=F,col.names=F,sep="\t")
geneID2GO <- readMappings(file="~/Lab Notebook/Chapter3/Tag_seq_analysis/Ahya_GOmap.txt",sep="\t",IDsep=",")

# Create background gene set:
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

total.set.gene.names <- rownames(assay(dds))
total.set.gene.names.df <- data.frame(gene = total.set.gene.names)
# write.csv(total.set.gene.names.df, 
#           file = '~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/total.set.gene.names.csv', 
#           row.names = FALSE)

# Read in set of all genes and make 'background' gene set
all <- read.csv("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/total.set.gene.names.csv")
allgenes <- all$gene

# Load Yellow module genes
yellow.module.genes <- read.csv("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/MEyellow.genes.csv")

# Modify yellow module gene set for topGO analysis
ahyacinthus.yel.mod.IG <- factor(as.numeric(allgenes%in%yellow.module.genes$gene))
names(ahyacinthus.yel.mod.IG) <- allgenes

##BP
GOahyacinthus <- new("topGOdata",ontology="BP",allGenes=ahyacinthus.yel.mod.IG, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
ahyacinthusFisher <- getSigGroups(GOahyacinthus,test.stat)
ahyacinthus.res <- GenTable(GOahyacinthus,classic=ahyacinthusFisher,topNodes=length(ahyacinthusFisher@score),numChar=100)
ahyacinthus.filt <- ahyacinthus.res[ahyacinthus.res$classic<0.01 & ahyacinthus.res$Significant>=10,]
# write.csv(ahyacinthus.filt,paste("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/topGO/ahyacinthus.yellow.mod.genes.GO_BP.csv",sep=""))

##MF
GOahyacinthus <- new("topGOdata",ontology="MF",allGenes=ahyacinthus.yel.mod.IG, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
ahyacinthusFisher <- getSigGroups(GOahyacinthus,test.stat)
ahyacinthus.res <- GenTable(GOahyacinthus,classic=ahyacinthusFisher,topNodes=length(ahyacinthusFisher@score), numChar=100)
ahyacinthus.filt <- ahyacinthus.res[ahyacinthus.res$classic<0.01 & ahyacinthus.res$Significant>=10,]
# write.csv(ahyacinthus.filt,paste("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/topGO/ahyacinthus.yellow.mod.genes.GO_MF.csv",sep=""))

##CC
GOahyacinthus <- new("topGOdata",ontology="CC",allGenes=ahyacinthus.yel.mod.IG, nodeSize=10,
                    annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
ahyacinthusFisher <- getSigGroups(GOahyacinthus,test.stat)
ahyacinthus.res <- GenTable(GOahyacinthus,classic=ahyacinthusFisher,topNodes=length(ahyacinthusFisher@score),numChar=100)
ahyacinthus.filt <- ahyacinthus.res[ahyacinthus.res$classic<0.01 & ahyacinthus.res$Significant>=10,]
# write.csv(ahyacinthus.filt,paste("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/topGO/ahyacinthus.yellow.mod.genes.GO_CC.csv",sep=""))

## Salmon module
# Load Salmon module genes
salmon.module.genes <- read.csv("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/MEsalmon.genes.csv")

# Modify salmon module gene set for topGO analysis
ahyacinthus.sal.mod.IG <- factor(as.numeric(allgenes%in%salmon.module.genes$gene))
names(ahyacinthus.sal.mod.IG) <- allgenes

##BP
GOahyacinthus <- new("topGOdata",ontology="BP",allGenes=ahyacinthus.sal.mod.IG, nodeSize=10,
                     annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
ahyacinthusFisher <- getSigGroups(GOahyacinthus,test.stat)
ahyacinthus.res <- GenTable(GOahyacinthus,classic=ahyacinthusFisher,topNodes=length(ahyacinthusFisher@score),numChar=100)
ahyacinthus.filt <- ahyacinthus.res[ahyacinthus.res$classic<0.01 & ahyacinthus.res$Significant>=10,]
# write.csv(ahyacinthus.filt,paste("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/topGO/ahyacinthus.salmon.mod.genes.GO_BP.csv",sep=""))

##MF
GOahyacinthus <- new("topGOdata",ontology="MF",allGenes=ahyacinthus.sal.mod.IG, nodeSize=10,
                     annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
ahyacinthusFisher <- getSigGroups(GOahyacinthus,test.stat)
ahyacinthus.res <- GenTable(GOahyacinthus,classic=ahyacinthusFisher,topNodes=length(ahyacinthusFisher@score), numChar=100)
ahyacinthus.filt <- ahyacinthus.res[ahyacinthus.res$classic<0.01 & ahyacinthus.res$Significant>=10,]
# write.csv(ahyacinthus.filt,paste("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/topGO/ahyacinthus.salmon.mod.genes.GO_MF.csv",sep=""))

##CC
GOahyacinthus <- new("topGOdata",ontology="CC",allGenes=ahyacinthus.sal.mod.IG, nodeSize=10,
                     annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
ahyacinthusFisher <- getSigGroups(GOahyacinthus,test.stat)
ahyacinthus.res <- GenTable(GOahyacinthus,classic=ahyacinthusFisher,topNodes=length(ahyacinthusFisher@score),numChar=100)
ahyacinthus.filt <- ahyacinthus.res[ahyacinthus.res$classic<0.01 & ahyacinthus.res$Significant>=10,]
# write.csv(ahyacinthus.filt,paste("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/topGO/ahyacinthus.salmon.mod.genes.GO_CC.csv",sep=""))


