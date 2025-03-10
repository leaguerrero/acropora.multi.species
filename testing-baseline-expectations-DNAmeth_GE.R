# Species standard deviation of normalized read counts and
# percent methylation 

# Expect higher DNA methylation in genes with lower expression variation and 
# and higher methylation in genes with higher expression

library(methylKit)
library(vegan)
library(aPCoA)
library(ggplot2)
library(plotly)
library(tidyr)
library(rtracklayer)
library(topGO)
library(ggdendro)
library(DESeq2)
library(dplyr)

# Load project meta data
my.dir <- "/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/"
meta.file.name <- "../meth_meta.csv"
meta.path <- file.path(my.dir, meta.file.name)
meta <- read.csv(meta.path, header = TRUE)
meta$ID <- as.factor(meta$ID)
meta$species <- as.factor(meta$species)

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


meth=methylKit::unite(filtered.myobj, 
                      min.per.group = 2L, # site found in at least 2 individual/species
                      destrand=TRUE)
pm <- percMethylation(meth, rowids = TRUE)


# Plot PCA Results
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


# Methylation object creation

# Acropora abrotanoides
file.list.abro <- list("C033_S22.rmdup_CpG.methylKit",
                       "C075_S29.rmdup_CpG.methylKit",
                       "C076_S28.rmdup_CpG.methylKit",
                       "C077_S30.rmdup_CpG.methylKit",
                       "C078_S26.rmdup_CpG.methylKit")

sample.id.abro=list("c033",
                    "c075",
                    "c076",
                    "c077",
                    "c078")

# Methylation object creation
setwd("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
setwd("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
myobj.abro=methRead(file.list.abro,                   
                    sample.id=sample.id.abro,       
                    assembly="Ahyacinthus.chrsV1",
                    treatment=c(0,0,0,0,0),
                    context="CpG",
                    dbtype = "tabix",
                    dbdir = "methylDB")
filtered.myobj.abro=filterByCoverage(myobj.abro,
                                     lo.count=10,
                                     lo.perc=NULL,
                                     hi.count=NULL,
                                     hi.perc=99.9)

# Acropora hyacinthus

file.list.hya <- list("C011_S25.rmdup_CpG.methylKit",
                      "C017_S23.rmdup_CpG.methylKit",
                      "C018_S27.rmdup_CpG.methylKit",
                      "C038_S31.rmdup_CpG.methylKit",
                      "C051_S24.rmdup_CpG.methylKit")

sample.id.hya=list("c011",
                   "c017",
                   "c018",
                   "c038",
                   "c051")

# Methylation object creation
setwd("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
myobj.hya=methRead(file.list.hya,                   
                   sample.id=sample.id.hya,       
                   assembly="Ahyacinthus.chrsV1",
                   treatment=c(0,0,0,0,0),
                   context="CpG",
                   dbtype = "tabix",
                   dbdir = "methylDB")
filtered.myobj.hya=filterByCoverage(myobj.hya,
                                    lo.count=10,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)

# Acropora lutkeni

file.list.lut <- list("C047_S36.rmdup_CpG.methylKit",
                      "C196_S35.rmdup_CpG.methylKit",
                      "C210_S41.rmdup_CpG.methylKit")

sample.id.lut=list("c047",
                   "c196",
                   "c210")

# Methylation object creation
setwd("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
myobj.lut=methRead(file.list.lut,                   
                   sample.id=sample.id.lut,       
                   assembly="Ahyacinthus.chrsV1",
                   treatment=c(0,0,0),
                   context="CpG",
                   dbtype = "tabix",
                   dbdir = "methylDB")
filtered.myobj.lut=filterByCoverage(myobj.lut,
                                    lo.count=10,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)

# Acropora pulchra

file.list.pul <- list("C004_S34.rmdup_CpG.methylKit",
                      "C007_S39.rmdup_CpG.methylKit",
                      "C013_S40.rmdup_CpG.methylKit",
                      "C016_S38.rmdup_CpG.methylKit",
                      "C009_S59.rmdup_CpG.methylKit")

sample.id.pul=list("c004",
                   "c007",
                   "c013",
                   "c016",
                   "c009")


# Methylation object creation
setwd("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
myobj.pul=methRead(file.list.pul,                   
                   sample.id=sample.id.pul,       
                   assembly="Ahyacinthus.chrsV1",
                   treatment=c(0,0,0,0,0),
                   context="CpG",
                   dbtype = "tabix",
                   dbdir = "methylDB")
filtered.myobj.pul=filterByCoverage(myobj.pul,
                                    lo.count=10,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)

# Acropora retusa

file.list.ret <- list("C030_S37.rmdup_CpG.methylKit",
                       "C032_S33.rmdup_CpG.methylKit",
                       "C040_S32.rmdup_CpG.methylKit",
                      "C094_S42.rmdup_CpG.methylKit")

sample.id.ret=list("c030",
                    "c032",
                    "c040",
                   "c094")

# Methylation object creation
setwd("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
myobj.ret=methRead(file.list.ret,                   
                    sample.id=sample.id.ret,       
                    assembly="Ahyacinthus.chrsV1",
                    treatment=c(0,0,0,0),
                    context="CpG",
                    dbtype = "tabix",
                    dbdir = "methylDB")
filtered.myobj.ret=filterByCoverage(myobj.ret,
                                     lo.count=10,
                                     lo.perc=NULL,
                                     hi.count=NULL,
                                     hi.perc=99.9)


# Acropora robusta

file.list.rob <- list("C090_S44.rmdup_CpG.methylKit",
                      "C123_S45.rmdup_CpG.methylKit",
                      "C167_S46.rmdup_CpG.methylKit",
                      "C211_S47.rmdup_CpG.methylKit")

sample.id.rob=list("c090",
                   "c123",
                   "c167",
                   "c211")

# Methylation object creation
setwd("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
myobj.rob=methRead(file.list.rob,                   
                   sample.id=sample.id.rob,       
                   assembly="Ahyacinthus.chrsV1",
                   treatment=c(0,0,0,0),
                   context="CpG",
                   dbtype = "tabix",
                   dbdir = "methylDB")
filtered.myobj.rob=filterByCoverage(myobj.rob,
                                    lo.count=10,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)

# Acropora tutuilensis

file.list.tut <- list("C014_S48.rmdup_CpG.methylKit",
                      "C024_S49.rmdup_CpG.methylKit",
                      "C066_S50.rmdup_CpG.methylKit",
                      "C071_S51.rmdup_CpG.methylKit",
                      "C122_S52.rmdup_CpG.methylKit",
                      "C137_S53.rmdup_CpG.methylKit",
                      "C164_S54.rmdup_CpG.methylKit")

sample.id.tut=list("c014",
                   "c024",
                   "c066",
                   "c071",
                   "c122",
                   "c137",
                   "c164")

# Methylation object creation
setwd("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
myobj.tut=methRead(file.list.tut,                   
                   sample.id=sample.id.tut,       
                   assembly="Ahyacinthus.chrsV1",
                   treatment=c(0,0,0,0,0,0,0),
                   context="CpG",
                   dbtype = "tabix",
                   dbdir = "methylDB")
filtered.myobj.tut=filterByCoverage(myobj.tut,
                                    lo.count=10,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)


# Acropora nasuta

file.list.nasuta <- list("C037_S55.rmdup_CpG.methylKit",
                          "C012_S56.rmdup_CpG.methylKit",
                       "C120_S58.rmdup_CpG.methylKit")

sample.id.nasuta=list("c037",
                    "c012",
                    "c120")

# Methylation object creation
setwd("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
myobj.nasuta=methRead(file.list.nasuta,                   
                    sample.id=sample.id.nasuta,       
                    assembly="Ahyacinthus.chrsV1",
                    treatment=c(0,0,0),
                    context="CpG",
                    dbtype = "tabix",
                    dbdir = "methylDB")
filtered.myobj.nasuta=filterByCoverage(myobj.nasuta,
                                     lo.count=10,
                                     lo.perc=NULL,
                                     hi.count=NULL,
                                     hi.perc=99.9)




# Methylation: Combine gene data with methylation object
gr.obj <- import("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/Ahyacinthuns.genes.gff") # This is the GFF file.
gr.df <- data.frame(gr.obj)
names(gr.df)
head(paste(gr.df$seqnames,gr.df$start,gr.df$end,sep="."))
# Methylation: Make a column with chr.start.end
gr.df$concat <- paste(gr.df$seqnames,gr.df$start,gr.df$end,sep=".")

# Methylation: Region counts for each species
genes.abro <- regionCounts(filtered.myobj.abro,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.hya <- regionCounts(filtered.myobj.hya,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.lut <- regionCounts(filtered.myobj.lut,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.pul <- regionCounts(filtered.myobj.pul,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.ret <- regionCounts(filtered.myobj.ret,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.rob <- regionCounts(filtered.myobj.rob,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.tut <- regionCounts(filtered.myobj.tut,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.nas <- regionCounts(filtered.myobj.nasuta,gr.obj,save.db=FALSE) # Gene counts for individuals

# Methylation: Make methylBase object for each species
genes.unite.abro <- methylKit::unite(genes.abro,
                                     destrand = FALSE, #Combine
                                     min.per.group = 1L, # in at least 1 colony per species
                                     save.db = F)
genes.unite.hya <- methylKit::unite(genes.hya,
                                    destrand = FALSE, #Combine
                                    min.per.group = 1L, # in at least 1 colony per species
                                    save.db = F)
genes.unite.lut <- methylKit::unite(genes.lut,
                                    destrand = FALSE, #Combine
                                    min.per.group = 1L, # in at least 1 colony per species
                                    save.db = F)
genes.unite.pul <- methylKit::unite(genes.pul,
                                    destrand = FALSE, #Combine
                                    min.per.group = 1L, # in at least 1 colony per species
                                    save.db = F)
genes.unite.ret <- methylKit::unite(genes.ret,
                                     destrand = FALSE, #Combine
                                     min.per.group = 1L, # in at least 1 colony per species
                                     save.db = F)
genes.unite.rob <- methylKit::unite(genes.rob,
                                    destrand = FALSE, #Combine
                                    min.per.group = 1L, # in at least 1 colony per species
                                    save.db = F)
genes.unite.tut <- methylKit::unite(genes.tut,
                                    destrand = FALSE, #Combine
                                    min.per.group = 1L, # in at least 1 colony per species
                                    save.db = F)
genes.unite.nas <- methylKit::unite(genes.nas,
                                     destrand = FALSE, #Combine
                                     min.per.group = 1L, # in at least 1 colony per species
                                     save.db = F)

# Methylation: Create methylation matrix per species
gene.pm.abro <- percMethylation(genes.unite.abro, rowids=TRUE)
gene.pm.hya <- percMethylation(genes.unite.hya, rowids=TRUE)
gene.pm.lut <- percMethylation(genes.unite.lut, rowids=TRUE)
gene.pm.pul <- percMethylation(genes.unite.pul, rowids=TRUE)
gene.pm.ret <- percMethylation(genes.unite.ret, rowids=TRUE)
gene.pm.rob <- percMethylation(genes.unite.rob, rowids=TRUE)
gene.pm.tut <- percMethylation(genes.unite.tut, rowids=TRUE)
gene.pm.nas <- percMethylation(genes.unite.nas, rowids=TRUE)

# Methylation: Making data frame for each species
gene.names.abro <- gr.df$ID[match(rownames(gene.pm.abro),gr.df$concat)]
rownames(gene.pm.abro) <- gene.names.abro
gene.pm.means.abro.df <- data.frame(gene = gene.names.abro, 
                                    pm.means = rowMeans(gene.pm.abro, na.rm = TRUE), 
                                    row.names = NULL)
names(gene.pm.means.abro.df)<- c("gene","pm.mean.abro")

gene.names.hya <- gr.df$ID[match(rownames(gene.pm.hya),gr.df$concat)]
rownames(gene.pm.hya) <- gene.names.hya
gene.pm.means.hya.df <- data.frame(gene = gene.names.hya, 
                                   pm.means = rowMeans(gene.pm.hya, na.rm = TRUE), 
                                   row.names = NULL)
names(gene.pm.means.hya.df)<- c("gene","pm.mean.hya")

gene.names.lut <- gr.df$ID[match(rownames(gene.pm.lut),gr.df$concat)]
rownames(gene.pm.lut) <- gene.names.lut
gene.pm.means.lut.df <- data.frame(gene = gene.names.lut, 
                                   pm.means = rowMeans(gene.pm.lut, na.rm = TRUE), 
                                   row.names = NULL)
names(gene.pm.means.lut.df)<- c("gene","pm.mean.lut")

gene.names.pul <- gr.df$ID[match(rownames(gene.pm.pul),gr.df$concat)]
rownames(gene.pm.pul) <- gene.names.pul
gene.pm.means.pul.df <- data.frame(gene = gene.names.pul, 
                                   pm.means = rowMeans(gene.pm.pul, na.rm = TRUE), 
                                   row.names = NULL)
names(gene.pm.means.pul.df)<- c("gene","pm.mean.pul")

gene.names.ret <- gr.df$ID[match(rownames(gene.pm.ret),gr.df$concat)]
rownames(gene.pm.ret) <- gene.names.ret
gene.pm.means.ret.df <- data.frame(gene = gene.names.ret, 
                                    pm.means = rowMeans(gene.pm.ret, na.rm = TRUE), 
                                    row.names = NULL)
names(gene.pm.means.ret.df)<- c("gene","pm.mean.ret")


gene.names.rob <- gr.df$ID[match(rownames(gene.pm.rob),gr.df$concat)]
rownames(gene.pm.rob) <- gene.names.rob
gene.pm.means.rob.df <- data.frame(gene = gene.names.rob, 
                                   pm.means = rowMeans(gene.pm.rob, na.rm = TRUE), 
                                   row.names = NULL)
names(gene.pm.means.rob.df)<- c("gene","pm.mean.rob")

gene.names.tut <- gr.df$ID[match(rownames(gene.pm.tut),gr.df$concat)]
rownames(gene.pm.tut) <- gene.names.tut
gene.pm.means.tut.df <- data.frame(gene = gene.names.tut, 
                                   pm.means = rowMeans(gene.pm.tut, na.rm = TRUE), 
                                   row.names = NULL)
names(gene.pm.means.tut.df)<- c("gene","pm.mean.tut")

gene.names.nas <- gr.df$ID[match(rownames(gene.pm.nas),gr.df$concat)]
rownames(gene.pm.nas) <- gene.names.nas
gene.pm.means.nas.df <- data.frame(gene = gene.names.nas, 
                                    pm.means = rowMeans(gene.pm.nas, na.rm = TRUE), 
                                    row.names = NULL)
names(gene.pm.means.nas.df)<- c("gene","pm.mean.nas")

# tag-seq data for each species
# File Set-up
my.dir <- "/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/HTseq_out"
my.metadata <- read.csv("/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/tagseq_meta_data.csv", header = TRUE)
my.files <- grep(".txt", list.files(my.dir), value=TRUE)

# Create file list and meta data files for each species
my.metadata.abro <- my.metadata[my.metadata$species == "Acropora.abrotanoides",]
prefixes.abro <- my.metadata.abro$seq.sample.id
abro.files <- my.files[grep(paste(prefixes.abro, collapse = "|"), my.files)]

my.metadata.hya <- my.metadata[my.metadata$species == "Acropora.hyacinthus",]
prefixes.hya <- my.metadata.hya$seq.sample.id
hya.files <- my.files[grep(paste(prefixes.hya, collapse = "|"), my.files)]

my.metadata.lut <- my.metadata[my.metadata$species == "Acropora.lutkeni",]
prefixes.lut <- my.metadata.lut$seq.sample.id
lut.files <- my.files[grep(paste(prefixes.lut, collapse = "|"), my.files)]

my.metadata.pul <- my.metadata[my.metadata$species == "Acropora.pulchra",]
prefixes.pul <- my.metadata.pul$seq.sample.id
pul.files <- my.files[grep(paste(prefixes.pul, collapse = "|"), my.files)]

my.metadata.ret <- my.metadata[my.metadata$species == "Acropora.retusa",]
prefixes.ret <- my.metadata.ret$seq.sample.id
ret.files <- my.files[grep(paste(prefixes.ret, collapse = "|"), my.files)]

my.metadata.rob <- my.metadata[my.metadata$species == "Acropora.robusta",]
prefixes.rob <- my.metadata.rob$seq.sample.id
rob.files <- my.files[grep(paste(prefixes.rob, collapse = "|"), my.files)]

my.metadata.tut <- my.metadata[my.metadata$species == "Acropora.tutuilensis",]
prefixes.tut <- my.metadata.tut$seq.sample.id
tut.files <- my.files[grep(paste(prefixes.tut, collapse = "|"), my.files)]

my.metadata.nas <- my.metadata[my.metadata$species == "Acropora.nasuta",]
prefixes.nas <- my.metadata.nas$seq.sample.id
nas.files <- my.files[grep(paste(prefixes.nas, collapse = "|"), my.files)]

# Create dds results data for each species
sampleNames.abro <- my.metadata$seq.sample.id[my.metadata$species == "Acropora.abrotanoides"]
sampleNames.hya <- my.metadata$seq.sample.id[my.metadata$species == "Acropora.hyacinthus"]
sampleNames.lut <- my.metadata$seq.sample.id[my.metadata$species == "Acropora.lutkeni"]
sampleNames.pul <- my.metadata$seq.sample.id[my.metadata$species == "Acropora.pulchra"]
sampleNames.ret <- my.metadata$seq.sample.id[my.metadata$species == "Acropora.retusa"]
sampleNames.rob <- my.metadata$seq.sample.id[my.metadata$species == "Acropora.robusta"]
sampleNames.tut <- my.metadata$seq.sample.id[my.metadata$species == "Acropora.tutuilensis"]
sampleNames.nas <- my.metadata$seq.sample.id[my.metadata$species == "Acropora.nasuta"]


my.sampleTable.abro <- data.frame(sampleName = sampleNames.abro, 
                             fileName = abro.files, 
                             condition = my.metadata.abro)

my.sampleTable.hya <- data.frame(sampleName = sampleNames.hya, 
                                  fileName = hya.files, 
                                  condition = my.metadata.hya)

my.sampleTable.lut <- data.frame(sampleName = sampleNames.lut, 
                                  fileName = lut.files, 
                                  condition = my.metadata.lut)

my.sampleTable.pul <- data.frame(sampleName = sampleNames.pul, 
                                  fileName = pul.files, 
                                  condition = my.metadata.pul)

my.sampleTable.ret <- data.frame(sampleName = sampleNames.ret, 
                                  fileName = ret.files, 
                                  condition = my.metadata.ret)

my.sampleTable.rob <- data.frame(sampleName = sampleNames.rob, 
                                  fileName = rob.files, 
                                  condition = my.metadata.rob)

my.sampleTable.tut <- data.frame(sampleName = sampleNames.tut, 
                                  fileName = tut.files, 
                                  condition = my.metadata.tut)

my.sampleTable.nas <- data.frame(sampleName = sampleNames.nas, 
                                  fileName = nas.files, 
                                  condition = my.metadata.nas)

# Convert variables to factors
factorVars <- c("condition.species", "condition.colony.id", "condition.treatment")

my.sampleTable.abro[, factorVars] <- lapply(my.sampleTable.abro[, factorVars], factor)
colnames(my.sampleTable.abro)

my.sampleTable.abro[, factorVars] <- lapply(my.sampleTable.abro[, factorVars], factor)
colnames(my.sampleTable.abro)

my.sampleTable.hya[, factorVars] <- lapply(my.sampleTable.hya[, factorVars], factor)
colnames(my.sampleTable.hya)

my.sampleTable.lut[, factorVars] <- lapply(my.sampleTable.lut[, factorVars], factor)
colnames(my.sampleTable.lut)

my.sampleTable.pul[, factorVars] <- lapply(my.sampleTable.pul[, factorVars], factor)
colnames(my.sampleTable.pul)

my.sampleTable.ret[, factorVars] <- lapply(my.sampleTable.ret[, factorVars], factor)
colnames(my.sampleTable.ret)

my.sampleTable.rob[, factorVars] <- lapply(my.sampleTable.rob[, factorVars], factor)
colnames(my.sampleTable.rob)


my.sampleTable.tut[, factorVars] <- lapply(my.sampleTable.tut[, factorVars], factor)
colnames(my.sampleTable.tut)

my.sampleTable.nas[, factorVars] <- lapply(my.sampleTable.nas[, factorVars], factor)
colnames(my.sampleTable.nas)


# Create DESeqDataSet Object
ddsHTSeq.abro <- DESeqDataSetFromHTSeqCount(
  sampleTable = my.sampleTable.abro,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep.abro <- rowSums(counts(ddsHTSeq.abro)) >= 10
ddsHTSeq.abro <- ddsHTSeq.abro[keep.abro,]


# Differential Expression Analysis
dds.abro <- DESeq(ddsHTSeq.abro) 
resultsNames(dds.abro)


# Create DESeqDataSet Object
ddsHTSeq.hya <- DESeqDataSetFromHTSeqCount(
  sampleTable = my.sampleTable.hya,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep.hya <- rowSums(counts(ddsHTSeq.hya)) >= 10
ddsHTSeq.hya <- ddsHTSeq.hya[keep.hya,]

# Differential Expression Analysis
dds.hya <- DESeq(ddsHTSeq.hya) 
resultsNames(dds.hya)

# Create DESeqDataSet Object
ddsHTSeq.lut <- DESeqDataSetFromHTSeqCount(
  sampleTable = my.sampleTable.lut,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep.lut <- rowSums(counts(ddsHTSeq.lut)) >= 10
ddsHTSeq.lut <- ddsHTSeq.lut[keep.lut,]


# Differential Expression Analysis
dds.lut <- DESeq(ddsHTSeq.lut) 
resultsNames(dds.lut)

# Create DESeqDataSet Object
ddsHTSeq.pul <- DESeqDataSetFromHTSeqCount(
  sampleTable = my.sampleTable.pul,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep.pul <- rowSums(counts(ddsHTSeq.pul)) >= 10
ddsHTSeq.pul <- ddsHTSeq.pul[keep.pul,]


# Differential Expression Analysis
dds.pul <- DESeq(ddsHTSeq.pul) 
resultsNames(dds.pul)

# Create DESeqDataSet Object
ddsHTSeq.ret <- DESeqDataSetFromHTSeqCount(
  sampleTable = my.sampleTable.ret,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep.ret <- rowSums(counts(ddsHTSeq.ret)) >= 10
ddsHTSeq.ret <- ddsHTSeq.ret[keep.ret,]


# Differential Expression Analysis
dds.ret <- DESeq(ddsHTSeq.ret) 
resultsNames(dds.ret)

# Create DESeqDataSet Object
ddsHTSeq.rob <- DESeqDataSetFromHTSeqCount(
  sampleTable = my.sampleTable.rob,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep.rob <- rowSums(counts(ddsHTSeq.rob)) >= 10
ddsHTSeq.rob <- ddsHTSeq.rob[keep.rob,]


# Differential Expression Analysis
dds.rob <- DESeq(ddsHTSeq.rob) 
resultsNames(dds.rob)

# Create DESeqDataSet Object
ddsHTSeq.tut <- DESeqDataSetFromHTSeqCount(
  sampleTable = my.sampleTable.tut,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep.tut <- rowSums(counts(ddsHTSeq.tut)) >= 10
ddsHTSeq.tut <- ddsHTSeq.tut[keep.tut,]


# Differential Expression Analysis
dds.tut <- DESeq(ddsHTSeq.tut) 
resultsNames(dds.tut)

# Create DESeqDataSet Object
ddsHTSeq.nas <- DESeqDataSetFromHTSeqCount(
  sampleTable = my.sampleTable.nas,
  directory = my.dir,
  design = ~ condition.treatment)
# Pre-filter
keep.nas <- rowSums(counts(ddsHTSeq.nas)) >= 10
ddsHTSeq.nas <- ddsHTSeq.nas[keep.nas,]


# Differential Expression Analysis
dds.nas <- DESeq(ddsHTSeq.nas) 
resultsNames(dds.nas)

heat.stress.results.abro <- results(dds.abro, 
                                  contrast = c("condition.treatment", 
                                               "heat", 
                                               "control"))
heat.stress.results.abro <- as.data.frame(heat.stress.results.abro[complete.cases(heat.stress.results.abro), ])
heat.stress.results.abro$gene <- row.names(heat.stress.results.abro)
heat.stress.results.abro$log10_baseMean <- log10(heat.stress.results.abro$baseMean)
row.names(heat.stress.results.abro) <- NULL
heat.stress.results.abro <- heat.stress.results.abro[, c("gene", setdiff(names(heat.stress.results.abro), "gene"))]


heat.stress.results.hya <- results(dds.hya, 
                                    contrast = c("condition.treatment", 
                                                 "heat", 
                                                 "control"))
heat.stress.results.hya <- as.data.frame(heat.stress.results.hya[complete.cases(heat.stress.results.hya), ])
heat.stress.results.hya$gene <- row.names(heat.stress.results.hya)
heat.stress.results.hya$log10_baseMean <- log10(heat.stress.results.hya$baseMean)
row.names(heat.stress.results.hya) <- NULL
heat.stress.results.hya <- heat.stress.results.hya[, c("gene", setdiff(names(heat.stress.results.hya), "gene"))]



heat.stress.results.lut <- results(dds.lut, 
                                    contrast = c("condition.treatment", 
                                                 "heat", 
                                                 "control"))
heat.stress.results.lut <- as.data.frame(heat.stress.results.lut[complete.cases(heat.stress.results.lut), ])
heat.stress.results.lut$gene <- row.names(heat.stress.results.lut)
heat.stress.results.lut$log10_baseMean <- log10(heat.stress.results.lut$baseMean)
row.names(heat.stress.results.lut) <- NULL
heat.stress.results.lut <- heat.stress.results.lut[, c("gene", setdiff(names(heat.stress.results.lut), "gene"))]


heat.stress.results.pul <- results(dds.pul, 
                                    contrast = c("condition.treatment", 
                                                 "heat", 
                                                 "control"))
heat.stress.results.pul <- as.data.frame(heat.stress.results.pul[complete.cases(heat.stress.results.pul), ])
heat.stress.results.pul$gene <- row.names(heat.stress.results.pul)
heat.stress.results.pul$log10_baseMean <- log10(heat.stress.results.pul$baseMean)
row.names(heat.stress.results.pul) <- NULL
heat.stress.results.pul <- heat.stress.results.pul[, c("gene", setdiff(names(heat.stress.results.pul), "gene"))]


heat.stress.results.ret <- results(dds.ret, 
                                    contrast = c("condition.treatment", 
                                                 "heat", 
                                                 "control"))
heat.stress.results.ret <- as.data.frame(heat.stress.results.ret[complete.cases(heat.stress.results.ret), ])
heat.stress.results.ret$gene <- row.names(heat.stress.results.ret)
heat.stress.results.ret$log10_baseMean <- log10(heat.stress.results.ret$baseMean)
row.names(heat.stress.results.ret) <- NULL
heat.stress.results.ret <- heat.stress.results.ret[, c("gene", setdiff(names(heat.stress.results.ret), "gene"))]


heat.stress.results.rob <- results(dds.rob, 
                                    contrast = c("condition.treatment", 
                                                 "heat", 
                                                 "control"))
heat.stress.results.rob <- as.data.frame(heat.stress.results.rob[complete.cases(heat.stress.results.rob), ])
heat.stress.results.rob$gene <- row.names(heat.stress.results.rob)
heat.stress.results.rob$log10_baseMean <- log10(heat.stress.results.rob$baseMean)
row.names(heat.stress.results.rob) <- NULL
heat.stress.results.rob <- heat.stress.results.rob[, c("gene", setdiff(names(heat.stress.results.rob), "gene"))]


heat.stress.results.tut <- results(dds.tut, 
                                    contrast = c("condition.treatment", 
                                                 "heat", 
                                                 "control"))
heat.stress.results.tut <- as.data.frame(heat.stress.results.tut[complete.cases(heat.stress.results.tut), ])
heat.stress.results.tut$gene <- row.names(heat.stress.results.tut)
heat.stress.results.tut$log10_baseMean <- log10(heat.stress.results.tut$baseMean)
row.names(heat.stress.results.tut) <- NULL
heat.stress.results.tut <- heat.stress.results.tut[, c("gene", setdiff(names(heat.stress.results.tut), "gene"))]

heat.stress.results.nas <- results(dds.nas, 
                                    contrast = c("condition.treatment", 
                                                 "heat", 
                                                 "control"))
heat.stress.results.nas <- as.data.frame(heat.stress.results.nas[complete.cases(heat.stress.results.nas), ])
heat.stress.results.nas$gene <- row.names(heat.stress.results.nas)
heat.stress.results.nas$log10_baseMean <- log10(heat.stress.results.nas$baseMean)
row.names(heat.stress.results.nas) <- NULL
heat.stress.results.nas <- heat.stress.results.nas[, c("gene", setdiff(names(heat.stress.results.nas), "gene"))]



# Report: Number of transcripts
nrow(heat.stress.results.abro)#5452 
nrow(heat.stress.results.hya) #9175
nrow(heat.stress.results.lut) #11092
nrow(heat.stress.results.pul) #2487
nrow(heat.stress.results.ret) #11112
nrow(heat.stress.results.rob) #10265
nrow(heat.stress.results.tut) #12849
nrow(heat.stress.results.nas) #5455

# Merging baseMeans and percent methylation
exp.and.perc.meth.abro <- inner_join(heat.stress.results.abro,gene.pm.means.abro.df, by = 'gene')
dim(exp.and.perc.meth.abro) 

# Gene expression and methylation: linear regression of the relationship
mod.ge.abro<-lm(log10_baseMean ~ pm.mean.abro, data = exp.and.perc.meth.abro)
summary(mod.ge.abro) # significant relationship between the percent methylation mean and the log10 base means

plot_l10bMs_methylation_abro <- ggplot(exp.and.perc.meth.abro, aes(x=pm.mean.abro, y=log10_baseMean)) + 
  geom_point(alpha = 0.3, color = "#3E6FCB") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle("Acropora abrotanoides") +
  xlab("Average % methylation") + 
  ylab("log10 transformed expression means") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)  
plot_l10bMs_methylation_abro 


plot_hs_l2fcSD_methylation_abro <- ggplot(exp.and.perc.meth.abro, aes(x=pm.mean.abro, y=lfcSE)) + 
  geom_point(alpha = 0.3, color = "#3E6FCB") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
 # ggtitle("Relationship between heat stress Log2FC SD and percent methylation\nA. abrotanoides") +
  xlab("Average % methylation") + 
  ylab("heat stress Log2FC SD")  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6 )
plot_hs_l2fcSD_methylation_abro

# Acropora hyacinthus
exp.and.perc.meth.hya <- inner_join(heat.stress.results.hya,gene.pm.means.hya.df, by = 'gene')
dim(exp.and.perc.meth.hya) 

# Gene expression and methylation: linear regression of the relationship
mod.ge.hya<-lm(log10_baseMean ~ pm.mean.hya, data = exp.and.perc.meth.hya)
summary(mod.ge.hya) # significant relationship between the percent methylation mean and the log10 base means

plot_l10bMs_methylation_hya <- ggplot(exp.and.perc.meth.hya, aes(x=pm.mean.hya, y=log10_baseMean)) + 
  geom_point(alpha = 0.3, color = "#3a86ff") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle("Acropora hyacinthus") +
  xlab("Average % methylation") + 
  ylab("log10 transformed expression means") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)  
plot_l10bMs_methylation_hya 


plot_hs_l2fcSD_methylation_hya <- ggplot(exp.and.perc.meth.hya, aes(x=pm.mean.hya, y=lfcSE)) + 
  geom_point(alpha = 0.3, color = "#3a86ff") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  #ggtitle("Relationship between heat stress Log2FC SD and percent methylation\n(A. hyacinthus)") +
  xlab("Average % methylation") + 
  ylab("heat stress Log2FC SD")  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6 )
plot_hs_l2fcSD_methylation_hya


# Acropora lutkeni
exp.and.perc.meth.lut <- inner_join(heat.stress.results.lut,gene.pm.means.lut.df, by = 'gene')
dim(exp.and.perc.meth.lut) 

# Gene expression and methylation: linear regression of the relationship
mod.ge.lut<-lm(log10_baseMean ~ pm.mean.lut, data = exp.and.perc.meth.lut)
summary(mod.ge.lut) # significant relationship between the percent methylation mean and the log10 base means

plot_l10bMs_methylation_lut <- ggplot(exp.and.perc.meth.lut, aes(x=pm.mean.lut, y=log10_baseMean)) + 
  geom_point(alpha = 0.3, color = "#FF9A17") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle("Acropora lutkeni") +
  xlab("Average % methylation") + 
  ylab("log10 transformed expression means") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)  
plot_l10bMs_methylation_lut 


plot_hs_l2fcSD_methylation_lut <- ggplot(exp.and.perc.meth.lut, aes(x=pm.mean.lut, y=lfcSE)) + 
  geom_point(alpha = 0.3, color = "#FF9A17") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
#  ggtitle("Relationship between heat stress Log2FC SD and percent methylation\n(A. lutkeni)") +
  xlab("Average % methylation") + 
  ylab("heat stress Log2FC SD")  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6 )
plot_hs_l2fcSD_methylation_lut


# Acropora robusta
exp.and.perc.meth.rob <- inner_join(heat.stress.results.rob,gene.pm.means.rob.df, by = 'gene')
dim(exp.and.perc.meth.rob) 

# Gene expression and methylation: linear regression of the relationship
mod.ge.rob<-lm(log10_baseMean ~ pm.mean.rob, data = exp.and.perc.meth.rob)
summary(mod.ge.rob) # significant relationship between the percent methylation mean and the log10 base means

plot_l10bMs_methylation_rob <- ggplot(exp.and.perc.meth.rob, aes(x=pm.mean.rob, y=log10_baseMean)) + 
  geom_point(alpha = 0.3, color = "#8338ec") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle("A. robusta") +
  xlab("Average % methylation") + 
  ylab("log10 transformed expression means") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)  
plot_l10bMs_methylation_rob 


plot_hs_l2fcSD_methylation_rob <- ggplot(exp.and.perc.meth.rob, aes(x=pm.mean.rob, y=lfcSE)) + 
  geom_point(alpha = 0.3, color = "#8338ec") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  #ggtitle("Relationship between heat stress Log2FC SD and percent methylation\n(A. robusta)") +
  xlab("Average % methylation") + 
  ylab("heat stress Log2FC SD")  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6 )
plot_hs_l2fcSD_methylation_rob

# Acropora pulchra
exp.and.perc.meth.pul <- inner_join(heat.stress.results.pul,gene.pm.means.pul.df, by = 'gene')
dim(exp.and.perc.meth.pul) 

# Gene expression and methylation: linear regression of the relationship
mod.ge.pul<-lm(log10_baseMean ~ pm.mean.pul, data = exp.and.perc.meth.pul)
summary(mod.ge.pul) # significant relationship between the percent methylation mean and the log10 base means

plot_l10bMs_methylation_pul <- ggplot(exp.and.perc.meth.pul, aes(x=pm.mean.pul, y=log10_baseMean)) + 
  geom_point(alpha = 0.3, color = "#fb5607") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle("Acropora pulchra") +
  xlab("Average % methylation") + 
  ylab("log10 transformed expression means") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)  
plot_l10bMs_methylation_pul 


plot_hs_l2fcSD_methylation_pul <- ggplot(exp.and.perc.meth.pul, aes(x=pm.mean.pul, y=lfcSE)) + 
  geom_point(alpha = 0.3, color = "#fb5607") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
#  ggtitle("Relationship between heat stress Log2FC SD and percent methylation\n(A. pulchra)") +
  xlab("Average % methylation") + 
  ylab("heat stress Log2FC SD")  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6 )
plot_hs_l2fcSD_methylation_pul

# Acropora retusa
exp.and.perc.meth.ret <- inner_join(heat.stress.results.ret,gene.pm.means.ret.df, by = 'gene')
dim(exp.and.perc.meth.ret) 

# Gene expression and methylation: linear regression of the relationship

mod.ge.ret<-lm(log10_baseMean ~ pm.mean.ret, data = exp.and.perc.meth.ret)
summary(mod.ge.ret) # significant relationship between the percent methylation mean and the log10 base means

plot_l10bMs_methylation_ret <- ggplot(exp.and.perc.meth.ret, aes(x=pm.mean.ret, y=log10_baseMean)) + 
  geom_point(alpha = 0.3, color = "#ff006e") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle("A. retusa") +
  xlab("Average % methylation") + 
  ylab("log10 transformed expression means") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)  
plot_l10bMs_methylation_ret 


plot_hs_l2fcSD_methylation_ret <- ggplot(exp.and.perc.meth.ret, aes(x=pm.mean.ret, y=lfcSE)) + 
  geom_point(alpha = 0.3, color = "#ff006e") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  #ggtitle("Relationship between heat stress Log2FC SD and percent methylation\n(A. retusa)") +
  xlab("Average % methylation") + 
  ylab("heat stress Log2FC SD")  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6 )
plot_hs_l2fcSD_methylation_ret

# Acropora nasuta
exp.and.perc.meth.nas <- inner_join(heat.stress.results.nas,gene.pm.means.nas.df, by = 'gene')
dim(exp.and.perc.meth.nas) 

# Gene expression and methylation: linear regression of the relationship


mod.ge.nas<-lm(log10_baseMean ~ pm.mean.nas, data = exp.and.perc.meth.nas)
summary(mod.ge.nas) # significant relationship between the percent methylation mean and the log10 base means

plot_l10bMs_methylation_nas <- ggplot(exp.and.perc.meth.nas, aes(x=pm.mean.nas, y=log10_baseMean)) + 
  geom_point(alpha = 0.3, color = "#c11cad") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle("Acropora nasuta") +
  xlab("Average % methylation") + 
  ylab("log10 transformed expression means") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)  
plot_l10bMs_methylation_nas 


plot_hs_l2fcSD_methylation_nas <- ggplot(exp.and.perc.meth.nas, aes(x=pm.mean.nas, y=lfcSE)) + 
  geom_point(alpha = 0.3, color = "#c11cad") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
 # ggtitle("Relationship between heat stress Log2FC SD and percent methylation\n(A. nasuta)") +
  xlab("Average % methylation") + 
  ylab("heat stress Log2FC SD")  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6 )
plot_hs_l2fcSD_methylation_nas

# Acropora tutuilensis
exp.and.perc.meth.tut <- inner_join(heat.stress.results.tut,gene.pm.means.tut.df, by = 'gene')
dim(exp.and.perc.meth.tut) 

# Gene expression and methylation: linear regression of the relationship

mod.ge.tut<-lm(log10_baseMean ~ pm.mean.tut, data = exp.and.perc.meth.tut)
summary(mod.ge.tut) # significant relationship between the percent methylation mean and the log10 base means

plot_l10bMs_methylation_tut <- ggplot(exp.and.perc.meth.tut, aes(x=pm.mean.tut, y=log10_baseMean)) + 
  geom_point(alpha = 0.3, color = "#F7C11E") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle("Acropora tutuilensis") +
  xlab("Average % methylation") + 
  ylab("log10 transformed expression means") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)  
plot_l10bMs_methylation_tut 


plot_hs_l2fcSD_methylation_tut <- ggplot(exp.and.perc.meth.tut, aes(x=pm.mean.tut, y=lfcSE)) + 
  geom_point(alpha = 0.3, color = "#F7C11E") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  #ggtitle("Relationship between heat stress Log2FC SD and percent methylation\n(A. tutuilensis)") +
  xlab("Average % methylation") + 
  ylab("heat stress Log2FC SD")  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6 )
plot_hs_l2fcSD_methylation_tut

# All species together?
genes.all.species <- regionCounts(filtered.myobj,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.unite.all.species <- methylKit::unite(genes.all.species,
                                     destrand = FALSE, #Combine
                                     min.per.group = 1L, # in at least 1 colony per species
                                     save.db = F)

gene.pm.all.species <- percMethylation(genes.unite.all.species, rowids=TRUE)
gene.names<- gr.df$ID[match(rownames(gene.pm.all.species),gr.df$concat)]
rownames(gene.pm.all.species) <- gene.names
gene.pm.means.all.species.df <- data.frame(gene = gene.names, 
                                    pm.means = rowMeans(gene.pm.all.species, na.rm = TRUE), 
                                    row.names = NULL)

# Gene expression all species
# File Set-up
my.dir <- "/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/HTseq_out"
my.metadata <- read.csv("/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/tagseq_meta_data.csv", header = TRUE)
my.files <- grep(".txt", list.files(my.dir), value=TRUE)

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
heat.stress.results.all.species <- results(dds, 
                                    contrast = c("condition.treatment", 
                                                 "heat", 
                                                 "control"))
heat.stress.results.all.species <- as.data.frame(heat.stress.results.all.species[complete.cases(heat.stress.results.all.species), ])
heat.stress.results.all.species$gene <- row.names(heat.stress.results.all.species)
heat.stress.results.all.species$log10_baseMean <- log10(heat.stress.results.all.species$baseMean)
row.names(heat.stress.results.all.species) <- NULL
heat.stress.results.all.species <- heat.stress.results.all.species[, c("gene", setdiff(names(heat.stress.results.all.species), "gene"))]

# Merging baseMeans and percent methylation
exp.and.perc.meth.all.species <- inner_join(heat.stress.results.all.species,gene.pm.means.all.species.df, by = 'gene')
dim(exp.and.perc.meth.all.species) 

# Gene expression and methylation: linear regression of the relationship
mod.ge.all.species<-lm(log10_baseMean ~ pm.means, data = exp.and.perc.meth.all.species)
summary(mod.ge.all.species) # significant relationship between the percent methylation mean and the log10 base means

plot_l10bMs_methylation_all.species <- ggplot(exp.and.perc.meth.all.species, aes(x=pm.means, y=log10_baseMean)) + 
  geom_point(alpha = 0.3, color = "#76877d") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle("All species") +
  xlab("Average % methylation") + 
  ylab("log10 transformed expression means") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)  
  plot_l10bMs_methylation_all.species 


plot_hs_l2fcSD_methylation_all.species <- ggplot(exp.and.perc.meth.all.species, aes(x=pm.means, y=lfcSE)) + 
  geom_point(alpha = 0.3, color = "#76877d") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
 # ggtitle("Relationship between heat stress\nLog2FC SD and percent methylation\n(All species)") +
  xlab("Average % methylation") + 
  ylab("heat stress Log2FC SD")  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)
plot_hs_l2fcSD_methylation_all.species

library(patchwork)
plot_l10bMs_methylation_all.species + plot_hs_l2fcSD_methylation_all.species

(plot_l10bMs_methylation_tut + plot_hs_l2fcSD_methylation_tut)/
(plot_l10bMs_methylation_lut + plot_hs_l2fcSD_methylation_lut)  
(plot_l10bMs_methylation_pul + plot_hs_l2fcSD_methylation_pul)/
  (plot_l10bMs_methylation_ret + plot_hs_l2fcSD_methylation_ret) 
(plot_l10bMs_methylation_nas + plot_hs_l2fcSD_methylation_nas)/
  (plot_l10bMs_methylation_rob + plot_hs_l2fcSD_methylation_rob) 
(plot_l10bMs_methylation_hya + plot_hs_l2fcSD_methylation_hya)/
  (plot_l10bMs_methylation_abro + plot_hs_l2fcSD_methylation_abro)


### Checking the % methylation of genes in yellow and salmon modules
# Acropora tutuilensis
exp.and.perc.meth.tut <- inner_join(heat.stress.results.tut,gene.pm.means.tut.df, by = 'gene')
dim(exp.and.perc.meth.tut) 

# Gene expression and methylation: linear regression of the relationship

mod.ge.tut<-lm(log10_baseMean ~ pm.mean.tut, data = exp.and.perc.meth.tut)
summary(mod.ge.tut) # significant relationship between the percent methylation mean and the log10 base means

plot_l10bMs_methylation_tut <- ggplot(exp.and.perc.meth.tut, aes(x=pm.mean.tut, y=log10_baseMean)) + 
  geom_point(alpha = 0.3, color = "#F7C11E") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle("Acropora tutuilensis") +
  xlab("Average % methylation") + 
  ylab("log10 transformed expression means") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)  
plot_l10bMs_methylation_tut

# Now make a data frame of the genes in the yellow module
dim(yellow.tut)
exp.and.perc.meth.tut.salmon <- exp.and.perc.meth.tut %>%
  filter(gene %in% MEsalmon.genes)
exp.and.perc.meth.tut.yellow <- exp.and.perc.meth.tut %>%
  filter(gene %in% MEyellow.genes)

plot_l10bMs_methylation_tut_salmon <- ggplot(exp.and.perc.meth.tut.salmon , aes(x=pm.mean.tut, y=log10_baseMean)) + 
  geom_point(alpha = 0.3, color = "#F7C11E") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle("Acropora tutuilensis (Salmon module genes)") +
  xlab("Average % methylation") + 
  ylab("log10 transformed expression means") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)  
plot_l10bMs_methylation_tut_salmon

plot_l10bMs_methylation_tut_yellow <- ggplot(exp.and.perc.meth.tut.yellow , aes(x=pm.mean.tut, y=log10_baseMean)) + 
  geom_point(alpha = 0.3, color = "#F7C11E") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle("Acropora tutuilensis (yellow module genes)") +
  xlab("Average % methylation") + 
  ylab("log10 transformed expression means") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)  
plot_l10bMs_methylation_tut_yellow


# Acropora abrotanoides
salmon.genes <- read.csv("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/MEsalmon.genes.csv")
yellow.genes <- read.csv("~/Lab Notebook/Chapter3/Tag_seq_analysis/Data_out/MEyellow.genes.csv")
exp.and.perc.meth.abro.salmon <- exp.and.perc.meth.abro %>%
  filter(gene %in% salmon.genes$x)
exp.and.perc.meth.abro.yellow <- exp.and.perc.meth.abro %>%
  filter(gene %in% yellow.genes$x)



plot_l10bMs_methylation_abro_salmon <- ggplot(exp.and.perc.meth.abro.salmon, aes(x=pm.mean.abro, y=log10_baseMean)) + 
  geom_point(alpha = 0.3, color = "#3E6FCB") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle("Acropora abrotanoides") +
  xlab("Average % methylation") + 
  ylab("log10 transformed expression means") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)  
plot_l10bMs_methylation_abro_salmon 

plot_l10bMs_methylation_abro_yellow <- ggplot(exp.and.perc.meth.abro.yellow, aes(x=pm.mean.abro, y=log10_baseMean)) + 
  geom_point(alpha = 0.3, color = "#3E6FCB") + geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_pubr() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle("Acropora abrotanoides") +
  xlab("Average % methylation") + 
  ylab("log10 transformed expression means") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)  
plot_l10bMs_methylation_abro_yellow

## Percent methylation in salmon and yellow gene modules
exp.and.perc.meth.abro$species <- "Acropora abrotanoides"
exp.and.perc.meth.abro$pm.mean <- exp.and.perc.meth.abro$pm.mean.abro
exp.and.perc.meth.abro <- exp.and.perc.meth.abro[,c(-9)]


exp.and.perc.meth.hya$species <- "Acropora hyacinthus"
exp.and.perc.meth.hya$pm.mean <- exp.and.perc.meth.hya$pm.mean.hya
exp.and.perc.meth.hya <- exp.and.perc.meth.hya[,c(-9)]

exp.and.perc.meth.lut$species <- "Acropora lutkeni"
exp.and.perc.meth.lut$pm.mean <- exp.and.perc.meth.lut$pm.mean.lut
exp.and.perc.meth.lut <- exp.and.perc.meth.lut[,c(-9)]

exp.and.perc.meth.nas$species <- "Acropora nasuta"
exp.and.perc.meth.nas$pm.mean <- exp.and.perc.meth.nas$pm.mean.nas
exp.and.perc.meth.nas <- exp.and.perc.meth.nas[,c(-9)]

exp.and.perc.meth.pul$species <- "Acropora pulchra"
exp.and.perc.meth.pul$pm.mean <- exp.and.perc.meth.pul$pm.mean.pul
exp.and.perc.meth.pul <- exp.and.perc.meth.pul[,c(-9)]

exp.and.perc.meth.ret$species <- "Acropora retusa"
exp.and.perc.meth.ret$pm.mean <- exp.and.perc.meth.ret$pm.mean.ret
exp.and.perc.meth.ret <- exp.and.perc.meth.ret[,c(-9)]

exp.and.perc.meth.rob$species <- "Acropora robusta"
exp.and.perc.meth.rob$pm.mean <- exp.and.perc.meth.rob$pm.mean.rob
exp.and.perc.meth.rob <- exp.and.perc.meth.rob[,c(-9)]

exp.and.perc.meth.tut$species <- "Acropora tutuilensis"
exp.and.perc.meth.tut$pm.mean <- exp.and.perc.meth.tut$pm.mean.tut
exp.and.perc.meth.tut <- exp.and.perc.meth.tut[,c(-9)]

all.species.pm.mean <- rbind(exp.and.perc.meth.tut,
                                        exp.and.perc.meth.rob,
                                        exp.and.perc.meth.ret,
                                        exp.and.perc.meth.pul,
                                        exp.and.perc.meth.nas,
                                        exp.and.perc.meth.lut,
                                        exp.and.perc.meth.hya,
                                        exp.and.perc.meth.abro)


all.species.pm.mean.yellow <- all.species.pm.mean[all.species.pm.mean$gene %in% yellow.genes$x, ]

all.species.pm.mean.yellow$species <- factor(all.species.pm.mean.yellow$species, levels = species.order)
yellow.module.log.pm <- ggplot(all.species.pm.mean.yellow, aes(x = species, y = log(pm.mean), fill = species )) +
  # # clouds
  introdataviz::geom_flat_violin(trim=TRUE,
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
  ggtitle("Percent methylation\nof genes in yellow module")
yellow.module.log.pm

yellow.module.pm <- ggplot(all.species.pm.mean.yellow, aes(x = species, y = pm.mean, fill = species)) +
  # # clouds
  introdataviz::geom_flat_violin(trim=TRUE,
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
  ggtitle("Percent methylation\nof genes in yellow module")
yellow.module.pm


# Salmon 

all.species.pm.mean.salmon <- all.species.pm.mean[all.species.pm.mean$gene %in% salmon.genes$x, ]

all.species.pm.mean.salmon$species <- factor(all.species.pm.mean.salmon$species, levels = species.order)
salmon.module.log.pm <- ggplot(all.species.pm.mean.salmon, aes(x = species, y = log(pm.mean), fill = species )) +
  # # clouds
  introdataviz::geom_flat_violin(trim=TRUE,
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
  ggtitle("Percent methylation\nof genes in salmon module")
salmon.module.log.pm

salmon.module.pm <- ggplot(all.species.pm.mean.salmon, aes(x = species, y = pm.mean, fill = species)) +
  # # clouds
  introdataviz::geom_flat_violin(trim=TRUE,
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
  ggtitle("Percent methylation\nof genes in salmon module")
salmon.module.pm

library(patchwork)

yellow.module.pm + salmon.module.pm + yellow.module.log.pm + salmon.module.log.pm


## Coorelation between the pm.mean and the basemean of the yellow and salmon genes?
model <- lm(log10_baseMean ~ pm.mean, data = all.species.pm.mean.yellow)

# Get the summary of the model
model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- signif(model_summary$coefficients[2, 4], digits=3)
r_squared <- model_summary$r.squared

# Print the results
cat("P-value:", p_value, "\n")
cat("R-squared:", r_squared, "\n")

# Create the plot

basemean.pm.meth.yellow <- ggplot(all.species.pm.mean.yellow, aes(x = pm.mean, y = log10_baseMean)) +
  geom_point(aes(color = species), size = 3, alpha = 0.7) +  # Adjust the size of the points and set color within aes
  labs(x = "pm.mean", y = "log(basemean)") + 
  scale_color_manual(values = official.species.colors) + # Label axes and legend
  geom_smooth(method = "lm", se = TRUE, color = "black") + 
  ggtitle("Percent methylation and overall expression\nof genes in Yellow Module") +
  annotate("text", x = max(all.species.pm.mean.yellow$pm.mean), 
           y = max(all.species.pm.mean.yellow$log10_baseMean), 
           label = paste("p-value =", p_value, "R-squared:", round(r_squared, 4)), 
           hjust = 1, vjust = 1, size = 3)  +
  theme_classic()   # Use a minimal theme

basemean.pm.meth.yellow


model <- lm(abs(log2FoldChange) ~ pm.mean, data = all.species.pm.mean.yellow)

# Get the summary of the model
model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- signif(model_summary$coefficients[2, 4], digits=3)
r_squared <- model_summary$r.squared

# Print the results
cat("P-value:", p_value, "\n")
cat("R-squared:", r_squared, "\n")

# Create the plot

l2fc.pm.meth.yellow <- ggplot(all.species.pm.mean.yellow, aes(x = pm.mean, y = abs(log2FoldChange))) +
  geom_point(aes(color = species), size = 3, alpha = 0.7) +  # Adjust the size of the points and set color within aes
  labs(x = "pm.mean", y = "abs(log2FoldChange)") + 
  scale_color_manual(values = official.species.colors) + # Label axes and legend
  geom_smooth(method = "lm", se = TRUE, color = "black") + 
  ggtitle("Percent methylation and abs l2FC\nof genes in Yellow Module") +
  annotate("text", x = max(all.species.pm.mean.yellow$pm.mean), 
           y = max(all.species.pm.mean.yellow$log2FoldChange), 
           label = paste("p-value =", p_value, "R-squared:", round(r_squared, 4)), 
           hjust = 1, vjust = 1, size = 3)  +
  theme_classic()   # Use a minimal theme
l2fc.pm.meth.yellow


## Coorelation between the pm.mean and the basemean of the salmon and salmon genes?
model <- lm(log10_baseMean ~ pm.mean, data = all.species.pm.mean.salmon)

# Get the summary of the model
model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- signif(model_summary$coefficients[2, 4], digits=3)
r_squared <- model_summary$r.squared

# Print the results
cat("P-value:", p_value, "\n")
cat("R-squared:", r_squared, "\n")

# Create the plot

basemean.pm.meth.salmon <- ggplot(all.species.pm.mean.salmon, aes(x = pm.mean, y = log10_baseMean)) +
  geom_point(aes(color = species), size = 3, alpha = 0.7) +  # Adjust the size of the points and set color within aes
  labs(x = "pm.mean", y = "log(basemean)") + 
  scale_color_manual(values = official.species.colors) + # Label axes and legend
  geom_smooth(method = "lm", se = TRUE, color = "black") + 
  ggtitle("Percent methylation and overall expression\nof genes in salmon Module") +
  annotate("text", x = max(all.species.pm.mean.salmon$pm.mean), 
           y = max(all.species.pm.mean.salmon$log10_baseMean), 
           label = paste("p-value =", p_value, "R-squared:", round(r_squared, 4)), 
           hjust = 1, vjust = 1, size = 3)  +
  theme_classic()   # Use a minimal theme

basemean.pm.meth.salmon


model <- lm(abs(log2FoldChange) ~ pm.mean, data = all.species.pm.mean.salmon)

# Get the summary of the model
model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- signif(model_summary$coefficients[2, 4], digits=3)
r_squared <- model_summary$r.squared

# Print the results
cat("P-value:", p_value, "\n")
cat("R-squared:", r_squared, "\n")

# Create the plot

l2fc.pm.meth.salmon <- ggplot(all.species.pm.mean.salmon, aes(x = pm.mean, y = abs(log2FoldChange))) +
  geom_point(aes(color = species), size = 3, alpha = 0.7) +  # Adjust the size of the points and set color within aes
  labs(x = "pm.mean", y = "abs(log2FoldChange)") + 
  scale_color_manual(values = official.species.colors) + # Label axes and legend
  geom_smooth(method = "lm", se = TRUE, color = "black") + 
  ggtitle("Percent methylation and abs l2FC\nof genes in salmon Module") +
  annotate("text", x = max(all.species.pm.mean.salmon$pm.mean), 
           y = max(all.species.pm.mean.salmon$log2FoldChange), 
           label = paste("p-value =", p_value, "R-squared:", round(r_squared, 4)), 
           hjust = 1, vjust = 1, size = 3)  +
  theme_classic()   # Use a minimal theme
l2fc.pm.meth.salmon


basemean.pm.meth.yellow + l2fc.pm.meth.yellow +basemean.pm.meth.salmon+ l2fc.pm.meth.salmon

## All genes in all species
model <- lm(log10_baseMean ~ pm.mean, data = all.species.pm.mean)

# Get the summary of the model
model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- signif(model_summary$coefficients[2, 4], digits=3)
r_squared <- model_summary$r.squared

# Print the results
cat("P-value:", p_value, "\n")
cat("R-squared:", r_squared, "\n")

# Create the plot
all.species.pm.mean$species <- factor(all.species.pm.mean$species, levels = species.order)

basemean.pm.meth<- ggplot(all.species.pm.mean, aes(x = pm.mean, y = log10_baseMean)) +
  geom_point(aes(color = species), size = 3, alpha = 0.4, show.legend = FALSE) +  # Adjust the size of the points and set color within aes
  labs(x = "pm.mean", y = "log(basemean)") + 
  scale_color_manual(values = official.species.colors) + # Label axes and legend
  geom_smooth(method = "lm", se = TRUE, color = "black") + 
  ggtitle("Percent methylation and overall expression\nof genes in all genes") +
  annotate("text", x = max(all.species.pm.mean$pm.mean), 
           y = max(all.species.pm.mean$log10_baseMean), 
           label = paste("p-value =", p_value, "R-squared:", round(r_squared, 4)), 
           hjust = 1, vjust = 1, size = 3)  +
  theme_classic()   # Use a minimal theme

basemean.pm.meth



model <- lm(lfcSE ~ pm.mean, data = all.species.pm.mean)

# Get the summary of the model
model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- signif(model_summary$coefficients[2, 4], digits=3)
r_squared <- model_summary$r.squared

# Print the results
cat("P-value:", p_value, "\n")
cat("R-squared:", r_squared, "\n")


hs_l2fcSD_methylation_all.species<- ggplot(all.species.pm.mean, aes(x = pm.mean, y = lfcSE)) +
  geom_point(aes(color = species), size = 3, alpha = 0.4) +  # Adjust the size of the points and set color within aes
  labs(x = "pm.mean", y = "heat stress Log2FC SD") + 
  scale_color_manual(values = official.species.colors) + # Label axes and legend
  geom_smooth(method = "lm", se = TRUE, color = "black") + 
  ggtitle("Relationship between heat stress\nLog2FC SD and percent methylation\n(All species)") +
  annotate("text", x = max(all.species.pm.mean$pm.mean), 
           y = max(all.species.pm.mean$log10_baseMean), 
           label = paste("p-value =", p_value, "R-squared:", round(r_squared, 4)), 
           hjust = 1, vjust = 1, size = 3)  +
  theme_classic()   # Use a minimal theme

hs_l2fcSD_methylation_all.species

basemean.pm.meth + hs_l2fcSD_methylation_all.species



model <- lm( abs(log2FoldChange) ~ pm.mean, data = all.species.pm.mean)

# Get the summary of the model
model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- signif(model_summary$coefficients[2, 4], digits=3)
r_squared <- model_summary$r.squared

# Print the results
cat("P-value:", p_value, "\n")
cat("R-squared:", r_squared, "\n")


hs_l2fc_methylation_all.species<- ggplot(all.species.pm.mean, aes(x = pm.mean, y = abs(log2FoldChange))) +
  geom_point(aes(color = species), size = 3, alpha = 0.4) +  # Adjust the size of the points and set color within aes
  labs(x = "pm.mean", y = "abs heat stress Log2FC") + 
  scale_color_manual(values = official.species.colors) + # Label axes and legend
  geom_smooth(method = "lm", se = TRUE, color = "black") + 
  ggtitle("Relationship between heat stress\nLog2FC and percent methylation\n(All species)") +
  annotate("text", x = max(all.species.pm.mean$pm.mean), 
           y = max(abs(all.species.pm.mean$log2FoldChange)), 
           label = paste("p-value =", p_value, "R-squared:", round(r_squared, 4)), 
           hjust = 1, vjust = 1, size = 3)  +
  theme_classic()   # Use a minimal theme

hs_l2fc_methylation_all.species



# plot_hs_l2fcSD_methylation_all.species <- ggplot(exp.and.perc.meth.all.species, aes(x=pm.means, y=lfcSE)) + 
#   geom_point(alpha = 0.3, color = "#76877d") + geom_smooth(method='lm', formula= y~x, color = "black") +
#   theme_pubr() +
#   # ggtitle("Relationship between heat stress\nLog2FC SD and percent methylation\n(All species)") +
#   xlab("Average % methylation") + 
#   ylab("heat stress Log2FC SD")  +
#   stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 6)
# plot_hs_l2fcSD_methylation_all.species
