# Multispecies Methylation
# Libraries
library(methylKit)
library(vegan)
library(aPCoA)
library(ggplot2)
library(plotly)
library(tidyr)
library(rtracklayer)
library(topGO)
library(ggdendro)

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


sample.ids <- sapply(filtered.myobj, function(x) x@sample.id)
num.records.list <- sapply(filtered.myobj, function(x) x@num.records)
num.records.df <- data.frame(colony.id = sample.ids, num.records = num.records.list)

num.records.df$species <- meta[match(num.records.df$colony.id, meta$ID),3]
species_color_palette <- c("#ffbe0b", # Acropora tutuilensis
                           "#fd8a09", # Acropora lutkeni
                           "#fb5607", # Acropora retusa
                           "#ff006e", # Acropora pulchra
                           "#c11cad", # Acrpora nasuta
                           "#8338ec", # Acropora robusta
                           "#3E6FCB", # Acropora hyacinthus
                           "#3a86ff") # Acropora abrotanoides 

# ridge plot CpG methylation calls for each colony
ggplot(num.records.df, aes(x = num.records, y = species, fill = species)) +
  geom_density_ridges(scale = 1) +
  theme_classic() +
  labs(title = "Distribution of Number of CpG calls by Species",
       x = "Number of CpG methylation calls",
       y = "Species") +
  theme(legend.position = "none") +
  scale_fill_manual(values = species_color_palette )
# Add Gr.obj
# Methylation: Combine gene data with methylation object
gr.obj <- import("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/Ahyacinthuns.genes.gff") # This is the GFF file.
gr.df <- data.frame(gr.obj)
names(gr.df)
# Methylation: Make a column with chr.start.end
gr.df$concat <- paste(gr.df$seqnames,gr.df$start,gr.df$end,sep=".")

# Methylation: Region counts for each species
genes.meth.obj <- regionCounts(filtered.myobj,gr.obj,save.db=FALSE)

genes.unite <- methylKit::unite(genes.meth.obj,
                                     destrand = FALSE, #Combine
                                     min.per.group = 2L, # in at least 2 colony per species
                                     save.db = F)

genes.pm <- percMethylation(genes.unite, rowids = TRUE)

gene.names <- gr.df$ID[match(rownames(genes.pm),gr.df$concat)]
rownames(genes.pm) <- gene.names
gene.pm.means <- data.frame(gene = gene.names, 
                                   pm.means = rowMeans(genes.pm, na.rm = TRUE), 
                                   row.names = NULL)

PCASamples(meth, screeplot=TRUE)
pc=PCASamples(meth,obj.return = TRUE, adj.lim=c(1,1))
str(pc)
# Isolate PC1 and PC2 from PC result returned from methylKit
pc.df <- as.data.frame(pc$x[,1:2])

# Add column for species and colony
colony.id <- as.data.frame(rownames(pc$x))
pc.df$colony.id <- colony.id
pc.df$species <- meta$species
names(pc.df) <- c("PC1","PC2","colony.id", "species")




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
pc.df$species <- factor(pc.df$species, levels = species.order)





species.PCA.plot <- ggplot(pc.df,aes(x=pc.df[,1],y=pc.df[,2], color = species)) +
  geom_point(shape = 15, size = 5) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = official.species.colors ) +
  xlab("PCA1") + 
  ylab("PCA2") +
  labs(title="PCA of CpG methylation", color = "Species") 
#ggplotly(species.PCA.plot)
species.PCA.plot
# Not sure if the clustering is an artifact of SNP variation between species
# of true methylation differences. Note: Site level CpG Methylation.


# Gene-level DNA methylation matrix for each species:
# Make methylation object for each species

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
                  "C016_S38.rmdup_CpG.methylKit")

sample.id.pul=list("c004",
                   "c007",
                   "c013",
                   "c016")


# Methylation object creation
setwd("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
myobj.pul=methRead(file.list.pul,                   
                   sample.id=sample.id.pul,       
                   assembly="Ahyacinthus.chrsV1",
                   treatment=c(0,0,0,0),
                   context="CpG",
                   dbtype = "tabix",
                   dbdir = "methylDB")
filtered.myobj.pul=filterByCoverage(myobj.pul,
                                    lo.count=10,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)

# Acropora retusa (a)

file.list.reta <- list("C030_S37.rmdup_CpG.methylKit",
                  "C032_S33.rmdup_CpG.methylKit",
                  "C040_S32.rmdup_CpG.methylKit")

sample.id.reta=list("c030",
                   "c032",
                   "c040")

# Methylation object creation
setwd("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
myobj.reta=methRead(file.list.reta,                   
                   sample.id=sample.id.reta,       
                   assembly="Ahyacinthus.chrsV1",
                   treatment=c(0,0,0),
                   context="CpG",
                   dbtype = "tabix",
                   dbdir = "methylDB")
filtered.myobj.reta=filterByCoverage(myobj.reta,
                                    lo.count=10,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)

# Acropora retusa (b)

file.list.retb <- list("C094_S42.rmdup_CpG.methylKit",
                       "C121_S43.rmdup_CpG.methylKit")

sample.id.retb=list("c094",
                    "c121")

# Methylation object creation
setwd("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
myobj.retb=methRead(file.list.retb,                   
                    sample.id=sample.id.retb,       
                    assembly="Ahyacinthus.chrsV1",
                    treatment=c(0,0),
                    context="CpG",
                    dbtype = "tabix",
                    dbdir = "methylDB")
filtered.myobj.retb=filterByCoverage(myobj.retb,
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


# Acropora spp 2

file.list.spp2 <- list("C012_S56.rmdup_CpG.methylKit",
                  "C025_S57.rmdup_CpG.methylKit",
                  "C120_S58.rmdup_CpG.methylKit")
                  
sample.id.spp2=list("c012",
                   "c025",
                   "c120")

# Methylation object creation
setwd("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
myobj.spp2=methRead(file.list.spp2,                   
                   sample.id=sample.id.spp2,       
                   assembly="Ahyacinthus.chrsV1",
                   treatment=c(0,0,0),
                   context="CpG",
                   dbtype = "tabix",
                   dbdir = "methylDB")
filtered.myobj.spp2=filterByCoverage(myobj.spp2,
                                    lo.count=10,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)




# Methylation: Combine gene data with methylation object
gr.obj <- import("/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/Ahyacinthuns.genes.gff") # This is the GFF file.
gr.df <- data.frame(gr.obj)
names(gr.df)
# Methylation: Make a column with chr.start.end
gr.df$concat <- paste(gr.df$seqnames,gr.df$start,gr.df$end,sep=".")

# Methylation: Region counts for each species
genes.abro <- regionCounts(filtered.myobj.abro,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.hya <- regionCounts(filtered.myobj.hya,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.lut <- regionCounts(filtered.myobj.lut,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.pul <- regionCounts(filtered.myobj.pul,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.reta <- regionCounts(filtered.myobj.reta,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.retb <- regionCounts(filtered.myobj.retb,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.rob <- regionCounts(filtered.myobj.rob,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.tut <- regionCounts(filtered.myobj.tut,gr.obj,save.db=FALSE) # Gene counts for individuals
genes.spp2 <- regionCounts(filtered.myobj.spp2,gr.obj,save.db=FALSE) # Gene counts for individuals

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
genes.unite.reta <- methylKit::unite(genes.reta,
                                destrand = FALSE, #Combine
                                min.per.group = 1L, # in at least 1 colony per species
                                save.db = F)
genes.unite.retb <- methylKit::unite(genes.retb,
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
genes.unite.spp2 <- methylKit::unite(genes.spp2,
                                destrand = FALSE, #Combine
                                min.per.group = 1L, # in at least 1 colony per species
                                save.db = F)

# Methylation: Create methylation matrix per species
gene.pm.abro <- percMethylation(genes.unite.abro, rowids=TRUE)
gene.pm.hya <- percMethylation(genes.unite.hya, rowids=TRUE)
gene.pm.lut <- percMethylation(genes.unite.lut, rowids=TRUE)
gene.pm.pul <- percMethylation(genes.unite.pul, rowids=TRUE)
gene.pm.reta <- percMethylation(genes.unite.reta, rowids=TRUE)
gene.pm.retb <- percMethylation(genes.unite.retb, rowids=TRUE)
gene.pm.rob <- percMethylation(genes.unite.rob, rowids=TRUE)
gene.pm.tut <- percMethylation(genes.unite.tut, rowids=TRUE)
gene.pm.spp2 <- percMethylation(genes.unite.spp2, rowids=TRUE)

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

gene.names.reta <- gr.df$ID[match(rownames(gene.pm.reta),gr.df$concat)]
rownames(gene.pm.reta) <- gene.names.reta
gene.pm.means.reta.df <- data.frame(gene = gene.names.reta, 
                               pm.means = rowMeans(gene.pm.reta, na.rm = TRUE), 
                               row.names = NULL)
names(gene.pm.means.reta.df)<- c("gene","pm.mean.reta")

gene.names.retb <- gr.df$ID[match(rownames(gene.pm.retb),gr.df$concat)]
rownames(gene.pm.retb) <- gene.names.retb
gene.pm.means.retb.df <- data.frame(gene = gene.names.retb, 
                               pm.means = rowMeans(gene.pm.retb, na.rm = TRUE), 
                               row.names = NULL)
names(gene.pm.means.retb.df)<- c("gene","pm.mean.retb")

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

gene.names.spp2 <- gr.df$ID[match(rownames(gene.pm.spp2),gr.df$concat)]
rownames(gene.pm.spp2) <- gene.names.spp2
gene.pm.means.spp2.df <- data.frame(gene = gene.names.spp2, 
                               pm.means = rowMeans(gene.pm.spp2, na.rm = TRUE), 
                               row.names = NULL)
names(gene.pm.means.spp2.df)<- c("gene","pm.mean.spp2")

# Methylation: Combine everything
gene.pm.species<- full_join(gene.pm.means.abro.df, gene.pm.means.hya.df, by = "gene") %>%
  full_join(., gene.pm.means.lut.df, by = "gene") %>%
  full_join(., gene.pm.means.pul.df, by = "gene") %>%
  full_join(., gene.pm.means.reta.df, by = "gene") %>%
  full_join(., gene.pm.means.retb.df, by = "gene") %>%
  full_join(., gene.pm.means.rob.df, by = "gene") %>%
  full_join(., gene.pm.means.tut.df, by = "gene") %>%
  full_join(., gene.pm.means.spp2.df, by = "gene")

colSums(is.na(gene.pm.species))

# write.csv(gene.pm.species, file = "/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/gene.pm.species.prelim.csv", row.names = FALSE)

# Maybe looking at CHH and CHG methylation will provide some insight?
# Files for methylKit sorted in order to group the species
setwd("~/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
file.list.CHH <- list("C033_S22.rmdup_CHH.methylKit",
                  "C075_S29.rmdup_CHH.methylKit",
                  "C076_S28.rmdup_CHH.methylKit",
                  "C077_S30.rmdup_CHH.methylKit",
                  "C078_S26.rmdup_CHH.methylKit",
                  "C011_S25.rmdup_CHH.methylKit",
                  "C017_S23.rmdup_CHH.methylKit",
                  "C018_S27.rmdup_CHH.methylKit",
                  "C038_S31.rmdup_CHH.methylKit",
                  "C051_S24.rmdup_CHH.methylKit",
                  "C047_S36.rmdup_CHH.methylKit",
                  "C196_S35.rmdup_CHH.methylKit",
                  "C210_S41.rmdup_CHH.methylKit",
                  "C004_S34.rmdup_CHH.methylKit",
                  "C007_S39.rmdup_CHH.methylKit",
                  "C013_S40.rmdup_CHH.methylKit",
                  "C016_S38.rmdup_CHH.methylKit",
                  "C030_S37.rmdup_CHH.methylKit",
                  "C032_S33.rmdup_CHH.methylKit",
                  "C040_S32.rmdup_CHH.methylKit",
                  "C094_S42.rmdup_CHH.methylKit",
                  "C090_S44.rmdup_CHH.methylKit",
                  "C123_S45.rmdup_CHH.methylKit",
                  "C167_S46.rmdup_CHH.methylKit",
                  "C211_S47.rmdup_CHH.methylKit",
                  "C014_S48.rmdup_CHH.methylKit",
                  "C024_S49.rmdup_CHH.methylKit",
                  "C066_S50.rmdup_CHH.methylKit",
                  "C071_S51.rmdup_CHH.methylKit",
                  "C122_S52.rmdup_CHH.methylKit",
                  "C137_S53.rmdup_CHH.methylKit",
                  "C164_S54.rmdup_CHH.methylKit",
                  "C037_S55.rmdup_CHH.methylKit",
                  "C012_S56.rmdup_CHH.methylKit",
                  "C120_S58.rmdup_CHH.methylKit",
                  "C009_S59.rmdup_CHH.methylKit")

file.list.CHG <- list("C033_S22.rmdup_CHG.methylKit",
                      "C075_S29.rmdup_CHG.methylKit",
                      "C076_S28.rmdup_CHG.methylKit",
                      "C077_S30.rmdup_CHG.methylKit",
                      "C078_S26.rmdup_CHG.methylKit",
                      "C011_S25.rmdup_CHG.methylKit",
                      "C017_S23.rmdup_CHG.methylKit",
                      "C018_S27.rmdup_CHG.methylKit",
                      "C038_S31.rmdup_CHG.methylKit",
                      "C051_S24.rmdup_CHG.methylKit",
                      "C047_S36.rmdup_CHG.methylKit",
                      "C196_S35.rmdup_CHG.methylKit",
                      "C210_S41.rmdup_CHG.methylKit",
                      "C004_S34.rmdup_CHG.methylKit",
                      "C007_S39.rmdup_CHG.methylKit",
                      "C013_S40.rmdup_CHG.methylKit",
                      "C016_S38.rmdup_CHG.methylKit",
                      "C030_S37.rmdup_CHG.methylKit",
                      "C032_S33.rmdup_CHG.methylKit",
                      "C040_S32.rmdup_CHG.methylKit",
                      "C094_S42.rmdup_CHG.methylKit",
                      "C090_S44.rmdup_CHG.methylKit",
                      "C123_S45.rmdup_CHG.methylKit",
                      "C167_S46.rmdup_CHG.methylKit",
                      "C211_S47.rmdup_CHG.methylKit",
                      "C014_S48.rmdup_CHG.methylKit",
                      "C024_S49.rmdup_CHG.methylKit",
                      "C066_S50.rmdup_CHG.methylKit",
                      "C071_S51.rmdup_CHG.methylKit",
                      "C122_S52.rmdup_CHG.methylKit",
                      "C137_S53.rmdup_CHG.methylKit",
                      "C164_S54.rmdup_CHG.methylKit",
                      "C037_S55.rmdup_CHG.methylKit",
                      "C012_S56.rmdup_CHG.methylKit",
                      "C120_S58.rmdup_CHG.methylKit",
                      "C009_S59.rmdup_CHG.methylKit")


myobj.CHH=methRead(file.list.CHH,                   
               sample.id=sample.id,       
               assembly="Ahyacinthus.chrsV1",
               treatment=c(0,0,0,0,0,1,1,1,1,1,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,6,6,6,7,7,7,3),
               context="CpG",
               dbtype = "tabix",
               dbdir = "methylDB")
filtered.myobj.CHH=filterByCoverage(myobj.CHH,
                                lo.count=10,
                                lo.perc=NULL,
                                hi.count=NULL,
                                hi.perc=99.9)


sample.ids.CHH <- sapply(filtered.myobj.CHH, function(x) x@sample.id)
num.records.list.CHH <- sapply(filtered.myobj.CHH, function(x) x@num.records)
num.records.df.CHH <- data.frame(colony.id = sample.ids.CHH, num.records = num.records.list.CHH)

num.records.df.CHH$species <- meta[match(num.records.df$colony.id, meta$ID),3]


meth.CHH=methylKit::unite(filtered.myobj.CHH, 
                          min.per.group = 2L,
                          destrand=TRUE)
pm.CHH <- percMethylation(meth.CHH)
PCASamples(meth.CHH, screeplot=TRUE)
pc.CHH=PCASamples(meth.CHH,obj.return = TRUE, adj.lim=c(1,1))
str(pc.CHH)
# Isolate PC1 and PC2 from PC result returned from methylKit
pc.CHH.df <- as.data.frame(pc.CHH$x[,1:2])

# Add column for species and colony
colony.id <- as.data.frame(rownames(pc.CHH$x))
pc.CHH.df$colony.id <- colony.id
pc.CHH.df$species <- meta$species
names(pc.CHH.df) <- c("PC1","PC2","colony.id", "species")

# Plot PCA Results
# species_color_palette.1 <- c("#98F5E1", # Acropora abrotanoides
#                              "#8EECF5", # Acropora hyacinthus
#                              "#FDE4CF", # Acropora lutkeni
#                              "#F1C0E8", # Acropora pulchra
#                              "#A3C4F3", # Acropora retusa
#                              "#90DBF4", # Acropora robusta
#                              "#FBF8CC", # Acropora tutuilensis
#                              "#CFBAF0") # Acrpora spp

pc.CHH.df$species <- factor(pc.CHH.df$species, levels = species.order)
species.PCA.plot.CHH <- ggplot(pc.CHH.df,aes(x=pc.CHH.df[,1],y=pc.CHH.df[,2], color = species)) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = official.species.colors) +
  xlab("PCA1") + 
  ylab("PCA2") +
  labs(title="PCA of DNA methylation (CHH context)", color = "Species") +
  geom_point(size = 5, shape = 17) 
#ggplotly(species.PCA.plot.CHH)
species.PCA.plot.CHH


## CHG Methylation
file.list.CHG <- list("C033_S22.rmdup_CHG.methylKit",
                      "C075_S29.rmdup_CHG.methylKit",
                      "C076_S28.rmdup_CHG.methylKit",
                      "C077_S30.rmdup_CHG.methylKit",
                      "C078_S26.rmdup_CHG.methylKit",
                      "C011_S25.rmdup_CHG.methylKit",
                      "C017_S23.rmdup_CHG.methylKit",
                      "C018_S27.rmdup_CHG.methylKit",
                      "C038_S31.rmdup_CHG.methylKit",
                      "C051_S24.rmdup_CHG.methylKit",
                      "C047_S36.rmdup_CHG.methylKit",
                      "C196_S35.rmdup_CHG.methylKit",
                      "C210_S41.rmdup_CHG.methylKit",
                      "C004_S34.rmdup_CHG.methylKit",
                      "C007_S39.rmdup_CHG.methylKit",
                      "C013_S40.rmdup_CHG.methylKit",
                      "C016_S38.rmdup_CHG.methylKit",
                      "C030_S37.rmdup_CHG.methylKit",
                      "C032_S33.rmdup_CHG.methylKit",
                      "C040_S32.rmdup_CHG.methylKit",
                      "C094_S42.rmdup_CHG.methylKit",
                      "C090_S44.rmdup_CHG.methylKit",
                      "C123_S45.rmdup_CHG.methylKit",
                      "C167_S46.rmdup_CHG.methylKit",
                      "C211_S47.rmdup_CHG.methylKit",
                      "C014_S48.rmdup_CHG.methylKit",
                      "C024_S49.rmdup_CHG.methylKit",
                      "C066_S50.rmdup_CHG.methylKit",
                      "C071_S51.rmdup_CHG.methylKit",
                      "C122_S52.rmdup_CHG.methylKit",
                      "C137_S53.rmdup_CHG.methylKit",
                      "C164_S54.rmdup_CHG.methylKit",
                      "C037_S55.rmdup_CHG.methylKit",
                      "C012_S56.rmdup_CHG.methylKit",
                      "C120_S58.rmdup_CHG.methylKit",
                      "C009_S59.rmdup_CHG.methylKit")


myobj.CHG=methRead(file.list.CHG,                   
                   sample.id=sample.id,       
                   assembly="Ahyacinthus.chrsV1",
                   treatment=c(0,0,0,0,0,1,1,1,1,1,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,6,6,6,7,7,7,3),
                   context="CpG",
                   dbtype = "tabix",
                   dbdir = "methylDB")
filtered.myobj.CHG=filterByCoverage(myobj.CHG,
                                    lo.count=10,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)

meth.CHG=methylKit::unite(filtered.myobj.CHG, 
                          min.per.group = 2L,
                          destrand=TRUE) # update this unite function
pm.CHG <- percMethylation(meth.CHG)
PCASamples(meth.CHG, screeplot=TRUE)
pc.CHG=PCASamples(meth.CHG,obj.return = TRUE, adj.lim=c(1,1))
str(pc.CHG)


# Isolate PC1 and PC2 from PC result returned from methylKit
pc.CHG.df <- as.data.frame(pc.CHG$x[,1:2])

# Add column for species and colony
colony.id <- as.data.frame(rownames(pc.CHG$x))
pc.CHG.df$colony.id <- colony.id
pc.CHG.df$species <- meta$species
names(pc.CHG.df) <- c("PC1","PC2","colony.id", "species")

# Plot PCA Results


pc.CHG.df$species <- factor(pc.CHG.df$species, levels = species.order)

species.PCA.plot.CHG <- ggplot(pc.CHG.df,aes(x=pc.CHG.df[,1],y=pc.CHG.df[,2], color = species)) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = official.species.colors) +
  xlab("PCA1") + 
  ylab("PCA2") +
  labs(title="PCA of DNA methylation (CHG context)", color = "Species") +
  geom_point(size = 5, shape = 18) 
#ggplotly(species.PCA.plot.CHG)
species.PCA.plot.CHG

## There is a small number of calls for each context for 
## in the variant corrected files. What about for the 
## non variant corrected methylation calls? I expect
## there will be more calls, but now sure what the difference is.

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

meth=methylKit::unite(filtered.myobj, destrand=TRUE)
pm <- percMethylation(meth)
PCASamples(meth, screeplot=TRUE)
pc=PCASamples(meth,obj.return = TRUE, adj.lim=c(1,1))
str(pc)
# Isolate PC1 and PC2 from PC result returned from methylKit
pc.df <- as.data.frame(pc$x[,1:2])

# Add column for species and colony
colony.id <- as.data.frame(rownames(pc$x))
pc.df$colony.id <- colony.id
pc.df$species <- meta$species
names(pc.df) <- c("PC1","PC2","colony.id", "species")




# Plot PCA Results
species_color_palette.1 <- c("#98F5E1", # Acropora abrotanoides
                             "#8EECF5", # Acropora hyacinthus
                             "#FDE4CF", # Acropora lutkeni
                             "#F1C0E8", # Acropora pulchra
                             "#A3C4F3", # Acropora retusa
                             "#90DBF4", # Acropora robusta
                             "#FBF8CC", # Acropora tutuilensis
                             "#CFBAF0") # Acrpora spp


species.PCA.plot <- ggplot(pc.df,aes(x=pc.df[,1],y=pc.df[,2], color = species)) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = species_color_palette.1) +
  xlab("PCA1") + 
  ylab("PCA2") +
  labs(title="PCA of DNA methylation", color = "Species") +
  geom_point(size = 5) 
ggplotly(species.PCA.plot)

meth.CHG.dendro <- clusterSamples(meth.CHG, dist="correlation", method="ward", plot=FALSE)
ggdendrogram(meth.CHG.dendro, rotate = FALSE, size = 2)


## CHH and CHG Methylation PCA (variant filtered) with no A. tutuilensis:
# Load project meta data
my.dir <- "/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/"
meta.file.name <- "../meth_meta.csv"
meta.path <- file.path(my.dir, meta.file.name)
meta <- read.csv(meta.path, header = TRUE)
meta$ID <- as.factor(meta$ID)
meta$species <- as.factor(meta$species)

meta.no.tutu <- meta[-c(26:32),]


setwd("~/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
file.list.CHH.no.tutu <- list("C033_S22.rmdup_CHH.methylKit",
                      "C075_S29.rmdup_CHH.methylKit",
                      "C076_S28.rmdup_CHH.methylKit",
                      "C077_S30.rmdup_CHH.methylKit",
                      "C078_S26.rmdup_CHH.methylKit",
                      "C011_S25.rmdup_CHH.methylKit",
                      "C017_S23.rmdup_CHH.methylKit",
                      "C018_S27.rmdup_CHH.methylKit",
                      "C038_S31.rmdup_CHH.methylKit",
                      "C051_S24.rmdup_CHH.methylKit",
                      "C047_S36.rmdup_CHH.methylKit",
                      "C196_S35.rmdup_CHH.methylKit",
                      "C210_S41.rmdup_CHH.methylKit",
                      "C004_S34.rmdup_CHH.methylKit",
                      "C007_S39.rmdup_CHH.methylKit",
                      "C013_S40.rmdup_CHH.methylKit",
                      "C016_S38.rmdup_CHH.methylKit",
                      "C030_S37.rmdup_CHH.methylKit",
                      "C032_S33.rmdup_CHH.methylKit",
                      "C040_S32.rmdup_CHH.methylKit",
                      "C094_S42.rmdup_CHH.methylKit",
                      "C090_S44.rmdup_CHH.methylKit",
                      "C123_S45.rmdup_CHH.methylKit",
                      "C167_S46.rmdup_CHH.methylKit",
                      "C211_S47.rmdup_CHH.methylKit",
                      #"C014_S48.rmdup_CHH.methylKit",
                      #"C024_S49.rmdup_CHH.methylKit",
                      #"C066_S50.rmdup_CHH.methylKit",
                      #"C071_S51.rmdup_CHH.methylKit",
                      #"C122_S52.rmdup_CHH.methylKit",
                      #"C137_S53.rmdup_CHH.methylKit",
                      #"C164_S54.rmdup_CHH.methylKit",
                      "C037_S55.rmdup_CHH.methylKit",
                      "C012_S56.rmdup_CHH.methylKit",
                      "C120_S58.rmdup_CHH.methylKit",
                      "C009_S59.rmdup_CHH.methylKit")

file.list.CHG.no.tutu <- list("C033_S22.rmdup_CHG.methylKit",
                      "C075_S29.rmdup_CHG.methylKit",
                      "C076_S28.rmdup_CHG.methylKit",
                      "C077_S30.rmdup_CHG.methylKit",
                      "C078_S26.rmdup_CHG.methylKit",
                      "C011_S25.rmdup_CHG.methylKit",
                      "C017_S23.rmdup_CHG.methylKit",
                      "C018_S27.rmdup_CHG.methylKit",
                      "C038_S31.rmdup_CHG.methylKit",
                      "C051_S24.rmdup_CHG.methylKit",
                      "C047_S36.rmdup_CHG.methylKit",
                      "C196_S35.rmdup_CHG.methylKit",
                      "C210_S41.rmdup_CHG.methylKit",
                      "C004_S34.rmdup_CHG.methylKit",
                      "C007_S39.rmdup_CHG.methylKit",
                      "C013_S40.rmdup_CHG.methylKit",
                      "C016_S38.rmdup_CHG.methylKit",
                      "C030_S37.rmdup_CHG.methylKit",
                      "C032_S33.rmdup_CHG.methylKit",
                      "C040_S32.rmdup_CHG.methylKit",
                      "C094_S42.rmdup_CHG.methylKit",
                      "C090_S44.rmdup_CHG.methylKit",
                      "C123_S45.rmdup_CHG.methylKit",
                      "C167_S46.rmdup_CHG.methylKit",
                      "C211_S47.rmdup_CHG.methylKit",
                      # "C014_S48.rmdup_CHG.methylKit",
                      # "C024_S49.rmdup_CHG.methylKit",
                      # "C066_S50.rmdup_CHG.methylKit",
                      # "C071_S51.rmdup_CHG.methylKit",
                      # "C122_S52.rmdup_CHG.methylKit",
                      # "C137_S53.rmdup_CHG.methylKit",
                      # "C164_S54.rmdup_CHG.methylKit",
                      "C037_S55.rmdup_CHG.methylKit",
                      "C012_S56.rmdup_CHG.methylKit",
                      "C120_S58.rmdup_CHG.methylKit",
                      "C009_S59.rmdup_CHG.methylKit")
# Samples
sample.id.no.tutu=list("c033",
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
               # "c014",
               # "c024",
               # "c066",
               # "c071",
               # "c122",
               # "c137",
               # "c164",
               "c037",
               "c012",
               "c120",
               "c009")

myobj.CHH.no.tutu=methRead(file.list.CHH.no.tutu,                   
                   sample.id=sample.id.no.tutu,       
                   assembly="Ahyacinthus.chrsV1",
                   treatment=c(0,0,0,0,0,1,1,1,1,1,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,7,7,7,3),
                   context="CpG",
                   dbtype = "tabix",
                   dbdir = "methylDB")
filtered.myobj.CHH.no.tutu=filterByCoverage(myobj.CHH.no.tutu,
                                    lo.count=10,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)

meth.CHH.no.tutu=methylKit::unite(filtered.myobj.CHH.no.tutu, destrand=TRUE)
pc.CHH.no.tutu=PCASamples(meth.CHH.no.tutu,obj.return = TRUE, adj.lim=c(1,1))

## CHG no tutu
myobj.CHG.no.tutu=methRead(file.list.CHG.no.tutu,                   
                           sample.id=sample.id.no.tutu,       
                           assembly="Ahyacinthus.chrsV1",
                           treatment=c(0,0,0,0,0,1,1,1,1,1,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,7,7,7,3),
                           context="CpG",
                           dbtype = "tabix",
                           dbdir = "methylDB")
filtered.myobj.CHG.no.tutu=filterByCoverage(myobj.CHG.no.tutu,
                                            lo.count=10,
                                            lo.perc=NULL,
                                            hi.count=NULL,
                                            hi.perc=99.9)

meth.CHG.no.tutu=methylKit::unite(filtered.myobj.CHG.no.tutu, destrand=TRUE)
pc.CHG.no.tutu=PCASamples(meth.CHG.no.tutu,obj.return = TRUE, adj.lim=c(1,1))


# Plot no tutu CHH and CHG methylation PCA:
# Isolate PC1 and PC2 from PC result returned from methylKit
pc.df.CHH.no.tutu <- as.data.frame(pc.CHH.no.tutu$x[,1:2])
pc.df.CHG.no.tutu <- as.data.frame(pc.CHG.no.tutu$x[,1:2])


# Add column for species and colony
colony.id <- as.data.frame(rownames(pc.CHH.no.tutu$x))
pc.df.CHH.no.tutu$colony.id <- colony.id
pc.df.CHH.no.tutu$species <- meta.no.tutu$species
names(pc.df.CHH.no.tutu) <- c("PC1","PC2","colony.id", "species")

colony.id <- as.data.frame(rownames(pc.CHG.no.tutu$x))
pc.df.CHG.no.tutu$colony.id <- colony.id
pc.df.CHG.no.tutu$species <- meta.no.tutu$species
names(pc.df.CHG.no.tutu) <- c("PC1","PC2","colony.id", "species")

# Plot PCA Results
official.species.colors.no.tutu <- official.species.colors[-1] # remove tutu?
pc.df.CHG.no.tutu$species <- factor(pc.df.CHG.no.tutu$species, levels = species.order)
pc.df.CHH.no.tutu$species <- factor(pc.df.CHH.no.tutu$species, levels = species.order)

species.PCA.plot.CHH.no.tutu <- ggplot(pc.df.CHH.no.tutu,aes(x=pc.df.CHH.no.tutu[,1],y=pc.df.CHH.no.tutu[,2], color = species)) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = official.species.colors.no.tutu) +
  xlab("PCA1") + 
  ylab("PCA2") +
  labs(title="PCA of DNA CHH methylation\n (A. tutuilensis removed)", color = "Species") +
  geom_point(size = 5, shape = 17) 
#ggplotly(species.PCA.plot.CHH.no.tutu)
species.PCA.plot.CHH.no.tutu 

species.PCA.plot.CHG.no.tutu <- ggplot(pc.df.CHG.no.tutu,aes(x=pc.df.CHG.no.tutu[,1],y=pc.df.CHG.no.tutu[,2], color = species)) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = official.species.colors.no.tutu) +
  xlab("PCA1") + 
  ylab("PCA2") +
  labs(title="PCA of DNA CHG methylation\n (A. tutuilensis removed)", color = "Species") +
  geom_point(size = 5, shape = 18) 
#ggplotly(species.PCA.plot.CHG.no.tutu)
species.PCA.plot.CHG.no.tutu 


#### PCA: NO ABRO
# Load project meta data
my.dir <- "/Users/tillandsia/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/"
meta.file.name <- "../meth_meta.csv"
meta.path <- file.path(my.dir, meta.file.name)
meta <- read.csv(meta.path, header = TRUE)
meta$ID <- as.factor(meta$ID)
meta$species <- as.factor(meta$species)

meta.no.abro <- meta[-c(1:5),]


setwd("~/Lab Notebook/Chapter3/WGBS_analysis/15_methyldackel_out_exc_vars/")
file.list.no.abro <- list( #"C033_S22.rmdup_CpG.methylKit",
                  #"C075_S29.rmdup_CpG.methylKit",
                  #"C076_S28.rmdup_CpG.methylKit",
                 # "C077_S30.rmdup_CpG.methylKit",
                  #"C078_S26.rmdup_CpG.methylKit",
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
sample.id.no.abro=list(#"c033",
               #"c075",
               #"c076",
              # "c077",
              # "c078",
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
myobj.no.abro=methRead(file.list.no.abro,                   
               sample.id=sample.id.no.abro,       
               assembly="Ahyacinthus.chrsV1",
               treatment=c(0,0,0,0,0,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,5,5,5,6,6,6,2),
               context="CpG",
               dbtype = "tabix",
               dbdir = "methylDB")
filtered.myobj.no.abro=filterByCoverage(myobj.no.abro,
                                lo.count=10,
                                lo.perc=NULL,
                                hi.count=NULL,
                                hi.perc=99.9)


meth.no.abro=methylKit::unite(filtered.myobj.no.abro, 
                      min.per.group = 2L, # site found in at least 2 individual/species
                      destrand=TRUE)

PCASamples(meth.no.abro, screeplot=TRUE)
pc <- methylKit::PCASamples(meth.no.abro,obj.return=TRUE, adj.lim=c(1,1))
str(pc)
# Isolate PC1 and PC2 from PC result returned from methylKit
pc.df <- as.data.frame(pc$x[,1:2])

# Add column for species and colony
colony.id <- as.data.frame(rownames(pc$x))
pc.df$colony.id <- colony.id
pc.df$species <- meta.no.abro$species
names(pc.df) <- c("PC1","PC2","colony.id", "species")




# Plot PCA Results
official.species.colors <- c("Acropora tutuilensis" = "#F7C11E",
                             "Acropora lutkeni" = "#FF9A17",
                             "Acropora pulchra" = "#fb5607",
                             "Acropora retusa" = "#ff006e",
                             "Acropora nasuta" = "#c11cad",
                             "Acropora robusta" = "#8338ec",
                             "Acropora hyacinthus" = "#3a86ff")

species.order <- c("Acropora tutuilensis", "Acropora lutkeni",
                   "Acropora pulchra", "Acropora retusa","Acropora nasuta",
                   "Acropora robusta","Acropora hyacinthus")
pc.df$species <- factor(pc.df$species, levels = species.order)





species.PCA.plot <- ggplot(pc.df,aes(x=pc.df[,1],y=pc.df[,2], color = species)) +
  geom_point(shape = 15, size = 5) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = official.species.colors ) +
  xlab("PCA1") + 
  ylab("PCA2") +
  labs(title="PCA of CpG methylation", color = "Species") +
  theme(legend.position = "none")
#ggplotly(species.PCA.plot)

species.PCA.plot + PCA_plot_hs_species_no_abro 

## No abro gene expression
library(DESeq2)
library(tidyverse)

# File Set-up
my.dir <- "/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/HTseq_out"
my.files <- grep(".txt", list.files(my.dir), value=TRUE)
my.metadata <- read.csv("/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/tagseq_meta_data.csv", header = TRUE)
my.metadata.no.abro <- my.metadata[-c(1:10),]
# Create sample table
sampleNames <- my.metadata.no.abro$seq.sample.id
my.sampleTable <- data.frame(sampleName = sampleNames, 
                             fileName = my.files[-(1:10)], 
                             condition = my.metadata.no.abro)
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
deseq_PCA_heat_stress_species$species <- my.metadata.no.abro$species

# official.species.colors <- c("Acropora tutuilensis" = "#F7C11E",
#                              "Acropora lutkeni" = "#FF9A17",
#                              "Acropora pulchra" = "#fb5607",
#                              "Acropora retusa" = "#ff006e",
#                              "Acropora sp" = "#c11cad",
#                              "Acropora robusta" = "#8338ec",
#                              "Acropora hyacinthus" = "#3a86ff",
#                              "Acropora abrotanoides" = "#3E6FCB")
# 
# species.order <- c("Acropora tutuilensis", "Acropora lutkeni",
#                    "Acropora pulchra", "Acropora retusa","Acropora sp","Acropora robusta","Acropora hyacinthus","Acropora abrotanoides")
# color.data.ed50$species <- factor(color.data.ed50$species, levels = species.order)

official.species.colors <- c("Acropora.tutuilensis" = "#F7C11E",
                             "Acropora.lutkeni" = "#FF9A17",
                             "Acropora.pulchra" = "#fb5607",
                             "Acropora.retusa" = "#ff006e",
                             "Acropora.unknown" = "#c11cad",
                             "Acropora.robusta" = "#8338ec",
                             "Acropora.hyacinthus" = "#3a86ff")

species.order <- c("Acropora.tutuilensis", "Acropora.lutkeni",
                   "Acropora.pulchra", "Acropora.retusa","Acropora.unknown",
                   "Acropora.robusta","Acropora.hyacinthus")
deseq_PCA_heat_stress_species$species <- factor(deseq_PCA_heat_stress_species$species, levels = species.order)


# Plot PCA Heat Stress
legend.names <- c("Acropora tutuilensis", "Acropora lutkeni",
                  "Acropora pulchra", "Acropora retusa","Acropora sp","Acropora robusta","Acropora hyacinthus") 
shape.names <- c('29 °C','36 °C') 
PCA_plot_hs_species <- ggplot(deseq_PCA_heat_stress_species,aes(x=PC1,y=PC2, color=species, shape=condition.treatment)) +
  theme_classic(base_size = 15) +
  labs(title="PCA of gene expression colored by species") +
  geom_point(size = 5, alpha = 0.7) +
  scale_color_manual(values = official.species.colors, labels = legend.names) +  
  #scale_shape_manual(labels = shape.names) +
  labs(color = "Species", shape = "Treatment") 



PCA_plot_hs_species


### Percent methylation of genes in yellow module across species
gene.pm.means.abro.df$species <- "Acropora abrotanoides"
gene.pm.means.hya.df$species <- "Acropora hyacinthus"
gene.pm.means.lut.df$species <- "Acropora lutkeni"
gene.pm.means.rob.df$species <- "Acropora robusta"
gene.pm.means.pul.df$species <- "Acropora pulchra"
gene.pm.means.ret.df$species <- "Acropora retusa"
gene.pm.means.unk.df$species <- "Acropora unknown"
gene.pm.means.tut.df$species <- "Acropora tutuilensis"

### Percent methylation of yellow module genes:
gene.pm.means.abro.df.yellow <- gene.pm.means.abro.df$gene %in% MEyellow.genes
sum(gene.pm.means.abro.df.yellow) #1514
gene.pm.means.abro.df.salmon <- gene.pm.means.abro.df$gene %in% MEsalmon.genes
sum(gene.pm.means.abro.df.salmon) #264
gene.pm.means.abro.df$me.yellow <- gene.pm.means.abro.df.yellow 
gene.pm.means.abro.df$me.salmon <- gene.pm.means.abro.df.salmon
gene.pm.means.abro.df <- gene.pm.means.abro.df %>%
  dplyr::mutate(no.module = !me.yellow & !me.salmon)


gene.pm.means.hya.df.yellow <- gene.pm.means.hya.df$gene %in% MEyellow.genes
sum(gene.pm.means.hya.df.yellow) #2103
gene.pm.means.hya.df.salmon <- gene.pm.means.hya.df$gene %in% MEsalmon.genes
sum(gene.pm.means.hya.df.salmon) #363
gene.pm.means.hya.df$me.yellow <- gene.pm.means.hya.df.yellow 
gene.pm.means.hya.df$me.salmon <- gene.pm.means.hya.df.salmon
gene.pm.means.hya.df <- gene.pm.means.hya.df %>%
  dplyr::mutate(no.module = !me.yellow & !me.salmon)

gene.pm.means.lut.df.yellow <- gene.pm.means.lut.df$gene %in% MEyellow.genes
sum(gene.pm.means.lut.df.yellow) #1809
gene.pm.means.lut.df.salmon <- gene.pm.means.lut.df$gene %in% MEsalmon.genes
sum(gene.pm.means.lut.df.salmon) #303
gene.pm.means.lut.df$me.yellow <- gene.pm.means.lut.df.yellow 
gene.pm.means.lut.df$me.salmon <- gene.pm.means.lut.df.salmon
gene.pm.means.lut.df <- gene.pm.means.lut.df %>%
  dplyr::mutate(no.module = !me.yellow & !me.salmon)

gene.pm.means.rob.df.yellow <- gene.pm.means.rob.df$gene %in% MEyellow.genes
sum(gene.pm.means.rob.df.yellow) #1519
gene.pm.means.rob.df.salmon <- gene.pm.means.rob.df$gene %in% MEsalmon.genes
sum(gene.pm.means.rob.df.salmon) #247
gene.pm.means.rob.df$me.yellow <- gene.pm.means.rob.df.yellow 
gene.pm.means.rob.df$me.salmon <- gene.pm.means.rob.df.salmon
gene.pm.means.rob.df <- gene.pm.means.rob.df %>%
  dplyr::mutate(no.module = !me.yellow & !me.salmon)

gene.pm.means.pul.df.yellow <- gene.pm.means.pul.df$gene %in% MEyellow.genes
sum(gene.pm.means.pul.df.yellow) #1950
gene.pm.means.pul.df.salmon <- gene.pm.means.pul.df$gene %in% MEsalmon.genes
sum(gene.pm.means.pul.df.salmon) #333
gene.pm.means.pul.df$me.yellow <- gene.pm.means.pul.df.yellow 
gene.pm.means.pul.df$me.salmon <- gene.pm.means.pul.df.salmon
gene.pm.means.pul.df <- gene.pm.means.pul.df %>%
  dplyr::mutate(no.module = !me.yellow & !me.salmon)

gene.pm.means.ret.df.yellow <- gene.pm.means.ret.df$gene %in% MEyellow.genes
sum(gene.pm.means.ret.df.yellow) #1759
gene.pm.means.ret.df.salmon <- gene.pm.means.ret.df$gene %in% MEsalmon.genes
sum(gene.pm.means.ret.df.salmon) #287
gene.pm.means.ret.df$me.yellow <- gene.pm.means.ret.df.yellow 
gene.pm.means.ret.df$me.salmon <- gene.pm.means.ret.df.salmon
gene.pm.means.ret.df <- gene.pm.means.ret.df %>%
  dplyr::mutate(no.module = !me.yellow & !me.salmon)

gene.pm.means.unk.df.yellow <- gene.pm.means.unk.df$gene %in% MEyellow.genes
sum(gene.pm.means.unk.df.yellow) #1625
gene.pm.means.unk.df.salmon <- gene.pm.means.unk.df$gene %in% MEsalmon.genes
sum(gene.pm.means.unk.df.salmon) #269
gene.pm.means.unk.df$me.yellow <- gene.pm.means.unk.df.yellow 
gene.pm.means.unk.df$me.salmon <- gene.pm.means.unk.df.salmon
gene.pm.means.unk.df <- gene.pm.means.unk.df %>%
  dplyr::mutate(no.module = !me.yellow & !me.salmon)

gene.pm.means.tut.df.yellow <- gene.pm.means.tut.df$gene %in% MEyellow.genes
sum(gene.pm.means.tut.df.yellow) #1769
gene.pm.means.tut.df.salmon <- gene.pm.means.tut.df$gene %in% MEsalmon.genes
sum(gene.pm.means.tut.df.salmon) #300
gene.pm.means.tut.df$me.yellow <- gene.pm.means.tut.df.yellow 
gene.pm.means.tut.df$me.salmon <- gene.pm.means.tut.df.salmon
gene.pm.means.tut.df <- gene.pm.means.tut.df %>%
  dplyr::mutate(no.module = !me.yellow & !me.salmon)

# Combine data frames for plotting salmon gene pm
pm.salmon.abro <- gene.pm.means.abro.df[gene.pm.means.abro.df$me.salmon,]
pm.salmon.abro <- pm.salmon.abro %>%
  rename(pm.mean = pm.mean.abro)
pm.salmon.hya <- gene.pm.means.hya.df[gene.pm.means.hya.df$me.salmon,]
pm.salmon.hya <- pm.salmon.hya %>%
  rename(pm.mean = pm.mean.hya)
pm.salmon.pul <- gene.pm.means.pul.df[gene.pm.means.pul.df$me.salmon,]
pm.salmon.pul <- pm.salmon.pul %>%
  rename(pm.mean = pm.mean.pul)
pm.salmon.lut <- gene.pm.means.lut.df[gene.pm.means.lut.df$me.salmon,]
pm.salmon.lut <- pm.salmon.lut %>%
  rename(pm.mean = pm.mean.lut)
pm.salmon.ret <- gene.pm.means.ret.df[gene.pm.means.ret.df$me.salmon,]
pm.salmon.ret <- pm.salmon.ret %>%
  rename(pm.mean = pm.mean.ret)
pm.salmon.rob <- gene.pm.means.rob.df[gene.pm.means.rob.df$me.salmon,]
pm.salmon.rob <- pm.salmon.rob %>%
  rename(pm.mean = pm.mean.rob)
pm.salmon.tut <- gene.pm.means.tut.df[gene.pm.means.tut.df$me.salmon,]
pm.salmon.tut <- pm.salmon.tut %>%
  rename(pm.mean = pm.mean.tut)
pm.salmon.unk <- gene.pm.means.unk.df[gene.pm.means.unk.df$me.salmon,]
pm.salmon.unk <- pm.salmon.unk %>%
  rename(pm.mean = pm.mean.unk)

salmon.df.list <- list(pm.salmon.abro, 
                       pm.salmon.hya, 
                       pm.salmon.pul,
                       pm.salmon.lut,
                       pm.salmon.ret,
                       pm.salmon.rob, 
                       pm.salmon.tut, 
                       pm.salmon.unk)
all.species.pm.salmon.df <- do.call(rbind, salmon.df.list)
official.species.colors <- c("Acropora tutuilensis" = "#F7C11E",
                             "Acropora lutkeni" = "#FF9A17",
                             "Acropora pulchra" = "#fb5607",
                             "Acropora retusa" = "#ff006e",
                             "Acropora unknown" = "#c11cad",
                             "Acropora robusta" = "#8338ec",
                             "Acropora hyacinthus" = "#3a86ff",
                             "Acropora abrotanoides" = "#3E6FCB")

species.order <- c("Acropora tutuilensis", "Acropora lutkeni",
                   "Acropora pulchra", "Acropora retusa","Acropora unknown",
                   "Acropora robusta","Acropora hyacinthus","Acropora abrotanoides")
all.species.pm.salmon.df$species <- factor(all.species.pm.salmon.df$species, levels = species.order)

# Plot
rain_height <- .1

pm.plot.all.samples.salmon <- ggplot(all.species.pm.salmon.df, aes(x = species, y = log(pm.mean), fill = species )) +
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
  ggtitle("Percent methylation\nof genes in Salmon Module")
pm.plot.all.samples.salmon

# Combine data frames for plotting yellow gene pm
pm.yellow.abro <- gene.pm.means.abro.df[gene.pm.means.abro.df$me.yellow,]
pm.yellow.abro <- pm.yellow.abro %>%
  rename(pm.mean = pm.mean.abro)
pm.yellow.hya <- gene.pm.means.hya.df[gene.pm.means.hya.df$me.yellow,]
pm.yellow.hya <- pm.yellow.hya %>%
  rename(pm.mean = pm.mean.hya)
pm.yellow.pul <- gene.pm.means.pul.df[gene.pm.means.pul.df$me.yellow,]
pm.yellow.pul <- pm.yellow.pul %>%
  rename(pm.mean = pm.mean.pul)
pm.yellow.lut <- gene.pm.means.lut.df[gene.pm.means.lut.df$me.yellow,]
pm.yellow.lut <- pm.yellow.lut %>%
  rename(pm.mean = pm.mean.lut)
pm.yellow.ret <- gene.pm.means.ret.df[gene.pm.means.ret.df$me.yellow,]
pm.yellow.ret <- pm.yellow.ret %>%
  rename(pm.mean = pm.mean.ret)
pm.yellow.rob <- gene.pm.means.rob.df[gene.pm.means.rob.df$me.yellow,]
pm.yellow.rob <- pm.yellow.rob %>%
  rename(pm.mean = pm.mean.rob)
pm.yellow.tut <- gene.pm.means.tut.df[gene.pm.means.tut.df$me.yellow,]
pm.yellow.tut <- pm.yellow.tut %>%
  rename(pm.mean = pm.mean.tut)
pm.yellow.unk <- gene.pm.means.unk.df[gene.pm.means.unk.df$me.yellow,]
pm.yellow.unk <- pm.yellow.unk %>%
  rename(pm.mean = pm.mean.unk)

yellow.df.list <- list(pm.yellow.abro, 
                       pm.yellow.hya, 
                       pm.yellow.pul,
                       pm.yellow.lut,
                       pm.yellow.ret,
                       pm.yellow.rob, 
                       pm.yellow.tut, 
                       pm.yellow.unk)
all.species.pm.yellow.df <- do.call(rbind, yellow.df.list)
official.species.colors <- c("Acropora tutuilensis" = "#F7C11E",
                             "Acropora lutkeni" = "#FF9A17",
                             "Acropora pulchra" = "#fb5607",
                             "Acropora retusa" = "#ff006e",
                             "Acropora unknown" = "#c11cad",
                             "Acropora robusta" = "#8338ec",
                             "Acropora hyacinthus" = "#3a86ff",
                             "Acropora abrotanoides" = "#3E6FCB")

species.order <- c("Acropora tutuilensis", "Acropora lutkeni",
                   "Acropora pulchra", "Acropora retusa","Acropora unknown",
                   "Acropora robusta","Acropora hyacinthus","Acropora abrotanoides")
all.species.pm.yellow.df$species <- factor(all.species.pm.yellow.df$species, levels = species.order)

# Plot
rain_height <- .1
epsilon <- 1e-6
all.species.pm.yellow.df$perc.meth.proportion <- all.species.pm.yellow.df$pm.mean/ 100
all.species.pm.yellow.df$perc.meth.logit <- log((all.species.pm.yellow.df$perc.meth.proportion + epsilon) / (1 - all.species.pm.yellow.df$perc.meth.proportion + epsilon))

pm.plot.all.samples.yellow <- ggplot(all.species.pm.yellow.df, aes(x = species, y = scale(pm.mean), fill = species )) +
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
  ggtitle("Percent methylation\nof genes in Yellow Module")
pm.plot.all.samples.yellow

pm.plot.all.samples.yellow + pm.plot.all.samples.salmon



##### Caluculating differential methylation across 8 species
pooled.meth=pool(meth,sample.ids=c("abro","hya", "lut","pul", "ret", "rob", "tut", "nas" ))
dm.pooledf=calculateDiffMeth(pooled.meth)

#### Genes instead of sites
class(genes.unite)
pooled.meth.genes=pool(genes.unite,sample.ids=c("abro","hya", "lut","pul", "ret", "rob", "tut", "nas" ))
dm.genes.pooledf=calculateDiffMeth(pooled.meth.genes, adjust = "fdr")


# Methylation: Attach gene names to the data
diff.meth.species.df  <- getData(dm.genes.pooledf)
diff.meth.species.df$concat <-paste(diff.meth.species.df$chr,diff.meth.species.df$start,diff.meth.species.df$end,sep=".")
diff.meth.species.df <- diff.meth.species.df %>%
  dplyr::select(concat, everything())

gene.names.dm.species <- gr.df$ID[match((diff.meth.species.df$concat),gr.df$concat)]

rownames(diff.meth.species.df) <- gene.names.dm.species

diff.meth.species.df <- diff.meth.species.df %>%
  dplyr::select(-concat, - chr, -start, -end)

diff.meth.species.df <- rownames_to_column(diff.meth.species.df, var = "rowname")
colnames(diff.meth.species.df)[1] <- "gene"

hist(diff.meth.species.df$qvalue, breaks = 100)
hist(diff.meth.species.df$meth.diff, breaks = 100)
qvalue.filt.diff.meth.species.df <- diff.meth.species.df %>%
  dplyr::filter(qvalue < 0.05). # 2860 genes match this criteria

# How many of these genes are in the salmon or yellow modules?


