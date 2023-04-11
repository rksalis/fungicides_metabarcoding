## DADA2 ITS pipeline via Bioconductor DADA2 Pipeline 
# R.K. Salis 20.10.2022

## Install DADA2 ##
install.packages("Rcpp")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("dada2")
BiocManager::install("ShortRead")
BiocManager::install("Biostrings")
BiocManager::install("phyloseq")
library(dada2); packageVersion("dada2") #1.14.1
library(ShortRead); packageVersion("ShortRead") #1.44.3
library(Biostrings); packageVersion("Biostrings") #2.54.0
library(seqinr); packageVersion("seqinr")
library(phyloseq); packageVersion("phyloseq") #1.30.0
library(ggplot2); packageVersion("ggplot2") #3.3.0
library(DESeq2); packageVersion("Biostrings")
library(plyr); packageVersion("plyr")
library(tidyverse) ; packageVersion("tidyverse") # 1.3.1
library(vegan) ; packageVersion("vegan") # 2.5.4
library(dendextend) ; packageVersion("dendextend") # 1.10.0
library(viridis) ; packageVersion("viridis")


############### MiSeq run 1 ###############
setwd("/fwd_only/run1/output")
path = "/fwd_only/run1/raw"
list.files(path)
f.names = as.vector(list.files(path, pattern = ".1.fastq.gz",full.names = F)) # extract the fwd seq file names as a vector
fnFs = sort(list.files(path, pattern = ".1.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "..?.fastq.gz"), `[`, 1)
FWD <- "GTGARTCATCGAATCTTTG"  ## fITS7 forward primer sequence 19 bp
REV <- "TCCTCCGCTTATTGATATGC"  ## ITS4 reverse primer sequence 20 bp
#verify the presence and orientation of the primers in the data
allOrients <- function(primer) { # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients
# “pre-filter” the sequences just to remove those with ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
filterAndTrim(fnFs, fnFs.filtN,  maxN = 0, multithread = TRUE)
# count the number of times the primers appear in the forward and reverse read, while considering all possible primer orientations.
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]))
#Remove Primers - with cutadapt
cutadapt <- "/miniconda3/bin/cutadapt" 
system2(cutadapt, args = "--version") 
path.cut <- "/WP1_ITS/fwd_only/run1/cutadapt"
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags,  "-n", 2, 
                             "-o", fnFs.cut[i], 
                             fnFs.filtN[i])) 
}
#read in the names of the cutadapt-ed FASTQ files and get the matched lists of forward and reverse fastq files.
list.files(path.cut)
cutFs <- sort(list.files(path.cut, pattern = ".1.fastq.gz", full.names = TRUE))
# Extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), "..?.fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
## Plot quality scores from .fastq files ##
plotQualityProfile(fnFs[1:2])
## Filter and trim (quality control) ##
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
#Filtering paraments: maxN=0 (DADA2 requires sequences contain no Ns), truncQ = 2,  rm.phix = TRUE and maxEE=2 (this means a maximum number of “expected errors” is allowed in a read)
out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = 2, 
                     truncQ = 2, minLen = 50, truncLen=150, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  
head(out)
saveRDS(out, "out.rds")
# File parsing
filtpath <- "/WP1_ITS/fwd_only/run1/cutadapt/filtered" #directory containing filtered fastqs
list.files(filtpath)
sample.names <- sapply(strsplit(basename(filtFs), "..?.fastq.gz"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
names(filtFs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Infer sequence variants
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filtFs[[sam]])
  dds[[sam]] <- dada(derep, err=errF, multithread=TRUE)
}
## Generate sequence table ##
seqtab <- makeSequenceTable(dds)
saveRDS(seqtab, "seqtab.rds")
dim(seqtab)
table(nchar(getSequences(seqtab)))

############### MiSeq run 2 ###############
setwd("/WP1_ITS/fwd_only/run2/output")
path = "/WP1_ITS/fwd_only/run2/raw"
list.files(path)
f.names = as.vector(list.files(path, pattern = ".1.fastq.gz",full.names = F)) # extract the fwd seq file names as a vector
fnFs = sort(list.files(path, pattern = ".1.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "..?.fastq.gz"), `[`, 1)
FWD <- "GTGARTCATCGAATCTTTG"  ## fITS7 forward primer sequence 19 bp
REV <- "TCCTCCGCTTATTGATATGC"  ## ITS4 reverse primer sequence 20 bp
#verify the presence and orientation of the primers in the data
allOrients <- function(primer) { # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
# “pre-filter” the sequences just to remove those with ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
filterAndTrim(fnFs, fnFs.filtN,  maxN = 0, multithread = TRUE)
# count the number of times the primers appear in the forward and reverse read, while considering all possible primer orientations.
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]))
#Remove Primers - with cutadapt
cutadapt <- "/miniconda3/bin/cutadapt" 
system2(cutadapt, args = "--version") 
path.cut <- "/WP1_ITS/fwd_only/run2/cutadapt"
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags,  "-n", 2, 
                             "-o", fnFs.cut[i], 
                             fnFs.filtN[i])) 
}
#read in the names of the cutadapt-ed FASTQ files and get the matched lists of forward and reverse fastq files.
list.files(path.cut)
cutFs <- sort(list.files(path.cut, pattern = ".1.fastq.gz", full.names = TRUE))
# Extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), "..?.fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
## Plot quality scores from .fastq files ##
plotQualityProfile(fnFs[1:2])
## Filter and trim (quality control) ##
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
#Filtering paraments: maxN=0 (DADA2 requires sequences contain no Ns), truncQ = 2,  rm.phix = TRUE and maxEE=2 (this means a maximum number of “expected errors” is allowed in a read)
out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = 2, 
                     truncQ = 2, minLen = 50, truncLen=150, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  
head(out)
saveRDS(out, "out.rds")
# File parsing
filtpath <- "/WP1_ITS/fwd_only/run2/cutadapt/filtered" #directory containing filtered fastqs
list.files(filtpath)
sample.names <- sapply(strsplit(basename(filtFs), "..?.fastq.gz"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
names(filtFs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Infer sequence variants
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filtFs[[sam]])
  dds[[sam]] <- dada(derep, err=errF, multithread=TRUE)
}
## Generate sequence table ##
seqtab <- makeSequenceTable(dds)
saveRDS(seqtab, "seqtab.rds")
dim(seqtab)
table(nchar(getSequences(seqtab)))

########################  Merge runs ##########################
setwd("/WP1_analyses/bioinformatics/ITS")
st1 <- readRDS("/WP1_ITS/fwd_only/run1/output/seqtab.rds")
st2 <- readRDS("/WP1_ITS/fwd_only/run2/output/seqtab.rds")
st.all <- mergeSequenceTables(st1, st2)
saveRDS(st.all, "seqtab.all.rds")
dim(st.all)
table(nchar(getSequences(st.all)))

## Remove chimeras ##
seqtab.chi <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.chi)/sum(st.all)
saveRDS(seqtab.chi, "seqtab.chi.rds", version=2)
#Inspect distribution of sequence lengths:
table(nchar(getSequences(seqtab.chi)))
## Track sequence loss ##
out1 <- readRDS("/WP1_ITS/fwd_only/run1/output/out.rds")
out2 <- readRDS("/WP1_ITS/fwd_only/run2/output/out.rds")
out.all <- mergeSequenceTables(out1, out2)
seqtab.all <- readRDS("seqtab.all.rds")
seqtab.chi <- readRDS("seqtab.chi.rds")
getN = function(x) sum(getUniques(x))
track = cbind(out.all, rowSums(seqtab.all), rowSums(seqtab.chi))
colnames(track) = c("input", "filtered", "denoised", "no chimeras")
head(track)
saveRDS(track, "track.rds")
write.csv(track, "track.csv")

## assignTaxonomy - UNITE - without custom database - 
# General fasta (Includes singletons set as RefS)
unite.taxa1 = assignTaxonomy(seqtab.chi, "sh_general_release_dynamic_04.02.2020.fasta")
saveRDS(unite.taxa1, "/fulldataset/unite.taxa1.rds")
## assignTaxonomy - UNITE - including long read custom database - 
unite.taxa1_LR = assignTaxonomy(seqtab.chi, "sh_general_release_dynamic_04.02.2020_LR.fasta")
saveRDS(unite.taxa1_LR, "/fulldataset/unite.taxa1_LR.rds")

#view the taxonomic assignments
unite.taxa.print1 = unite.taxa1
rownames(unite.taxa.print1) = NULL
head(unite.taxa.print1)
unite.taxa.print1_LR = unite.taxa1_LR
rownames(unite.taxa.print1_LR) = NULL
head(unite.taxa.print1_LR)

#create fasta file
seqnum <- paste0("ASV", seq(ncol(seqtab.chi)))
uniqueSeqs <- as.list(colnames(seqtab.chi))
write.fasta(uniqueSeqs, seqnum, "/fulldataset/WP1_ITSfulldataset.fasta")

##Construct Phyloseq Object 
unite.taxa1 <- readRDS("/fulldataset/unite.taxa1.rds")
unite.taxa1_LR <- readRDS("/fulldataset/unite.taxa1_LR.rds")
seqtab.chi <- readRDS("/fulldataset/seqtab.chi.rds")

## Upload metadata ##
sd.DFG.ITS1 = read.csv("/fulldataset/mapping_DFG1_ITS2.csv")
head(sd.DFG.ITS1)
###Rename samples for metadata
sd.DFG.ITS <- sd.DFG.ITS1[, -1]
row.names(sd.DFG.ITS)<- sd.DFG.ITS1$Sample_ID #sample data row names must align with dada2 rowname outputs
sd.DFG.ITS = as.data.frame(sd.DFG.ITS)
head(sd.DFG.ITS)

#unite.taxa1
DFG.ITS.ex1 = phyloseq(tax_table(unite.taxa1), otu_table(seqtab.chi, taxa_are_rows = FALSE), sample_data(sd.DFG.ITS))
DFG.ITS.ex1
###Rename sequence variants
a.vec = as.vector(1:4020)  #number should reflect your total ASVs
a.nam = cbind("ASV", a.vec)
a.nam = as.data.frame(a.nam)
asv.names = paste0(a.nam$V1, a.nam$a.vec)
asv.names = as.data.frame(asv.names)
head(asv.names)
# apply ASV names to sequence table
taxa_names(DFG.ITS.ex1) = asv.names$asv.names
taxa_names(DFG.ITS.ex1)
DFG.ITS.ex1
saveRDS(DFG.ITS.ex1, '/fulldataset/DFG.ITS.ex1.rds')
# write files for full dataset, no taxonomic filtering or sample filtering (includes PCR negative control and field controls etc)
write.csv(otu_table(DFG.ITS.ex1), file = "/fulldataset/asv_table_WP1_ITSfulldataset1.csv")
write.csv(tax_table(DFG.ITS.ex1), file = "/fulldataset/tax_table_WP1_ITSfulldataset1.csv")

#unite.taxa1_LR with long reads
DFG.ITS.ex1_LR = phyloseq(tax_table(unite.taxa1_LR), otu_table(seqtab.chi, taxa_are_rows = FALSE), sample_data(sd.DFG.ITS))
DFG.ITS.ex1_LR
###Rename sequence variants
a.vec = as.vector(1:4020)  #number should reflect your total ASVs
a.nam = cbind("ASV", a.vec)
a.nam = as.data.frame(a.nam)
asv.names = paste0(a.nam$V1, a.nam$a.vec)
asv.names = as.data.frame(asv.names)
head(asv.names)
# apply ASV names to sequence table
taxa_names(DFG.ITS.ex1_LR) = asv.names$asv.names
taxa_names(DFG.ITS.ex1_LR)
DFG.ITS.ex1_LR
saveRDS(DFG.ITS.ex1_LR, 'fulldataset/DFG.ITS.ex1_LR.rds')
# write files for full dataset, no taxonomic filtering or sample filtering (includes PCR negative control and field controls etc)
write.csv(otu_table(DFG.ITS.ex1_LR), file = "/fulldataset/asv_table_WP1_ITSfulldataset1_LR.csv")
write.csv(tax_table(DFG.ITS.ex1_LR), file = "/fulldataset/tax_table_WP1_ITSfulldataset1_LR.csv")
write.csv(unite.taxa1_LR, file = "/fulldataset/unite.taxa1_LR.csv")

#Include rescued taxonomy from phylogenetic trees
unite.taxa1_LR_rescued1 = read.csv("/fulldataset/unite.taxa1_LR_rescued.csv")
unite.taxa1_LR_rescued <- unite.taxa1_LR_rescued1[, -1]
row.names(unite.taxa1_LR_rescued)<- unite.taxa1_LR_rescued1$X #sample data row names must align with dada2 rowname outputs
unite.taxa1_LR_rescued <- as.matrix(unite.taxa1_LR_rescued)

DFG.ITS.ex1_LR_res = phyloseq(tax_table(unite.taxa1_LR_rescued), otu_table(seqtab.chi, taxa_are_rows = FALSE), sample_data(sd.DFG.ITS))
DFG.ITS.ex1_LR_res


###Rename sequence variants
a.vec = as.vector(1:4020)  #number should reflect your total ASVs
a.nam = cbind("ASV", a.vec)
a.nam = as.data.frame(a.nam)
asv.names = paste0(a.nam$V1, a.nam$a.vec)
asv.names = as.data.frame(asv.names)
head(asv.names)
# apply ASV names to sequence table
taxa_names(DFG.ITS.ex1_LR_res) = asv.names$asv.names
taxa_names(DFG.ITS.ex1_LR_res)
DFG.ITS.ex1_LR_res
saveRDS(DFG.ITS.ex1_LR_res, '/fulldataset/DFG.ITS.ex1_LR_res.rds')
# write files for full dataset, no taxonomic filtering or sample filtering (includes PCR negative control and field controls etc)
write.csv(otu_table(DFG.ITS.ex1_LR_res), file = "/fulldataset/asv_table_WP1_ITSfulldataset1_LR_res.csv")
write.csv(tax_table(DFG.ITS.ex1_LR_res), file = "/fulldataset/tax_table_WP1_ITSfulldataset1_LR_res.csv")

DFG.ITS.ex.full<- DFG.ITS.ex1_LR_res
DFG.ITS.ex.full
#remove PCR negative control
DFG.ITS.ex.full.neg = subset_samples(DFG.ITS.ex.full, fungicide_treatment == "Negative")
DFG.ITS.ex.full.neg 
DFG.ITS.ex.full.neg <- filter_taxa(DFG.ITS.ex.full.neg, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
DFG.ITS.ex.full.neg
write.csv(otu_table(DFG.ITS.ex.full.neg), file = "/taxfilt/neg_asv_table_DFG.ITS.ex.full.neg.csv")
write.csv(tax_table(DFG.ITS.ex.full.neg), file = "/taxfilt/neg_tax_table_DFG.ITS.ex.full.neg.csv")
DFG.ITS.ex.full.noneg = subset_samples(DFG.ITS.ex.full, fungicide_treatment != "Negative")
DFG.ITS.ex.full.noneg 
###taxonomic filtering
rank_names(DFG.ITS.ex.full.noneg)
#create a table of read counts for each Phylum present in the dataset
table(tax_table(DFG.ITS.ex.full.noneg)[, "Phylum"], exclude = NULL)
#remove taxa for which Phylum is NA
#DFG.ITS.ex.full.noneg0 <- subset_taxa(DFG.ITS.ex.full.noneg, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
#table(tax_table(DFG.ITS.ex.full.noneg0)[, "Phylum"], exclude = NULL)
DFG.ITS.ex.full.noneg0 <-DFG.ITS.ex.full.noneg
write.csv(otu_table(DFG.ITS.ex.full.noneg0), file = "/taxfilt/asv_table_DFG.ITS.ex.full.noneg0.csv")
write.csv(tax_table(DFG.ITS.ex.full.noneg0), file = "/taxfilt/tax_table_DFG.ITS.ex.full.noneg0.csv")
write.csv(sample_data(DFG.ITS.ex.full.noneg0), file = "/taxfilt/sample_data_DFG.ITS.ex.full.noneg0.csv")
saveRDS(DFG.ITS.ex.full.noneg0, "/taxfilt/DFG.ITS.ex.full.noneg0.rds")

#for AH only - remove asvs for which genus is NA
table(tax_table(DFG.ITS.ex.full.noneg0)[, "Genus"], exclude = NULL)
DFG.ITS.ex.full.noneg0g <- subset_taxa(DFG.ITS.ex.full.noneg0, !is.na(Genus) & !Genus %in% c("", "uncharacterized"))
table(tax_table(DFG.ITS.ex.full.noneg0g)[, "Genus"], exclude = NULL)
DFG.ITS.ex.full.noneg0g
write.csv(otu_table(DFG.ITS.ex.full.noneg0g), file = "/taxfilt/asv_table_DFG.ITS.ex.full.noneg0g.csv")
write.csv(tax_table(DFG.ITS.ex.full.noneg0g), file = "/taxfilt/tax_table_DFG.ITS.ex.full.noneg0g.csv")
write.csv(sample_data(DFG.ITS.ex.full.noneg0g), file = "/taxfilt/sample_data_DFG.ITS.ex.full.noneg0g.csv")
saveRDS(DFG.ITS.ex.full.noneg0g, "/taxfilt/DFG.ITS.ex.full.noneg0g.rds")

#subset for paper 1: fungicide (TU-1) negative control, without shredders and remove field controls
ITS.nonegF.field  = subset_samples(DFG.ITS.ex.full.noneg0, fungicide_treatment == "field")
ITS.nonegF.field  <- filter_taxa(ITS.nonegF.field, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
ITS.nonegF.field
write.csv(otu_table(ITS.nonegF.field), file = "/taxfilt/asv_fieldsamples_ITS.nonegF.field.csv")
write.csv(tax_table(ITS.nonegF.field), file = "/taxfilt/tax_table_fieldsamples_ITS.nonegF.field.csv")
ITS.nofield  = subset_samples(DFG.ITS.ex.full.noneg0, fungicide_treatment != "field")
ITS.nofield <- filter_taxa(ITS.nofield, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
ITS.nofield
ITS.Fto  = subset_samples(ITS.nofield, shredder == "w/o")
ITS.Fto
ITS.Fto <- filter_taxa(ITS.Fto, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
ITS.Fto
ITS.Fto_CT1T2  = subset_samples(ITS.Fto, fungicide_treatment != "TU-3")
ITS.Fto_CT1T2
ITS.Fto_CT1  = subset_samples(ITS.Fto_CT1T2, fungicide_treatment != "TU-2")
ITS.Fto_CT1  <- filter_taxa(ITS.Fto_CT1, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
ITS.Fto_CT1
write.csv(otu_table(ITS.Fto_CT1), file = "/taxfilt/asv_table_ITS.Fto_CT1.csv")
write.csv(tax_table(ITS.Fto_CT1), file = "/taxfilt/tax_table_ITS.Fto_CT1.csv")
write.csv(sample_data(ITS.Fto_CT1), file = "/taxfilt/sample_data_ITS.Fto_CT1.csv")
saveRDS(ITS.Fto_CT1, "/taxfilt/ITS.Fto_CT1.rds")

#summarise to species level
ITS.Fto_CT1 <- readRDS("/taxfilt/ITS.Fto_CT1.rds")
tax.clean <- data.frame(tax_table(ITS.Fto_CT1))
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}
tax_table(ITS.Fto_CT1) <- as.matrix(tax.clean)
tax_table(ITS.Fto_CT1)
ITS.Fto_CT1_sp <- tax_glom(ITS.Fto_CT1, taxrank=rank_names(ITS.Fto_CT1)[7])
ITS.Fto_CT1_sp
write.csv(otu_table(ITS.Fto_CT1_sp), file = "/taxfilt/asv_table_ITS.Fto_CT1_sp.csv")
write.csv(tax_table(ITS.Fto_CT1_sp), file = "/taxfilt/tax_table_ITS.Fto_CT1_sp.csv")
write.csv(sample_data(ITS.Fto_CT1_sp), file = "/taxfilt/sample_data_ITS.Fto_CT1_sp.csv")
saveRDS(ITS.Fto_CT1_sp, "/taxfilt/ITS.Fto_CT1_sp.rds")

#normalise 
ITS.Fto_CT1_deseq <- phyloseq_to_deseq2(ITS.Fto_CT1, ~fungicide_treatment)
deseq_counts_vst <- varianceStabilizingTransformation(ITS.Fto_CT1_deseq, blind = TRUE)
vst_trans_count_tab <- assay(deseq_counts_vst)# extract transformed table
#sample data
#calculate Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")
plot(euc_clust)
euc_dend <- as.dendrogram(euc_clust, hang=0.1)
plot(euc_dend, ylab="VST Euc. dist.")
#construct transformed phyloseq
tax_table_ITS.Fto_CT1_1 = read.csv("/taxfilt/tax_table_ITS.Fto_CT1.csv")
tax_table_ITS.Fto_CT1 <- tax_table_ITS.Fto_CT1_1[, -1]
row.names(tax_table_ITS.Fto_CT1)<- tax_table_ITS.Fto_CT1_1$X #sample data row names must align with dada2 rowname outputs
tax_table_ITS.Fto_CT1 <- as.matrix(tax_table_ITS.Fto_CT1)
ITS.Fto_CT1_vst <- phyloseq(tax_table(tax_table_ITS.Fto_CT1), otu_table(vst_trans_count_tab, taxa_are_rows=T), sample_data(sd.DFG.ITS))#transformed
ITS.Fto_CT1_vst
df_ITS.Fto_CT1_vst <- otu_table(ITS.Fto_CT1_vst)
df_ITS.Fto_CT1_vst[df_ITS.Fto_CT1_vst < 0] <- 0 # Set negative values to 0
df_ITS.Fto_CT1_vst   
write.csv(df_ITS.Fto_CT1_vst, file = "/taxfilt/asv_table_ITS.Fto_CT1_vst.csv")
write.csv(tax_table(ITS.Fto_CT1_vst), file = "/taxfilt/tax_table_ITS.Fto_CT1_vst.csv")
write.csv(sample_data(ITS.Fto_CT1_vst), file = "/taxfilt/sample_data_ITS.Fto_CT1_vst.csv")
saveRDS(ITS.Fto_CT1_vst, "/taxfilt/ITS.Fto_CT1_vst.rds")


###calculate number of ASVs at each level for UNITE, UNITE+LR and UNITE+LR+res
#unite+LR+res
tax_table_ITS.Fto_CT1_unite_LR_res <- read.csv("/taxfilt/tax_table_ITS.Fto_CT1.csv")
head(tax_table_ITS.Fto_CT1_unite_LR_res)
count(tax_table_ITS.Fto_CT1_unite_LR_res, vars = "Kingdom")
count(tax_table_ITS.Fto_CT1_unite_LR_res, vars = "Phylum")
count(tax_table_ITS.Fto_CT1_unite_LR_res, vars = "Class")
count(tax_table_ITS.Fto_CT1_unite_LR_res, vars = "Order")
count(tax_table_ITS.Fto_CT1_unite_LR_res, vars = "Family")
count(tax_table_ITS.Fto_CT1_unite_LR_res, vars = "Genus")
count(tax_table_ITS.Fto_CT1_unite_LR_res, vars = "Species")

#unite
DFG.ITS.ex1 <- readRDS("/fulldataset/DFG.ITS.ex1.rds")
DFG.ITS.ex.full_unite <- DFG.ITS.ex1
DFG.ITS.ex.full.neg_unite = subset_samples(DFG.ITS.ex.full_unite, fungicide_treatment == "Negative")
DFG.ITS.ex.full.neg_unite <- filter_taxa(DFG.ITS.ex.full.neg_unite, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
DFG.ITS.ex.full.noneg_unite = subset_samples(DFG.ITS.ex.full_unite, fungicide_treatment != "Negative")
rank_names(DFG.ITS.ex.full.noneg_unite)

table(tax_table(DFG.ITS.ex.full.noneg_unite)[, "Phylum"], exclude = NULL)
ITS.nonegF.field_unite  = subset_samples(DFG.ITS.ex.full.noneg_unite, fungicide_treatment == "field")
ITS.nonegF.field_unite  <- filter_taxa(ITS.nonegF.field_unite, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
ITS.nofield_unite  = subset_samples(DFG.ITS.ex.full.noneg_unite, fungicide_treatment != "field")
ITS.nofield_unite <- filter_taxa(ITS.nofield_unite, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
ITS.Fto_unite  = subset_samples(ITS.nofield_unite, shredder == "w/o")
ITS.Fto_unite <- filter_taxa(ITS.Fto_unite, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
ITS.Fto_CT1T2_unite  = subset_samples(ITS.Fto_unite, fungicide_treatment != "TU-3")
ITS.Fto_CT1_unite  = subset_samples(ITS.Fto_CT1T2_unite, fungicide_treatment != "TU-2")
ITS.Fto_CT1_unite  <- filter_taxa(ITS.Fto_CT1_unite, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
ITS.Fto_CT1_unite
write.csv(otu_table(ITS.Fto_CT1_unite), file = "/taxfilt/asv_table_ITS.Fto_CT1_unite.csv")
write.csv(tax_table(ITS.Fto_CT1_unite), file = "/taxfilt/tax_table_ITS.Fto_CT1_unite.csv")
write.csv(sample_data(ITS.Fto_CT1_unite), file = "/taxfilt/sample_data_ITS.Fto_CT1_unite.csv")
saveRDS(ITS.Fto_CT1_unite, "/taxfilt/ITS.Fto_CT1_unite.rds")
tax_table_ITS.Fto_CT1_unite <- read.csv("/taxfilt/tax_table_ITS.Fto_CT1_unite.csv")
head(tax_table_ITS.Fto_CT1_unite)
table(tax_table(ITS.Fto_CT1_unite)[, "Kingdom"], exclude = NULL)
count(tax_table_ITS.Fto_CT1_unite, vars = "Kingdom")
count(tax_table_ITS.Fto_CT1_unite, vars = "Phylum")
count(tax_table_ITS.Fto_CT1_unite, vars = "Class")
count(tax_table_ITS.Fto_CT1_unite, vars = "Order")
count(tax_table_ITS.Fto_CT1_unite, vars = "Family")
count(tax_table_ITS.Fto_CT1_unite, vars = "Genus")
count(tax_table_ITS.Fto_CT1_unite, vars = "Species")

##unite+LR
DFG.ITS.ex1_LR <- readRDS("/fulldataset/DFG.ITS.ex1_LR.rds")
DFG.ITS.ex.full_unite_LR <- DFG.ITS.ex1_LR
DFG.ITS.ex.full.neg_unite_LR = subset_samples(DFG.ITS.ex.full_unite_LR, fungicide_treatment == "Negative")
DFG.ITS.ex.full.neg_unite_LR <- filter_taxa(DFG.ITS.ex.full.neg_unite_LR, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
DFG.ITS.ex.full.noneg_unite_LR = subset_samples(DFG.ITS.ex.full_unite_LR, fungicide_treatment != "Negative")
rank_names(DFG.ITS.ex.full.noneg_unite_LR)
table(tax_table(DFG.ITS.ex.full.noneg_unite_LR)[, "Phylum"], exclude = NULL)
#DFG.ITS.ex.full.noneg0_unite_LR <- subset_taxa(DFG.ITS.ex.full.noneg_unite_LR, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
#table(tax_table(DFG.ITS.ex.full.noneg0_unite_LR)[, "Phylum"], exclude = NULL)
DFG.ITS.ex.full.noneg0_unite_LR <- subset_taxa(DFG.ITS.ex.full.noneg_unite_LR, (Phylum!="p__unidentified") | is.na(Phylum))
ITS.nonegF.field_unite_LR  = subset_samples(DFG.ITS.ex.full.noneg_unite_LR, fungicide_treatment == "field")
ITS.nonegF.field_unite_LR  <- filter_taxa(ITS.nonegF.field_unite_LR, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
ITS.nofield_unite_LR  = subset_samples(DFG.ITS.ex.full.noneg_unite_LR, fungicide_treatment != "field")
ITS.nofield_unite_LR <- filter_taxa(ITS.nofield_unite_LR, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
ITS.Fto_unite_LR  = subset_samples(ITS.nofield_unite_LR, shredder == "w/o")
ITS.Fto_unite_LR <- filter_taxa(ITS.Fto_unite_LR, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
ITS.Fto_CT1T2_unite_LR  = subset_samples(ITS.Fto_unite_LR, fungicide_treatment != "TU-3")
ITS.Fto_CT1_unite_LR  = subset_samples(ITS.Fto_CT1T2_unite_LR, fungicide_treatment != "TU-2")
ITS.Fto_CT1_unite_LR  <- filter_taxa(ITS.Fto_CT1_unite_LR, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
ITS.Fto_CT1_unite_LR
write.csv(otu_table(ITS.Fto_CT1_unite_LR), file = "/taxfilt/asv_table_ITS.Fto_CT1_unite_LR.csv")
write.csv(tax_table(ITS.Fto_CT1_unite_LR), file = "/taxfilt/tax_table_ITS.Fto_CT1_unite_LR.csv")
write.csv(sample_data(ITS.Fto_CT1_unite_LR), file = "/taxfilt/sample_data_ITS.Fto_CT1_unite_LR.csv")

saveRDS(ITS.Fto_CT1_unite_LR, "/taxfilt/ITS.Fto_CT1_unite_LR.rds")
tax_table_ITS.Fto_CT1_unite_LR <- read.csv("/taxfilt/tax_table_ITS.Fto_CT1_unite_LR.csv")
head(tax_table_ITS.Fto_CT1_unite_LR)
count(tax_table_ITS.Fto_CT1_unite_LR, vars = "Kingdom")
count(tax_table_ITS.Fto_CT1_unite_LR, vars = "Phylum")
count(tax_table_ITS.Fto_CT1_unite_LR, vars = "Class")
count(tax_table_ITS.Fto_CT1_unite_LR, vars = "Order")
count(tax_table_ITS.Fto_CT1_unite_LR, vars = "Family")
count(tax_table_ITS.Fto_CT1_unite_LR, vars = "Genus")
count(tax_table_ITS.Fto_CT1_unite_LR, vars = "Species")

##Subset aquatic hyphomyetes
tax_table_DFG.ITS.ex.full.noneg0_AH1 = read.csv("/taxfilt/tax_table_DFG.ITS.ex.full.noneg0g_AH.csv")
tax_table_DFG.ITS.ex.full.noneg0_AH <- tax_table_DFG.ITS.ex.full.noneg0_AH1[, -1]
row.names(tax_table_DFG.ITS.ex.full.noneg0_AH)<- tax_table_DFG.ITS.ex.full.noneg0_AH1$X #sample data row names must align with dada2 rowname outputs
tax_table_DFG.ITS.ex.full.noneg0_AH <- as.matrix(tax_table_DFG.ITS.ex.full.noneg0_AH)
DFG.ITS.ex.full.noneg0_AH = phyloseq(tax_table(tax_table_DFG.ITS.ex.full.noneg0_AH), otu_table(seqtab.chi, taxa_are_rows = FALSE), sample_data(sd.DFG.ITS))
DFG.ITS.ex.full.noneg0_AH
saveRDS(DFG.ITS.ex.full.noneg0_AH, '/taxfilt/DFG.ITS.ex.full.noneg0_AH.rds')
#subset AH for paper 1
ITS.nofield_AH  = subset_samples(DFG.ITS.ex.full.noneg0_AH, fungicide_treatment != "field")
ITS.nofield_AH <- filter_taxa(ITS.nofield_AH, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
saveRDS(ITS.nofield_AH, "taxfilt/ITS.nofield_AH.rds")
ITS.Fto_AH  = subset_samples(ITS.nofield_AH, shredder == "w/o")
ITS.Fto_AH <- filter_taxa(ITS.Fto_AH, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
ITS.Fto_CT1T2_AH  = subset_samples(ITS.Fto_AH, fungicide_treatment != "TU-3")
ITS.Fto_CT1_AH  = subset_samples(ITS.Fto_CT1T2_AH, fungicide_treatment != "TU-2")
ITS.Fto_CT1_AH  <- filter_taxa(ITS.Fto_CT1_AH, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
ITS.Fto_CT1_AH
write.csv(otu_table(ITS.Fto_CT1_AH), file = "/taxfilt/asv_table_ITS.Fto_CT1_AH.csv")
write.csv(tax_table(ITS.Fto_CT1_AH), file = "/taxfilt/tax_table_ITS.Fto_CT1_AH.csv")
write.csv(sample_data(ITS.Fto_CT1_AH), file = "/taxfilt/sample_data_ITS.Fto_CT1_AH.csv")
saveRDS(ITS.Fto_CT1_AH, "/taxfilt/ITS.Fto_CT1_AH.rds")

#Summarise aquatic hyphomyetes asvs to species level
tax.clean <- data.frame(tax_table(ITS.Fto_CT1_AH))
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}

tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

tax_table(ITS.Fto_CT1_AH) <- as.matrix(tax.clean)
tax_table(ITS.Fto_CT1_AH)
ITS.Fto_CT1_AHsp <- tax_glom(ITS.Fto_CT1_AH, taxrank=rank_names(ITS.Fto_CT1_AH)[7])
ITS.Fto_CT1_AHsp
write.csv(otu_table(ITS.Fto_CT1_AHsp), file = "/taxfilt/asv_table_ITS.Fto_CT1_AHsp.csv")
write.csv(tax_table(ITS.Fto_CT1_AHsp), file = "/taxfilt/tax_table_ITS.Fto_CT1_AHsp.csv")
write.csv(sample_data(ITS.Fto_CT1_AHsp), file = "/taxfilt/sample_data_ITS.Fto_CT1_AHsp.csv")
saveRDS(ITS.Fto_CT1_AHsp, "/taxfilt/ITS.Fto_CT1_AHsp.rds")


# Transform data to realtive abundance - proportions 
ITS.Fto_CT1.ra <- transform_sample_counts(ITS.Fto_CT1, function(otu) otu/sum(otu))
write.csv(otu_table(ITS.Fto_CT1.ra), file = "taxfilt/asv_table_ITS.Fto_CT1.ra.csv")
write.csv(tax_table(ITS.Fto_CT1.ra), file = "taxfilt/tax_table_ITS.Fto_CT1.ra.csv")
write.csv(sample_data(ITS.Fto_CT1.ra), file = "taxfilt/sample_data_ITS.Fto_CT1.ra.csv")
saveRDS(ITS.Fto_CT1.ra, "taxfilt/ITS.Fto_CT1.ra.rds")
ITS.Fto_CT1_sp.ra <- transform_sample_counts(ITS.Fto_CT1_sp, function(otu) otu/sum(otu))
write.csv(otu_table(ITS.Fto_CT1_sp.ra), file = "taxfilt/asv_table_ITS.Fto_CT1_sp.ra.csv")
write.csv(tax_table(ITS.Fto_CT1_sp.ra), file = "taxfilt/tax_table_ITS.Fto_CT1_sp.ra.csv")
write.csv(sample_data(ITS.Fto_CT1_sp.ra), file = "taxfilt/sample_data_ITS.Fto_CT1_sp.ra.csv")
saveRDS(ITS.Fto_CT1_sp.ra, "taxfilt/ITS.Fto_CT1_sp.ra.rds")
ITS.Fto_CT1_AH.ra <- transform_sample_counts(ITS.Fto_CT1_AH, function(otu) otu/sum(otu))
write.csv(otu_table(ITS.Fto_CT1_AH.ra), file = "taxfilt/asv_table_ITS.Fto_CT1_AH.ra.csv")
write.csv(tax_table(ITS.Fto_CT1_AH.ra), file = "taxfilt/tax_table_ITS.Fto_CT1_AH.ra.csv")
write.csv(sample_data(ITS.Fto_CT1_AH.ra), file = "taxfilt/sample_data_ITS.Fto_CT1_AH.ra.csv")
saveRDS(ITS.Fto_CT1_AH.ra, "taxfilt/ITS.Fto_CT1_AH.ra.rds")
ITS.Fto_CT1_AHsp.ra <- transform_sample_counts(ITS.Fto_CT1_AHsp, function(otu) otu/sum(otu))
write.csv(otu_table(ITS.Fto_CT1_AHsp.ra), file = "taxfilt/asv_table_ITS.Fto_CT1_AHsp.ra.csv")
write.csv(tax_table(ITS.Fto_CT1_AHsp.ra), file = "taxfilt/tax_table_ITS.Fto_CT1_AHsp.ra.csv")
write.csv(sample_data(ITS.Fto_CT1_AHsp.ra), file = "taxfilt/sample_data_ITS.Fto_CT1_AHsp.ra.csv")
saveRDS(ITS.Fto_CT1_AHsp.ra, "taxfilt/ITS.Fto_CT1_AHsp.ra.rds")