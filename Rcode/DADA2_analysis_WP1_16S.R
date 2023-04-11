## Bioconductor Pipeline ## DFG1 - 16S

## Install DADA2 ##
## The first stage is to install DADA2 if it is not already available on your computer.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
BiocManager::install("phyloseq")
BiocManager::install("coseq")
library(dada2)
library(phyloseq)
library(ggplot2)
library(coseq)


setwd("/WP1_16S/output")

## Read files ##
path = "/WP1_16S/raw"
# check the file path is correct using list.files(path).
list.files(path)

# we take this opportunity to extract the file names as a
# vector, which will be used later in the script
f.names = as.vector(list.files(path, pattern = "_1.fastq.bz2", 
                               full.names = F))
r.names = as.vector(list.files(path, pattern = "_2.fastq.bz2", 
                               full.names = F))

# pattern identifies all files in the path ending in
# '_1.fastq'. This orders and separates forward and
# reverse reads.
fnFs = sort(list.files(path, pattern = "_1.fastq.bz2", full.names = TRUE))
fnRs = sort(list.files(path, pattern = "_2.fastq.bz2", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_lib"), `[`, 1)

## Plot quality scores from .fastq files ##
# illustration of plot outputs (first 2 only)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

## Filter and trim (quality control) ##
filtF_path = file.path(path, "filtF")  #place filtered forward files in 'filtered' subdirectory
filtR_path = file.path(path, "filtR")  #place filtered reverse files in 'filtered' subdirectory
filtFs = file.path(filtF_path, paste0(f.names, "_F_filt.fastq.gz"))
filtRs = file.path(filtR_path, paste0(r.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(240, 
                                                              220), trimLeft = c(19, 20), maxN = 0, maxEE = c(2, 2), truncQ = 2, 
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE)

head(out)  #check the results

###dada2 dereplicate and merge paired-end reads ##for big data - by processing samples independently
# File parsing
filtpathF <- "/WP1_16S/output/filtF" #directory containing filtered forward fastqs
filtpathR <- "/WP1_16S/output/filtR" ##directory containing filtered reverse fastqs
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_lib"), `[`, 1) # 
sample.namesR <- sapply(strsplit(basename(filtRs), "_lib"), `[`, 1) # 
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)

## Generate sequence table ##
## Construct a sequence table from the merged sequence object. Viewing this table will reveal a number of amplicon lengths which are outside the expected range of your region, these are likely to be erroneous
seqtab = makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

## Plot amplicon length distribution ##
## Plot amplicon length frequency to visualise the distribution of amplicon sizes, informing the amplicon length range to conserve versus the length range to discard.
ex = as.data.frame(table(nchar(getSequences(seqtab))))
ex.plot = ggplot(ex, aes(Var1, Freq)) + geom_bar(stat = "identity", 
                                                 aes(fill = Var1)) + ylab("Frequency") + xlab("Merged Sequence Length") + 
  scale_x_discrete(breaks = c(220, 240, 250, 260, 270, 300, 
                              330, 360)) + theme_bw() + theme(panel.grid.major = element_blank(), 
                                                              panel.grid.minor = element_blank()) + theme(legend.position = "none")
ex.plot

## Excise amplicons of desired length ##
## In this case, only amplicons between 250 and 256bp in length are selected. Amplicons outside of this range are discarded.
seqtab.ex = seqtab[, nchar(colnames(seqtab)) %in% seq(250, 256)]
table(nchar(getSequences(seqtab.ex)))

## Remove chimeras ##
## Note that if you lose a large number of sequence reads at this stage, it is likely that primer sequences still remain in your reads. These should have been removed manually at the trim and truncate stage.
seqtab.ex.chi = removeBimeraDenovo(seqtab.ex, method = "consensus", 
                                   multithread = T, verbose = T)
sum(seqtab.ex.chi)/sum(seqtab.ex)
saveRDS(seqtab.ex.chi, "/WP1_16S/output/seqtab.ex.chi.rds")

## Track sequence loss ##
## Track the number of sequence reads per sample remaining after each stage of processing. This provides an opportunity to identify points in the pipeline where large proportions of sequences were lost, which may be worth revisiting and tweaking. Commonly, 15-20% of reads are lost over the course of the pipeline.
getN = function(x) sum(getUniques(x))
track = cbind(out, sapply(mergers, getN), 
              rowSums(seqtab.ex), rowSums(seqtab.ex.chi))
colnames(track) = c("input", "filtered", "merged", 
                    "tabled", "no chim")
rownames(track) = f.names
head(track)
saveRDS(track, "/WP1_16S/output/track.rds")
write.csv(track, "/WP1_16S/output/track.csv")


## Assign taxonomy ##
silva.taxa = assignTaxonomy(seqtab.ex.chi, "/WP1_16S/output/silva_nr_v132_train_set.fa.gz", 
                            multithread = T)
saveRDS(silva.taxa, "/WP1_16S/output/silva.taxa.rds")
#species level assignments
silva.taxa.sp <- addSpecies(silva.taxa, "/WP1_16S/output/silva_species_assignment_v132.fa.gz") 
saveRDS(silva.taxa.sp, "/WP1_16S/output/silva.taxa.sp.rds")

#view the taxonomic assignments
silva.taxa.print = silva.taxa # Removing sequence rownames for display only
rownames(silva.taxa.print) = NULL
head(silva.taxa.print)

#create fasta file
library(seqinr)
seqnum <- paste0("ASV", seq(ncol(seqtab.ex.chi)))
uniqueSeqs <- as.list(colnames(seqtab.ex.chi))
write.fasta(uniqueSeqs, seqnum, "/WP1_16S/output/WP1_16Sfulldataset.fasta")

#create phyloseq object
seqtab.chi.16S <- readRDS("/WP1_16S/output/seqtab.ex.chi.rds")
silva.taxa.16S <- readRDS("/WP1_16S/output/silva.taxa.sp.rds")
B16S.ex = phyloseq(tax_table(silva.taxa.16S), otu_table(seqtab.chi.16S, taxa_are_rows = FALSE))
B16S.ex
rownames(otu_table(B16S.ex))
rownames(tax_table(B16S.ex))
## Upload metadata ##
sd.B16S1 = read.csv("/WP1_16S/output/mapping_DFG1_16S.csv")
head(sd.B16S1)
###Rename samples for metadata
sd.B16S <- sd.B16S1[, -1]
row.names(sd.B16S)<- sd.B16S1$Sample #sample data row names must align with dada2 rowname outputs
sd.B16S = as.data.frame(sd.B16S)
head(sd.B16S)
B16S.ex = phyloseq(tax_table(silva.taxa.16S), otu_table(seqtab.chi.16S, taxa_are_rows = FALSE), sample_data(sd.B16S))
B16S.ex
###Rename sequence variants
a.vec = as.vector(1:14383)  #number should reflect your total ASVs
a.nam = cbind("ASV", a.vec)
a.nam = as.data.frame(a.nam)
asv.names = paste0(a.nam$V1, a.nam$a.vec)
asv.names = as.data.frame(asv.names)
head(asv.names)
# apply ASV names to sequence table
taxa_names(B16S.ex) = asv.names$asv.names
taxa_names(B16S.ex)
B16S.ex

# write files for full dataset, no taxonomic filtering or sample filtering (includes PCR negative control and field controls etc)
write.csv(otu_table(B16S.ex), file = "/WP1_16S/output/asv_table_WP1_16Sfulldataset.csv")
write.csv(tax_table(B16S.ex), file = "/WP1_16S/output/tax_table_WP1_16Sfulldataset.csv")


#### IMPORT tree
#B16S.ex.full <- merge_phyloseq(B16S.ex, nwk.tree)
B16S.ex.full <- B16S.ex
B16S.ex.full
#remove PCR negative control
B16S.ex.full.neg = subset_samples(B16S.ex.full, fungicide_treatment == "Negative")
B16S.ex.full.neg 
B16S.ex.full.neg <- filter_taxa(B16S.ex.full.neg, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
B16S.ex.full.neg 
write.csv(otu_table(B16S.ex.full.neg), file = "/WP1_16S/output/taxfilt/neg_asv_table_B16S.ex.full.neg.csv")
write.csv(tax_table(B16S.ex.full.neg), file = "/WP1_16S/output/taxfilt/neg_tax_table_DB16S.ex.full.neg.csv")
B16S.ex.full.noneg = subset_samples(B16S.ex.full, fungicide_treatment != "Negative")
B16S.ex.full.noneg <- filter_taxa(B16S.ex.full.noneg, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
B16S.ex.full.noneg
###taxonomic filtering
rank_names(B16S.ex.full.noneg)
table(tax_table(B16S.ex.full.noneg)[, "Kingdom"], exclude = NULL)
B16S.ex.full.noneg
#just filter bacteria
table(tax_table(B16S.ex.full.noneg)[, "Kingdom"], exclude = NULL)
B16S.ex.full.nonegB <- subset_taxa(B16S.ex.full.noneg, Kingdom %in% "Bacteria")
B16S.ex.full.nonegB
B16S.ex.full.noneg_other <- subset_taxa(B16S.ex.full.noneg, !Kingdom %in% "Bacteria")
B16S.ex.full.noneg_other
table(tax_table(B16S.ex.full.nonegB)[, "Kingdom"], exclude = NULL)
table(tax_table(B16S.ex.full.noneg_other)[, "Kingdom"], exclude = NULL)
#create a table of read counts for each Phylum present in the dataset
table(tax_table(B16S.ex.full.noneg)[, "Phylum"], exclude = NULL)
#remove taxa for which Phylum is NA
B16S.ex.full.noneg0 <- subset_taxa(B16S.ex.full.noneg, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
table(tax_table(B16S.ex.full.noneg0)[, "Phylum"], exclude = NULL)
B16S.ex.full.noneg0
#remove Chloroplasts (at Order level) and Mitochondria (Family)
B16S.ex.full.noneg0NC <- subset_taxa(B16S.ex.full.noneg0, (Order!="Chloroplast") | is.na(Order))
B16S.ex.full.noneg0NC
B16S.ex.full.noneg0NCM <- subset_taxa(B16S.ex.full.noneg0NC, (Family!="Mitochondria") | is.na(Family))
B16S.ex.full.noneg0NCM
#subset field controls -
B16S.ex.full.noneg0NCM.field  = subset_samples(B16S.ex.full.noneg0NCM, fungicide_treatment == "field")
B16S.ex.full.noneg0NCM.field 
B16S.ex.full.noneg0NCM.field  <- filter_taxa(B16S.ex.full.noneg0NCM.field, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
B16S.ex.full.noneg0NCM.field
write.csv(otu_table(B16S.ex.full.noneg0NCM.field), file = "/WP1_16S/output/taxfilt/asv_fieldsamples_16S.noneg0NCM.field.csv")
write.csv(tax_table(B16S.ex.full.noneg0NCM), file = "/WP1_16S/output/taxfilt/tax_table_fieldsamples_16S.noneg0NCM.field.csv")
B16S.nofield  = subset_samples(B16S.ex.full.noneg0NCM, fungicide_treatment != "field")
B16S.Fto  = subset_samples(B16S.nofield, shredder == "w/o")
B16S.Fto <- filter_taxa(B16S.Fto, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
B16S.Fto
write.csv(otu_table(B16S.Fto), file = "/WP1_16S/output/taxfilt/asv_table_B16S.Fto.csv")
write.csv(tax_table(B16S.Fto), file = "/WP1_16S/output/taxfilt/tax_table_B16S.Fto.csv")
write.csv(sample_data(B16S.Fto), file = "/WP1_16S/output/taxfilt/sample_data_B16S.Fto.csv")
saveRDS(B16S.Fto, "/WP1_16S/output/taxfilt/B16S.Fto.rds")
B16S.Fto_CT1T2  = subset_samples(B16S.Fto, fungicide_treatment != "TU-3")
B16S.Fto_CT1  = subset_samples(B16S.Fto_CT1T2, fungicide_treatment != "TU-2")
B16S.Fto_CT1  <- filter_taxa(B16S.Fto_CT1, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
B16S.Fto_CT1
saveRDS(B16S.Fto_CT1, "/WP1_16S/output/taxfilt/B16S.Fto_CT1.rds")
write.csv(otu_table(B16S.Fto_CT1), file = "/WP1_16S/output/taxfilt/asv_table_B16S.Fto_CT1.csv")
write.csv(tax_table(B16S.Fto_CT1), file = "/WP1_16S/output/taxfilt/tax_table_B16S.Fto_CT1.csv")
write.csv(sample_data(B16S.Fto_CT1), file = "/WP1_16S/output/taxfilt/sample_data_B16S.Fto_CT1.csv")
B16S.Fto_CT1 <- readRDS("taxfilt/B16S.Fto_CT1.rds")

# Transform data to proportions 
B16S.Fto_CT1.ra <- transform_sample_counts(B16S.Fto_CT1, function(otu) otu/sum(otu))
write.csv(otu_table(B16S.Fto_CT1.ra), file = "/WP1_16S/output/taxfilt/asv_table_B16S.Fto_CT1.ra.csv")
write.csv(tax_table(B16S.Fto_CT1.ra), file = "/WP1_16S/output/taxfilt/tax_table_B16S.Fto_CT1.ra.csv")
write.csv(sample_data(B16S.Fto_CT1.ra), file = "/WP1_16S/output/taxfilt/sample_data_B16S.Fto_CT1.ra.csv")
saveRDS(B16S.Fto_CT1.ra, "/WP1_16S/output/taxfilt/B16S.Fto_CT1.ra.rds")

