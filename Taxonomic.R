#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="Path where R1 and R2 files available", metavar="character"),
  make_option(c("-f", "--forward"), type="character", default=NULL, 
              help="Prefix of R1 read that is common in all Samples", metavar="character"),
  make_option(c("-r", "--reverse"), type="character", default=NULL, 
              help="Prefix of R2 read that is common in all Samples", metavar="character"),
  make_option(c("-s", "--silva"), type="character", default=NULL, 
              help="Path to the Silva databse for Taxonomy assignment", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$path, opt$forward, opt$reverse, opt$silva)){
  print_help(opt_parser)
  stop("All four argument must be supplied (input file).n", call.=FALSE)
}
library(dada2)
path <- opt$path
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern=opt$forward, full.names = TRUE))
fnRs <- sort(list.files(path, pattern=opt$reverse, full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
#errF <- learnErrors(filtFs, multithread=TRUE)
#errR <- learnErrors(filtRs, multithread=TRUE)
#plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, selfConsist = TRUE)
dadaRs <- dada(filtRs, selfConsist = TRUE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, justConcatenate = TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
taxa <- assignTaxonomy(seqtab.nochim, opt$silva, multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa, "taxonomic Classification.csv")