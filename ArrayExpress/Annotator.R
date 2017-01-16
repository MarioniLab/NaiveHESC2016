########################################################################################
# Print out important bits from the metadata.

collected <- list()
for (sample in c("2383", "2384", "2677", "2678", "2739", "2740")){
    all.files <- read.table(file.path("..", "counts", paste0("genic_counts_", sample, ".tsv.gz")), 
                            nrows=1, stringsAsFactor=FALSE)
    prefixes <- as.character(all.files[-c(1:2)])

    snames <- sub(".*_(3600STDY[0-9]+)_.*", "\\1", prefixes)
    if (sample %in% c("2383", "2384")) { batch <- 1 }
    else { batch <- 2 }
    if (sample %in% c("2383", "2677", "2678")) { phenotype <- "naive" }
    else { phenotype <- "primed" }

    prelim <- data.frame(Sample=snames, Run=sample, Batch=batch, Phenotype=phenotype)
    collected[[sample]] <- rbind(data.frame(prelim, File=paste0(prefixes, "_R1.fastq.gz")),
                                 data.frame(prelim, File=paste0(prefixes, "_R2.fastq.gz")))
}

collected <- do.call(rbind, collected)
table(table(collected$Sample)) # should all be 4; two technical replicates for each paired-end library.

########################################################################################
# Getting all the file names.

relink <- "make_links.sh"
write(file=relink, c("set -e", "set -u", "mkdir fastq"), ncol=1)

fpath <- "/run/user/1753941046/gvfs/smb-share:server=jmlab-data,share=jmlab/group_folders/lun01/Reik/primed_HSCs"
all.md5 <- list()
for (x in unique(collected$Run)) {
    curpath <- file.path(fpath, paste0("real-", x, "_20161024"), "fastq")
    chosen <- list.files(curpath, pattern="fastq.gz$")
    write(file=relink, paste0("ln -s ", file.path(curpath, chosen), " ", file.path("fastq", chosen)), append=TRUE, ncol=1)

    curmd5 <- read.table(file.path(curpath, "md5.all"), stringsAsFactors=FALSE)
    curmd5 <- curmd5[match(chosen, curmd5[,2]),] 
    all.md5[[x]]  <- curmd5
}                                                                

all.md5 <- do.call(rbind, all.md5)
m <- match(collected$File, all.md5[,2])
any(is.na(m))
any(duplicated(m))
collected$MD5 <- all.md5[m,1]

########################################################################################
# Constructing the sdrf.tsv file.

output <- list()
output[["Source Name"]] <- collected$Sample
output[["Characteristics[organism]"]] <- "Homo sapiens"
output[["Characteristics[cell line]"]] <- "H9"
output[["Material Type"]] <- "RNA"
output[[paste0(rep(c("Protocol REF", "Performer"), 6), collapse="\t")]] <- paste0(c("Obtaining H9 cells", "Ferdinand von Meyenn",
                                                                                    "Culturing H9 cells", "Ferdinand von Meyenn",
                                                                                    "Reverse transcription", "Ferdinand von Meyenn",
                                                                                    "Extracting RNA", "Ferdinand von Meyenn",
                                                                                    "Creating libraries","Ferdinand von Meyenn"
                                                                                    ), collapse="\t")
output[["Extract Name"]] <- collected$Sample
output[["Comment[LIBRARY_LAYOUT]"]] <- "PAIRED"
output[["Comment[LIBRARY_SELECTION]"]] <- "Oligo-dT"
output[["Comment[LIBRARY_SOURCE]"]] <- "TRANSCRIPTOMIC"
output[["Comment[LIBRARY_STRAND]"]] <- "not applicable"
output[["Comment[LIBRARY_STRATEGY]"]] <- "RNA-seq"
output[["Comment[NOMINAL_LENGTH]"]] <- 295
output[["Comment[NOMINAL_SDEV]"]] <- 25
output[["Comment[ORIENTATION]"]] <- "5'-3'-3'-5'"
output[["Protocol REF\tPerformer"]] <- "Sequencing libraries\tFerdinand von Meyenn"
output[["Assay Name"]] <- collected$Sample
output[["Technology Type"]] <- "sequencing assay"
output[["Comment[experiment batch]"]] <- collected$Batch
output[["Comment[sequencing run]"]] <- collected$Run
output[["Array Data File"]] <- collected$File
output[["Protocol REF"]] <- "Assigning reads to genes"
output[["Derived Array Data File"]] <- paste0("genic_counts_", collected$Run, ".tsv")
output[["Comment[MD5]"]] <- collected$MD5
output[["Factor Value[phenotype]"]] <- collected$Phenotype

output$check.names <- FALSE
sdrf <- do.call(data.frame, output)
write.table(file="sdrf.tsv", sdrf, row.names=FALSE, sep="\t", quote=FALSE)
