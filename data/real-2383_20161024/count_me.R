
anno.files <- file.path("../../genomes/annotation/", c("hg38.gtf", "ERCC92.gtf"))
bam.files <- list.files("bam", full=TRUE, pattern="bam$")
stat.file <- "all_qual.tsv"
ispet <- TRUE
# NOT STRAND-SPECIFIC.

# Files of interest.
bam.files
anno.files

if (!exists("ispet")) { 
    ispet <- FALSE
} 
if (!exists("strandspec")) {
    strandspec <- 0
}
if (!exists("minq")) { 
    minq <- 10
}
if (!exists("additional")) {
    additional <- list()
}

ispet
strandspec
minq
additional

if (length(anno.files)==1L) {
    file.symlink(anno.files, "temp.gtf")
} else {
    system(paste(c("cat", anno.files, "> temp.gtf"), collapse=" "))
}

# Running featureCounts.

additional$minMQS <- minq
additional$isPairedEnd <- ispet
additional$strandSpecific <- strandspec
require(Rsubread)
out <- do.call(featureCounts, c(list(files=bam.files, annot.ext="temp.gtf", isGTFAnnotationFile=TRUE, nthreads=4), additional))

# Saving counts to file, with gene names.
colnames(out$counts) <- sub("\\.bam$", "", basename(bam.files))
final <- data.frame(GeneID=rownames(out$counts), Length=out$annotation$Length, out$counts, check.names=FALSE)
write.table(file="genic_counts.tsv", final, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Augmenting the stats.
my.stats <- as.data.frame(t(out$stat[,-1]))
colnames(my.stats) <- out$stat[,1]
rownames(my.stats) <- colnames(out$counts)
write.table(file="my_qual.tsv", my.stats, col.names=NA, quote=FALSE, sep="\t")

# Saving the session information.
unlink("temp.gtf")
sessionInfo()

