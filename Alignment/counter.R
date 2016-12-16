 # Files of interest.

bam.files <- list.files(path="./bam", full.names=TRUE, pattern = '.bam$')

anno.files <- c("/lustre/jmlab/resources/annotation/processed/hg38.gtf", "/lustre/jmlab/resources/annotation/original/ERCC92.gtf")

if (!exists("ispet")) { 
    ispet <- TRUE 
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
final <- data.frame(GeneID=rownames(out$counts), Length=out$annotation$Length, out$counts)
write.table(file="genic_counts.tsv", final, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Augmenting the stats.
my.stats <- as.data.frame(t(out$stat[,-1]))
colnames(my.stats) <- out$stat[,1]
rownames(my.stats) <- colnames(out$counts)
write.table(file="my_qual.tsv", my.stats, col.names=NA, quote=FALSE, sep="\t")

# Saving the session information.
unlink("temp.gtf")
sessionInfo()

