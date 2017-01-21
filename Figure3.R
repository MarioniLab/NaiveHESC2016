source("central.R")

go.primed <- read.table(file.path("results-primed", "go_primed.tsv"), header=TRUE, sep="\t")
go.naive <- read.table(file.path("results-primed", "go_naive.tsv"), header=TRUE, sep="\t")
kegg.primed <- read.table(file.path("results-primed", "kegg_primed.tsv"), header=TRUE, sep="\t", quote="")
kegg.naive <- read.table(file.path("results-primed", "kegg_naive.tsv"), header=TRUE, sep="\t", quote="")

n.sets <- cbind(Primed=c(sum(go.primed$FDR <= 0.05), sum(kegg.primed$FDR <= 0.05)),
                Naive=c(sum(go.naive$FDR <= 0.05), sum(kegg.naive$FDR <= 0.05)))
colors <- c("grey20", "grey80")

# Figure 3A
pdf(file=file.path(figdir, "3a.pdf"), width=5)
xout <- barplot(n.sets, beside=TRUE, ylab="Number of significant gene sets", cex.axis=1.2, 
                cex.names=1.4, cex.lab=1.4, ylim=c(0, max(n.sets)*1.2), col=colors)
text(xout, n.sets, n.sets, pos=3, cex=1.2)
legend("topright", fill=colors, legend=c("GO", "KEGG"), cex=1.2)
dev.off()

# Figure 3B
var.naive <- read.table(file.path("results-naive", "var.tsv"), header=TRUE)
var.primed <- read.table(file.path("results-primed", "var.tsv"), header=TRUE)

library(scater)
sce <- readRDS("sce_all.rds")
stopifnot(identical(rownames(sce), rownames(var.naive)))
stopifnot(identical(rownames(sce), rownames(var.primed)))
dbio <- var.primed$bio - var.naive$bio
dbio[! ( (var.naive$bio > 0 | var.primed$bio > 0) & var.primed$bio!=var.naive$bio)] <- NA

go.interest <- c("GO:0008380", # RNA splicing
                 "GO:0007049", # Cell cycle
                 "GO:0006412", # Translation
                 "GO:0006406", # mRNA export from nucleus
                 "GO:0000209") # protein polyubiquitination
library(org.Hs.eg.db)
go.sets <- select(org.Hs.eg.db, keys=go.interest, keytype="GOALL", column="ENTREZID")
go.sets <- split(go.sets$ENTREZID, go.sets$GOALL)

kegg.interest <- c("04120", # Ubiquitin mediated proteolysis
                   "03013", # RNA transport 
                   "03040", # Spliceosome   
                   "03060", # Protein export
                   "04130") # SNARE interactions in vesicular transport
kegg.sets <- select(org.Hs.eg.db, keys=kegg.interest, keytype="PATH", column="ENTREZID")
kegg.sets <- split(kegg.sets$ENTREZID, kegg.sets$PATH)

to.plot <- list()
all.sets <- c(go.sets, kegg.sets)
for (term in names(all.sets)) {
    current <- fData(sce)$entrezgene %in% all.sets[[term]]
    curbio <- dbio[current]
    to.plot[[term]] <- curbio[!is.na(curbio)]
}

names(to.plot) <- c("RNA splicing", "Cell cycle", "Translation", "mRNA export\nfrom nucleus", "Polyubiquitination",
                    "Ubiquitin-mediated\nproteolysis", "RNA transport", "Spliceosome", "Protein export", "SNARE interactions\nin vesicular transport")

pdf(file.path(figdir, "3b.pdf"), width=10)
par(mar=c(11.5, 4.5, 4.1, 2.1))
boxplot(to.plot, las=2, ylim=c(-3, 4), col=terrain.colors(length(to.plot)), 
            ylab=expression(Delta~"biological component (primed - naive)"), 
            cex.axis=1.2, cex.lab=1.3, cex.names=1.2, pch=16, cex=0.5)
abline(h=0, col="red", lwd=2, lty=2)

par(xpd=TRUE)
coords <- par()$usr
liney <- coords[4] + (coords[4]-coords[3])*0.05
extra <- 0.3
segments(1 - extra, liney, length(go.sets) + extra, liney, lwd=2)
segments(length(go.sets)+1 - extra, liney, length(all.sets) + extra, liney, lwd=2)
mtext("GO", at=(1+length(go.sets))/2, line=1.2, cex =1.4)
mtext("KEGG", at=(1+length(kegg.sets))/2 + length(go.sets), line=1.2, cex=1.4)
dev.off()

### Supplements

# A
naive.nhvg <- nrow(read.table(file.path("results-naive", "hvg.tsv"), header=TRUE))
primed.nhvg <- nrow(read.table(file.path("results-primed", "hvg.tsv"), header=TRUE))
naive.genenum <- sum(var.naive$mean > 0)
primed.genenum <- sum(var.primed$mean > 0)

pdf(file=file.path(figdir, "s3a.pdf"), width=7, height=7)
par(mar=c(5.1, 5.1, 4.1, 2.1), mfrow=c(1,2))
gcounts <- c(Primed=primed.genenum, Naive=naive.genenum)
out <- barplot(gcounts/1e3, col = c(primed.col, naive.col), ylim=c(0, 16), width=0.5,
        ylab=expression("Number of expressed genes ("*10^3*")"), cex.axis=1.2, cex.names=1.4, cex.lab=1.4)
text(out, gcounts/1e3, gcounts, pos=3, cex=1.2)

gcounts <- c(Primed=primed.nhvg, Naive=naive.nhvg)
out <- barplot(gcounts/1e3, col = c(primed.col, naive.col), ylim=c(0, 7), width=0.5,
        ylab=expression("Number of HVGs ("*10^3*")"), cex.axis=1.2, cex.names=1.4, cex.lab=1.4)
text(out, gcounts/1e3, gcounts, pos=3, cex=1.2)
dev.off()

# B 
pdf(file=file.path(figdir, "s3b.pdf"), width=7, height=7)
par(mar=c(5.1, 5.1, 4.1, 2.1), mfrow=c(1,2))
gcounts <- c(Primed=sum(sce$phenotype=="primed"), Naive=sum(sce$phenotype=="naive"))
out <- barplot(gcounts, col = c(primed.col, naive.col), ylim=c(0, 450), width=0.5,
        ylab="Number of cells", cex.axis=1.2, cex.names=1.4, cex.lab=1.4)
text(out, gcounts, gcounts, pos=3, cex=1.2)

by.pheno <- list(Primed=sizeFactors(sce)[sce$phenotype=="primed"],
                 Naive=sizeFactors(sce)[sce$phenotype=="naive"])
boxplot(by.pheno, col=c(primed.col, naive.col), ylab="Size factor per cell", 
        log="y", cex.axis=1.2, cex.names=1.4, cex.lab=1.4, range=0)
dev.off()

# C
primed.phases <- table(sce$phase[sce$phenotype=="primed"])
naive.phases <- table(sce$phase[sce$phenotype=="naive"])
pdf(file=file.path(figdir, "s3c.pdf"), width=8, height=6)
barplot(rbind(primed.phases, naive.phases), beside = TRUE, cex.axis=1.2, cex.names=1.5, cex.lab=1.5,
            col = c(primed.col, naive.col), ylab="Number of cells")
legend("topright", legend = c("Primed", "Naive"), fill = c(primed.col, naive.col), cex=1.5)
dev.off()
