source("central.R")

hvg.naive <- read.table(file.path("results-naive", "hvg.tsv"), header=TRUE)
hvg.primed <- read.table(file.path("results-primed", "hvg.tsv"), header=TRUE)
var.naive <- read.table(file.path("results-naive", "var.tsv"), header=TRUE)
var.primed <- read.table(file.path("results-primed", "var.tsv"), header=TRUE)
sce_naive <- readRDS("sce_naive.rds")
sce_primed <- readRDS("sce_primed.rds")

# Figure 3A
pdf(file=file.path(figdir, "3a.pdf"), width = 4)
par(mar = c(5.1, 5.5, 4.1, 2.1), las = 1)
hvg.num <- c(Primed=nrow(hvg.primed), Naive=nrow(hvg.naive))
out <- barplot(hvg.num/1e3, ylab=expression("Number of highly variable genes ("*10^3*")"),
        ylim = c(0, 7), col = c(primed.col, naive.col), beside = TRUE,
        width = 0.7, space = c(0.2), cex.axis = 1.2, cex.lab = 1.3, cex.names = 1.3)
text(x = out, y = hvg.num/1e3, label = hvg.num, cex = 1.2, pos = 3)
dev.off()

# Figure 3B
library(org.Hs.eg.db)
anno <- select(org.Hs.eg.db, key="GO:0000278", keytype="GOALL", column="ENSEMBL")
library(scater)
stopifnot(identical(rownames(sce_naive), rownames(var.naive)))
is.cycling <- fData(sce_naive)$ensembl %in% anno$ENSEMBL

pdf(file=file.path(figdir, "3b.pdf"))
plot(var.naive$bio, var.primed$bio, xlim=c(0, 12), ylim=c(0, 12), pch=16, cex=0.5, col="grey80",
     xlab="Biological component (naive)", ylab="Biological component (primed)", cex.axis=1.2, cex.lab=1.4)
hvgs <- intersect(rownames(hvg.naive), rownames(hvg.primed))
is.shared <- rownames(var.naive) %in% hvgs
points(var.naive$bio[is.shared], var.primed$bio[is.shared], col="black", pch=16, cex=0.8)
abline(0, 1, col="red", lwd=2)

points(var.naive$bio[is.cycling & is.shared], var.primed$bio[is.cycling & is.shared], col="coral", pch=16, cex=1)
legend(12, 0, xjust=1, yjust=0, legend=c("Not shared", "Shared HVG", "Shared + cell cycle"), 
       col=c("grey80", "black", "coral"), pch=16, cex=1.2)
dev.off()    
    
# Figure 3C
anno <- select(org.Hs.eg.db, key="GO:0007420", keytype="GOALL", column="ENSEMBL")
is.brain <- fData(sce_naive)$ensembl %in% anno$ENSEMBL
anno <- select(org.Hs.eg.db, key="GO:0042384", keytype="GOALL", column="ENSEMBL")
is.cilia <- fData(sce_naive)$ensembl %in% anno$ENSEMBL

unique.primed <- hvg.primed[!rownames(hvg.primed) %in% rownames(hvg.naive),]
unique.naive <- hvg.naive[!rownames(hvg.naive) %in% rownames(hvg.primed),]
op <- seq_len(nrow(unique.primed))/nrow(unique.primed) * 100
on <- seq_len(nrow(unique.naive))/nrow(unique.naive) * 100
yp <- sort(unique.primed$bio)
yn <- sort(unique.naive$bio)

pdf(file=file.path(figdir, "3c.pdf"))
par(mar = c(5.1, 5.1, 4.1, 2.1))
plot(on, yn, ylab="Biological component", col=naive.col, xlab="Relative rank (%)", type="l", lwd=3, cex.axis=1.2, cex.lab=1.4)
lines(op, yp, col=primed.col, lwd=3)

is.primed.brain <- is.brain[match(rownames(unique.primed), rownames(var.primed))]
bonus <- 0.1
primed.sub.col <- "royalblue3"
naive.sub.col <- "goldenrod3"
points(op[is.primed.brain], yp[is.primed.brain]+bonus, pch=25, bg=primed.sub.col)
is.naive.brain <- is.brain[match(rownames(unique.naive), rownames(var.naive))]
points(on[is.naive.brain], yn[is.naive.brain]+bonus, pch=25, bg=naive.sub.col)

is.primed.cilia <- is.cilia[match(rownames(unique.primed), rownames(var.primed))]
points(op[is.primed.cilia], yp[is.primed.cilia], pch=21, bg=primed.sub.col)
is.naive.cilia <- is.cilia[match(rownames(unique.naive), rownames(var.naive))]
points(on[is.naive.cilia], yn[is.naive.cilia], pch=21, bg=naive.sub.col)

legend(0, max(yp, yn), 
       legend=c(sprintf("Primed-only HVG (%i)", length(op)),
                sprintf("Primed-only + brain development"),
                sprintf("Primed-only + cilia assembly"),
                sprintf("Naive-only HVG (%i)", length(on)),
                sprintf("Naive-only + brain development"),
                sprintf("Naive-only + cilia assembly")),
       col=c(primed.col, "black", "black", naive.col, "black", "black"),
       pt.bg=c(primed.col, primed.sub.col, primed.sub.col, naive.col, naive.sub.col, naive.sub.col),
       pch=c(NA, 25, 21, NA, 25, 21),
       lwd=c(2,NA,NA,2,NA,NA), cex=1.2)
dev.off()

### Supplements

# A
naive.genenum <- sum(var.naive$mean > 0)
primed.genenum <- sum(var.primed$mean > 0)
sce_all <- readRDS("sce_all.rds")
sce_primed <- readRDS("sce_primed.rds")

pdf(file=file.path(figdir, "s3a.pdf"), width=8, height=5)
par(mar=c(5.1, 5.1, 4.1, 2.1), mfrow=c(1,3))
gcounts <- c(Primed=primed.genenum, Naive=naive.genenum)
out <- barplot(gcounts/1e3, col = c(primed.col, naive.col), ylim=c(0, 16), width=0.5,
        ylab=expression("Number of expressed genes  ("*10^3*")"), cex.axis=1.2, cex.names=1.4, cex.lab=1.4)
text(out, gcounts/1e3, gcounts, pos=3, cex=1.2)

gcounts <- c(Primed=unname(ncol(sce_primed)), Naive=unname(ncol(sce_naive)))
out <- barplot(gcounts, col = c(primed.col, naive.col), ylim=c(0, 450), width=0.5,
        ylab="Number of cells", cex.axis=1.2, cex.names=1.4, cex.lab=1.4)
text(out, gcounts, gcounts, pos=3, cex=1.2)

by.pheno <- list(Primed=sizeFactors(sce_all)[sce_all$phenotype=="primed"],
                 Naive=sizeFactors(sce_all)[sce_all$phenotype=="naive"])
boxplot(by.pheno, col=c(primed.col, naive.col), ylab="Size factor per cell", log="y", cex.axis=1.2, cex.names=1.4, cex.lab=1.4, range=0)
dev.off()

# B
de.res <- read.table(file.path("results-overall", "de.tsv"), header=TRUE)
de.res <- de.res[rownames(var.naive),]

x <- -de.res$logFC
xrange <- range(x[!is.na(x)]) # Due to the spikes, when subsetting by row names above.
y <- var.primed$bio - var.naive$bio
yrange <- range(y[!is.na(y)])

pdf(file=file.path(figdir, "s3b.pdf"), width=8, height=8)
par(mfrow=c(2,2), mar=c(4.1, 4.1, 3.1, 1.1))
for (mode in 1:4) {
    if (mode==1) {
        X <- x
        Y <- y
        main <- "All genes"
    } else if (mode==2) {
        hvgs <- intersect(rownames(hvg.naive), rownames(hvg.primed))
        is.shared <- rownames(var.naive) %in% hvgs
        X <- x[is.shared]
        Y <- y[is.shared]
        main <- "Shared HVGs"
    } else if (mode==3) {
        hvgs <- setdiff(rownames(hvg.primed), rownames(hvg.naive))
        is.primed <- rownames(var.naive) %in% hvgs
        X <- x[is.primed]
        Y <- y[is.primed]
        main <- "Primed-only HVGs"
    } else if (mode==4) {
        hvgs <- setdiff(rownames(hvg.naive), rownames(hvg.primed))
        is.naive <- rownames(var.naive) %in% hvgs
        X <- x[is.naive]
        Y <- y[is.naive]
        main <- "Naive-only HVGs"
    }
    smoothScatter(X, Y, main=main, xlim=xrange, ylim=yrange,
                  ylab=expression(Delta~"bio. component (primed - naive)"),
                  xlab=expression(Log[2]~"fold change (primed/naive)")) 
}
dev.off()

# C

pdf(file=file.path(figdir, "s3c.pdf"), width=12, height=6)
par(mfrow=c(1,2), mar=c(5.1, 4.1, 2.1, 2.1))
plot(var.primed$mean, var.primed$total, col="grey80", pch=16, ylim=c(0, 15), xlim=c(0, 16),
     xlab=expression("Average"~log[2]~"count"), ylab="Variance", cex.axis=1.2, cex.lab=1.4)
is.hvg <- rownames(var.primed) %in% rownames(hvg.primed)
points(var.primed$mean[is.hvg], var.primed$total[is.hvg], col=primed.col, pch=16)
o <- order(var.primed$mean)
lines(var.primed$mean[o], var.primed$tech[o], col="red", lwd=2)
legend(16, 15, xjust=1, yjust=1, 
       legend=c("Not HVG",
                "Primed HVG", 
                "Technical"), col=c("grey80", primed.col, "red"),
       pch=c(16, 16, NA), lwd=c(NA, NA, 2))

plot(var.naive$mean, var.naive$total, col="grey80", pch=16, ylim=c(0, 15), xlim=c(0, 16),
     xlab=expression("Average"~log[2]~"count"), ylab="Variance",  cex.axis=1.2, cex.lab=1.4)
is.hvg <- rownames(var.naive) %in% rownames(hvg.naive)
points(var.naive$mean[is.hvg], var.naive$total[is.hvg], col=naive.col, pch=16)
o <- order(var.naive$mean)
lines(var.naive$mean[o], var.naive$tech[o], col="red", lwd=2)
legend(16, 15, xjust=1, yjust=1, 
       legend=c("Not HVG",
                "Naive HVG", 
                "Technical"), col=c("grey80", naive.col, "red"),
       pch=c(16, 16, NA), lwd=c(NA, NA, 2))
dev.off()

# D

primed.phases <- table(sce_primed$phase)
naive.phases <- table(sce_naive$phase)
pdf(file=file.path(figdir, "s3d.pdf"), width=8, height=6)
barplot(rbind(primed.phases, naive.phases), beside = TRUE, cex.axis=1.2, cex.names=1.5, cex.lab=1.5,
            col = c(primed.col, naive.col), ylab="Number of cells")
legend("topright", legend = c("Primed", "Naive"), fill = c(primed.col, naive.col), cex=1.5)
dev.off()
