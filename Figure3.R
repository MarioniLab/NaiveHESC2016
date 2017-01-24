source("central.R")

# Figure 3A
hvg.naive <- read.table(file.path("results-naive", "hvg.tsv"), header=TRUE)
hvg.primed <- read.table(file.path("results-primed", "hvg.tsv"), header=TRUE)

pdf(file=file.path(figdir, "3a.pdf"), width=5, useDingbats=FALSE)
ylimit <- c(0, 6)
nb <- hvg.naive$bio
pb <- hvg.primed$bio
boxplot(list(Naive=nb, Primed=pb), pch=16, cex=0.5,
        col=c(naive.col, primed.col), ylim=ylimit, 
        ylab="Biological component", cex.axis=1.2, cex.lab=1.4, cex.names=1.4)
text(1, min(nb), pos=1, nrow(hvg.naive), cex=1.2)
text(2, min(pb), pos=1, nrow(hvg.primed), cex=1.2)
dev.off()
### wilcox.test(hvg.naive$bio, Primed=hvg.primed$bio) # should give < 2e-16

# Figure 3B

library(scater)
sce_all <- readRDS("sce_all.rds")
library(org.Hs.eg.db)
cell.anno <- select(org.Hs.eg.db, keys="GO:0007049", keytype="GOALL", columns="ENTREZID")
in.cell.cycle <- fData(sce_all)$entrezgene %in% cell.anno$ENTREZID

var.naive <- read.table(file.path("results-naive", "var.tsv"), header=TRUE)
var.primed <- read.table(file.path("results-primed", "var.tsv"), header=TRUE)

cycle.naive <- var.naive[in.cell.cycle & rownames(var.naive) %in% rownames(hvg.naive),]
cycle.primed <- var.primed[in.cell.cycle & rownames(var.primed) %in% rownames(hvg.primed),]
rn <- rank(-cycle.naive$bio)
rp <- rank(-cycle.primed$bio)
is <- which(is.na(var.naive$FDR))
osn <- is[order(var.naive$bio[is], decreasing=TRUE)]
osp <- is[order(var.primed$bio[is], decreasing=TRUE)]

pdf(file=file.path(figdir, "3b.pdf"), width=7, useDingbats=FALSE)
ylimits <- c(0, 12)
xlimits <- c(0, 50)
plot(rn, cycle.naive$bio, pch=16, col=naive.col, ylim=ylimits, xlim=xlimits,
     xlab="Rank", ylab="Biological component", cex.axis=1.2, cex.lab=1.4)
points(rp, cycle.primed$bio, pch=16, col=primed.col)
o <- order(rn)
lines(rn[o], cycle.naive$bio[o], col=naive.col)
o <- order(rp)
lines(rp[o], cycle.primed$bio[o], col=primed.col)
points(var.naive$bio[osn], pch=1, col=naive.col, ylim=c(0, 12))
points(var.primed$bio[osp], pch=1, col=primed.col)

legend(xlimits[2], ylimits[2], xjust=1, yjust=1, lty=c(1, 1, 0, 0),
       legend=c("Cell cycle HVG (naive)", "Cell cycle HVG (primed)", 
                "Spike-in transcript (naive)", "Spike-in transcript (primed)"),
       col=c(naive.col, primed.col), pch=c(16, 16, 1, 1), cex=1.2)

#intersect(rownames(cycle.primed[rp <= 50,]), rownames(cycle.naive[rn <= 50,]))
to.annotate <- c("CDK1", "PLK1", "CCNA2", "CENPE", "NUF2")
xoffset <- c(5, 2, 3,  5, 10)  
yoffset <- c(-0.2, -2, 3, -1, 2)     
positions <- c(4, 1, 4, 1, 4)
mn <- match(to.annotate, rownames(cycle.naive))
hostx <- rn[mn]
extrax <- hostx + xoffset
hosty <- cycle.naive$bio[mn]
extray <- hosty + yoffset
text(extrax, extray, to.annotate, pos=positions)
segments(hostx, hosty, extrax, extray)
mp <- match(to.annotate, rownames(cycle.primed))
hostx2 <- rp[mp]
hosty2 <- cycle.primed$bio[mp]
segments(hostx2, hosty2, extrax, extray)
dev.off()

# Figure 3C

splice.anno <- select(org.Hs.eg.db, keys="03040", keytype="PATH", columns="ENTREZID")
in.splice <- fData(sce_all)$entrezgene %in% splice.anno$ENTREZID
splice.naive <- var.naive[in.splice & rownames(var.naive) %in% rownames(hvg.naive),]
splice.primed <- var.primed[in.splice & rownames(var.primed) %in% rownames(hvg.primed),]

yp <- density(splice.primed$bio)
yn <- density(splice.naive$bio)
yn$y <- yn$y * nrow(splice.naive)/nrow(splice.primed)

pdf(file.path(figdir, "3c.pdf"), width=7)
par(mar = c(5.1, 5.1, 4.1, 2.1), xpd=TRUE)
plot(yp$x, yp$y, xlab="Biological component", ylab="", yaxt="n", type="n", cex.axis=1.2, cex.lab=1.4, bty="n")
polygon(c(min(yp$x), yp$x, max(yp$x)), c(0, yp$y, 0), border=NA, col="lightblue")
lines(yp$x, yp$y, lwd=2, col=primed.col)
text(4, 0, sprintf("Primed\n(%i)", nrow(splice.primed)), pos=3, cex=1.2)
polygon(c(min(yn$x), yn$x, max(yn$x)), c(0, yn$y, 0), border=NA, col="gold1")
lines(yn$x, yn$y, lwd=2, col=naive.col)
text(0.9, 0, sprintf("Naive\n(%i)", nrow(splice.naive)), pos=3, cex=1.2)

#rownames(splice.primed)[order(splice.primed$bio, decreasing=TRUE)]
to.annotate <- c("SLU7", "NCBP1", "THOC2", "EFTUD2", "SRSF4")
xoffset <- c(0.5, 0.5, 0, 1, 1)  
yoffset <- c(0.05, 0.05, 0.05, 0.05, 0.05)
positions <- c(3, 3, 3, 4, 4)
mn <- match(to.annotate, rownames(splice.primed))
hostx <- splice.primed$bio[mn]
extrax <- hostx + xoffset
hosty <- approx(yp$x, yp$y, xout=hostx)$y
extray <- hosty + yoffset
text(extrax, extray, to.annotate, pos=positions)
segments(hostx, hosty, extrax, extray)
dev.off()

### Supplements

# A
naive.genenum <- sum(var.naive$mean > 0)
primed.genenum <- sum(var.primed$mean > 0)

pdf(file=file.path(figdir, "s3a.pdf"), width=8, height=5)
par(mar=c(5.1, 5.1, 4.1, 1.1), mfrow=c(1,3))
gcounts <- c(Primed=primed.genenum, Naive=naive.genenum)
out <- barplot(gcounts/1e3, col = c(primed.col, naive.col), ylim=c(0, 16), width=0.5,
        ylab=expression("Number of expressed genes ("*10^3*")"), cex.axis=1.2, cex.names=1.4, cex.lab=1.4)
text(out, gcounts/1e3, gcounts, pos=3, cex=1.2)

gcounts <- c(Primed=sum(sce_all$phenotype=="primed"), Naive=sum(sce_all$phenotype=="naive"))
out <- barplot(gcounts, col = c(primed.col, naive.col), ylim=c(0, 450), width=0.5,
        ylab="Number of cells", cex.axis=1.2, cex.names=1.4, cex.lab=1.4)
text(out, gcounts, gcounts, pos=3, cex=1.2)

by.pheno <- list(Primed=sizeFactors(sce_all)[sce_all$phenotype=="primed"],
                 Naive=sizeFactors(sce_all)[sce_all$phenotype=="naive"])
boxplot(by.pheno, col=c(primed.col, naive.col), ylab="Size factor per cell", 
        log="y", cex.axis=1.2, cex.names=1.4, cex.lab=1.4, range=0)
dev.off()

# B

pdf(file=file.path(figdir, "s3b.pdf"), width=10, height=6, useDingbats=FALSE)
par(mfrow=c(1,2))
for (mode in 1:2) { 
    if (mode==1) {
        nb <- sort(hvg.naive$bio)
        pb <- sort(hvg.primed$bio)
        ylimits <- c(0, 6)
        ylab <- "Biological component"
        xlab <- "by component"
    } else {
        nb <- sort(hvg.naive$total/hvg.naive$tech)
        pb <- sort(hvg.primed$total/hvg.primed$tech)
        ylimits <- c(1, 5)
        ylab <- "Total:technical ratio"
        xlab <- "by ratio"
    }    

    N <- 4
    nint <- round(seq(from=0, to=length(nb), length.out=N+1))
    pint <- round(seq(from=0, to=length(pb), length.out=N+1))
    to.plot <- list()
    for (i in seq_len(N)) {
        ncur <- (1+nint[i]):nint[i+1]
        pcur <- (1+pint[i]):pint[i+1]
        to.plot <- c(to.plot, list(nb[ncur], pb[pcur]))
    }

    xcoords <- rbind(seq_len(N)*2.5, seq_len(N)*2.5 + 1)
    boxplot(to.plot, at=xcoords, xaxt="n", xlab=sprintf("Quartile %s", xlab), ylab=ylab, ylim=ylimits,
            col=c(naive.col, primed.col), pch=16, cex=0.5, cex.axis=1.2, cex.lab=1.4)
    axis(1, at=colMeans(xcoords), label=seq_len(N), cex.axis=1.2)    
    legend(xcoords[1]-0.5, ylimits[2], fill=c(naive.col, primed.col), legend=c("Naive", "Primed"), cex=1.2)
}
dev.off()

# C
primed.phases <- table(sce_all$phase[sce_all$phenotype=="primed"])
naive.phases <- table(sce_all$phase[sce_all$phenotype=="naive"])
pdf(file=file.path(figdir, "s3c.pdf"), width=8, height=6)
barplot(rbind(primed.phases, naive.phases), beside = TRUE, cex.axis=1.2, cex.names=1.5, cex.lab=1.5,
            col = c(primed.col, naive.col), ylab="Number of cells")
legend("topright", legend = c("Primed", "Naive"), fill = c(primed.col, naive.col), cex=1.5)
dev.off()
