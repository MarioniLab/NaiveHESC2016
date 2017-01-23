source("central.R")

# Figure 3A
hvg.naive <- read.table(file.path("results-naive", "hvg.tsv"), header=TRUE)
hvg.primed <- read.table(file.path("results-primed", "hvg.tsv"), header=TRUE)

pdf(file=file.path(figdir, "3a.pdf"), width=5)
ylimit <- c(0, 6)
nb <- hvg.naive$bio
pb <- hvg.primed$bio
boxplot(list(Naive=nb, Primed=pb), pch=16, cex=0.5,
        col=c(naive.col, primed.col), ylim=ylimit, main="Top 1000 HVGs",
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
in.both <- rownames(sce_all) %in% intersect(rownames(hvg.naive), rownames(hvg.primed))

var.naive <- read.table(file.path("results-naive", "var.tsv"), header=TRUE)
var.primed <- read.table(file.path("results-primed", "var.tsv"), header=TRUE)

pdf(file=file.path(figdir, "3b.pdf"), width=7)
x <- var.naive$bio[in.cell.cycle]
y <- var.primed$bio[in.cell.cycle]
limits <- c(0, max(c(x, y)))
plot(x, y, xlim=limits, ylim=limits, pch=16, col="grey", cex.axis=1.2, cex.lab=1.4,
     xlab="Biological component (naive)", ylab="Biological component (primed)")
x <- var.naive$bio[in.cell.cycle & in.both]
y <- var.primed$bio[in.cell.cycle & in.both]
points(x, y, pch=16, col="black")
abline(0,1,col="red", lwd=2, lty=2)

legend(limits[1], limits[2], legend=c("Not shared HVG", "Shared HVG"), col=c("grey", "black"), pch=16, cex=1.2)
text(9, 10,  labels=sum(x < y), col=primed.col, cex=1.2)
text(10, 9,  labels=sum(x > y), col=naive.col, cex=1.2)
dev.off()

# Figure 3C

go.res <- read.table(file.path("results-primed", "go_diff_hvg.tsv"), sep="\t", row.names=1, header=TRUE, comment="", quote="")
kegg.res <- read.table(file.path("results-primed", "kegg_diff_hvg.tsv"), header=TRUE, sep="\t", quote="", row.names=1, 
                        colClasses=c("character", "character", rep("numeric", 3)))

go.interest <- c("GO:0008380", # RNA splicing
                 "GO:0042254", # ribosome biogenesis
                 "GO:0006412", # translation
                 "GO:0043604") # amide biosynthetic process
to.plot <- list()
extra.data <- list()
for (term in go.interest) {
    to.plot[[term]] <- c(go.res[term,"Prop.primed"], go.res[term, "Prop.naive"])
    extra.data[[term]] <- c(go.res[term,"Num.primed"], go.res[term, "Num.naive"])
}

kegg.interest <- c("03040", # Spliceosome
                   "03013", # RNA transport 
                   "03050") # Proteasome
for (term in kegg.interest) {
    to.plot[[term]] <- c(kegg.res[term,"Prop.primed"], kegg.res[term, "Prop.naive"])
    extra.data[[term]] <- c(kegg.res[term,"Num.primed"], kegg.res[term, "Num.naive"])
}

to.plot <- do.call(cbind, to.plot) * 100
colnames(to.plot) <- c("RNA splicing", "Ribosome biogenesis", "Translation", "Amide biosynthesis",
                       "Spliceosome", "RNA transport", "Proteosome")
colnames(to.plot) <- sprintf("%s\n(%s)", colnames(to.plot), c(go.interest, kegg.interest))
to.plot <- to.plot[,ncol(to.plot):1]
extra.data <- do.call(cbind, extra.data)
extra.data <- extra.data[,ncol(extra.data):1]

pdf(file.path(figdir, "3c.pdf"), width=7)
par(mar=c(5.1, 11, 4.1, 2.1))
xout <- barplot(to.plot, beside=TRUE, col=c(primed.col, naive.col), las=1, horiz=TRUE,
                xlab="Percentage of HVGs", xlim=c(0, 5.5), cex.axis=1.2, cex.lab=1.4, cex.names=1.2)
text(to.plot-0.1, xout, extra.data, pos=4) 
ucoords <- par()$usr
legend(ucoords[2] - (ucoords[2] - ucoords[1])*0.05, xjust=1,
       ucoords[3] + (ucoords[4] - ucoords[3])*0.05, yjust=0,
       legend=c("Naive", "Primed"), fill=c(naive.col, primed.col), cex=1.2)
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

pdf(file=file.path(figdir, "s3b.pdf"), width=10, height=6)
par(mfrow=c(1,2))
for (mode in 1:2) { 
    if (mode==1) {
        nb <- sort(hvg.naive$bio)
        pb <- sort(hvg.primed$bio)
        ylimits <- c(0, 6)
        ylab <- "Biological component"
    } else {
        nb <- sort(hvg.naive$bio/hvg.naive$tech)
        pb <- sort(hvg.primed$bio/hvg.primed$tech)
        ylimits <- c(0, 4)
        ylab <- "Biological:technical ratio"
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
    boxplot(to.plot, at=xcoords, xaxt="n", xlab="Quartile", ylab=ylab, ylim=ylimits,
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
