source("central.R")

hvg.naive <- read.table(file.path("results-naive", "hvg.tsv"), header=TRUE)
hvg.primed <- read.table(file.path("results-primed", "hvg.tsv"), header=TRUE)
var.naive <- read.table(file.path("results-naive", "var.tsv"), header=TRUE)
var.primed <- read.table(file.path("results-primed", "var.tsv"), header=TRUE)
de.res <- read.table(file.path("results-overall", "de.tsv"), header=TRUE)
de.res <- de.res[rownames(var.naive),]

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
sce_naive <- readRDS("sce_naive.rds")
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
legend(12, 0, xjust=1, yjust=0, legend=c("Not shared", "Shared HVG", "Shared, cell cycle"), 
       col=c("grey80", "black", "coral"), pch=16, cex=1.2)
dev.off()    
    
# Figure 2C
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
plot(op, yp, ylab="Biological component", col=primed.col, xlab="Relative rank (%)", type="l", lwd=3, cex.axis=1.2, cex.lab=1.4)
lines(on, yn, col=naive.col, lwd=3)

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
                sprintf("Primed-only, brain development"),
                sprintf("Primed-only, cilia assembly"),
                sprintf("Naive-only HVG (%i)", length(on)),
                sprintf("Naive-only, brain development"),
                sprintf("Naive-only, cilia assembly")),
       col=c(primed.col, "black", "black", naive.col, "black", "black"),
       pt.bg=c(primed.col, primed.sub.col, primed.sub.col, naive.col, naive.sub.col, naive.sub.col),
       pch=c(NA, 25, 21, NA, 25, 21),
       lwd=c(2,NA,NA,2,NA,NA), cex=1.2)
dev.off()


# Supplements

# Plotting all of the HVGs together in one pile.
x <- -de.res$logFC
y <- var.primed$bio - var.naive$bio
plot(x, y, pch=16, cex=0.5,
     ylab="Difference in variance (primed - naive)",
     xlab=expression(Log[2]~"fold change (primed/naive)")) 

hvgs <- c(rownames(hvg.naive), rownames(hvg.primed))
cols <- rep(c(naive.col, primed.col), c(nrow(hvg.naive), nrow(hvg.primed)))
set.seed(100)
s <- sample(length(cols))
hvgs <- hvgs[s]
cols <- cols[s]
m <- match(hvgs, rownames(var.naive))
points(x[m], y[m], col=cols, pch=16, cex=0.5)


sce_primed <- readRDS(file = paste0(objectpath, "new_sce_primed_object"))
sce_naive <- readRDS(file = paste0(objectpath, "new_sce_naive_object"))

# A

par(mfrow=c(1,1))
primed.phases <- table(sce_primed$phase)
naive.phases <- table(sce_naive$phase)

pdf(file="./figures/sup.figure3a.pdf")
barplot(rbind(primed.phases, naive.phases), beside = TRUE, col = c(primed.col, naive.col), ylab="Number of cells")
legend("top", legend = c("Primed", "Naive"), fill = c(primed.col, naive.col), bty="n", cex=1.5)
dev.off()

# B
naive.means <- log(rowMeans(counts(sce_naive)))
primed.means <- log(rowMeans(counts(sce_primed)))

naive.genenum <- table(naive.means > 1)["TRUE"]
primed.genenum <- table(primed.means > 1)["TRUE"]
shared <- length(intersect(rownames(sce_naive[naive.means >1,]), rownames(sce_primed[primed.means >1,])))

num.matrix <- cbind(primed.genenum, naive.genenum, shared)
colnames(num.matrix) <- c("Primed", "Naive", "Shared")

pdf(file="./figures/sup.figure3b.pdf")
par(mfrow=c(1,1), mar = c(5.1, 5.5, 4.1, 2.1), las = 1)
plot <- barplot(num.matrix, col = c(primed.col, naive.col, "#CC79A7"), beside = TRUE, names.arg = c("Primed", "Naive", "Shared"), 
                ylab="Number of expressed genes", xlim = c(0,3), width = 0.3, space = c(0,.5),
                cex.axis = .75, cex.lab = 1.3, cex.names = 1., ylim=c(0,14700))
text(x = plot, y = num.matrix, label = num.matrix, cex = 0.75, pos = 3)
dev.off()

# C

pdf(file="./figures/sup.figure3c.pdf")
par(mfrow=c(1,2), mar=c(5.1, 4.1, 4.1, 0.1))
hist(primed.means, xlim = c(0,11), breaks = 100, col=primed.col, xlab="log mean", main="Primed")
hist(naive.means, xlim = c(0,11), breaks = 100, col=naive.col, xlab="log mean", main="Naive")
dev.off()

# D

naive.size <- sce_naive$size_factor
primed.size <- sce_primed$size_factor

pdf(file="./figures/sup.figure3d.pdf")
par(mfrow=c(1,2), mar=c(5.1, 4.1, 4.1, 1), las=1)
hist(primed.size, xlim = c(0,8), breaks = 50, col=primed.col, xlab="Size factor", main="Primed")
hist(naive.size, xlim = c(0,8), breaks = 50, col=naive.col, xlab="Size factor", main="Naive", ylab = "")
dev.off()

# E
var.out_naive <- readRDS(file = paste0(objectpath, "var.naive"))
var.out_primed <-  readRDS(file = paste0(objectpath, "var.primed"))

pdf(file="./figures/sup.figure3e.pdf")

logFCs <- (rowMeans(exprs(sce_primed))) - (rowMeans(exprs(sce_naive)))
biovar <- var.out_primed$bio - var.out_naive$bio
colours <- rep("black", nrow(sce_naive))
colours[which(rownames(sce_naive) %in% rownames(hvg.out_naive))] <- naive.col
colours[which(rownames(sce_naive) %in% rownames(hvg.out_primed))] <- primed.col
colours[which(rownames(sce_naive) %in% shared_genes)] <- "#CC79A7"
plot <- cbind(logFCs, biovar, colours)
#plot <- plot[order(plot[,3], decreasing = T),]
plot(plot[,"logFCs"],plot[,"biovar"], pch=16, col = plot[,"colours"], cex=0.6, ylab = expression(paste(Delta, " Biological variance")), xlab="mean [logFC]")
legend("bottomright", legend = c("HVG [primed]", "HVG [naive]", "HVG [shared]"), 
       col = c(primed.col, naive.col, "#CC79A7"), bty="o", pch=16, box.col = "white", bg="white")

dev.off()
