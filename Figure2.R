source("central.R")

# Loading in the expression values.
library(scater)
sce <- readRDS("sce_all.rds")
pops <- read.table(file.path("results-naive", "groups.tsv"), stringsAsFactors=FALSE, header=TRUE)
m <- match(colnames(sce), pops$Cell)
pops <- pops[m,]

###############################
# Figure 2A

# Loading in the DE marker genes.
de.naive <- read.table(file.path("results-naive", "markers_trans_vs_naive.tsv"), header=TRUE, nrows=50)
de.primed <- read.table(file.path("results-naive", "markers_trans_vs_primed.tsv"), header=TRUE, nrows=50)
o <- order(-de.naive$logFC)
top.naive <- rownames(de.naive)[o]
o <- order(de.primed$logFC)
top.primed <- rownames(de.primed)[o]

heat.naive.vals <- exprs(sce)[top.naive,]
heat.naive.vals <- heat.naive.vals - rowMeans(heat.naive.vals)
heat.primed.vals <- exprs(sce)[top.primed,]
heat.primed.vals <- heat.primed.vals - rowMeans(heat.primed.vals)

colramp <- colorRampPalette(c("navy", "white", "orangered"))(50)
bound <- 5
colbreaks <- seq(-bound, bound, length.out=length(colramp)+1)
heat.naive.vals[heat.naive.vals >= bound] <- bound 
heat.naive.vals[heat.naive.vals <= -bound] <- -bound 
heat.naive.vals[heat.naive.vals >= bound] <- bound 
heat.naive.vals[heat.naive.vals <= -bound] <- -bound 

pdf(file=file.path(figdir, "2a.pdf"), width = 10, height = 10)
layout(rbind(c(13,10,11,12,13), c(8,1,2,3,7),c(9,4,5,6,7)), widths=c(0.2, 1, 1, 1, 0.2), heights=c(0.1, 1, 1))

# Adding heatmaps.
par(mar=c(0.5, 0.5, 0.5, 0.5))
image(t(heat.naive.vals[,pops$Type=="naive"]), axes=FALSE, col=colramp, breaks=colbreaks)
image(t(heat.naive.vals[,pops$Type=="transition"]), axes=FALSE, col=colramp, breaks=colbreaks)
image(t(heat.naive.vals[,pops$Type=="primed"]), axes=FALSE, col=colramp, breaks=colbreaks)

image(t(heat.primed.vals[,pops$Type=="naive"]), axes=FALSE, col=colramp, breaks=colbreaks)
image(t(heat.primed.vals[,pops$Type=="transition"]), axes=FALSE, col=colramp, breaks=colbreaks)
image(t(heat.primed.vals[,pops$Type=="primed"]), axes=FALSE, col=colramp, breaks=colbreaks)

# Adding the colorbar.
par(mar=c(1.1, 0.5, 1.1, 2))
plot(c(0, 10), range(colbreaks)*2, type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
axis(4, at=(-bound):bound, las=1)
for (i in seq_along(colramp)){ 
    rect(0,colbreaks[i],10,colbreaks[i+1], col=colramp[i], border=colramp[i])
}

# Adding references to the type of gene.
par(mar=c(0.5, 0.5, 0.5, 0.2))
plot(0,0,type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', xlim=c(0, 1), ylim=c(0, 1))
abline(v=1, col="grey30", lwd=5)
text(0.5, 0.5, "DE genes (transition versus naive)", srt=90, cex=1.5)

plot(0,0,type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', xlim=c(0, 1), ylim=c(0, 1))
abline(v=1, col="grey70", lwd=5)
text(0.5, 0.5, "DE genes (transition versus primed)", srt=90, cex=1.5)

# Adding references to the type of cells.
par(mar=c(0.2, 0.5, 0.5, 0.5))
plot(0,0,type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', xlim=c(0, 1), ylim=c(0, 1))
abline(h=0, col=naive.col, lwd=5)
text(0.5, 0.5, sprintf("Naive cells (%i)", sum(pops$Type=="naive")), cex=1.5)

plot(0,0,type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', xlim=c(0, 1), ylim=c(0, 1))
abline(h=0, col=trans.col, lwd=5)
text(0.5, 0.5, sprintf("Transition cells (%i)", sum(pops$Type=="transition")), cex=1.5)

plot(0,0,type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', xlim=c(0, 1), ylim=c(0, 1))
abline(h=0, col=primed.col, lwd=5)
text(0.5, 0.5, sprintf("Primed cells (%i)", sum(pops$Type=="primed")), cex=1.5)
dev.off()

###############################
# Figure 2B
pdf(file=file.path(figdir, "2b.pdf"), width = 6, height = 10, useDingbats=FALSE)
markers <- c("KLF4", "KLF17", "DPPA3", "TFCP2L1", "NANOG")

par(las=1, mar = c(2.1, 4.5, 0.1, 0.5), mfrow=c(3,1))
for (ptype in c("naive", "transition", "primed")) {
    object <- sce[,pops$Type==ptype]
    ugly.plot <- plotExpression(object, features = markers, col = "phenotype")
    plot_data <- ggplot_build(ugly.plot)
  
    pcol <- switch(ptype, naive=naive.col, primed=primed.col, transition=trans.col)      
    plot(plot_data$data[[1]]$x, plot_data$data[[1]]$y, col = pcol, 
         ylab=bquote(.(stuff)~log[2]*"-expression", list(stuff=paste0(toupper(substring(ptype, 1, 1)), substring(ptype, 2)))),
         pch=16, xlab = "" , cex.lab = 1.5, xaxt = "n", cex=2, ylim=c(0,11), cex.axis=1.2)
    if (ptype == "primed"){
        axis(1, at=seq_along(markers), labels = levels(plot_data$plot$data$Feature), cex.axis = 1.5)
    } else {
        axis(1, at=seq_along(markers), labels=character(length(markers)))
    }
}

dev.off()

###############################
# Figure S2A
pcs <- readRDS(file.path("results-overall", "pcs.rds"))

all.colors <- character(ncol(sce))
all.colors[pops$Type=="naive"] <- naive.col
all.colors[pops$Type=="transition"] <- trans.col
all.colors[pops$Type=="primed"] <- primed.col

pdf(file=file.path(figdir, "s2a.pdf"), width=9, height=8, useDingbats=FALSE)
layout(cbind(1,2), width=c(5, 1))
par(mar=c(5.1, 4.2, 4.1, 1.1))
var1 <- round(attr(pcs, "percentVar")[1]*100)
var2 <- round(attr(pcs, "percentVar")[2]*100)
plot(pcs[,1], pcs[,2], col = all.colors,
     pch=16, xlab = paste0('Component 1: ', var1, '%'), ylab = paste0('Component 2: ', var2, '%') ,
     cex.lab = 1.5 , cex.axis=1.2)
par(mar=c(5.1, 0.1, 4.1, 0.1))
plot.new()
legend("topleft", legend=c("Naive","Transition", "Primed"), col=c(naive.col, trans.col, primed.col), pch=16, cex=1.4, bty='n')
dev.off()

###############################
# Figure S2B

imprinted <- c("ENSG00000140443", # IGF1R
               "ENSG00000128739", # SNRPN
               "ENSG00000053438", # NNAT
               "ENSG00000101294", # HM13
               "ENSG00000198300", # PEG3
               "ENSG00000214548") # MEG3

heat.imprint.vals <- exprs(sce)[match(imprinted, rowData(sce)$ensembl_gene_id),]
heat.imprint.vals <- heat.imprint.vals - rowMeans(heat.imprint.vals)

pdf(file=file.path(figdir, "s2b.pdf"), width = 10, height = 5)
layout(rbind(c(9,6,7,8,9), c(5,1,2,3,4)), widths=c(0.2, 1, 1, 1, 0.2), heights=c(0.1, 1, 1))
par(mar=c(0.5, 0.5, 0.5, 0.5))
image(t(heat.imprint.vals[,pops$Type=="naive"]), axes=FALSE, col=colramp, breaks=colbreaks)
image(t(heat.imprint.vals[,pops$Type=="transition"]), axes=FALSE, col=colramp, breaks=colbreaks)
image(t(heat.imprint.vals[,pops$Type=="primed"]), axes=FALSE, col=colramp, breaks=colbreaks)

# Adding the colorbar.
par(mar=c(1.1, 0.5, 1.1, 2))
plot(c(0, 10), range(colbreaks)*2, type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
axis(4, at=(-bound):bound, las=1)
for (i in seq_along(colramp)){ 
    rect(0,colbreaks[i],10,colbreaks[i+1], col=colramp[i], border=colramp[i])
}

# Adding references to the type of gene.
par(mar=c(0.5, 0.5, 0.5, 0.2))
plot(0,0,type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', xlim=c(0, 1), ylim=c(0, 1))
coords <- par()$usr
ypts <- seq(coords[3], coords[4], length.out=length(imprinted)*2+1)
text(x=0.5, y=ypts[seq_along(imprinted)*2], rownames(heat.imprint.vals), cex=1.2)

# Adding references to the type of cells.
par(mar=c(0.2, 0.5, 0.5, 0.5))
plot(0,0,type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', xlim=c(0, 1), ylim=c(0, 1))
abline(h=0, col=naive.col, lwd=5)
text(0.5, 0.5, sprintf("Naive cells (%i)", sum(pops$Type=="naive")), cex=1.5)

plot(0,0,type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', xlim=c(0, 1), ylim=c(0, 1))
abline(h=0, col=trans.col, lwd=5)
text(0.5, 0.5, sprintf("Transition cells (%i)", sum(pops$Type=="transition")), cex=1.5)

plot(0,0,type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', xlim=c(0, 1), ylim=c(0, 1))
abline(h=0, col=primed.col, lwd=5)
text(0.5, 0.5, sprintf("Primed cells (%i)", sum(pops$Type=="primed")), cex=1.5)
dev.off()

# Figure S2C
pdf(file=file.path(figdir, "s2c.pdf"), width = 6, height = 10, useDingbats=FALSE)
markers <- c("ABCG2", "CLDN4", "VGLL1", "GATA2", "GATA3", "ERP27")

par(las=1, mar = c(2.1, 4.5, 0.1, 0.5), mfrow=c(3,1))
for (ptype in c("naive", "transition", "primed")) {
  object <- sce[,pops$Type==ptype]
  ugly.plot <- plotExpression(object, features = markers, col = "phenotype")
  plot_data <- ggplot_build(ugly.plot)
  
  pcol <- switch(ptype, naive=naive.col, primed=primed.col, transition=trans.col)      
  plot(plot_data$data[[1]]$x, plot_data$data[[1]]$y, col = pcol, 
       ylab=bquote(.(stuff)~log[2]*"-expression", list(stuff=paste0(toupper(substring(ptype, 1, 1)), substring(ptype, 2)))),
       pch=16, xlab = "" , cex.lab = 1.5, xaxt = "n", cex=2, ylim=c(0,11), cex.axis=1.2)
  if (ptype == "primed"){
    axis(1, at=seq_along(markers), labels = levels(plot_data$plot$data$Feature), cex.axis = 1.5)
  } else {
    axis(1, at=seq_along(markers), labels=character(length(markers)))
  }
}

dev.off()

# Figure S2D
topgo <- read.table(file=file.path("results-naive", "go_unique_trans.tsv"), sep="\t", header = TRUE) 
topgo <- topgo[1:10,]
pdf(file=file.path(figdir, "s2d.pdf"), width = 6, height = 10, useDingbats=FALSE)
plot <- barplot(topgo$P.DE, horiz = TRUE, col='white', ylab = 'GO Term', xlab='Logged P-Val')
x<-0.5*topgo$P.DE
text(x[8:10], plot[8:10],topgo$Term[8:10])
text((x*2)[1:7], plot[1:7], topgo$Term[1:7], pos =4)
dev.off()


