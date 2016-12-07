# Figure 1

# B
sce <- readRDS(file = "~/Documents/vMeyenn/paper_files/objects/sce_new_object")
chosen <- readRDS(file = "~/Documents/vMeyenn/paper_files/objects/chosen1")
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))

pdf(file="~/Documents/vMeyenn/paper_files/figures/figure1b.pdf", width=10)

ugly.plot <- plotPCA(sce, exprs_values="exprs", colour_by="KLF4", 
                     feature_set = chosen) + fontsize + theme(legend.title = element_text(colour="black", 
                    size=15), legend.text = element_text(size=10), legend.key.size = unit(12, "mm"))
plot_data <- ggplot_build(ugly.plot)

par(mfrow=c(1,1), mar = c(5.1, 5.1, 4.1, 6.1), las = 1)
plot(plot_data$data[[2]]$x, plot_data$data[[2]]$y, col = plot_data$data[[2]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab = plot_data$plot$labels$y, cex.lab = 1.5 , cex=1)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.95, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)

dev.off()

# C
primed.col <- readRDS(file = "~/Documents/vMeyenn/paper_files/objects/primed.col")
naive.col <- readRDS(file = "~/Documents/vMeyenn/paper_files/objects/naive.col")

res <- readRDS(file = "~/Documents/vMeyenn/paper_files/objects/res_all")
is.de <- readRDS(file = "~/Documents/vMeyenn/paper_files/objects/is.de_all")


main.col <- "grey10"

naive.genes <- c("KLF4", "KLF17", "DPPA3", "DNMT3L", "DPPA5")
primed.genes <- c("DUSP6", "THY1")

pdf(file="~/Documents/vMeyenn/paper_files/figures/figure1c.pdf")

 par(mar = c(5.1, 5.1, 4.1, 2.1), las=1)
plot(log2(rowMeans(counts(sce)[is.de==0,])), res$table$logFC[is.de==0], cex=0.3, ylim=c(-22.5, 22.5), pch = 16, 
     ylab="logFC", xlab="Mean Average", cex.lab=1.5, col=main.col)
  points(log2(rowMeans(counts(sce)[is.de==1,])), res$table$logFC[is.de==1], cex=0.3, col=primed.col, pch = 16)
  points(log2(rowMeans(counts(sce)[is.de==-1,])), res$table$logFC[is.de==-1], cex=0.3, col=naive.col, pch = 16)
  points(log2(rowMeans(counts(sce)[c(naive.genes, primed.genes),])), res$table[c(naive.genes, primed.genes),]$logFC, cex=0.8)
  legend("topright", legend=c("Upregulated in Primed Cells", "Upregulated in Naive Cells"), col=c(primed.col, naive.col), pch=16, bty = "n")
  
  text(log2(rowMeans(counts(sce)[naive.genes[c(-2,-4)],]))+1.05, (res$table[naive.genes[c(-2,-4)],]$logFC), 
       labels = naive.genes[c(-2,-4)], col=naive.col)
  
  text(log2(mean(counts(sce)["DUSP6",]))+1, (res$table["DUSP6",]$logFC), labels = "DUSP6", col=primed.col)
  text(log2(mean(counts(sce)["THY1",])), (res$table["THY1",]$logFC)-1, labels = "THY1", col=primed.col)
  text(log2(mean(counts(sce)["DNMT3L",])), (res$table["DNMT3L",]$logFC)+1, labels = "DNMT3L", col=naive.col)
  text(log2(mean(counts(sce)["KLF17",]))-0.5, (res$table["KLF17",]$logFC)-1, labels = "KLF17", col=naive.col)

dev.off()


