# Create figure directory and specifiy path of necessary R-objects
dir.create(file.path("./figures"), showWarnings = F)
objectpath <- file.path("./objects/")

# Load colours
primed.col <- readRDS(file = paste0(objectpath, "primed.col"))
naive.col <- readRDS(file = paste0(objectpath, "naive.col"))
trans.col <- readRDS(file = paste0(objectpath, "trans.col"))

# Figure 2

# A

library("pheatmap")
library("gridExtra")
library("grid")

re.pheno <- readRDS(file = paste0(objectpath, "pheno"))

# lower third
shared.vals <- readRDS(file = paste0(objectpath, "heat.vals"))
shared.vals[shared.vals < -5] <- -5
shared.vals[shared.vals > 5] <- 5
shared.vals <- shared.vals[which(duplicated(rownames(shared.vals))),]

genes.shared <- rep("shared",nrow(shared.vals))
genes.shared <- factor(genes.shared, levels = "shared")
genes.shared <- as.data.frame(genes.shared)
colnames(genes.shared) <- "Gene Type"
rownames(genes.shared) <- rownames(shared.vals)

# upper two thirs
heat.vals <- readRDS(file = paste0(objectpath, "heat.vals"))
heat.vals[heat.vals < -5] <- -5
heat.vals[heat.vals > 5] <- 5

naive <- re.pheno[which(re.pheno == "naive")]
primed <- re.pheno[which(re.pheno == "primed")]
trans <- re.pheno[which(re.pheno == "transition")]

naive <- as.data.frame(naive)
rownames(naive) <- colnames(heat.vals[,which(re.pheno == "naive")])
primed <- as.data.frame(primed)
rownames(primed) <- colnames(heat.vals[,which(re.pheno == "primed")])
trans <- as.data.frame(trans)
rownames(trans) <- colnames(heat.vals[,which(re.pheno == "transition")])

gene.type <- c(rep("primed", nrow(heat.vals)/2), rep("naive", nrow(heat.vals)/2))
gene.type <- gene.type[-which(rownames(heat.vals) %in% rownames(shared.vals))]
gene.type <- factor(gene.type, levels = c( "primed", "naive"))
gene.type <- as.data.frame(gene.type)
gene.type <- as.data.frame(gene.type)

heat.vals <- heat.vals[-which(rownames(heat.vals) %in% rownames(shared.vals)),]
colnames(gene.type) <- "Transition vs"
rownames(gene.type) <- rownames(heat.vals)


# Create three new heatmaps with only the stuff additionally to this one

pdf(file="./figures/figure2a.pdf", onefile = FALSE, width = 10, height = 15)
# naive vs trans and primed vs trans
p1 <- pheatmap(heat.vals[,which(re.pheno == "naive")], color = colorRampPalette(c("navy", "white", "orangered"))(50), legend=T, 
               annotation_row = gene.type, cluster_rows = FALSE, gaps_row = length(which(gene.type=="naive")),
               show_colnames = FALSE, cluster_cols = FALSE, fontsize = 4, annotation_col =naive, show_rownames = FALSE,
               annotation_colors = list("naive" = c("naive" = naive.col), "Transition vs" = c(primed = "grey20", naive ="grey60")),
               annotation_legend = FALSE, annotation_names_row = FALSE)


p2 <- pheatmap(heat.vals[,which(re.pheno == "primed")], color = colorRampPalette(c("navy", "white", "orangered"))(50), legend=T, 
               annotation_row = gene.type, cluster_rows = FALSE, gaps_row = length(which(gene.type=="naive")), 
               show_colnames = FALSE, cluster_cols = FALSE, fontsize = 4, annotation_col = primed, show_rownames = FALSE,
               annotation_colors = list("primed" = c("primed" = primed.col), "Transition vs" = c(primed = "grey20", naive ="grey60")),
               annotation_legend = FALSE, annotation_names_row = FALSE)

p3 <-  pheatmap(heat.vals[,which(re.pheno == "transition")], color = colorRampPalette(c("navy", "white", "orangered"))(50), legend=T, 
                annotation_row = gene.type, cluster_rows = FALSE, gaps_row = length(which(gene.type=="naive")),
                show_colnames = FALSE, cluster_cols = FALSE, fontsize = 4, annotation_col = trans, show_rownames = FALSE,
                annotation_colors = list("trans" = c("transition" = trans.col), "Transition vs" = c(primed = "grey20", naive ="grey60")),
                annotation_legend = FALSE, annotation_names_row = FALSE)
# shared genes
p4 <- pheatmap(shared.vals[,which(re.pheno == "naive")], color = colorRampPalette(c("navy", "white", "orangered"))(50), legend=T, 
               annotation_row = genes.shared, cluster_rows = FALSE, fontsize_row = 11,
               show_colnames = FALSE, cluster_cols = FALSE, fontsize = 4, annotation_col = naive, show_rownames = TRUE,
               annotation_colors = list("naive" = c("naive" = naive.col), "Gene Type" = c(shared = "grey80")),
               annotation_legend = FALSE, annotation_names_row = FALSE) 

p5 <- pheatmap(shared.vals[,which(re.pheno == "primed")], color = colorRampPalette(c("navy", "white", "orangered"))(50), legend=T, 
               annotation_row = genes.shared, cluster_rows = FALSE, 
               show_colnames = FALSE, cluster_cols = FALSE, fontsize = 4, annotation_col = primed, show_rownames = FALSE,
               annotation_colors = list("primed" = c("primed" = primed.col), "Gene Type" = c(shared = "grey80")),
               annotation_legend = FALSE, annotation_names_row = FALSE)

p6 <- pheatmap(shared.vals[,which(re.pheno == "transition")], color = colorRampPalette(c("navy", "white", "orangered"))(50), legend=T, 
               annotation_row = genes.shared, cluster_rows = FALSE, 
               show_colnames = FALSE, cluster_cols = FALSE, fontsize = 4, annotation_col = trans, show_rownames = FALSE,
               annotation_colors = list("trans" = c("transition" = trans.col), "Gene Type" = c(shared = "grey80")),
               annotation_legend = FALSE, annotation_names_row = FALSE)

plot1.grob <- p1$gtable$grob[[1]] #naive
xlab1.grob <- p1$gtable$grob[[2]] 
pheno.grob <- p1$gtable$grob[[3]]  
ylab.grob <- p1$gtable$grob[[4]]  
legend.grob <- p1$gtable$grob[[5]]  

plot2.grob <- p2$gtable$grob[[1]] #primed
xlab2.grob <- p2$gtable$grob[[2]]  

plot3.grob <- p3$gtable$grob[[1]] #transition
xlab3.grob<- p3$gtable$grob[[2]]  
empty.grob <- nullGrob()

ylab2.grob <- p5$gtable$grob[[4]]  
plot4.grob <- p4$gtable$grob[[1]] #naive
plot5.grob <- p5$gtable$grob[[1]] #primed
plot6.grob <- p6$gtable$grob[[1]] #transition

grid.mat <- matrix(1,40,41)  #empty matrix
grid.mat[1,4:37] <- 2      #xtext 
grid.mat[3:27,1:2] <- 3    #ytext
grid.mat[4:38,40:41] <- 4  #legend
grid.mat[4:27,3] <- 5      #ylab
#plots
grid.mat[4:27,4:14] <- 6   #naive
grid.mat[4:27,16:26] <- 7  #trans
grid.mat[4:27,28:38] <- 8  #primed
#xlabs
grid.mat[3,4:14] <- 9      #naive
grid.mat[3,16:26] <- 10    #trans
grid.mat[3,28:38] <- 11    #primed
grid.mat[2, 4:5] <- 12     #number of naive
grid.mat[2, 16:17] <- 13   #number of transition
grid.mat[2, 28:29] <- 14   #number of primed

# shared genes
grid.mat[29:38,4:14] <- 15 #naive plot
grid.mat[29:38,16:26] <- 16 #trans plot
grid.mat[29:38,28:38] <- 17 #primed plot
grid.mat[29:38,3] <- 18    #ylab
grid.mat[29:38,1:2] <- 19  #ytext
grid.mat[29:38,39:41] <- 20  #rownames

shared.rownames <- p4$gtable$grob[[2]]
xtext <- textGrob("Phenotype", gp = gpar(cex=1.25))
ytext <- textGrob("DE genes", rot=90, gp = gpar(cex=1.5))
ytext.shared <- textGrob("Shared DE genes", rot=90, gp = gpar(cex=1.5))
xnaive <- textGrob(nrow(naive))
xprimed <- textGrob(nrow(primed))
xtrans <- textGrob(nrow(trans))

grid.arrange(grobs = list(empty.grob, xtext, ytext, legend.grob, ylab.grob, plot1.grob, 
                          plot3.grob, plot2.grob, xlab1.grob, xlab3.grob, xlab2.grob, 
                          xnaive, xtrans, xprimed, plot4.grob, plot6.grob, plot5.grob, ylab2.grob, 
                          ytext.shared, shared.rownames),
             layout_matrix = grid.mat) 

dev.off()


# Legend 2A
pdf(file="./figures/figure2a_legend.pdf", width = 10, height = 5)
plot.new()
par(mar=c(1.1,1.1,1.1,1.1))
legend("left", legend = c("Naive", "Primed", "Transition"), fill=c(naive.col, primed.col, trans.col), 
       title = "Population", bty="n", horiz = T)
legend("right", legend = c("Trans vs Naive", "Trans vs Primed", "Shared"), fill=c("grey20", "grey60", "grey80"), title = "DE Genes", bty="n", horiz = T)

dev.off()

# B
library("edgeR")
res_all <- readRDS(file = paste0(objectpath, "res_all"))
top.genes <- topTags(res_all, n=Inf)$table
top.genes <- head(rownames(top.genes[top.genes$PValue<0.05,]), n=10)
sce_trans <- readRDS(file = "~/Documents/vMeyenn/paper_files/objects/new_sce_trans")
chosen <- readRDS(file = "~/Documents/vMeyenn/paper_files/objects/chosen1")
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
sce_trans$Phenotype <- sce_trans$phenotype

pdf(file="./figures/figure2b.pdf", width = 10, height = 15)

par(las=1, mar = c(3.1, 2.1, 1.1, 0.5))
layout(matrix(c(5,5,5,1,2,3,4,4,4), ncol=3), width = c(1,9,2))
for (ptype in c("naive", "transition", "primed")) {
  object <- sce_trans[,pData(sce_trans)$phenotype==ptype]
  ugly.plot <- plotExpression(object, features = c("KLF4", "KLF17", "DPPA3", "TFCP2L1", "NANOG"), col = "phenotype")
  plot_data <- ggplot_build(ugly.plot)

  if (ptype == "primed"){plot_data$data[[1]]$colour[] <- primed.col
  } else if (ptype == "naive") {plot_data$data[[1]]$colour[] <- naive.col
  } else if (ptype == "transition") {plot_data$data[[1]]$colour[] <- trans.col}

  plot(plot_data$data[[1]]$x, plot_data$data[[1]]$y, col = plot_data$data[[1]]$colour,
     pch=16, xlab = "" , ylab ="", cex.lab = 1.5, xaxt = "n", cex=2, ylim=c(0,11))
  if (ptype == "primed"){
  axis(1, at=1:5, labels = levels(plot_data$plot$data$Feature), cex.axis = 2)}
}
plot.new()
legend("center", title = "Population",  legend = c("Naive", "Transition", "Primed"), 
       pch=16, cex=1.75, col=c(naive.col, trans.col, primed.col), bty='n')
par(las = 0)
mtext("Log expression", line = -2, side=2, padj=2, outer=TRUE, cex = 1.5)

dev.off()

# C
pdf(file="./figures/figure2c.pdf", width = 10)

ugly.plot <- plotPCA(sce_trans, exprs_values="exprs", colour_by="Phenotype", feature_set = chosen) + 
  fontsize + theme(legend.title = element_text(size=15), legend.text = element_text(size=12))
plot_data <- ggplot_build(ugly.plot)

plot_data$data[[2]]$fill[which(plot_data$data[[2]]$fill=="#FF9E4A")] <- primed.col # Primed
plot_data$data[[2]]$fill[which(plot_data$data[[2]]$fill=="#729ECE")] <- naive.col # Naive
plot_data$data[[2]]$fill[which(plot_data$data[[2]]$fill=="#67BF5C")] <- trans.col # Transition

par(mar = c(5.1, 5.1, 4.1, 6.1), las = 1)
plot(plot_data$data[[2]]$x, plot_data$data[[2]]$y, col = plot_data$data[[2]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab = plot_data$plot$labels$y, cex.lab = 1.5)
par(xpd=TRUE)
legend(23.7,10.37, title = "Population",  legend = c("Primed", "Naive", "Transition"), pch=16, cex=1.25, col=c(primed.col, naive.col, trans.col), bty='n')

dev.off()
