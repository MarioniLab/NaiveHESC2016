source("central.R")


# (A) should be three PCA (or t-SNE) plots for the naive condition. Each PCA plot should be constructed only with lineage markers; 
#the PCA plots are the same, except that they are coloured by different lineage markers, all of which should show no structure. 
#(B) should be the same as (A) except for the primed condition. 
#(C) should contain two heatmaps for the naive condition, for showing significant marker/reader and marker/writer correlations as discussed above. 
#(D) should contain the same type of heatmaps for the primed condition, showing the exact same genes.

# Figure 3A
library(scater)
library(ggplot2)
library(grid)

tsne1 <- readRDS("naive_tsne1.Rds")
tsne2 <- readRDS("naive_tsne2.Rds")
tsne3 <- readRDS("naive_tsne3.Rds")

pdf(file=file.path(figdir, "3a.pdf"), height=5, width=10, useDingbats=FALSE)

par(mfrow=c(1,3), mar = c(5.1, 5.1, 4.1, 6.1), las = 1)
plot_data <- ggplot_build(tsne1)
plot(plot_data$data[[2]]$x, plot_data$data[[2]]$y, col = plot_data$data[[2]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab = plot_data$plot$labels$y, cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.295, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)

par(mar = c(5.1, 4.1, 4.1, 6.1), las = 1)
plot_data <- ggplot_build(tsne2)
plot(plot_data$data[[2]]$x, plot_data$data[[2]]$y, main = 'naive', col = plot_data$data[[2]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab ='', cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.62, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)

par(mar = c(5.1, 4.1, 4.1, 6.1), las = 1)
plot_data <- ggplot_build(tsne3)
plot(plot_data$data[[2]]$x, plot_data$data[[2]]$y, col = plot_data$data[[2]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab = '', cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.96, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)

dev.off()

# Figure 3B
tsne1 <- readRDS("primed_tsne1.Rds")
tsne2 <- readRDS("primed_tsne2.Rds")
tsne3 <- readRDS("primed_tsne3.Rds")

pdf(file=file.path(figdir, "3b.pdf"), height=5, width=10, useDingbats=FALSE)

par(mfrow=c(1,3), mar = c(5.1, 5.1, 4.1, 6.1), las = 1)
plot_data <- ggplot_build(tsne1)
plot(plot_data$data[[2]]$x, plot_data$data[[2]]$y, col = plot_data$data[[2]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab = plot_data$plot$labels$y, cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.295, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)

par(mar = c(5.1, 4.1, 4.1, 6.1), las = 1)
plot_data <- ggplot_build(tsne2)
plot(plot_data$data[[2]]$x, plot_data$data[[2]]$y, main = 'primed', col = plot_data$data[[2]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab ='', cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.62, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)

par(mar = c(5.1, 4.1, 4.1, 6.1), las = 1)
plot_data <- ggplot_build(tsne3)
plot(plot_data$data[[2]]$x, plot_data$data[[2]]$y, col = plot_data$data[[2]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab ='', cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.96, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)

dev.off()

# 3C

library(pheatmap)
library(gridExtra)

naive.epi.heat <- readRDS('naive.epi.Rds')

heat.dists <- dist(current_mat1)
heat.tree <- hclust(heat.dists)

heat.dist2 <- dist(t(current_mat1))
heat.tree2 <- hclust(heat.dist2)

naive.epi.heat[which(naive.epi.heat>0.5, arr.ind = TRUE)] <- 0.5
naive.epi.heat[which(naive.epi.heat< (-0.5), arr.ind = TRUE)] <- -0.5

pdf(file.path(figdir, "3c.pdf"), width=5, height = 5, onefile=FALSE)
pheatmap(naive.epi.heat, breaks=seq(-0.5, 0.5, length.out=101), main = 'naive', color = colorRampPalette(c("navy", "white", "orangered"))(101),
         treeheight_row = 0, treeheight_col = 0, cluster_rows = heat.tree, cluster_cols = heat.tree2)
dev.off()


# 3D

primed.epi.heat <- readRDS('primed.epi.Rds')

primed.epi.heat[which(primed.epi.heat>0.5, arr.ind = TRUE)] <- 0.5
primed.epi.heat[which(primed.epi.heat< (-0.5), arr.ind = TRUE)] <- -0.5   ##to save each plot into a list. note the [[4]]

pdf(file.path(figdir, "3d.pdf"), width=5, height = 5, onefile=FALSE)
pheatmap(primed.epi.heat, breaks=seq(-0.5, 0.5, length.out=101), main = 'primed', color = colorRampPalette(c("navy", "white", "orangered"))(101),
         treeheight_row = 0, treeheight_col = 0, cluster_rows = heat.tree, cluster_cols = heat.tree2)
dev.off()

### Supplements

# Figure S3A

naive.plot <- readRDS("naive_plot.Rds")
primed.plot <- readRDS("primed_plot.rds")

pdf(file=file.path(figdir, "s3a.pdf"), width=10, useDingbats=FALSE)
par(mfrow=c(1,2), mar = c(5.1, 5.1, 4.1, 4.1), las = 1)
plot_data <- ggplot_build(naive.plot)
plot(plot_data$data[[2]]$x, plot_data$data[[2]]$y, main = 'naive', col = plot_data$data[[2]]$fill, 
           pch=16, xlab = plot_data$plot$labels$x , ylab = plot_data$plot$labels$y, cex.lab = 1.5 , cex=1.5)
par(mar = c(5.1, 4.1, 4.1, 5.1), las = 1)
plot_data <- ggplot_build(primed.plot)
plot(plot_data$data[[2]]$x, plot_data$data[[2]]$y, main = 'primed', col = plot_data$data[[2]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab = plot_data$plot$labels$y, cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.95, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)
dev.off()


# Figure S3B
colour.code <- readRDS('fdr.corr.Rds')

colour.scheme <- c('grey', 'black', naive.col, primed.col)
colour.code <- colour.scheme[colour.code]

pdf(file.path(figdir, "s3b.pdf"), width=10, height = 10, onefile=FALSE)
plot(current_mat1, current_mat2, xlim = c(-0.5, 0.5) , ylim = c(-0.5, 0.5), xlab = 'Naive correlations', ylab='Primed correlations', pch=16, col = colour.code)
legend('topleft', legend = c('Total: 625', 'Sig in both: 6', 'Sig in naive: 214', 'Sig in primed: 1'), pch = 16, col = c('grey', 'black', naive.col, primed.col), bty='n', cex=1.5)
dev.off()