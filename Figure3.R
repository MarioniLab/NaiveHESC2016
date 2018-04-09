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
     pch=16, xlab = plot_data$plot$labels$x , ylab = plot_data$plot$labels$y, cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.62, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)

par(mar = c(5.1, 4.1, 4.1, 6.1), las = 1)
plot_data <- ggplot_build(tsne3)
plot(plot_data$data[[2]]$x, plot_data$data[[2]]$y, col = plot_data$data[[2]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab = plot_data$plot$labels$y, cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.96, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)

dev.off()


# Figure 3C
library(pheatmap)
library(gridExtra)

naive.read.heat <- readRDS('naive.read.Rds')
naive.mod.heat <- readRDS('naive.mod.Rds')

naive.mod.heat[which(naive.mod.heat>0.5, arr.ind = TRUE)] <- 0.5
naive.mod.heat[which(naive.mod.heat< (-0.5), arr.ind = TRUE)] <- -0.5
selection <- tail(order(apply(naive.mod.heat^2,1, mean, na.rm = TRUE)), n=25)
naive.mod.heat <- naive.mod.heat[selection,]

naive.read.heat[which(naive.read.heat>0.5, arr.ind = TRUE)] <- 0.5
naive.read.heat[which(naive.read.heat< (-0.5), arr.ind = TRUE)] <- -0.5   ##to save each plot into a list. note the [[4]]
selection <- tail(order(apply(naive.read.heat^2,1, mean, na.rm = TRUE)), n=25)
naive.read.heat <- naive.read.heat[selection,]

pdf(file.path(figdir, "3c1.pdf"), width=5, height = 5, onefile=FALSE)
pheatmap(naive.mod.heat, breaks=seq(-0.5, 0.5, length.out=101), color = colorRampPalette(c("navy", "white", "orangered"))(101),
         treeheight_row = 0, treeheight_col = 0)
dev.off()
pdf(file.path(figdir, "3c2.pdf"), width=5, height = 5, onefile=FALSE)
pheatmap(naive.read.heat, breaks=seq(-0.5, 0.5, length.out=101), color = colorRampPalette(c("navy", "white", "orangered"))(101),
         treeheight_row = 0, treeheight_col = 0)
dev.off()

# Figure 3D

primed.mod.heat <- readRDS('primed.mod.Rds')
primed.read.heat <- readRDS('primed.read.Rds')

primed.mod.heat[which(primed.mod.heat>0.5, arr.ind = TRUE)] <- 0.5
primed.mod.heat[which(primed.mod.heat< (-0.5), arr.ind = TRUE)] <- -0.5
selection <- tail(order(apply(primed.mod.heat^2,1, mean, na.rm = TRUE)), n=25)
primed.mod.heat <- primed.mod.heat[selection,]

primed.read.heat[which(primed.read.heat>0.5, arr.ind = TRUE)] <- 0.5
primed.read.heat[which(primed.read.heat< (-0.5), arr.ind = TRUE)] <- -0.5
selection <- tail(order(apply(primed.read.heat^2,1, mean, na.rm = TRUE)), n=25)
primed.read.heat <- primed.read.heat[selection,]

pdf(file.path(figdir, "3d1.pdf"), width=5, height = 5, onefile=FALSE)
pheatmap(primed.mod.heat, breaks=seq(-0.5, 0.5, length.out=101), color = colorRampPalette(c("navy", "white", "orangered"))(101),
         treeheight_row = 0, treeheight_col = 0)
dev.off()
pdf(file.path(figdir, "3d2.pdf"), width=5, height = 5, onefile=FALSE)
pheatmap(primed.read.heat, breaks=seq(-0.5, 0.5, length.out=101), color = colorRampPalette(c("navy", "white", "orangered"))(101),
         treeheight_row = 0, treeheight_col = 0)
dev.off()

### Supplements

# Figure 3A

naive.plot <- readRDS("naive_plot.Rds")
primed.plot <- readRDS("primed_plot.rds")

pdf(file=file.path(figdir, "s3b.pdf"), width=10, useDingbats=FALSE)
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
