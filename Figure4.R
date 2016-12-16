# Create figure directory and specifiy path of necessary R-objects
dir.create(file.path("./figures"), showWarnings = F)
objectpath <- file.path("./objects/")

# Figure 4a
library("viridis")
stages <- c("E3", "E4", "E4.late","E5.early", "E5", "E6", "E7")
nums <- readRDS(file = paste0(objectpath, "nums.petro"))
nums <- nums*100

pdf(file = "./figures/figure4a.pdf", width = 25, height = 12)

par(bty='n', mar = c(4.1, 4.1, 5.1, 0.2),  las = 1)
layout(matrix(1:15,ncol=5, byrow = TRUE), width = c(1,5,5,5,5), heights = c(5,5,1))
plot.new()

for (stage in stages){
  main <- stage
  if (stage == "E5"){plot.new()}
  if (stage == "E4.late"){main <- "E4 [primed]"}
  if (stage == "E5.early"){main <- "E5 [early]"}
  current_nums <- subset(nums, rownames(nums)==stage)
  cell.state <- paste0(round(current_nums[,1]), "|", round(current_nums[,2]))
  state.num <- table(cell.state)
  
  num.per.cell <- pmin(6 , state.num[cell.state])
  
  col.scale <- viridis(6)
  cell.col <- col.scale[num.per.cell]
  
  plot(current_nums, xlim = c(0,80), ylim = c(0,80), cex.axis = 2.5, yaxt="n", cex.lab = 2,
       ylab="", main = main, pch=19, cex=2, xlab = "", xaxt="n",  
       cex.main = 3, col=cell.col)
  axis(side = 1, lwd.ticks = 0, cex.axis=2)
  if (stage == "E3") {axis(side = 2, lwd.ticks = 0, cex.axis=2)}
  if (stage == "E5") {axis(side = 2, lwd.ticks = 0, cex.axis=2)}

  abline(0,1, col = "indianred2", lty=1, cex=2)
}

mtext("Proportion of naive genes [%]", side=1, padj=-2, outer=TRUE, cex = 2.5)

par(las = 0)
mtext("  Proportion of primed genes [%]", line = -2, side=2, padj=2, outer=TRUE, cex = 2.5)

legend_image <- as.raster(matrix(rev(col.scale), ncol=1))
plot(c(0,7),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=2.65, y = 1, labels = "Number of cells", cex = 3)
text(x=2.2, y = seq(0,0.8,l=6), labels = seq(1,6,l=6), cex = 3)
rasterImage(legend_image, 1, 0, 1.8, 0.8, interpolate=FALSE)

dev.off()

# Figure 4 B

stages <- readRDS(file = paste0(objectpath, "stages.mouse"))
nums <- readRDS(file =  paste0(objectpath, "nums.mouse"))
nums <- nums*100


pdf(file = "./figures/figure4b.pdf", width = 25, height = 6)

par(bty='n', mar = c(8.1, 4.1, 3.1, 0.2),  las = 1)
layout(matrix(1:7,ncol=7), width = c(1,5,5,5,5,5,2))
plot.new()

for (stage in stages){
  current_nums <- subset(nums, rownames(nums)==stage)
  cell.state <- paste0(round(current_nums[,1],1), "|", round(current_nums[,2], 1))
  state.num <- table(cell.state)
  
  num.per.cell <- pmin(10, state.num[cell.state])
  
  col.scale <- viridis(10)
  cell.col <- col.scale[num.per.cell]
    
  plot(current_nums, xlim = c(0,31), ylim = c(0,31), cex.axis = 1.5, yaxt="n", cex.lab = 2,
       ylab="", main = stage, pch=19, cex=2, xlab = "",  
       cex.main = 3, col=cell.col, xaxt="n")
  axis(side = 1, lwd.ticks = 0, cex.axis=2)
  if (stage == "E3.5") {axis(side = 2, lwd.ticks = 0, cex.axis=2)}
  abline(0,1, col = "indianred2", lty=1, cex=2)
}

mtext("Proportion of naive genes [%]", side=1, padj=-2, outer=TRUE, cex = 1.75)

par(las = 0)
mtext("  Proportion of primed genes [%]", line = -2, side=2, padj=2, outer=TRUE, cex = 1.75)

legend_image <- as.raster(matrix(rev(col.scale), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=1.2, y = 1, labels = "Number", cex = 2.5)
text(x=1.2, y = 0.925, labels = "of cells", cex = 2.5)
text(x=1.5, y = seq(0,0.8,l=4), labels = seq(1,10,l=4), cex = 2)
rasterImage(legend_image, 0, 0, 0.9,0.8, interpolate=FALSE)

dev.off()

# Figure 4 C

stages <- readRDS(file =  paste0(objectpath, "stages.monkey"))
nums <- readRDS(file =  paste0(objectpath, "nums.monkey"))
nums <- nums*100

pdf(file = "./figures/figure4c.pdf", width = 25, height = 12)
par(bty='n', mar = c(2.1, 4.1, 5.1, 0.2),  las = 1)
layout(matrix(c(1:11, 6, 12:17),ncol=6, byrow = TRUE), width = c(1,5,5,5,5,2), heights = c(5,5,1))
plot.new()

for (stage in stages){
  
  if (stage == "E13"){
    legend_image <- as.raster(matrix(rev(col.scale), ncol=1))
    plot(c(0,7),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', xaxt="n", yaxt="n")
    text(x=2.65, y = 0.975, labels = "Number", cex = 3)
    text(x=2.65, y = 0.925, labels = "of cells", cex = 3)
    text(x=4.5, y = seq(0.2,0.8,l=5), labels = seq(1,5,l=5), cex = 3.5)
    rasterImage(legend_image, 0, 0.2, 2.5, 0.8, interpolate = F)
    
    plot.new()
    }
  main <- stage
  main <- sub("0", "", main)
  current_nums <- subset(nums, rownames(nums)==stage)
  cell.state <- paste0(round(current_nums[,1],1), "|", round(current_nums[,2], 1))
  state.num <- table(cell.state)
  
  num.per.cell <- pmin(5 , state.num[cell.state])
  
  col.scale <- viridis(5)
  cell.col <- col.scale[num.per.cell]
  
  plot(current_nums, xlim = c(0,25), ylim = c(0,25), cex.axis = 2.5, yaxt="n", cex.lab = 2,
       ylab="", main = main, pch=19, cex=2, xlab = "", xaxt="n",  
       cex.main = 3, col=cell.col)
  axis(side = 1, lwd.ticks = 0, cex.axis=2)
  if (stage == "E06") {axis(side = 2, lwd.ticks = 0, cex.axis=2)}
  if (stage == "E13") {axis(side = 2, lwd.ticks = 0, cex.axis=2)}
  abline(0,1, col = "indianred2", lty=1, cex=2)
}

mtext("Proportion of naive genes [%]", side=1, padj=-2, outer=TRUE, cex = 2.5)

par(las = 0)
mtext("  Proportion of primed genes [%]", line = -2, side=2, padj=2, outer=TRUE, cex = 2.5)

dev.off()

# 4 Supplements

all_markers <- readRDS(file =  paste0(objectpath, "markers.map"))
naive_markers <- readRDS(file =  paste0(objectpath, "n.markers.map"))
primed_markers <- readRDS(file =  paste0(objectpath, "p.markers.map"))
sce_trans <- readRDS(file =  paste0(objectpath, "new_sce_trans"))

pdf(file = "./figures/sup.figure4.pdf", width = 25, height = 6)

par(bty='n', mar = c(8.1, 4.1, 3.1, 0.2),  las = 1)
layout(matrix(1:5,ncol=5), width = c(1,5,5,5,2))
plot.new()

for (object in c(sce_trans[,pData(sce_trans)$phenotype=="primed"], 
                 sce_trans[,pData(sce_trans)$phenotype=="naive"], 
                 sce_trans[,pData(sce_trans)$phenotype=="transition"])) {
  exprs_marked <- counts(object)[rownames(all_markers),]
  
  naive_num <- list()
  primed_num <- list()
  
  for (cell in colnames(counts(object))){ 
    naive_genes <- exprs_marked[which(counts(object)[rownames(naive_markers),cell] > threshold), cell]
    primed_genes <- exprs_marked[which(counts(object)[rownames(primed_markers),cell] > threshold), cell]
    
    naive_num[cell] <- length(naive_genes)/nrow(naive_markers)
    primed_num[cell] <- length(primed_genes)/nrow(primed_markers)
  } 
  if (all(object$phenotype=="primed")){
    name <- "Primed"
  } else if(all(object$phenotype=="naive")){
    name <- "Naive"
  } else if(all(object$phenotype=="transition")){
    name <- "Transition"
  } 
  
  naive_num <- matrix(matrix(unlist(naive_num)) * 100)
  primed_num <- matrix(matrix(unlist(primed_num)) * 100)
  
  cell.state <- paste0(round(naive_num), "|", round(primed_num))
  state.num <- table(cell.state)
  
  num.per.cell <- pmin(10, state.num[cell.state])
  
  col.scale <- viridis(10)
  cell.col <- col.scale[num.per.cell]
  
  plot(naive_num, primed_num, xlim = c(0,75), ylim = c(0,75), cex.axis = 1.5, yaxt=yaxt, cex.lab = 2,
       ylab="", main = name, pch=19, cex=2, xlab = "", xaxt="n", yaxt="n",  
       cex.main = 3, col=cell.col)
  axis(side = 1, lwd.ticks = 0, cex.axis=2)
  if (all(object$phenotype=="primed")) {axis(side = 2, lwd.ticks = 0, cex.axis=2)}
  abline(0,1, col = "indianred2", lty=1, cex=2)
}

mtext("Proportion of naive genes [%]", side=1, padj=-2, outer=TRUE, cex = 1.75)

par(las = 0)
mtext("  Proportion of primed genes [%]", line = -2, side=2, padj=2, outer=TRUE, cex = 1.75)

legend_image <- as.raster(matrix(rev(col.scale), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=0.725, y = 1, labels = "Number", cex = 2.5)
text(x=0.725, y = 0.925, labels = "of cells", cex = 2.5)
text(x=1.5, y = seq(0,0.8,l=4), labels = seq(1,10,l=4), cex = 2)
rasterImage(legend_image, 0, 0, 0.8,0.8, interpolate=FALSE)

dev.off()
