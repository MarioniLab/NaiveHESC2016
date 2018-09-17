source("Figures/central.R")
remapdir <- "analysis/results-remapping"

# Make mapping plots.
mapper <- function(x, y, colors, width=0.01, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, ...) {
    x.id <- ceiling(x/width)
    y.id <- ceiling(y/width)
    ids <- paste0(x.id, "|", y.id)
    num <- table(ids)    
    x2 <- as.integer(sub("\\|.*", "", names(num)))
    y2 <- as.integer(sub(".*\\|", "", names(num)))

    if (is.null(xlim)) xlim <- c(0, max(x))
    if (is.null(ylim)) ylim <- c(0, max(y))
    if (is.null(xlab)) xlab <- "Proportion of naive markers"
    if (is.null(ylab)) ylab <- "Proportion of primed markers"
    plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)

    num <- pmin(num, length(colors))
    col <- colors[num]
    adj <- width * 0.1
    rect((x2-1)*width + adj,
         (y2-1)*width + adj, 
         x2*width - adj, 
         y2*width -adj, col=col, border=NA)
    
    abline(0, 1, col="red", lwd=2, lty=2)
}

# Make other plots
plotTheRest <- function(colors, skiprow=TRUE) {
    par(mar=c(2.1, 0.1, 4.1, 1.5)) # Colorbar
    plot(c(0, 20), c(0, 15), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
    for (i in seq_along(colors)){ 
        rect(5,i,10,i+1, col=colors[i], border=colors[i])
        text(10, i+0.5, ifelse(i==length(colors), paste0(i, "+"), i), pos=4, cex=2)
    }
    text(10, length(colors)+1, "# cells", cex=2.5, pos=3)
    
    if (skiprow) plot.new()
    
    par(mar=c(2.1, 1.5, 4.1, 0.1), xpd=TRUE) 
    plot(0, 0, bty="n", type="n", xaxt="n", yaxt="n", xlab="", ylab="", cex.lab=1.5)
    text(0,0, "Proportion of primed markers", srt=90, cex=3)
    
    par(mar=c(2.1, 2.1, 0.1, 2.1), xpd=TRUE) 
    plot(0, 0, bty="n", type="n", xaxt="n", yaxt="n", xlab="", ylab="", cex.lab=1.5)
    text(0,0, "Proportion of naive markers", cex=3)
}

library(viridis)



# Figure 4A

blah <- read.table(file.path(remapdir, "mouse_coords.tsv"), header=TRUE)
colors <- viridis(5)

mat <- matrix(1:6, nrow=2, byrow=TRUE)
n <- length(mat)
final.mat <- cbind(n+3, mat, c(n+1, n+2))
final.mat <- rbind(final.mat, n+4)
lwidths <- c(0.2, rep(1, ncol(mat)), 0.3)
lheights <- c(rep(1, nrow(mat)), 0.2)

pdf(file.path(figdir, "4a.pdf"), width=sum(lwidths)*5, height=sum(lheights)*5)
par(mar=c(2.1, 2.1, 4.1, 2.1))
layout(final.mat, width=lwidths, height=lheights)

i <- 1L
for (stage in c("E3.5", "E4.5", "E5.5", "E6.5", "E6.75")){ 
    current <- blah$stage==stage
    x <- blah$naive
    y <- blah$primed
    mapper(x[current], y[current], colors, xlim=c(0, 0.3), ylim=c(0, 0.3), xaxt="n", yaxt="n",
           main=stage, cex.main=2.5, cex.axis=1.4) 
    axis(2, cex.axis=1.7, labels=(i==mat[1,1] || i==mat[2,1]))
    axis(1, cex.axis=1.7, labels=(i %in% mat[2,] || i %in% mat[,3]))
    i <- i+1L
}

plot.new()
plotTheRest(colors)
dev.off()

# Figure 4B

blah <- read.table(file.path(remapdir, "monkey_coords.tsv"), header=TRUE)
colors <- viridis(5)

mat <- matrix(1:8, nrow=2, byrow=TRUE)
n <- length(mat)
final.mat <- cbind(n+3, mat, c(n+1, n+2))
final.mat <- rbind(final.mat, n+4)
lwidths <- c(0.2, rep(1, ncol(mat)), 0.3)
lheights <- c(rep(1, nrow(mat)), 0.2)

pdf(file.path(figdir, "4b.pdf"), width=sum(lwidths)*5, height=sum(lheights)*5)
par(mar=c(2.1, 2.1, 4.1, 2.1))
layout(final.mat, width=lwidths, height=lheights)

i <- 1L
for (stage in c("E06", "E07", "E08", "E09", "E13", "E14", "E16", "E17")) {
    current <- blah$stage==stage
    x <- blah$naive
    y <- blah$primed
    mapper(x[current], y[current], colors, xlim=c(0, 0.5), ylim=c(0, 0.5), xaxt="n", yaxt="n",
           main=sub("E0", "E", stage), cex.main=2.5, cex.axis=1.4, width=0.02) 
    axis(2, cex.axis=1.7, labels=(i==mat[1,1] || i==mat[2,1]))
    axis(1, cex.axis=1.7, labels=(i %in% mat[2,]))
    i <- i + 1L
}

plotTheRest(colors)
dev.off()

# Figure 4C

blah <- read.table(file.path(remapdir, "petro_coords.tsv"), header=TRUE)
colors <- viridis(8)

mat <- matrix(1:8, nrow=2, byrow=TRUE)
n <- length(mat)
final.mat <- cbind(n+3, mat, c(n+1, n+2))
final.mat <- rbind(final.mat, n+4)
lwidths <- c(0.2, rep(1, ncol(mat)), 0.3)
lheights <- c(rep(1, nrow(mat)), 0.2)

pdf(file.path(figdir, "4c.pdf"), width=sum(lwidths)*5, height=sum(lheights)*5)
par(mar=c(2.1, 2.1, 4.1, 2.1))
layout(final.mat, width=lwidths, height=lheights)

i <- 1L
for (stage in c("E3", "E4", "E4.late","E5.early", "E5", "E6", "E7")) {
  current <- blah$stage==stage
  x <- blah$naive
  y <- blah$primed
  mapper(x[current], y[current], colors, xlim=c(0, 0.8), ylim=c(0, 0.8), xaxt="n", yaxt="n",
         main=sub("\\.", " ", stage), cex.main=2.5, width=0.02) 
  axis(2, cex.axis=1.7, labels=(i==mat[1,1] || i==mat[2,1]))
  axis(1, cex.axis=1.7, labels=(i %in% mat[2,] || i %in% mat[,4]))
  i <- i+1L
}

plot.new()
plotTheRest(colors)
dev.off()

# 4 Supplements

blah <- read.table(file.path(remapdir, "own_coords.tsv"), header=TRUE)
colors <- viridis(10)

mat <- matrix(1:3, nrow=1, byrow=TRUE)
n <- length(mat)
final.mat <- cbind(n+2, mat, n+1)
final.mat <- rbind(final.mat, n+3)
lwidths <- c(0.2, rep(1, ncol(mat)), 0.3)
lheights <- c(rep(1, nrow(mat)), 0.2)

pdf(file.path(figdir, "s4a.pdf"), width=sum(lwidths)*5, height=sum(lheights)*5)
par(mar=c(2.1, 2.1, 4.1, 2.1))
layout(final.mat, width=lwidths, height=lheights)

i <- 1L
for (type in c("naive", "transition", "primed")) {
    current <- blah$type==type
    x <- blah$naive
    y <- blah$primed
    mapper(x[current], y[current], colors, xlim=c(0, 0.6), ylim=c(0, 0.6), yaxt='n', 
           main=paste0(toupper(substring(type, 1, 1)), substring(type, 2)),
           cex.main=2.5, cex.axis=1.4, width=0.02) 
    axis(2, cex.axis=1.7, labels=(i==mat[1,1]))
    i <- i + 1L
}

plotTheRest(colors, skiprow=FALSE)
dev.off()

# Naive/transition markers
# Adjust functions for transition 
mapper_trans <- function(x, y, colors, width=0.01, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, ...) {
  x.id <- ceiling(x/width)
  y.id <- ceiling(y/width)
  ids <- paste0(x.id, "|", y.id)
  num <- table(ids)    
  x2 <- as.integer(sub("\\|.*", "", names(num)))
  y2 <- as.integer(sub(".*\\|", "", names(num)))
  
  if (is.null(xlim)) xlim <- c(0, max(x))
  if (is.null(ylim)) ylim <- c(0, max(y))
  if (is.null(xlab)) xlab <- "Proportion of transition markers"
  if (is.null(ylab)) ylab <- "Proportion of naive markers"
  plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  
  num <- pmin(num, length(colors))
  col <- colors[num]
  adj <- width * 0.1
  rect((x2-1)*width + adj,
       (y2-1)*width + adj, 
       x2*width - adj, 
       y2*width -adj, col=col, border=NA)
  
  abline(0, 1, col="red", lwd=2, lty=2)
}

# Make other plots
plotTheRest_trans <- function(colors, skiprow=TRUE) {
  par(mar=c(2.1, 0.1, 4.1, 1.5)) # Colorbar
  plot(c(0, 20), c(0, 15), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
  for (i in seq_along(colors)){ 
    rect(5,i,10,i+1, col=colors[i], border=colors[i])
    text(10, i+0.5, ifelse(i==length(colors), paste0(i, "+"), i), pos=4, cex=2)
  }
  text(10, length(colors)+1, "# cells", cex=2.5, pos=3)
  
  if (skiprow) plot.new()
  
  par(mar=c(2.1, 1.5, 4.1, 0.1), xpd=TRUE) 
  plot(0, 0, bty="n", type="n", xaxt="n", yaxt="n", xlab="", ylab="", cex.lab=1.5)
  text(0,0, "Proportion of naive markers", srt=90, cex=3)
  
  par(mar=c(2.1, 2.1, 0.1, 2.1), xpd=TRUE) 
  plot(0, 0, bty="n", type="n", xaxt="n", yaxt="n", xlab="", ylab="", cex.lab=1.5)
  text(0,0, "Proportion of transition markers", cex=3)
}

#S4B

blah <- read.table(file.path(remapdir, "petro_coords_trans.tsv"), header=TRUE)
colors <- viridis(8)

mat <- matrix(1:8, nrow=2, byrow=TRUE)
n <- length(mat)
final.mat <- cbind(n+3, mat, c(n+1, n+2))
final.mat <- rbind(final.mat, n+4)
lwidths <- c(0.2, rep(1, ncol(mat)), 0.3)
lheights <- c(rep(1, nrow(mat)), 0.2)

pdf(file.path(figdir, "s4b.pdf"), width=sum(lwidths)*5, height=sum(lheights)*5)
par(mar=c(2.1, 2.1, 4.1, 2.1))
layout(final.mat, width=lwidths, height=lheights)

i <- 1L
for (stage in c("E3", "E4", "E4.late","E5.early", "E5", "E6", "E7")) {
  current <- blah$stage==stage
  x <- blah$transition
  y <- blah$naive
  mapper_trans(x[current], y[current], colors, xlim=c(0, 0.8), ylim=c(0, 0.8), xaxt="n", yaxt="n",
               main=sub("\\.", " ", stage), cex.main=2.5, width=0.02) 
  axis(2, cex.axis=1.7, labels=(i==mat[1,1] || i==mat[2,1]))
  axis(1, cex.axis=1.7, labels=(i %in% mat[2,] || i %in% mat[,4]))
  i <- i+1L
}

plot.new()
plotTheRest_trans(colors)
dev.off()


