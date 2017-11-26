#  Functions for calculating the variance and testing for differential variability.

library(edgeR)
library(BiocParallel)

getBioDisp <- function(counts, design, is.spike, sf.cell=1, sf.spike=1, plot=FALSE, niters=20, BPPARAM=SerialParam()) {
    # Getting the dispersion estimates.
    ncells <- ncol(counts)
    sf.cell <- rep(sf.cell, length.out=ncells)
    sf.spike <- rep(sf.spike, length.out=ncells)
    original <- .getDispersion(counts, is.spike, design=design, sf.cell=sf.cell, sf.spike=sf.spike, plot=plot)
    if (niters <= 0L) {
        return(original)
    }

    # Bootstrapping and repeating the estimation.
    bootstraps <- bplapply(seq_len(niters), FUN=function(x, counts, is.spike, design, sf.cell, sf.spike) {
        ncells <- ncol(counts)
        chosen <- sample(ncells, ncells, replace=TRUE)
        .getDispersion(counts[,chosen,drop=FALSE], is.spike, design=design[chosen,,drop=FALSE], 
                       sf.cell=sf.cell[chosen], sf.spike=sf.spike[chosen], plot=FALSE)$bio
    }, counts=counts, is.spike=is.spike, design=design, sf.cell=sf.cell, sf.spike=sf.spike, BPPARAM=BPPARAM)
  
    # Estimating the precision of the dispersion from the boostraps.
    collected <- do.call(rbind, bootstraps) 
    est.var <- apply(collected, 2, var)
    original$precision <- 1/pmax(est.var, 1e-4)
    return(original)
}

.getDispersion <- function(current, is.spike, design, sf.cell, sf.spike, plot=FALSE, do.scale=FALSE) {
    sf.cell <- sf.cell/mean(sf.cell)
    sf.spike <- sf.spike/mean(sf.spike)
    cell.counts <- current[!is.spike,]
    spike.counts <- current[is.spike,]

    # Estimating for the spike-ins.
    o.spike <- makeCompressedMatrix(log(sf.spike), dim(spike.counts), byrow=TRUE)
    y.spike <- estimateDisp(spike.counts, design, offset=o.spike, prior.df=0, trend.method="none", min.row.sum=0)$tagwise.dispersion

    # Calculating the average abundance for the spike-ins.
    mean.lib <- mean(colSums(current))
    ave.spike <- aveLogCPM(spike.counts, lib.size=mean.lib*sf.spike)
 
    # Figuring out if the trend is caused by the scaling factor, assuming the trend in the NB dispersions is negligible.
    if (do.scale) {
        Y <- log2(y.spike[keep])
        A <- ave.spike[keep]
        scale <- max(1, 2^median(Y+A) - 1)
        cat(scale)
    
        spike.counts <- spike.counts/scale
        cell.counts <- cell.counts/scale
    
        # Re-estimating the spike-in parameters.
        ave.spike <- aveLogCPM(spike.counts, lib.size=mean.lib*sf.spike)
        y.spike <- estimateDisp(spike.counts, design, offset=o.spike, prior.df=0, trend.method="none", min.row.sum=0)$tagwise.dispersion
    }

    # Estimating for the cells.
    o.cell <- makeCompressedMatrix(log(sf.cell), dim(cell.counts), byrow=TRUE)
    ave.cell <- aveLogCPM(cell.counts, lib.size=mean.lib*sf.cell)
    y.cell <- estimateDisp(cell.counts, design, offset=o.cell, prior.df=0, trend.method="none", min.row.sum=0)$tagwise.dispersion   

    # Fitting a trend and making sure it looks okay.
    used <- ave.spike > aveLogCPM(5, mean.lib)
    A.spike <- ave.spike[used]
    D.spike <- log(y.spike[used])
    fit <- smooth.spline(D.spike ~ A.spike, df=4)
    if (plot) {
        plot(ave.cell, y.cell, xlab="Average log CPM", ylab="Tagwise dispersion", log="y")
        points(ave.spike, y.spike, col="red", pch=16)
        curve(exp(predict(fit, x)$y), add=TRUE, col="dodgerblue")
    }

    # Returning statistics for all rows, including spike-ins.
    total.mean <- total.disp <- numeric(nrow(current))
    total.mean[!is.spike] <- ave.cell
    total.mean[is.spike] <- ave.spike
    total.disp[!is.spike] <- y.cell
    total.disp[is.spike] <- y.spike
    inferred.tech <- exp(predict(fit, total.mean)$y)
    return(list(total=total.disp, bio=total.disp - inferred.tech, tech=inferred.tech, mean=total.mean))
}


