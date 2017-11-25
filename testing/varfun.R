#  Functions for calculating the variance and testing for differential variability.

library(edgeR)
library(limma)

VARFUN <- function(counts, is.spike, sf.cell=1, sf.spike=1, plot=FALSE) {
    bio.mat <- mean.mat <- weight.mat <- matrix(0, sum(!is.spike), length(counts))
    for (x in seq_along(counts)) { 
        current <- counts[[x]]
        design <- matrix(1, ncol(current), 1)

        # Estimating for the cells.
        cell.counts <- current[!is.spike,]
        o.cell <- makeCompressedMatrix(log(sf.cell), dim(cell.counts), byrow=TRUE)
        y.cell <- estimateDisp(cell.counts, design, offset=o.cell, prior.df=0, trend.method="none", min.row.sum=0)$tagwise.dispersion

        # Estimating for the spike-ins.
        spike.counts <- current[is.spike,]
        o.spike <- makeCompressedMatrix(log(sf.spike), dim(spike.counts), byrow=TRUE)
        y.spike <- estimateDisp(spike.counts, design, offset=o.spike, prior.df=0, trend.method="none", min.row.sum=0)$tagwise.dispersion

        # Calculating the average abundances.
        mean.lib <- mean(colSums(current))
        ave.cell <- aveLogCPM(cell.counts, lib.size=mean.lib*sf.cell)
        ave.spike <- aveLogCPM(spike.counts, lib.size=mean.lib*sf.spike)

        # Fitting a trend and making sure it looks okay.
        D.spike <- log(y.spike)
        fit <- smooth.spline(D.spike ~ ave.spike, df=4)
        if (plot) {
            plot(ave.cell, y.cell, xlab="Average log CPM", ylab="Tagwise dispersion")
            points(ave.spike, y.spike, col="red", pch=16)
            curve(exp(predict(fit, data.frame(ave.spike=x))$y[,1]), add=TRUE, col="dodgerblue")
        }

        inferred.tech <- exp(predict(fit, data.frame(ave.spike=ave.cell))$y[,1])
        bio.mat[,x] <- y.cell - inferred.tech
        mean.mat[,x] <- ave.cell
        weight.mat[,x] <- nrow(design) - ncol(design) # proportional to the precision of the variance estimate.
    }
    return(list(bio=bio.mat, mean=mean.mat, weights=weight.mat))
}

