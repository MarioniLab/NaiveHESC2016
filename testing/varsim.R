############################################################################
# Simulation parameters

nspikes <- 100
ngenes <- 10000
spike.means <- 2^runif(nspikes, 1, 10)
tech.disp <- function(mu) { 100/mu + 0.5 }
cell.means <- 2^runif(ngenes, 1, 10)
bio.disp <- runif(ngenes, 0.1, 1)
is.spike <- rep(c(TRUE, FALSE), c(nspikes, ngenes)) # First set are spike-ins.

GENERATE_COUNTS <- function(ncells=200, disp.mult=1, mean.mult=1, scale=1) {
     spike.means <- spike.means * mean.mult
     spike.disp <- tech.disp(spike.means) * disp.mult
     spike.data <- matrix(rnbinom(nspikes*ncells, mu=spike.means, size=1/spike.disp), ncol=ncells)
     
     cell.means <- cell.means * mean.mult
     cell.disp <- tech.disp(cell.means) * disp.mult + bio.disp
     cell.data <- matrix(rnbinom(ngenes*ncells, mu=cell.means, size=1/cell.disp), ncol=ncells)
     return(rbind(spike.data, cell.data)*scale)
}

#  Testing for differential variability.

source("varfun.R")
library(limma)

VARFUN <- function(counts, is.spike, plot=FALSE) {
    bio.mat <- mean.mat <- weight.mat <- matrix(0, sum(!is.spike), length(counts))
    for (x in seq_along(counts)) { 
        current <- counts[[x]]
        design <- matrix(1, ncol(current), 1)
        output <- getBioDisp(current, is.spike=is.spike, design=design, plot=plot, BPPARAM=MulticoreParam(3), niters=0)

        bio.mat[,x] <- output$bio[!is.spike]
        mean.mat[,x] <- output$mean[!is.spike]
        weight.mat[,x] <- ncol(current) - 1L # residual d.f.
    }
    return(list(bio=bio.mat, mean=mean.mat, weights=weight.mat))
}

TESTFUN <- function(bio.mat, mean.mat, weight.mat, design, coef=ncol(design)) {
    keep <- rowSums(is.na(bio.mat))==0L
    fit <- lmFit(bio.mat[keep,], design=design, weights=weight.mat[keep,])
    fit$Amean <- rowMeans(mean.mat)[keep]
    fit <- eBayes(fit, robust=TRUE, trend=TRUE)
    return(fit$p.value[,coef])    
}

############################################################################
# Execution.

for (nbatches in c(2, 5)) { 
    for (scenario in 1:5) { 
        grouping <- gl(2,nbatches)
        design <- model.matrix(~grouping)

        # Setting up simulation parameters.   
        ncells1 <- ncells2 <- 200
        disp1 <- disp2 <- 1
        mean1 <- mean2 <- 1
        scale1 <- scale2 <- 1
    
        if (scenario==1) {
            method <- "Simple"
        } else if (scenario==2) {
            method <- "MeanShift"
            mean2 <- 5
        } else if (scenario==3) {
            method <- "DispShift"
            disp2 <- 5
        } else if (scenario==4) {
            method <- "DiffCells"
            ncells2 <- 1000
        } else if (scenario==5) {
            method <- "Scaling"
            scale2 <- 5
        }

        # Spawning simulation points.
        ncells <- c(ncells1, ncells2)[grouping]
        disp.mult <- c(disp1, disp2)[grouping]
        mean.mult <- c(mean1, mean2)[grouping]
        scaling <- c(scale1, scale2)[grouping]

        # Running through simulations 
        png(sprintf("diagnostics_%s_%i.png", method, nbatches), width=12, height=5*nbatches, units="in", res=150, pointsize=15)
        par(mfrow=c(nbatches,2), mar=c(4.1, 4.1, 0.1, 0.1))   
   
        niters <- 10
        collected.p <- vector("list", niters)
        for (it in seq_len(niters)) { 
            output <- vector('list', length(disp.mult))
            for (x in seq_along(grouping)) {
                output[[x]] <- GENERATE_COUNTS(ncells=ncells[x], 
                                               disp.mult=disp.mult[x], 
                                               mean.mult=mean.mult[x], 
                                               scale=scaling[x])
            }
            
            mats.out <- VARFUN(output, is.spike, plot=(it==1))
            pvals <- TESTFUN(mats.out$bio, mats.out$mean, mats.out$weights, design)
            collected.p[[it]] <- pvals
        }
        dev.off()
    
        # Making a quantile-log-ratio plot
        pdf(sprintf("results_%s_%i.pdf", method, nbatches))
        lower.threshold <- 10/ngenes
        ybounds <- sum(unlist(collected.p) < lower.threshold)/sum(lengths(collected.p)) / lower.threshold
        ylim <- c(pmin(0.5, ybounds/2), pmax(2, ybounds*2))

        plot(1,1, type="n", xlim=c(lower.threshold, 1), log="xy", ylim=ylim,
             ylab="Observed/expected", xlab="Specified type II error rate",
             main=method)
        for (it in seq_along(collected.p)) {
            pvals <- sort(collected.p[[it]])
            lines(pvals, (seq_along(pvals)/ngenes)/pvals, col="black")
        }
        dev.off()            
    }
}
