library(edgeR)

# With respect to the number of cells, the mean, or the dispersion.
ntries <- 10000
ncells <- c(10, 20, 50, 100, 200, 500)
disps <- c(0.05, 0.2, 0.5, 2, 5)
means <- c(0.1, 1, 10, 100, 1000)
output <- array(dim=c(length(ncells), length(disps), length(means)))

for (i in seq_along(ncells)) { 
    N <- ncells[i]
    for (j in seq_along(disps)) { 
        D <- disps[j]
        for (k in seq_along(means)) {
            M <- means[k]

            counts <- matrix(rnbinom(ntries*N, mu=M, size=1/D), ncol=N)
            design <- matrix(1, N, 1)
            y <- estimateDisp(counts, design, prior.df=0, trend="none", min.row.sum=0)$tagwise
            output[i,j,k] <- var(y)     
        }
    }   
}

# Making plots with respect to the number of cells.
library(viridis)

pdf("precision_by_cells.pdf")
colors <- viridis(length(means))
for (j in seq_along(disps)) {
    plot(1,1, type="n", xlab="Number of cells", ylab="Variance of dispersion", log="xy", 
         xlim=range(ncells), ylim=range(output), main=paste("Dispersion is", disps[j]))
    for (k in seq_along(means)) {
        current <- output[,j,k]
        lines(ncells, current, col=colors[k], lwd=2)
    }
}
dev.off()

pdf("precision_by_mean.pdf")
colors <- viridis(length(disps))
for (i in seq_along(ncells)) {
    plot(1,1, type="n", xlab="Mean", ylab="Variance of dispersion", log="xy", 
         xlim=range(means), ylim=range(output), main=paste("Number of cells is", ncells[i]))
    for (j in seq_along(disps)) {
        current <- output[i,j,]
        lines(means, current, col=colors[j], lwd=2)
    }
}
dev.off()

pdf("precision_by_disp.pdf")
colors <- viridis(length(means))
for (i in seq_along(ncells)) {
    plot(1,1, type="n", xlab="Dispersion", ylab="Variance of dispersion", log="xy", 
         xlim=range(disps), ylim=range(output), main=paste("Number of cells is", ncells[i]))
    for (k in seq_along(means)) {
        current <- output[i,,k]
        lines(disps, current, col=colors[k], lwd=2)
    }
}
dev.off()
