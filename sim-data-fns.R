# helper functions

sim_data <- function (nind, nloci) {
    indivs <- paste("I",1:nind,sep="")
    snps <- paste("SNP",1:nloci,sep="")
    locs <- cbind( lat=runif(nind), lon=runif(nind) )
    distmat <- sqrt( outer( locs[,1], locs[,1], "-" )^2 + outer( locs[,2], locs[,2], "-" )^2 )
    nugget <- 0.2
    covmat <- exp(-distmat*2) + nugget * diag(nind)
    untf.freqs <- matrix( rnorm(nind*nloci), nrow=nloci ) %*% chol(covmat)
    mean.freqs <- 5*rnorm(nloci)
    freqs <- 1/(1+exp(-untf.freqs+mean.freqs[row(untf.freqs)]))
    genotypes <- matrix( rbinom(nloci*nind,size=2,prob=freqs), nrow=nloci )
    bases <- c("A","C","G","T")
    major <- sample(1:4,size=nloci,replace=TRUE)
    minor <- 1+(((major-1)+sample(1:3,size=nloci,replace=TRUE))%%4)
    alleles <- matrix( paste(
            bases[ cbind(major,minor)[cbind(as.vector(row(genotypes)),1+as.vector(genotypes>0)) ] ],
            bases[ cbind(major,minor)[cbind(as.vector(row(genotypes)),1+as.vector(genotypes>1)) ] ],
        sep="/" ), nrow=nloci )
    colnames(alleles) <- colnames(genotypes) <- rownames(locs) <- indivs
    rownames(alleles) <- rownames(genotypes) <- snps
    return( list( locs=locs, covmat=covmat, genotypes=genotypes, alleles=alleles ) )
}
