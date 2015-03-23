# helper functions

sim_data <- function (nind, nloci, do.alleles=TRUE,
        locs=cbind( lat=runif(nind), lon=runif(nind) ),
        env=ifelse( locs[,1] > 0.5, 1, 0 ),
        a0=0.7,aD=1,aE=0,a2=1,nugget=0.2) {
    indivs <- paste("I",1:nind,sep="")
    snps <- paste("SNP",1:nloci,sep="")
    D <- sqrt( outer( locs[,1], locs[,1], "-" )^2 + outer( locs[,2], locs[,2], "-" )^2 )
    E <- abs( outer( env, env, "-" ) )
    covmat <- exp(-sqrt(aD*D^2+aE*E^2)^a2)/a0 + nugget * diag(nind)
    untf.freqs <- matrix( rnorm(nind*nloci), nrow=nloci ) %*% chol(covmat)
    mean.freqs <- 5*rnorm(nloci)
    freqs <- 1/(1+exp(-untf.freqs+mean.freqs[row(untf.freqs)]))
    genotypes <- matrix( rbinom(nloci*nind,size=2,prob=freqs), nrow=nloci )
    bases <- c("A","C","G","T")
    major <- sample(1:4,size=nloci,replace=TRUE)
    minor <- 1+(((major-1)+sample(1:3,size=nloci,replace=TRUE))%%4)
    alleles <- if (do.alleles) { matrix( paste(
            bases[ cbind(major,minor)[cbind(as.vector(row(genotypes)),1+as.vector(genotypes>0)) ] ],
            bases[ cbind(major,minor)[cbind(as.vector(row(genotypes)),1+as.vector(genotypes>1)) ] ],
        sep="/" ), nrow=nloci ) } else { NULL }
    colnames(alleles) <- colnames(genotypes) <- rownames(locs) <- indivs
    rownames(alleles) <- rownames(genotypes) <- snps
    return( list( locs=locs, env=env, covmat=covmat, genotypes=genotypes, alleles=alleles, 
            params=list(a0=a0,aD=aD,aE=aE,a2=a2,nugget=nugget) ) )
}
