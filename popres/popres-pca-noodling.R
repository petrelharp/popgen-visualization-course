
###
# euro pca

source("indivinfo.R")

# read in the inner product matrix
cprod <- as.matrix(read.table("crossprod_chr8.tsv",row.names=1,check.names=FALSE))
cprod <- as.matrix(read.table("crossprod_all.tsv",row.names=1,check.names=FALSE))

# unweighted
cprod.eigen <- eigen(cprod)
stopifnot(all(rownames(cprod)==indivinfo$SUBJID))

plot.pca( cprod.eigen$vectors, pcs=1:3 )

# unweighted, normalized
nmat <- diag(nrow(cprod)) - matrix(1/nrow(cprod),nrow=nrow(cprod),ncol=ncol(cprod))

n.cprod.eigen <- eigen( (nmat %*% cprod %*% nmat) )

plot.pca( n.cprod.eigen$vectors, pcs=1:3 )


# weighted
weightings <- 1/pmax(10,sqrt( nsamples[indivinfo$COUNTRY_SELF] ))

w.cprod.eigen <- eigen( cprod * outer(weightings,weightings,"*") )

plot.pca( w.cprod.eigen$vectors, pcs=1:3 )

# weighted, normalized
nmat <- diag(nrow(cprod)) - matrix(1/nrow(cprod),nrow=nrow(cprod),ncol=ncol(cprod))

nw.cprod.eigen <- eigen( (nmat %*% cprod %*% nmat) * outer(weightings,weightings,"*") )

plot.pca( nw.cprod.eigen$vectors, pcs=1:3 )


########
## covariance
covmat <- as.matrix(read.table("crossprod_all-covariance.tsv",row.names=1,check.names=FALSE))
cov.eigen <- eigen(covmat)

plot.pca(cov.eigen$vectors, pcs=1:3)

# normalized
nmat <- diag(nrow(covmat)) - matrix(1/nrow(covmat),nrow=nrow(covmat),ncol=ncol(covmat))
ncov.eigen <- eigen(nmat%*%covmat%*%nmat)

plot.pca(ncov.eigen$vectors, pcs=1:3)

# weighted, normalized
wncov.eigen <- eigen( (nmat%*%covmat%*%nmat) * outer(weightings,weightings,"*") )

plot.pca(wncov.eigen$vectors, pcs=1:3)

# weighted, normalized-weighted
nwmat <- diag(nrow(covmat)) - weightings[col(covmat)]
wnwcov.eigen <- eigen(nmat%*%covmat%*%t(nmat) * outer(weightings,weightings,"*") )

plot.pca(wnwcov.eigen$vectors, pcs=1:3)


#######
## subset
sub.inds <- ( indivinfo$COUNTRY_SELF %in% c("Spain","France","Italy") )
sub.cprod <- cprod[sub.inds,sub.inds]
sub.nmat <- diag(nrow(sub.cprod)) - matrix(1/nrow(sub.cprod),nrow=nrow(sub.cprod),ncol=ncol(sub.cprod))
sub.eigen <- eigen( sub.nmat %*% sub.cprod %*% sub.nmat )

plot.pca( sub.eigen$vectors, pcs=1:3, info=indivinfo[sub.inds,] )

# or, cheap version:
cheap.eigen <- eigen( (nmat%*%covmat%*%nmat)[sub.inds,sub.inds] )

plot.pca( cheap.eigen$vectors, pcs=1:3, info=indivinfo[sub.inds,] )



