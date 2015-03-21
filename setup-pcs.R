

tfam <- read.table("tfam-info.txt",header=TRUE)
sample.info <- read.table("sample-info",header=TRUE,sep="\t")
omit.labels <- c("Australian","Canadian","Central_European","White_/_Caucasian","European","Japanese","West_Indian","Fijian","Trinidad_and_Tobago","South_African","Seychellois","Not_Specified","missing","Filipino","Grenadian","African_American")
sample.info$usethese <- with(sample.info, WORLD_LABEL %in% c("European","Southern_Asia","Western_Asia","Other") & ( ! ORIGIN_LABEL %in% omit.labels ) )

world.samps <- unlist( tapply( which(sample.info$usethese), sample.info$WORLD_LABEL[sample.info$usethese], function (x) { if (length(x)==0) { NULL } else { sample(x,min(length(x),200)) } } ) )
sample.info$use.world <- ( 1:nrow(sample.info) %in% world.samps )

write.table(subset(sample.info,use.world),file="world-info.tsv",row.names=FALSE, sep="\t")

cov.all <- read.table("cov_for_whole_genome.txt",header=TRUE)
cov.chr8 <- read.table("cov_data_for_chr8.txt",header=TRUE)

cov.all.world <- cov.all[sample.info$use.world,sample.info$use.world]
write.table(cov.all.world,file="cov_world_whole_genome.tsv",row.names=FALSE,sep="\t")

cov.chr8.world <- cov.chr8[sample.info$use.world,sample.info$use.world]
write.table(cov.chr8.world,file="cov_world_chr8.tsv",row.names=FALSE,sep="\t")


## Euro samples
# To trim down to snps in pairwise r2 no more than .25:
#     for chr in $(seq 22); do plink --file ../dbgap_chr/POPRES_chr$chr --indep-pairwise 50 5 0.25 --out POPRES_chr${chr}_prunesnps; done
# To turn a .ped file into what we want with plink, use  e.g.
#     for chr in $(seq 22); do plink --file ../dbgap_chr/POPRES_chr$chr --extract POPRES_chr${chr}_prunesnps.prune.in --recodeA --out POPRES_chr${chr}_pruned; done
#  which produces POPRES_chrXX_prunesnps.raw
#
# IN: peter@phoebe:~/projects/ibd/data/POPRES/clean_pca

txp <- function(x,y=NULL) {
    x[is.na(x)] <- 0
    if (!is.null(y)) { y[is.na(y)] <- 0 }
    tcrossprod(x,y)
}

indivinfo <- read.csv("indivinfo.csv",header=TRUE,stringsAsFactors=FALSE)
chr.files <- list.files(".",pattern="POPRES_chr.*_pruned.raw")

require(parallel)

chroms <- 1:22
# for (chrom in chroms) {
mclapply( chroms, function (chrom) {
    # inner products
        cat("chrom ", chrom, "\n")
        filename <- paste("POPRES_chr",chrom,"_pruned.raw",sep="")
        outbase <- paste("crossprod_chr",chrom,sep="")
        genotypes <- read.table(filename,header=TRUE)
        genotypes <- genotypes[match(indivinfo$SUBJID,genotypes$IID),]
        indivs <- genotypes[,1:6]
        genotypes <- as.matrix( genotypes[7:ncol(genotypes)] )
        # denominator
        nshared <- tcrossprod( !is.na(genotypes) )
        # G G^T  with missings removed
        iprod <- txp( genotypes )
        rownames(iprod) <- colnames(iprod) <- indivinfo$SUBJID
        rownames(nshared) <- colnames(nshared) <- indivinfo$SUBJID
        write.csv(iprod,file=paste(outbase,"-inner_prod.csv",sep=''),row.names=TRUE)
        write.csv(nshared,file=paste(outbase,"-nonmissing.csv",sep=''),row.names=TRUE)
        if (FALSE) {
            # svd
            iprod.eigen <- eigen(iprod/nshared)
            # normalization
            nmat <- diag(nrow(iprod)) - matrix(1/nrow(iprod),nrow=nrow(iprod),ncol=nrow(iprod))
            niprod.eigen <- eigen( nmat %*% (iprod/nshared) %*% nmat )
            write.csv(niprod.eigen$vectors[,1:10], file="niprod-eigen.csv",row.names=TRUE)
        }
    } , mc.cores=detectCores() )


##  combined

fnames <- function (chrom) { 
    outbase <- paste("crossprod_chr",chrom,sep="")
    c( ibase=paste(outbase,"-inner_prod.csv",sep=''), nbase=paste(outbase,"-nonmissing.csv",sep='') )
} 

chrom <- 1
odd.iprod <- even.iprod <- iprod <- read.csv(fnames(chrom)[1],row.names=1,check.names=FALSE)
even.iprod[] <- 0
odd.nshared <- even.nshared <- nshared <- read.csv(fnames(chrom)[2],row.names=1,check.names=FALSE)
even.nshared[] <- 0
for ( chrom in setdiff(chroms,1) ) {
    cat("chrom ", chrom, "\n")
    new.iprod <- read.csv(fnames(chrom)[1],row.names=1,check.names=FALSE)
    iprod <- iprod + new.iprod
    new.nshared <- read.csv(fnames(chrom)[2],row.names=1,check.names=FALSE)
    nshared <- nshared + new.nshared
    if ( (chrom%%2)==1 ) {
        odd.iprod <- odd.iprod + new.iprod
        odd.nshared <- odd.nshared + new.nshared
    } else {
    }
}
write.csv(iprod,file="crossprod_all-inner_prod.csv",row.names=TRUE)
write.csv(nshared,file="crossprod_all-nonmissing.csv",row.names=TRUE)
write.csv(iprod/nshared,file="crossprod_all.csv",row.names=TRUE)
write.csv(iprod,file="crossprod_odd-inner_prod.csv",row.names=TRUE)
write.csv(nshared,file="crossprod_odd-nonmissing.csv",row.names=TRUE)
write.csv(iprod/nshared,file="crossprod_odd",row.names=TRUE)
write.csv(iprod,file="crossprod_even-inner_prod.csv",row.names=TRUE)
write.csv(nshared,file="crossprod_even-nonmissing.csv",row.names=TRUE)
write.csv(iprod/nshared,file="crossprod_even.csv",row.names=TRUE)


# covariances
# for (chrom in chroms) {
mclapply( chroms, function (chrom) {
        cat("chrom ", chrom, "\n")
        filename <- paste("POPRES_chr",chrom,"_pruned.raw",sep="")
        outbase <- paste("crossprod_chr",chrom,sep="")
        genotypes <- read.table(filename,header=TRUE)
        genotypes <- genotypes[match(indivinfo$SUBJID,indivs$IID),]
        indivs <- genotypes[,1:6]
        genotypes <- as.matrix( genotypes[7:ncol(genotypes)] )
        covmat <- cov(t(genotypes),use="pairwise")
        write.csv(covmat,file=paste(outbase,"-covariance.csv",sep=''),row.names=TRUE)
    } , mc.cores=detectCores() )


# combined

cname <- function (chrom) { 
    outbase <- paste("crossprod_chr",chrom,sep="")
    paste(outbase,"-covariance.csv",sep='')
} 

chrom <- 1
covmat <- read.csv(cname(chrom),row.names=1,check.names=FALSE)
for ( chrom in setdiff(chroms,1) ) {
    cat("chrom ", chrom, "\n")
    covmat <- covmat + read.csv(cname(chrom),row.names=1,check.names=FALSE)
}
covmat <- covmat/length(chroms)
write.csv(covmat,file="crossprod_all-covariance.csv",row.names=TRUE)

# relabel data for class
indivinfo$ORIG_ID <- indivinfo$SUBJID
indivinfo$SUBJID <- sample(1:nrow(indivinfo))
write.table(indivinfo[,1:2],file="indivinfo.tsv",sep="\t",row.names=FALSE)
cprod <- as.matrix(read.csv("orig_data/crossprod_all.csv",row.names=1,check.names=FALSE))
rownames(cprod) <- colnames(cprod) <- indivinfo$SUBJID[match(indivinfo$ORIG_ID,rownames(cprod))]
write.table(cprod,file="crossprod_all.csv",sep="\t",row.names=TRUE)

cprod <- as.matrix(read.csv("orig_data/crossprod_chr8-inner_prod.csv",row.names=1,check.names=FALSE)) / as.matrix(read.csv("orig_data/crossprod_chr8-nonmissing.csv",row.names=1,check.names=FALSE))
rownames(cprod) <- colnames(cprod) <- indivinfo$SUBJID[match(indivinfo$ORIG_ID,rownames(cprod))]
write.table(cprod,file="crossprod_chr8.csv",sep="\t",row.names=TRUE)
