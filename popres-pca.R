
###
# euro pca

# relabel data for class
if (FALSE) {
    indivinfo$ORIG_ID <- indivinfo$SUBJID
    indivinfo$SUBJID <- sample(1:nrow(indivinfo))
    write.table(indivinfo[,1:2],file="indivinfo.tsv",sep="\t",row.names=FALSE)
    cprod <- as.matrix(read.csv("orig_data/crossprod_all.csv",row.names=1,check.names=FALSE))
    rownames(cprod) <- colnames(cprod) <- indivinfo$SUBJID[match(indivinfo$ORIG_ID,rownames(cprod))]
    write.table(cprod,file="crossprod_all.csv",sep="\t",row.names=TRUE)

    cprod <- as.matrix(read.csv("orig_data/crossprod_chr8-inner_prod.csv",row.names=1,check.names=FALSE)) / as.matrix(read.csv("orig_data/crossprod_chr8-nonmissing.csv",row.names=1,check.names=FALSE))
    rownames(cprod) <- colnames(cprod) <- indivinfo$SUBJID[match(indivinfo$ORIG_ID,rownames(cprod))]
    write.table(cprod,file="crossprod_chr8.csv",sep="\t",row.names=TRUE)
}

###
indivinfo <- read.csv("indivinfo.csv",header=TRUE)

countrycols <- rainbow( nlevels(indivinfo$COUNTRY_SELF) )
names(countrycols) <- levels(indivinfo$COUNTRY_SELF)
# Abbreviations
country.abbrev.list <- do.call( rbind, list( 
        c("Belgium", "BE"), 
        c("France", "FR"), 
        c("Austria", "AT"), 
        c("Bulgaria", "BG"), 
        c("Italy", "IT"), 
        c("Poland", "PL"), 
        c("Czech Republic", "CZ"), 
        c("Cyprus", "CY"), 
        c("Portugal", "PT"), 
        c("Denmark", "DK"), 
        c("Latvia", "LV"), 
        c("Romania", "RO"), 
        c("Germany", "DE"), 
        c("Lithuania", "LT"), 
        c("Slovenia", "SI"), 
        c("Estonia", "EE"), 
        c("Luxembourg", "LU"), 
        c("Slovakia", "SK"), 
        c("Ireland", "IE"), 
        c("Hungary", "HU"), 
        c("Finland", "FI"), 
        c("Greece", "EL"), 
        c("Malta", "MT"), 
        c("Sweden", "SE"), 
        c("Spain", "ES"), 
        c("Netherlands", "NL"), 
        c("United Kingdom", "UK"), 
        c("Iceland", "IS"), 
        c("Norway", "NO"), 
        c("Liechtenstein", "LI"), 
        c("Switzerland", "CH"), 
        c("Croatia", "HR"), 
        c("Montenegro", "ME"), 
        c("Turkey", "TR"), 
        c("Albania", "AL"), 
        c("Serbia", "RS"), 
        c("Bosnia and Herzegovina", "BA"), 
        c("Armenia", "AM"), 
        c("Belarus", "BY"), 
        c("Georgia", "GE"), 
        c("Azerbaijan", "AZ"), 
        c("Moldova", "MD"), 
        c("Ukraine", "UA"), 
        c("Algeria", "DZ"), 
        c("Lebanon", "LB"), 
        c("Syria", "SY"), 
        c("Russia", "RU"), 
        c("Swiss French", "CHf"), 
        c("Swiss German", "CHd"), 
        c("Yugoslavia", "YU"),
        c("Bosnia", "BO"),     
        c("Croatia", "CR"),    
        c("England", "EN"),    
        c("Kosovo", "KO"),
        c("Macedonia", "MA"),
        c("Montenegro", "MO"),
        c("Scotland", "SC"),
        c("Serbia", "SR")
    ) )
countryabbrevs <- country.abbrev.list[ match(levels(indivinfo$COUNTRY_SELF),country.abbrev.list[,1]), 2 ]
names(countryabbrevs) <- levels(indivinfo$COUNTRY_SELF)

plot.pca <- function (xy,pcs=2:1,clabs=c(names(nsamples)[nsamples>10]),info=indivinfo) {
    countries <- with(info, levels(COUNTRY_SELF)[as.numeric(COUNTRY_SELF)] )
    pairs( xy[,pcs], col=adjustcolor(countrycols[countries],.5), pch=20, cex=1, 
        panel=function(x,y,...) {
            points( x, y, ... )
            mean.pca <- cbind(
                    tapply( x, info$COUNTRY_SELF, mean ),
                    tapply( y, info$COUNTRY_SELF, mean ) 
                )[clabs,]
            points( mean.pca, col=adjustcolor(countrycols[rownames(mean.pca)],.45), pch=20, cex=5 )
            text( mean.pca, labels=countryabbrevs[rownames(mean.pca)], )
        } )
}

# read in the inner product matrix
cprod <- as.matrix(read.csv("crossprod_all.csv",row.names=1,check.names=FALSE))
cprod <- as.matrix(read.csv("crossprod_chr8-inner_prod.csv",row.names=1,check.names=FALSE)) / as.matrix(read.csv("crossprod_chr8-nonmissing.csv",row.names=1,check.names=FALSE))

# unweighted
cprod.eigen <- eigen(cprod)
stopifnot(all(rownames(cprod)==indivinfo$SUBJID))

plot.pca( cprod.eigen$vectors, pcs=1:3 )

# unweighted, normalized
nmat <- diag(nrow(cprod)) - matrix(1/nrow(cprod),nrow=nrow(cprod),ncol=ncol(cprod))

n.cprod.eigen <- eigen( (nmat %*% cprod %*% nmat) )

plot.pca( n.cprod.eigen$vectors, pcs=1:3 )


# weighted
nsamples <- table(indivinfo$COUNTRY_SELF)
weightings <- 1/pmax(10,sqrt( nsamples[indivinfo$COUNTRY_SELF] ))

w.cprod.eigen <- eigen( cprod * outer(weightings,weightings,"*") )

plot.pca( w.cprod.eigen$vectors, pcs=1:3 )

# weighted, normalized
nmat <- diag(nrow(cprod)) - matrix(1/nrow(cprod),nrow=nrow(cprod),ncol=ncol(cprod))

nw.cprod.eigen <- eigen( (nmat %*% cprod %*% nmat) * outer(weightings,weightings,"*") )

plot.pca( nw.cprod.eigen$vectors, pcs=1:3 )


## subset
sub.inds <- ( indivinfo$COUNTRY_SELF %in% c("Spain","France","Italy") )
sub.cprod <- cprod[sub.inds,sub.inds]
sub.nmat <- diag(nrow(sub.cprod)) - matrix(1/nrow(sub.cprod),nrow=nrow(sub.cprod),ncol=ncol(sub.cprod))
sub.eigen <- eigen( sub.nmat %*% sub.cprod %*% sub.nmat )

plot.pca( sub.eigen$vectors, pcs=1:3, info=indivinfo[sub.inds,] )
