###
indivinfo <- read.table("indivinfo.tsv",header=TRUE)
nsamples <- table(indivinfo$COUNTRY_SELF)

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

plot.pca <- function (xy,pcs=1:2,clabs=c(names(nsamples)[nsamples>10]),info=indivinfo) {
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

