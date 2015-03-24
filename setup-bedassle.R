#!/usr/bin/Rscript

source("sim-data-fns.R")
require(BEDASSLE)
source("bedassle-fixes.R")

basedir <- "bedassle-ex"

data(HGDP.bedassle.data)

thisdir <- file.path(basedir,"hgdp-ex")
dir.create(thisdir,showWarnings=FALSE)

mcmc1.message <- MCMC(
         counts = HGDP.bedassle.data$allele.counts,
         sample_sizes = HGDP.bedassle.data$sample.sizes,
         D = HGDP.bedassle.data$GeoDistance,
         E = HGDP.bedassle.data$EcoDistance,
         k = HGDP.bedassle.data$number.of.populations,
         loci = HGDP.bedassle.data$number.of.loci,
         delta = 0.001,
         aD_stp = 0.0018,
         aE_stp = 0.04,
         a2_stp = 0.0035,
         thetas_stp = 0.07,
         mu_stp = 0.17,
         ngen = 100,
         printfreq = 10,
         savefreq = 100,
         samplefreq = 5,
         prefix = file.path(thisdir,"HGDP_EX1_")
)

mcmc1 <- gsub(".*runs to ","",gsub(" [.]$","",gsub("*","1",mcmc1.message,fixed=TRUE)))
show(load(mcmc1))


pdf(file="mcmc1_diagnostics.pdf")
plot_all_trace(mcmc1)
plot_all_marginals(mcmc1)
plot_all_joint_marginals(mcmc1)
plot_all_acceptance_rates(mcmc1)
dev.off()

restart.file <- gsub("output","restart",mcmc1)
make.continuing.params(mcmc1,file=restart.file)

mcmc2.message <- MCMC(
         counts = HGDP.bedassle.data$allele.counts,
         sample_sizes = HGDP.bedassle.data$sample.sizes,
         D = HGDP.bedassle.data$GeoDistance,
         E = HGDP.bedassle.data$EcoDistance,
         k = HGDP.bedassle.data$number.of.populations,
         loci = HGDP.bedassle.data$number.of.loci,
         delta = 0.001,
         aD_stp = 0.0020,
         aE_stp = 0.08,
         a2_stp = 0.0035,
         thetas_stp = 0.02,
         mu_stp = 0.17,
         ngen = 100,
         printfreq = 10,
         savefreq = 100,
         samplefreq = 5,
         prefix = file.path(thisdir,"HGDP_EX2_"),
         continue = TRUE,
         continuing.params = restart.file
)

mcmc2 <- gsub(".*runs to ","",gsub(" [.]$","",gsub("*","1",mcmc2.message,fixed=TRUE)))
show(load(mcmc2))


pdf(file="mcmc2_diagnostics.pdf")
plot_all_trace(mcmc2)
plot_all_marginals(mcmc2)
plot_all_joint_marginals(mcmc2)
plot_all_acceptance_rates(mcmc2)
dev.off()


restart.file <- gsub("output","restart",mcmc2)
make.continuing.params(mcmc2,file=restart.file)

longrun.script <- paste("#!/usr/bin/Rscript
setwd('", dirname(file.path(getwd(),restart.file)), "')
library(BEDASSLE)
data(HGDP.bedassle.data)
mcmc.message <- MCMC(
         counts = HGDP.bedassle.data$allele.counts,
         sample_sizes = HGDP.bedassle.data$sample.sizes,
         D = HGDP.bedassle.data$GeoDistance,
         E = HGDP.bedassle.data$EcoDistance,
         k = HGDP.bedassle.data$number.of.populations,
         loci = HGDP.bedassle.data$number.of.loci,
         delta = 0.001,
         aD_stp = 0.0020,
         aE_stp = 0.10,
         a2_stp = 0.0035,
         thetas_stp = 0.02,
         mu_stp = 0.17,
         ngen = 2e6,
         printfreq = 1e4,
         savefreq = 2e5,
         samplefreq = 1000,
         prefix = 'HGDP_EX3_',
         continue = TRUE,
         continuing.params = '", basename(restart.file), "'
)
", sep='')

script.dir <- file.path(basedir, "scripts")
script.file <- file.path(script.dir,gsub("Robj$","R",basename(restart.file))) 
dir.create( script.dir, showWarnings=FALSE )
cat(longrun.script, file=script.file )


### sim data
if (FALSE) {

    nind <- 50; nloci <- 1000
    locs <- data.frame( lat=runif(nind), lon=runif(nind) )
    env <- ifelse( locs$lat>0.5, 1, 0 )
    sim <- sim_data(nind=nind,nloci=nloci,locs=locs,env=env,aE=2)

    D <- sqrt( outer( sim$locs[,1], sim$locs[,1], "-" )^2 + outer( sim$locs[,2], sim$locs[,2], "-" )^2 )
    E <- abs( outer( sim$env, sim$env, "-" ) )

    MCMC_run <- MCMC(	
            counts = t(sim$genotypes),
            sample_sizes = matrix(2,nrow=nind,ncol=nloci),
            D = D,  # geographic distances
            E = E,  # environmental distances
            k = nind, loci = nloci,  # dimensions of the data
            delta = 0.0001,  # a small, positive, number
            aD_stp = 0.01,   # step sizes for the MCMC
            aE_stp = 0.01,
            a2_stp = 0.025,
            thetas_stp = 0.2,
            mu_stp = 0.35,
            ngen = 100, 		# number of steps (2e6)
            printfreq = 10,		# print progress (1000)
            savefreq = 2e5,     # save out current state
            samplefreq = 5,		# record current state for posterior (2000)
            prefix = file.path(basedir,"discrete-sim_"),   # filename prefix
            continue=FALSE,
            continuing.params=NULL
        )
}
