#!/usr/bin/Rscript

source("sim-data-fns.R")

basedir <- "bedassle-ex"

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
        prefix = "CGW_1",   # filename prefix
        continue=FALSE,
        continuing.params=NULL
    )
