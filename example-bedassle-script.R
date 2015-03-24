#!/usr/bin/Rscript
library(BEDASSLE)
setwd("bedassle-ex/sim")
load("sim-data.Robj")
MCMC(   counts = t(sim$genotypes), sample_sizes = matrix(2,nrow=nind,ncol=nloci),
        D = D,  E = E,  k = nind, loci = nloci,  delta = 0.0001,
        aD_stp = 0.5, aE_stp = 0.5, a2_stp = 0.02, thetas_stp = 0.2, mu_stp = 0.25,
        ngen = 1e6,                # should take about 5.5 hours
        printfreq = 10000,         # no need to write out all the time
        savefreq = 1e5,            # save every 0.55 hours
        samplefreq = 250,          # as determined above
        prefix = "CGW_4_",   # NEW filename prefix
        continue=TRUE,                       # CONTINUE from previous run
        continuing.params="CGW_3_MCMC_final1.Robj")

