
Visualizing geographic structure and demographic history
========================================================

*Goals/skills:* 

1. Describe different causes and patterns of population structure. 
2. Use principal components analysis (PCA) to visualize population structure. 
3. Use BEDASSLE to fit models of isolation by distance and/or environment.

*Incidental skills:*

1. Interpreting the output of a Markov chain Monte Carlo (MCMC) algorithm.

*Prerequisites:*

1. Install R 
2. and the R package BEDASSLE: 
``` 
install.packages("BEDASSLE")
```


Outline
-------

1.   PCA and visualization
    -   Goal: low-dimensional summary.
    -   In practice: examples from the literature
    -   Words about PCA
        -   direction of maximum variation: weights on alleles
        -   decomposition of variance
        -   matrix decompositions
        -   spatial structure and Fourier modes
    -   Technical issues
        -   normalization
        -   weighting and sample sizes
        -   correlated markers
2.   Doing PCA: hands-on
    -   Three populations, simulated
        -   computing the covariance matrix
        -   PCA
    -   POPRES (as in http://www.ncbi.nlm.nih.gov/pubmed/18758442)
        -   covariance matrix provided
        -   downweighting populations
        -   subsets of the genome
3.   Continuous geography
    -   Nearby things are more similar than distant ones
        -   ... because they are more closely related.
        -   Migration, coalescence, covariance, and genetic distance.
    -   Isolation by distance/environment/ecology/etcetera
        -   resistance distance
    -   Phenomenological model: correlated allele frequencies
4.   Using BEDASSLE
    -   BEDASSLE's parameterization
    -   Lightning introduction to MCMC
        -   acceptance rates
        -   likelihood profile
        -   mixing and stationarity
    -   Simulated data

