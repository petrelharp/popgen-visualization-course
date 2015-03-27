Visualizing geographic structure and demographic history: (very) short course
======================================================================

This is the material for a short (3 hour) course I taught at the UCLA/La Kretz Workshop in Conservation Genomics, 22-27 March, 2015.

To see the slides directly, go to [the github-pages site](http://petrelharp.github.io/popgen-visualization-course/)


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

To get the data used in the examples, download and uncompress [this tarball](http://phoebe.usc.edu/CGW/CGW-data.tar.gz) to the folder you've put this git repository in.


In this repository
------------------

1. [The presentation.](slides.html)
2. [Example of running PCA on some data.](popres/popres-pca.html)
3. [Example of running BEDASSLE on some simulated data.](bedassle-ex/view-bedassle-sim.html)
4. Example of running BEDASSLE on two real datasets: [teosinte](bedassle-ex/view-bedassle-teosinte.html) and the [HGDP](bedassle-ex/view-bedassle-hgdp.html).

Each of the above html files are made from the corresponding Rmd files of the same name;
you can read the source code to see what happens behind the scenes when the R code isn't in the presentation.

Outline
-------

1.   PCA and visualization (45min)
    -   Goal: low-dimensional summary.
    -   In practice: examples from the literature
    -   Properties of PCA
        -   direction of maximum variation: weights on alleles
        -   decomposition of variance
        -   matrix decompositions
        -   spatial structure and Fourier modes
    -   Technical issues
        -   normalization
        -   weighting and sample sizes
        -   correlated markers
2.   Doing PCA: hands-on (45min)
    -   Three populations, simulated
        -   computing the covariance matrix
        -   PCA
    -   POPRES (as in http://www.ncbi.nlm.nih.gov/pubmed/18758442)
        -   covariance matrix provided
        -   downweighting populations
        -   subsets of the genome
3.   Continuous geography (45min)
    -   Nearby things are more similar than distant ones
        -   ... because they are more closely related.
        -   Migration, coalescence, covariance, and genetic distance.
    -   Isolation by distance/environment/ecology/etcetera
        -   resistance distance
    -   Phenomenological model: correlated allele frequencies
4.   Using BEDASSLE (45min)
    -   BEDASSLE's parameterization
    -   Lightning introduction to MCMC
        -   acceptance rates
        -   likelihood profile
        -   mixing and stationarity
    -   Simulated data



License and reuse
-----------------

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work, including the code and [presentation](slides.html), is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
