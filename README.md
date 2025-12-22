## Iterating marginalized Bayes maps for likelihood maximization with application to nonlinear panel models

This repository contains the data and `.Rnw` (R-code) needed to reproduce the article: "Iterating marginalized Bayes maps for likelihood maximization with application to nonlinear panel models" by Jesse Wheeler, Aaron J. Abkemeir, and Edward L. Ionides.

The sub-folder `data/` contains the output of fitted models and simulation studies.
Each of the results are saved as `.rds` files, and were computed on a high-performance computer (HPC).
In the case of the simulation studies involving the stochastic Gompterz population model, the calculations were completed using the aid of the `batchtools` package in `R`. 

All code that was used to generate the results is available upon request to the corresponding author, Jesse Wheeler:  jessewheeler@isu.edu.

### Measles data

The UK-measles dataset used in the study were originally made available by:

<p style="text-indent: -1.5em; margin-left: 1.5em;">
Korevaar H, Metcalf CJ, Grenfell BT. (2020) "Structure, space and size: competing drivers of variation in urban and rural measles transmission." J. R. Soc. Interface 17: 20200010. http://dx.doi.org/10.1098/rsif.2020.0010
</p>

The data can also be accessed in the `uk_measles` data object in the `panelPomp` package, which was used for building and fitting the models in this article.


<p style="text-indent: -1.5em; margin-left: 1.5em;">
Bretó, C., Wheeler, J., King, A. A., Ionides, E. L. 2025. “panelPomp: Analysis of Panel Data via Partially Observed Markov Processes in R”. 
The R Journal, 17, 180-199. 10.32614/RJ-2025-009.
</p>


