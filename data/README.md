## Description

All data that are used in this study are publically available, and readily accessible in the `panelPomp` `R` package, which is available on CRAN. 

The `.rds` files contained here are the output of computations and simulation studies described in the article.
Code used to create these files is available upon request to the study authors: jessewheeler@isu.edu. 

### Gompertz simulation study

The Gompertz simulation study is broken down into three files for ease of computing, and later combined in the `.Rnw` file in the parent folder. 

- `GompertzRL3.rds`
- `GompertzBigUnitRL3.rds`
- `GompertzHugeRL3.rds`

These three files represent distinct tests for the efficacy of PIF, MPIF, and a naive Kalman-Filter approach on various levels of unit-dimensions.

### Model fit to measles data

The `measles.rds` contains the results of fitting multiple versions of the measles models to the UK-measles data.

