FastWilcox and fasterWilcox are c++ functions intended to speed up the 
wilcox.test(x,y,alternative = "two.sided", paired = F,exact = F,correct = T) function in R:

To use them, Rcpp is required

The R script evaluates the performances of the two functions against the function gficf:::rcpp_parallel_WMU_test()
taken from https://github.com/dibbelab/gficf/blob/master/src/mann_whitney.cpp

The latter function serve the same purpose as the first two but is not standalone


FasterWilcox() improves performances even more than fastWilcox() at the cost of precision

Both functions return a pValue