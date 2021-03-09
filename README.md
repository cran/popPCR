# popPCR

R package for fitting droplet fluorescence populations of dPCR amplitude data using Expectation Maximization.

## Installation

Install from CRAN
```r
install.packages("popPCR")
```

Install from GitHub

```r
install.packages("devtools")
devtools::install_github("Zeroh729/popPCR")
```

## Usage
```
popPCR(x, dist = "t")
```

3 example datasets are available upon import
```r
library(popPCR)

hist(x_onePop, breaks = 100)     # dPCR sample w/ 1 population
hist(x_twoPop, breaks = 100)     # dPCR sample w/ 2 populations
hist(x_multiPop, breaks = 100)   # dPCR sample w/ >=3 populations
```
  
Case 1. One population sample
```r
result <- popPCR(x_onePop, dist = "t")
#        Populations detected : 1
#        Total droplets : 8000
#        Positive : 1 (0.01%)
#        Negative : 7999 (99.99%)
#
#        Target copies in sample          : 2.9414 ( 95% CI: [ -2.8237 , 8.7064 ] )
#        Mean target copies per partition : 1e-04 ( 95% CI: [ -1e-04 , 4e-04 ] )

# Increasing negProbThres makes negative classification stricter
result <- popPCR(x_onePop, dist = "t", negProbThres = 1e-4)  
#        Populations detected : 1
#        Total droplets : 8000
#        Positive : 691 (8.64%)
#        Negative : 7309 (91.36%)
#
#        Target copies in sample          : 2125.5312 ( 95% CI: [ 1966.9936 , 2284.0688 ] )
#        Mean target copies per partition : 0.0903 ( 95% CI: [ 0.0836 , 0.0971 ] )
```
  
Case 2. Two population sample
```r
result <- popPCR(x_twoPop, dist = "t")
#        Populations detected : 2
#        Total droplets : 10254
#        Positive : 8693 (84.78%)
#        Negative : 1561 (15.22%)
#
#        Target copies in sample          : 44290.3819 ( 95% CI: [ 43215.6408 , 45365.1231 ] )
#        Mean target copies per partition : 1.8823 ( 95% CI: [ 1.8367 , 1.928 ] )
```
  
Case 3. Multiple population sample
```r
result <- popPCR(x_multiPop, dist = "t", maxComponents = 4)
#        Populations detected : 4
#        Total droplets : 1814
#        Positive : 44 (2.43%)
#        Negative : 1252 (69.02%)
#        Rain (1) : 258 (14.22%)
#        Rain (2) : 260 (14.33%)
#
#        Target copies in sample          : 8724.5195 ( 95% CI: [ 7999.0578 , 9449.9812 ] )
#        Mean target copies per partition : 0.3708 ( 95% CI: [ 0.34 , 0.4016 ] )

# In the output above, we see 2 rain populations! Let's examine its plot.
plot(density(x_multiPop))

# We can see that Rain (1) is very close to the Negative population.
# Let's include droplets in Rain (1) in the negative droplet count.
nNegative <- result@dropletCount$neg + result@dropletCount$rain1
nTotal <- result@dropletCount$total

# Re-estimate concentration as follows
newEstimates <- calculateConc(nNegative, nTotal, volSamp = 20, volDrp = 0.85)
newEstimates
#    Output:
#       $lambda
#          lambda     lower     upper
#       0.1834247 0.1627763 0.2040731
#
#       $conc
#           conc    lower    upper 
#       4315.875 3830.031 4801.719 
```
  
Print results summary
```r
result <- popPCR(x_twoPop, dist = "t")
printSummaryFit(result)
#        Results of fitting a 2-component t mixture model
#
#        Negative Population
#        Mix prop. : 0.1522
#        Mu        : 2136.7435
#        Sigma     : 4126.8357
#        Dof       : 12.3562
#
#        Positive Population
#        Mix prop. : 0.8478
#        Mu        : 7580.1275
#        Sigma     : 42621.1894
#        Dof       : 2.415
```
  
Available `dist` values : `normal`, `skewed-normal`, `t`, and `skewed-t`  
Use `?popPCR` to view complete documentation.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[GNU GPLv3.0](https://choosealicense.com/licenses/gpl-3.0/)
