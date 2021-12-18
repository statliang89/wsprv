# wsprv: weighted selection probability for rare variant analysis
This package is developed to locate individual rare variants after a multiple-phenotype based group test finds significance. <br>
It takes a long time to compute weighted selection probabilities for individual rare variants when the bootstrap replication is large and a large number of rare variants exists.
We recommend users to use 'rep=100' to compute the weighted selection probability each rare variant.
In order to decrease the computation time in examples, we only use 'rep=5' to compute weighted selection probabilities. <br>  

# Installation
## "devtools" package is required if you don't have it.  
install.packages('devtools') <br>
For Windows users, you also need to install Rtools from https://cran.r-project.org/bin/windows/Rtools <br>

library(devtools) <br>
install_github("statliang89/wsprv") <br>
