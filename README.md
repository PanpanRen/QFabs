# Introduction to the qfabs R Package
The source code of the R package is provided in folder qfabs_package.

## Overview
R package "qfabs" provides a minorization-maximization forward and backward stagewise algorithm for high-dimensional quantile regression and high-dimensional integrative quantile regression. 
qfabs uses coordinate descent with a fixed step size which consists of both forward and backward steps. At each step, the first-order Taylor's expansion is used to reflect the main part of the increment. For given lambda2, a tuning parameter before the contrast penalty, the contrast penalty can be incorporated into the loss function.

## Installation
You can install the qfabs package in R using the following steps:
    install.packages("devtools")
    library(devtools)
    install_github("PanpanRen/qfabs/qfabs_package")

## Usage
The detail of the usage of the package can be found in https://github.com/PanpanRen/qfabs/qfabs_package/inst/qfabs-manual.pdf.

# Introduction to the simulation code
The simulation results can be reproduced by using the code in folder simulation. The detail description can be found in simulation/README.md.

# Introduction to the real datasets and real data code
The real datasets and real data code are provided in folder realdata. The detail description can be found in realdata/README.md.
