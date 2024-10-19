# Introduction to the qfabs R Package

## Overview
R package "qfabs" provides a minorization-maximization forward and backward stagewise algorithm for high-dimensional quantile regression and high-dimensional integrative quantile regression. 
qfabs uses coordinate descent with a fixed step size which consists of both forward and backward steps. At each step, the first-order Taylor's expansion is used to reflect the main part of the increment. For given lambda2, a tuning parameter before the contrast penalty, the contrast penalty can be incorporated into the loss function.

## Installation
You can install the qfabs package in R using the following steps:
    install.packages("devtools")
    library(devtools)
    install_github("PanpanRen/qfabs/qfabs_package")

## Usage
The detail of the usage of the package can be found in https://github.com/PanpanRen/qfabs/tree/main/qfabs_package/inst/qfabs-manual.pdf.


# Introduction to the file `names.md`
The file `names.md` includes the detailed explanation of the meanings of the response variable and the explanatory variables used in the real data analysis.


# Replicating the results of (Hebei, Shandong, Fujian) for tau=0.5 in Table 1 of the manuscript using real1.r

## Overview
The `real1.r` script is used to replicate the results of (Hebei, Shandong, Fujian) for tau=0.5 in Table 1 of the main paper.

## Usage
1. Ensure that you have the necessary R packages installed, such as 'qfabs' and 'readxl'.

2. Navigate to the directory where the `real1.r` script is located using the `setwd()` function if necessary.

3. Open your R terminal and run the script using the following commands:

    ```
    Rscript real1.r
    ```
    Or you can select the entire script and run it in the R environment.

4. After running the script, the output will be generated, typically printed directly to the terminal. Then, the estimated coefficients and the results of random split can be observed.


# Replicating the results of (Zhejiang, Anhui, Shaanxi) for tau=0.5 in Table 1 of the manuscript using real2.r

## Overview
The `real2.r` script is used to replicate the results of (Zhejiang, Anhui, Shaanxi) for tau=0.5 in Table 1 of the main paper.

## Usage
1. Ensure that you have the necessary R packages installed, such as 'qfabs' and 'readxl'.

2. Navigate to the directory where the `real2.r` script is located using the `setwd()` function if necessary.

3. Open your R terminal and run the script using the following commands:

    ```
    Rscript real2.r
    ```
    Or you can select the entire script and run it in the R environment.

4. After running the script, the output will be generated, typically printed directly to the terminal. Then, the estimated coefficients and the results of random split can be observed.

The results for tau=0.7 can be similarly obtained by changing q=0.5 to q=0.7 in the scripts `real1.r` and `real2.r`.
