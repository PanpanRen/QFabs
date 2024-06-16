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
The detail of the usage of the package can be found in https://github.com/PanpanRen/qfabs/qfabs_package/inst/qfabs-manual.pdf.


# Replicating Figure 2(a), Figure 2(b), Figure 3(a) and Figure 3(b) in the manuscript using code1.r

## Overview
The `code1.r` script is used to replicate Figure 2(a), Figure 2(b), Figure 3(a) and Figure 3(b) in the main paper.

## Usage
1. Ensure that you have the necessary R packages installed, such as 'qfabs', 'mnormt' and 'ggplot2'.

2. Navigate to the directory where the `code1.r` script is located using the `setwd()` function if necessary.

3. Open your R terminal and run the script using the following commands:

    ```
    Rscript code1.r
    ```
    Or you can select the entire script and run it in the R environment.

4. After running the script, the output will be generated, typically in the form of `png` and `eps` files.


# Replicating Figure 2(c), Figure 2(d), Figure 3(c) and Figure 3(d) using code2.r

## Overview
The `code2.r` script is used to replicate Figure 2(c), Figure 2(d), Figure 3(c) and Figure 3(d) in the main paper.

## Usage
1. Ensure that you have the necessary R packages installed, such as 'qfabs', 'mnormt' and 'ggplot2'.

2. Navigate to the directory where the `code2.r` script is located using the `setwd()` function if necessary.

3. Open your R terminal and run the script using the following commands:

    ```
    Rscript code2.r
    ```
    Or you can select the entire script and run it in the R environment.

4. After running the script, the output will be generated, typically in the form of `png` and `eps` files. 


