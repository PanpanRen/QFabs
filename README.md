# qfabs
A minorization-maximization forward and backward stagewise algorithm for high-dimensional quantile regression and high-dimensional integrative quantile regression.

qfabs uses coordinate descent with a fixed step size which consists of both forward and backward steps. At each step, the first-order Taylor's expansion is used to reflect the main part of the increment. For given lambda2, a tuning parameter before the contrast penalty, the contrast penalty can be incorporated into the loss function.

# Installation

    #install Rtools 3.5 (http://cran.r-project.org/bin/windows/Rtools)
    #install.packages("devtools")
    #install.packages("Rcpp")
    library(devtools)
    install_github("PanpanRen/qfabs")

# Usage

- [x] [qfabs-manual](https://github.com/PanpanRen/qfabs/inst/qfabs-manual.pdf) ------------ Details of the usage of the package.

# Example

    library(qfabs)
    library(mnormt)

    M = 2
    n1 = n2 = 20
    n = c(n1, n2)
    p = 50
    beta = matrix(0, p, M)
    index1 = 1:10
    index2 = 2:11
    beta[index1, 1] = runif(10, 0.2, 1.0)
    beta[index2, 2] = runif(10, 0.4, 1.4)
    distr = "gaussian"
    rho = 0.5

    dat = generator(n, p, beta, distr, rho)
    x = dat$x
    y = dat$y
    tau = 0.2
    Lambda2 = 0.001

    fit <- mqfabs(y, x, n, tau, Lambda2)



# References

Integrative analysis of regional differences in elder support with high-dimensional quantile regression. Manuscript.

# Development
The R-package is developed by Panpan Ren (ren08067146@163.com), Xu Liu, Xiao Zhang, Peng Zhan and Tingting Qiu.




