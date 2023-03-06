#' Generate samples for multiple datasets
#'
#' This function generates data including response and covariates for multiple datasets.
#' Data can be heavy-tailed.
#' Morever, there may be a high-correlated pattern among covariates.
#' @param n The sample sizes in multiple datasets. An integer vector in \eqn{R^{M}}, where \eqn{M} is the number of datasets.
#' @param p The number of covariates. The number of covariates in each dataset is the same.
#' @param beta A \eqn{p \times M} numeric matrix, which is the true coefficients for the \eqn{M} datasets.
#' @param distr The error's distribution, including "gaussian", "t3", "cauchy" and "logistic".
#' @param rho The strength of correlation among covariates.
#'
#' @return A list
#' \itemize{
#'   \item x - A \eqn{N \times p} design matrix, where \eqn{N} is the total sample size in \eqn{M} datasets.
#'   \item y - A length \eqn{N} vector of the response for \eqn{M} datasets.
#' }
#' @examples
#' library("mnormt")
#' M = 2
#' n1 = n2 = 20
#' n = c(n1, n2)
#' p = 50
#' beta = matrix(0, p, M)
#' index1 = 1:10
#' index2 = 2:11
#' beta[index1, 1] = runif(10, 0.2, 1.0)
#' beta[index2, 2] = runif(10, 0.4, 1.4)
#' distr = "t3"
#' rho = 0.8
#' dat = generator(n, p, beta, distr, rho)
#' x = dat$x
#' y = dat$y
#' @export

generator <- function(n, p, beta, distr = "gaussian", rho = 0.5){
    M = length(n)
    N = 0
    for(m in 1:M){
        N = N + n[m]
    }

    # x
    x = matrix(NA, N, p)
    x[1:n[1],] = rmnorm(n[1], mean = rep(0, p), varcov = outer(1:p, 1:p, FUN = function(x, y) rho^(abs(x-y))))
    if(m!=1){
        # x[1:n[1],] = rmnorm(n[1], mean = rep(0, p), varcov = outer(1:p, 1:p, FUN = function(x, y) rho^(abs(x - y))))
        for(m in 1:(M-1)){
            x[(n[m]+1):(n[m]+n[m+1]),] = rmnorm(n[m+1], mean = rep(0, p), varcov = outer(1:p, 1:p, FUN = function(x, y) rho^(abs(x - y))))
        }
    }

    # error
    error = vector("list", M)
    if(distr == "gaussian"){
        for(m in 1:M){
            error[[m]] = rnorm(n[m], 0, 1)
        }
    }
    else if(distr == "t3"){
        for(m in 1:M){
            error[[m]] = rt(n[m], 3)
        }
    }
    else if(distr == "cauchy"){
        for(m in 1:M){
            error[[m]] = rcauchy(n[m])
        }
    }
    else if(distr == "logistic"){
        for(m in 1:M){
            error[[m]] = rlogis(n[m])
        }
    }

    # y
    y = rep(NA, N)
    y[1:n[1]] = x[1:n[1],] %*% beta[,1] + error[[1]]
    if(m!=1){
        for(m in 1:(M-1)){
            index = (n[m]+1):(n[m]+n[m+1])
            y[index] = x[index,] %*% beta[,(m+1)] + error[[m+1]]
        }
    }

    val = list(x = x, y = y)

    return(val)
} 