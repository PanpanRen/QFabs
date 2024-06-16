# First install the following packages:
# install.packages("devtools")
# library(devtools)
# install_github("PanpanRen/qfabs/qfabs_package")
# install.packages("readxl")
library(qfabs)
library(readxl)

loss1 <- function(r, q) {
    c = ifelse(r<0, (q-1)*r, q*r)
    l = mean(c)
    l
}

split <- function(y1,y2,y3,x1,x2,x3,q,epsilon,delta,xi,max.iter,K = 50,M = 3){
    n1 = length(y1)
    n2 = length(y2)
    n3 = length(y3)
    d = ncol(x1)
    n1_train = ceiling(n1*2/3)
    n2_train = ceiling(n2*2/3)
    n3_train = ceiling(n3*2/3)
    n_train = c(n1_train, n2_train, n3_train)
    number1 = rep(NA, K)
    error1 = rep(NA, K)
    number2 = rep(NA, K)
    error2 = rep(NA, K)
    for(k in 1:K){
        cat("k=",k,"\n")
        set.seed(10000+k)
        index1 = sample(c(1:n1), size = n1_train)
        index2 = sample(c(1:n2), size = n2_train)
        index3 = sample(c(1:n3), size = n3_train)
        y1_train = y1[index1]
        y2_train = y2[index2]
        y3_train = y3[index3]
        x1_train = x1[index1,]
        x2_train = x2[index2,]
        x3_train = x3[index3,]
        y1_test = y1[-index1]
        y2_test = y2[-index2]
        y3_test = y3[-index3]
        x1_test = cbind(1,x1[-index1,])
        x2_test = cbind(1,x2[-index2,])
        x3_test = cbind(1,x3[-index3,])
        y_train = c(y1_train,y2_train,y3_train)
        x_train = rbind(x1_train,x2_train,x3_train)
        y_test = c(y1_test,y2_test,y3_test)
        x_test = rbind(x1_test,x2_test,x3_test)

        # Lasso
        result1 = mqfabs(y_train, x_train, n = n_train, tau = q, Lambda2 = 0, epsilon = epsilon, delta = delta, xi = xi, max.iter = max.iter, gamma = 0)
        betahat = matrix(result1$beta,d+1,M)
        number1[k] = d - sum(rowSums(betahat[-1,]!=0)==0)
        r11 = y1_test - x1_test%*%betahat[,1]
        r12 = y2_test - x2_test%*%betahat[,2]
        r13 = y3_test - x3_test%*%betahat[,3]
        error1[k] = loss1(r11,q)+loss1(r12,q)+loss1(r13,q)

        # CLasso
        result1 = mqfabs(y_train, x_train, n = n_train, tau = q, nlambda2 = 20, Lambda2_thre = 1e-4, epsilon = epsilon, delta = delta, xi = xi, max.iter = max.iter, gamma = 0)
        betahat1 = matrix(result1$beta, d+1, M)
        number2[k] = d - sum(rowSums(betahat1[-1,]!=0)==0)
        r21 = y1_test - x1_test%*%betahat1[,1]
        r22 = y2_test - x2_test%*%betahat1[,2]
        r23 = y3_test - x3_test%*%betahat1[,3]
        error2[k] = loss1(r21,q)+loss1(r22,q)+loss1(r23,q)
    }
    return(cbind(number1,error1,number2,error2))
}

data1 = read_excel("dataset/data_zhejiang.xlsx")
data2 = read_excel("dataset/data_anhui.xlsx")
data3 = read_excel("dataset/data_shaanxi.xlsx")

y1 = data1[,1]
y1 = as.vector(unlist(y1))

y2 = data2[,1]
y2 = as.vector(unlist(y2))

y3 = data3[,1]
y3 = as.vector(unlist(y3))

x1 = data1[,-1]
n1 = nrow(x1)
x1 = matrix(unlist(x1),nrow = n1)

x2 = data2[,-1]
n2 = nrow(x2)
x2 = matrix(unlist(x2),nrow = n2)

x3 = data3[,-1]
n3 = nrow(x3)
x3 = matrix(unlist(x3),nrow = n3)

y = c(y1,y2,y3)
x = rbind(x1,x2,x3)
M = 3
d = ncol(x)
q = 0.5
epsilon = 0.0001
delta = 1e-12 
xi = 1e-10 
max.iter = 500000

# Lasso
result1 = mqfabs(y, x, n = c(n1,n2,n3), tau = q, Lambda2 = 0, epsilon = epsilon, delta = delta, xi = xi, max.iter = max.iter, gamma = 0)
betahat = matrix(result1$beta,d+1,M)
cat("the estimated results of Lasso is: ", "\n")
print(betahat)

# CLasso
coer = c(abs(t(x1) %*% (y1-quantile(y1,q)))/n1,abs(t(x2) %*% (y2-quantile(y2,q)))/n2,abs(t(x3) %*% (y3-quantile(y3,q)))/n3)
lambda2.max = quantile(coer,0.9)
lambda2.min <- lambda2.max * 10^(-4) 
grid.n <- 1*10^2 
lambda2.grid <- exp(seq(log(lambda2.min), log(lambda2.max), length.out=grid.n))
result1 = mqfabs(y, x, n = c(n1,n2,n3), tau = q, Lambda2 = lambda2.grid, epsilon = epsilon, delta = delta, xi = xi, max.iter = max.iter, gamma = 0)
betahat1 = matrix(result1$beta, d+1, M)
cat("the estimated results of CLasso is: ", "\n")
print(betahat1)

# random split
output = split(y1,y2,y3,x1,x2,x3,q,epsilon,delta,xi,max.iter)
number1 = output[,1]
error1 = output[,2]
number2 = output[,3]
error2 = output[,4]
cat("the average number of variables selected using Lasso is: ", mean(number1,na.rm = TRUE), "\n", 
"the average number of variables selected using CLasso is: ", mean(number2, na.rm = TRUE), "\n", 
"the prediction error using Lasso is: ", mean(error1,na.rm=TRUE), "\n", 
"the prediction error using CLasso is: ", mean(error2,na.rm=TRUE), "\n", 
"the standard error of the prediction error using Lasso is: ", sd(error1,na.rm = TRUE), "\n",
"the standard error of the prediction error using CLasso is: ", sd(error2,na.rm = TRUE))