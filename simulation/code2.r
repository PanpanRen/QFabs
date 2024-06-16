# First install the following packages:
# install.packages("devtools")
# library(devtools)
# install_github("PanpanRen/qfabs/qfabs_package")
# install.packages("mnormt")
# install.packages("ggplot2")
library(qfabs)
library(mnormt)
library(ggplot2)

assig <- function(n_args) {
	# example: n_args = c(2, 3, 4)
	cargs <- vector("list", length(n_args))
	for(i in 1:length(n_args)) cargs[[i]] = 1:n_args[i]
	# expand.grid: Create a data frame from all combinations of the supplied vectors or factors
	t(expand.grid(cargs))
}

# generate data
generatedata <- function(n1, n2, n3, p, beta1, beta2, beta3, distr, rho){
    varcov = diag(p)
    x1 = rmnorm(n1, mean = rep(0,p), varcov = outer(1:p, 1:p, FUN = function(x, y) rho^(abs(x - y))))
    x2 = rmnorm(n2, mean = rep(0,p), varcov = outer(1:p, 1:p, FUN = function(x, y) rho^(abs(x - y))))
    x3 = rmnorm(n3, mean = rep(0,p), varcov = outer(1:p, 1:p, FUN = function(x, y) rho^(abs(x - y))))
    if (distr == 1) {
        epsilon1 <- rnorm(n1, 0, 1)
        epsilon2 <- rnorm(n2, 0, 1)
        epsilon3 <- rnorm(n3, 0, 1)
    }
    else if (distr == 2) {
        epsilon1 <- rt(n1, 3)
        epsilon2 <- rt(n2, 3)
        epsilon3 <- rt(n3, 3)
    }
    else if (distr == 3) {
        epsilon1 <- rcauchy(n1)
        epsilon2 <- rcauchy(n2)
        epsilon3 <- rcauchy(n3)
    }
    else {
        epsilon1 <- rlogis(n1)
        epsilon2 <- rlogis(n2)
        epsilon3 <- rlogis(n3)
    }
    y1 <- x1 %*% beta1 + epsilon1
    y2 <- x2 %*% beta2 + epsilon2
    y3 <- x3 %*% beta3 + epsilon3

    list(x=list(x1=x1,x2=x2, x3=x3), y=list(y1=y1,y2=y2,y3=y3))
}

main <- function(n1,n2,n3,d,q,distr) {
    beta1 <- rep(0, d)
    beta2 <- rep(0, d)
    beta3 <- rep(0, d)
    index1 = index2 = c(1:8,10:11)
    index3 = c(1:10)
    set.seed(1000)
    beta1[1:5] = runif(5,0.6,1.4)
    beta1[6:10] = runif(5,-1.4,-0.6)
    diff = runif(10,-0.2,0.2)
    beta2[index2] = beta1[index1] + diff
    diff = runif(10,-0.2,0.2)
    beta3[index3] = beta2[index2] + diff

    rho = 0.5
    M = 3
    epsilon = 0.01
    delta = 1e-12 
    xi = 1e-10 
    max.iter = 4000
    
    K = 1000    # replication times

    L1 = matrix(NA, K, 2)
    FP = matrix(NA, K, 2)
    TP = matrix(NA, K, 2)
    for(k in 1:K){
        if(k%%10==0){cat("k = ", k, "\n")}
        set.seed(10000+k)
        data <- generatedata(n1, n2, n3, d, beta1, beta2, beta3, distr, rho)
        x1 = data$x$x1
        x2 = data$x$x2
        x3 = data$x$x3
        x = rbind(x1,x2,x3)
        y1 = data$y$y1
        y2 = data$y$y2
        y3 = data$y$y3
        y = c(y1,y2,y3)

        # Lasso
        result1 = mqfabs(y, x, n = c(n1,n2,n3), tau = q, Lambda2 = 0, epsilon = epsilon, delta = delta, xi = xi, max.iter = max.iter)
        betahat0 = result1$beta
        betahat0 = matrix(betahat0, d+1, 3)
        L1[k,1] = sum(abs(betahat0[-1,1]-beta1))+sum(abs(betahat0[-1,2]-beta2))+sum(abs(betahat0[-1,3]-beta3))
        TP[k,1] = sum(betahat0[-1,1][index1]!=0)+sum(betahat0[-1,2][index2]!=0)+sum(betahat0[-1,3][index3]!=0)
        FP[k,1] = sum(betahat0[-1,1][-index1]!=0)+sum(betahat0[-1,2][-index2]!=0)+sum(betahat0[-1,3][-index3]!=0)
        
        # CLasso
        result1 = mqfabs(y, x, n = c(n1,n2,n3), tau = q, nlambda2 = 20, epsilon = epsilon, delta = delta, xi = xi, max.iter = max.iter)
        betahat1 = result1$beta
        betahat1 = matrix(betahat1, d+1, 3)
        L1[k,2] = sum(abs(betahat1[-1,1]-beta1))+sum(abs(betahat1[-1,2]-beta2))+sum(abs(betahat1[-1,3]-beta3))
        TP[k,2] = sum(betahat1[-1,1][index1]!=0)+sum(betahat1[-1,2][index2]!=0)+sum(betahat1[-1,3][index3]!=0)
        FP[k,2] = sum(betahat1[-1,1][-index1]!=0)+sum(betahat1[-1,2][-index2]!=0)+sum(betahat1[-1,3][-index3]!=0)
    }
    l1 = apply(L1, 2, mean, na.rm = T)
    tp = apply(TP, 2, mean, na.rm = T)
    fp = apply(FP, 2, mean, na.rm = T)
    output = cbind(l1,tp,fp)
    return(output)
}

# bar plot of three measures: L1, TP and FP
myplot <- function(data){
    bp2 <- ggplot(data, aes(x=p1, y=value1, fill=method1)) +
        
    geom_bar(position = "dodge", stat="identity") + 
    
    facet_grid(distr1 ~ measure1, scales = "fixed", labeller=label_parsed)+
    
    guides(fill=guide_legend(title=NULL, nrow = 1)) + 

    xlab(expression(p)) + 
    
    theme_bw() + 
    
    theme(plot.title = element_text(size=25,face = "bold",hjust=0.5),
        axis.title.x = element_text(size=25,face = "bold"),
        axis.text.x = element_text(size=25,face = "bold"),
        axis.text.y = element_text(size=25,face = "bold"),
        axis.title.y = element_blank(),
        legend.position="bottom",
        strip.text = element_text(size=25,face = "bold"),
        legend.title=element_text(size=25,face = "bold"),
        legend.text=element_text(size=25,face = "bold")) + 
    scale_fill_manual(values=c( "#d3ba68", "#d5695d"))
    
    bp2+scale_y_continuous(limits = c(0, 1))
}

exam <- function(){
    distrr = 1:4
    dd = c(500,1000)
    nn = c(100,200)
    qq = c(0.2,0.5)

    n_args = c(length(distrr), length(dd), length(nn), length(qq)) # 4*2*2*2 = 32
    jobs = assig(n_args)

    result = matrix(NA, ncol(jobs)*2, 3)
    for(number in 1:ncol(jobs)){
        id = jobs[,number]
        distr = distrr[id[1]]
        d = dd[id[2]]
        n1 = n2 = n3 = nn[id[3]]
        q = qq[id[4]]
        cat("n1=",n1,"d=",d,"q=",q,"distr=",distr,"\n")

        result[(2*number-1):(2*number),] = main(n1,n2,n3,d,q,distr)
    }
    colnames(result) = c("L1","TP","FP")
    rownames(result) = rep(c("Lasso", "CLasso"),32)
    print(result)

    max_ = apply(result, 2, max, na.rm = T)
    result[,1] = result[,1]/max_[1]
    result[,2] = result[,2]/max_[2]
    result[,3] = result[,3]/max_[3]

    P = c(500,1000)
    Mlist1 <- c("CLasso", "Lasso")
    measure = rep(c(expression(L[1]),"TP","FP"),each = 16)
    method = rep(rep(c("Lasso","CLasso"),4),3)
    p = rep(rep(c(500,1000),each = 8),3)
    distr = rep(rep(c("N(0,1)","t(3)","Cauchy(0,1)","Logistic(0,1)"),each = 2),6)

    p1      = factor(p, levels=as.character(P))
    measure1  = factor(measure, levels=as.character(c(expression(L[1]),"TP","FP")))
    method1  = factor(method, levels=Mlist1)
    distr1 = factor(distr,levels=c("N(0,1)","t(3)","Cauchy(0,1)","Logistic(0,1)"))

    # Figure 2(c) in the manuscript: n = 100, q = 0.2, scenario 2
    subset = result[1:(nrow(result)/4),]
    value = as.vector(subset)

    data1 = data.frame(
        value1    = value, 
        p1      = p1,
        measure1  = measure1,
        method1  = method1, 
        distr1 = distr1
    )

    plot1 = myplot(data1)
    
    filename = paste("n=100","rho=0.5","tau=0.2","s2", "bar",".png", sep="_")
    ggsave(filename, plot1, width=28, height=28, units='cm')
    
    filename = paste("n=100","rho=0.5","tau=0.2","s2", "bar",".eps", sep="_")
    ggsave(filename, plot1, width=28, height=28, units='cm')

    # Figure 2(d) in the manuscript: n = 200, q = 0.2, scenario 2
    subset = result[(nrow(result)/4+1):(nrow(result)/2),]
    value = as.vector(subset)

    data1 = data.frame(
        value1    = value, 
        p1      = p1,
        measure1  = measure1,
        method1  = method1, 
        distr1 = distr1
    )

    plot1 = myplot(data1)
    
    filename = paste("n=200","rho=0.5","tau=0.2","s2", "bar",".png", sep="_")
    ggsave(filename, plot1, width=28, height=28, units='cm')
    
    filename = paste("n=200","rho=0.5","tau=0.2","s2", "bar",".eps", sep="_")
    ggsave(filename, plot1, width=28, height=28, units='cm')

    # Figure 3(c) in the manuscript: n = 100, q = 0.5, scenario 2
    subset = result[(nrow(result)/2+1):(nrow(result)*3/4),]
    value = as.vector(subset)

    data1 = data.frame(
        value1    = value, 
        p1      = p1,
        measure1  = measure1,
        method1  = method1, 
        distr1 = distr1
    )

    plot1 = myplot(data1)
    
    filename = paste("n=100","rho=0.5","tau=0.5","s2", "bar",".png", sep="_")
    ggsave(filename, plot1, width=28, height=28, units='cm')
    
    filename = paste("n=100","rho=0.5","tau=0.5","s2", "bar",".eps", sep="_")
    ggsave(filename, plot1, width=28, height=28, units='cm')

    # Figure 3(d) in the manuscript: n = 200, q = 0.5, scenario 2
    subset = result[(nrow(result)*3/4+1):nrow(result),]
    value = as.vector(subset)

    data1 = data.frame(
        value1    = value, 
        p1      = p1,
        measure1  = measure1,
        method1  = method1, 
        distr1 = distr1
    )

    plot1 = myplot(data1)
    
    filename = paste("n=200","rho=0.5","tau=0.5","s2", "bar",".png", sep="_")
    ggsave(filename, plot1, width=28, height=28, units='cm')
    
    filename = paste("n=200","rho=0.5","tau=0.5","s2", "bar",".eps", sep="_")
    ggsave(filename, plot1, width=28, height=28, units='cm')
}

exam()