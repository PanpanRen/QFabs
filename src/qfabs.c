#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

void printArrayDouble(const double *arr, int n, int ncol){
	for (int i=0; i<n; i++){
		printf("%f   ",arr[i]);
		if(!((i+1)%ncol)) printf("\n");
	}
	printf("\n");
}

void printArrayDoubleInt(const int *arr, int n, int ncol){
	for (int i=0; i<n; i++){
		printf("%d   ",arr[i]);
		if(!((i+1)%ncol)) printf("\n");
	}
	printf("\n");
}

void printArrayDouble1(const double *arr, int n, int nrow, int ncol){
	if(nrow>n) fprintf(stderr, "nrow must not be greater than n!");
	for (int i=0; i<nrow; i++){
		for(int j=0; j<ncol; j++){
			printf("%f   ",arr[j*n+i]);
		}
		printf("\n");
	}
}

void preprocess(double *x, int *n, int M, int p, double *sdx){
    int i, j, N=n[M-1], m, tmp, tmp2, tmp3;
    double cor, var, *x_j, tmp1;

    for(m=0; m<M; m++){
        tmp2 = m*p;
        sdx[tmp2] = 1.0;
        if(m==0){tmp=0;}
        else {tmp=n[m-1];}
        tmp3 = n[m];
        for(j=1; j<p; j++){
            cor = 0.0;
            var = 0.0;
            // x_j = x + (tmp2+j)*N;
            x_j = x + j*N;
            for(i=tmp;i<tmp3;i++){
                tmp1 = x_j[i];
                cor += tmp1;
                var += tmp1*tmp1;
            }
            tmp1 = cor/(tmp3-tmp);
            sdx[tmp2+j] = sqrt((var-(tmp3-tmp)*tmp1*tmp1)/(tmp3-tmp-1));
        }
    }
}

int sgn(double x){
    double eps = 1e-8;
    if(x>eps){return 1;}
    else if(x<-eps){return -1;}
    else {return 0;}
}

double _loss(int *n, int M, int p, double q, double lambda2, double *r, double *beta){
    int i,m,tmp1,tmp2,tmp3,l,j,a;
    double tmp, tmp4, tmp5, loss1, loss2 = 0.0, penal;


    penal = 0.0;
    for(m=0;m<M;m++){
        tmp1 = m*p;
        if(m==0){tmp2=0;}
        else {tmp2=n[m-1];}
        tmp3 = n[m];
        loss1 = 0.0;
        for(i=tmp2;i<tmp3;i++){
            tmp = r[i];
            loss1 += tmp*(q-(tmp<0.0));
        }
        loss2 += loss1/(tmp3-tmp2); 
        // printf("%lf\n",r[n[m]-1]);
        // penal = 0.0;
        for(l=0;l<M;l++){
            tmp = 0.0;
            if(l!=m){
                tmp2 = l*p;
                for(j=1;j<p;j++){
                    tmp4 = beta[tmp1+j];
                    tmp5 = beta[tmp2+j];
                    a = 1;
                    if((sgn(tmp4)==1 && sgn(tmp5)==-1) || (sgn(tmp4)==-1 && sgn(tmp5)==1))  
                    {
                        a=0;
                    }
                    // a = sgn(tmp4)==sgn(tmp5) ? 1 : 0;
                    // if((beta[tmp3]==0.0 && beta[tmp4]!=0.0) || (beta[tmp3]!=0.0 && beta[tmp4]==0.0)){
                    //     a=1;
                    // }
                    tmp += a*(tmp4 - tmp5)*(tmp4-tmp5);
                }
            }
            penal += tmp;
        }
        // loss = loss + loss2 + lambda2*penal;
    }
    // printf("%lf\n",penal);
    loss2 += lambda2*penal;

    return loss2;
}

double _loss1(int n1, int n2, double q, double *r){
    int i;
    double tmp, loss = 0.0;

    for(i=n1;i<n2;i++){
        tmp = r[i];
        loss += tmp*(q-(tmp<0.0));
    }
    loss /= (n2 - n1); 

    return loss;
}

void der(double *x, double *sdx, int *n, int M, int p, double q, double lambda2, double *r, double delta, double *derivative, double *beta){
    int i,j,k,m,l,a,N=n[M-1],tmp2,tmp3,tmp4, temp1, temp2;
    double tmp,tmp1,temp,*x_j, *r1;
    r1 = (double*)malloc(sizeof(double) *N);

    tmp1 = 1.0-2.0*q;
    for(k=0;k<N;k++){
        tmp = r[k];
        r1[k] = tmp1 - tmp/(delta+fabs(tmp));
    }
    for(m=0;m<M;m++){
        tmp2 = m*p;
        if(m==0){tmp3=0;}
        else {tmp3=n[m-1];}
        tmp4 = n[m];
        for(j=0;j<p;j++){
            temp1 = tmp2 + j;
            temp = 0.0;
            if(j!=0){
                for(l=0;l<M;l++){
                    temp2 = l*p+j;
                    if(l!=m){
                        tmp = beta[temp1];
                        tmp1 = beta[temp2];
                        a = 1;
                        if((sgn(tmp)==1 && sgn(tmp1)==-1) || (sgn(tmp)==-1 && sgn(tmp1)==1))
                        {
                            a=0;
                        }
                        // a = sgn(beta[temp1])==sgn(beta[temp2]) ? 1 : 0;
                        // if((beta[temp1]==0.0 && beta[temp2]!=0.0) || (beta[temp1]!=0.0 && beta[temp2]==0.0)){
                        //     a=1;
                        // }
                        temp += a*(tmp - tmp1);
                    }
                    
                }
            }
            tmp = 0.0;
            // x_j = x + (tmp2+j)*N;
            x_j = x + j*N;
            for(i=tmp3;i<tmp4;i++){
                tmp += x_j[i]*r1[i];
            }
            tmp/=sdx[temp1];
            // printf("%lf\n",tmp);
            derivative[temp1] = tmp*0.5/(tmp4-tmp3) + 4*lambda2*temp;
        }
    }
    free(r1);
}

int estimate(double *y, double *x, double *sdx, int *active, int *n, int M, int p, double q, double lambda2, double epsilon, double delta, double xi, int max_iter, double *vbic, double *lambda1, double *betahat){
    int m, P=M*p, N=n[M-1], *num, j, row, column, count, i, nm, location,k, number, c, l, tmp1, tmp3, tmp5;
    double *r1, *derivative, *r2, loss1, value, *x_j, loss2, bic, temp, temp1, tmp, tmp2, tmp4, dactive;
    r1 = (double*)malloc(sizeof(double) *N);
    derivative = (double*)malloc(sizeof(double) *P);
    r2 = (double*)malloc(sizeof(double) *N);
    // betahat1 = (double*)calloc(max_iter*P, sizeof(double));
    // activeb = (int*)calloc(P, sizeof(int));
    num = (int*)calloc(M, sizeof(int));
    // printf("%d\n",M);
    preprocess(x, n, M, p, sdx);
    // printArrayDouble(sdx,P,20);

    for(m=0;m<M;m++){
        tmp1 = m*p;
        for(j=0;j<p;j++){
            betahat[(tmp1+j)] = 0.0;
        }
    }
    
    der(x, sdx, n, M, p, q, lambda2, y, delta, derivative, betahat);
    // printArrayDouble(derivative,P,20);
    tmp2 = 0.0;
    for(m=0;m<M;m++){
        tmp1 = m*p;
        for(j=0;j<p;j++){
            value = fabs(derivative[tmp1+j]);
            if(value>tmp2){
                tmp2 = value;
                row = j;
                column = m;
            }
        }
    }
    // printArrayDouble(y,n[M-1],20);

    loss1 = _loss(n, M, p, q, lambda2, y, betahat);
    // printArrayDouble(y,n[M-1],20);
    // printf("%lf\n",loss1);
    count = column*p+row;
    // printf("%d\n",count);
    active[0] = count;
    num[column] += 1;
    number = 1;
    value = epsilon;
    if(derivative[count] > 0.0){value *= -1.0;}
    // printf("%lf\n",value);
    betahat[count] = value;
    // betahat1[count*max_iter] = value;
    temp1 = value/sdx[count];
    x_j = x + row*N;
    for(m=0;m<M;m++){
        if(m==0){tmp1=0;}
        else{tmp1 = n[m-1];}
        if(m==column){
            for(i=tmp1;i<n[m];i++){
                r2[i] = y[i] - x_j[i]*temp1;
            }
        }
        else{
            for(i=tmp1;i<n[m];i++){
                r2[i] = y[i];
            }
        }
    }
    loss2 = _loss(n, M, p, q, lambda2, r2, betahat);
    lambda1[0] = (loss1-loss2)/epsilon;
    // printf("%lf\n",loss2);
    bic = 0.0;
    for(m=0;m<M;m++){
        if(m==0){tmp5=0;}
        else{tmp5 = n[m-1];}
        tmp1 = n[m];
        tmp3 = num[m];
        nm = tmp1 - tmp5;
        bic = bic + 2*nm*_loss1(tmp5, tmp1, q, r2) + tmp3*log(nm) + 2.0*lchoose(p-1, tmp3);
    }
    temp = bic;
    location = 0;
    vbic[0] = bic;

    for(k=0;k<max_iter;k++){
        c = 0;
        for(l=0;l<number;l++){
            if((active[l]%p)!=0){c++;break;}
        }
        if(c==0){
            der(x, sdx, n, M, p, q, lambda2, r2, delta, derivative, betahat);
            loss1 = _loss(n, M, p, q, lambda2, r2, betahat);
            tmp = 0.0;
            for(m=0;m<M;m++){
                tmp1 = m*p;
                for(j=0;j<p;j++){
                    value = fabs(derivative[tmp1+j]);
                    if(value>tmp){
                        tmp = value;
                        row = j;
                        column = m;
                    }
                }
            }
            count = column*p+row;
            if(row!=0){
                active[number] = count;
                num[column] += 1;
                number += 1;
            }
            value = epsilon;
            if(derivative[count] > 0.0){value *= -1.0;}
            betahat[count] += value;
            temp1 = value/sdx[count];
            x_j = x + row*N;
            for(m=0;m<M;m++){
                if(m==0){tmp1=0;}
                else{tmp1 = n[m-1];}
                if(m==column){
                    for(i=tmp1;i<n[m];i++){
                        r1[i] = r2[i] - x_j[i]*temp1;
                        r2[i] = r1[i];
                    }
                }
                else{
                    for(i=tmp1;i<n[m];i++){
                        r1[i] = r2[i];
                    }
                }
            }
            loss2 = _loss(n, M, p, q, lambda2, r1, betahat);
            tmp = (loss1-loss2)/epsilon;
            tmp2 = lambda1[k];
            lambda1[k+1]    = tmp < tmp2 ? tmp : tmp2;
        }
        else{
            der(x, sdx, n, M, p, q, lambda2, r2, delta, derivative, betahat);
            loss1 = _loss(n, M, p, q, lambda2, r2, betahat);
            tmp = 1e14;
            for(l=0;l<number;l++){
                if((active[l]%p)!=0){
                    m = active[l];
                    dactive = derivative[m];
                    if(betahat[m]>0.0){dactive *= -1.0;}
                    if(dactive<tmp){tmp = dactive; count = m;c=l;}
                }
            }
            value = epsilon;
            tmp = betahat[count];
            if(tmp>0.0){value*=-1.0;}
            betahat[count] = tmp + value;
            temp1 = value/sdx[count];
            row = count%p;
            column = (count - row)/p;
            x_j = x + row*N;
            for(m=0;m<M;m++){
                if(m==0){tmp1=0;}
                else{tmp1 = n[m-1];}
                if(m==column){
                    for(i=tmp1;i<n[m];i++){
                        r1[i] = r2[i] - x_j[i]*temp1;
                    }
                }
                else{
                    for(i=tmp1;i<n[m];i++){
                        r1[i] = r2[i];
                    }
                }
            }
            loss2 = _loss(n, M, p, q, lambda2, r1, betahat);
            tmp2 = lambda1[k];
            if(loss2-loss1-tmp2*epsilon < -xi){
                if(fabs(betahat[count]) < xi){
                    betahat[count] = 0.0;
                    for(l=c;l<number;l++){
                        active[l] = active[l+1];
                    }
                    number--;
                    num[column]--;
                }
                for(i=0;i<N;i++){
                    r2[i] = r1[i];
                }
                lambda1[k+1] = tmp2;
            }
            else{
                betahat[count] -= value;
                tmp4 = 0.0;
                for(m=0;m<M;m++){
                    tmp1 = m*p;
                    for(j=0;j<p;j++){
                        value = fabs(derivative[tmp1+j]);
                        if(value>tmp4){
                            tmp4 = value;
                            row = j;
                            column = m;
                        }
                    }
                }
                count = column*p+row;
                c = 0;
                for(l=0;l<number;l++){
                    if(count == active[l]){c++;break;}
                }
                if(c==0){
                    active[number] = count;
                    number +=1;
                    num[column] += 1;
                }

                value = epsilon;
                if(derivative[count]>0.0){value *= -1.0;}
                betahat[count] += value;
                temp1 = value/sdx[count];
                x_j = x + row*N;
                for(m=0;m<M;m++){
                    if(m==0){tmp1=0;}
                    else{tmp1 = n[m-1];}
                    if(m==column){
                        for(i=tmp1;i<n[m];i++){
                            r1[i] = r2[i] - x_j[i]*temp1;
                            r2[i] = r1[i];
                        }
                    }
                    else{
                        for(i=tmp1;i<n[m];i++){
                            r1[i] = r2[i];
                        }
                    }
                }
                loss2 = _loss(n, M, p, q, lambda2, r1, betahat);
                tmp = (loss1-loss2)/epsilon;
                lambda1[k+1]    = tmp < tmp2 ? tmp : tmp2;
            }
        }
        // printf("%lf\n",loss2);
        // for(j=0;j<P;j++){
        //     betahat1_j = betahat1 + j*max_iter;
        //     tmp2 = betahat1_j[k];
        //     if(tmp2!=0.0){
        //         betahat1_j[k+1] = tmp2;
        //     }
        // }
        // betahat1[count*max_iter+k+1] = betahat[count];
        for(l=0;l<number;l++){
            tmp1 = active[l];
            if((tmp1%p)==0){
                m = tmp1/p;
                num[m]--;
            }
        }
        bic = 0.0;
        for(m=0;m<M;m++){
            if(m==0){tmp1=0;}
            else{tmp1 = n[m-1];}
            tmp3 = n[m];
            tmp5 = num[m];
            nm = tmp3 - tmp1;
            bic = bic + 2*nm*_loss1(tmp1, tmp3, q, r2) + tmp5*log(nm) + 2.0*lchoose(p-1, tmp5);
        }
        // printf("%lf\n",bic);
        // printf("%d\n",num[M-1]);
        for(l=0;l<number;l++){
            tmp1 = active[l];
            if((tmp1%p)==0){
                m = tmp1/p;
                num[m]++;
            }
        }
        vbic[k+1] = bic;
        if(bic<temp){
            temp = bic;
            location = k;
        }
        if(lambda1[k+1] < xi){break;}
    }
    // for(j=0;j<P;j++){
    //     betahat1_j = betahat1 + j*max_iter;
    //     betahat[j] = betahat1_j[location+1];
    //     betahat[j] /= sdx[j];
    // }

    free(r1);
    free(derivative);
    free(r2);
    free(num);

    return(k);
}

int estimate1(double *y, double *x, double *sdx, int *active, int *n, int M, int p, double q, double lambda2, double epsilon, double delta, double xi, int max_iter, double *vbic, double *lambda1, double *betahat){
    int m, P=M*p, N=n[M-1], *num, j, row, column, count, i, nm, location,k, number, c, l, tmp1, tmp3, tmp5;
    double *r1, *derivative, *r2, *betahat1, loss1, value, *x_j, loss2, bic, temp, temp1, tmp, tmp2, tmp4, dactive, *betahat1_j;
    r1 = (double*)malloc(sizeof(double) *N);
    derivative = (double*)malloc(sizeof(double) *P);
    r2 = (double*)malloc(sizeof(double) *N);
    betahat1 = (double*)calloc(max_iter*P, sizeof(double));
    // activeb = (int*)calloc(P, sizeof(int));
    num = (int*)calloc(M, sizeof(int));

    preprocess(x, n, M, p, sdx);
    // printf("%d\n",M);
    // printArrayDouble(sdx,P,20);

    for(m=0;m<M;m++){
        tmp1 = m*p;
        for(j=0;j<p;j++){
            betahat[(tmp1+j)] = 0.0;
        }
    }
    
    der(x, sdx, n, M, p, q, lambda2, y, delta, derivative, betahat);
    // printArrayDouble(derivative,P,20);
    tmp2 = 0.0;
    for(m=0;m<M;m++){
        tmp1 = m*p;
        for(j=0;j<p;j++){
            value = fabs(derivative[tmp1+j]);
            if(value>tmp2){
                tmp2 = value;
                row = j;
                column = m;
            }
        }
    }
    // printArrayDouble(y,n[M-1],20);

    loss1 = _loss(n, M, p, q, lambda2, y, betahat);
    // printArrayDouble(y,n[M-1],20);
    // printf("%lf\n",loss1);
    count = column*p+row;
    // printf("%d\n",count);
    active[0] = count;
    num[column] += 1;
    number = 1;
    value = epsilon;
    if(derivative[count] > 0.0){value *= -1.0;}
    // printf("%lf\n",value);
    betahat[count] = value;
    betahat1[count*max_iter] = value;
    temp1 = value/sdx[count];
    x_j = x + row*N;
    for(m=0;m<M;m++){
        if(m==0){tmp1=0;}
        else{tmp1 = n[m-1];}
        if(m==column){
            for(i=tmp1;i<n[m];i++){
                r2[i] = y[i] - x_j[i]*temp1;
            }
        }
        else{
            for(i=tmp1;i<n[m];i++){
                r2[i] = y[i];
            }
        }
    }
    loss2 = _loss(n, M, p, q, lambda2, r2, betahat);
    lambda1[0] = (loss1-loss2)/epsilon;
    // printf("%lf\n",loss2);
    bic = 0.0;
    for(m=0;m<M;m++){
        if(m==0){tmp5=0;}
        else{tmp5 = n[m-1];}
        tmp1 = n[m];
        tmp3 = num[m];
        nm = tmp1 - tmp5;
        bic = bic + 2*nm*_loss1(tmp5, tmp1, q, r2) + tmp3*log(nm) + 2.0*lchoose(p-1, tmp3);
    }
    temp = bic;
    location = 0;
    vbic[0] = bic;

    for(k=0;k<max_iter;k++){
        c = 0;
        for(l=0;l<number;l++){
            if((active[l]%p)!=0){c++;break;}
        }
        if(c==0){
            der(x, sdx, n, M, p, q, lambda2, r2, delta, derivative, betahat);
            loss1 = _loss(n, M, p, q, lambda2, r2, betahat);
            tmp = 0.0;
            for(m=0;m<M;m++){
                tmp1 = m*p;
                for(j=0;j<p;j++){
                    value = fabs(derivative[tmp1+j]);
                    if(value>tmp){
                        tmp = value;
                        row = j;
                        column = m;
                    }
                }
            }
            count = column*p+row;
            if(row!=0){
                active[number] = count;
                num[column] += 1;
                number += 1;
            }
            value = epsilon;
            if(derivative[count] > 0.0){value *= -1.0;}
            betahat[count] += value;
            temp1 = value/sdx[count];
            x_j = x + row*N;
            for(m=0;m<M;m++){
                if(m==0){tmp1=0;}
                else{tmp1 = n[m-1];}
                if(m==column){
                    for(i=tmp1;i<n[m];i++){
                        r1[i] = r2[i] - x_j[i]*temp1;
                        r2[i] = r1[i];
                    }
                }
                else{
                    for(i=tmp1;i<n[m];i++){
                        r1[i] = r2[i];
                    }
                }
            }
            loss2 = _loss(n, M, p, q, lambda2, r1, betahat);
            tmp = (loss1-loss2)/epsilon;
            tmp2 = lambda1[k];
            lambda1[k+1]    = tmp < tmp2 ? tmp : tmp2;
        }
        else{
            der(x, sdx, n, M, p, q, lambda2, r2, delta, derivative, betahat);
            loss1 = _loss(n, M, p, q, lambda2, r2, betahat);
            tmp = 1e14;
            for(l=0;l<number;l++){
                if((active[l]%p)!=0){
                    m = active[l];
                    dactive = derivative[m];
                    if(betahat[m]>0.0){dactive *= -1.0;}
                    if(dactive<tmp){tmp = dactive; count = m;c=l;}
                }
            }
            value = epsilon;
            tmp = betahat[count];
            if(tmp>0.0){value*=-1.0;}
            betahat[count] = tmp + value;
            temp1 = value/sdx[count];
            row = count%p;
            column = (count - row)/p;
            x_j = x + row*N;
            for(m=0;m<M;m++){
                if(m==0){tmp1=0;}
                else{tmp1 = n[m-1];}
                if(m==column){
                    for(i=tmp1;i<n[m];i++){
                        r1[i] = r2[i] - x_j[i]*temp1;
                    }
                }
                else{
                    for(i=tmp1;i<n[m];i++){
                        r1[i] = r2[i];
                    }
                }
            }
            loss2 = _loss(n, M, p, q, lambda2, r1, betahat);
            tmp2 = lambda1[k];
            if(loss2-loss1-tmp2*epsilon < -xi){
                if(fabs(betahat[count]) < xi){
                    betahat[count] = 0.0;
                    for(l=c;l<number;l++){
                        active[l] = active[l+1];
                    }
                    number--;
                    num[column]--;
                }
                for(i=0;i<N;i++){
                    r2[i] = r1[i];
                }
                lambda1[k+1] = tmp2;
            }
            else{
                betahat[count] -= value;
                tmp4 = 0.0;
                for(m=0;m<M;m++){
                    tmp1 = m*p;
                    for(j=0;j<p;j++){
                        value = fabs(derivative[tmp1+j]);
                        if(value>tmp4){
                            tmp4 = value;
                            row = j;
                            column = m;
                        }
                    }
                }
                count = column*p+row;
                c = 0;
                for(l=0;l<number;l++){
                    if(count == active[l]){c++;break;}
                }
                if(c==0){
                    active[number] = count;
                    number +=1;
                    num[column] += 1;
                }

                value = epsilon;
                if(derivative[count]>0.0){value *= -1.0;}
                betahat[count] += value;
                temp1 = value/sdx[count];
                x_j = x + row*N;
                for(m=0;m<M;m++){
                    if(m==0){tmp1=0;}
                    else{tmp1 = n[m-1];}
                    if(m==column){
                        for(i=tmp1;i<n[m];i++){
                            r1[i] = r2[i] - x_j[i]*temp1;
                            r2[i] = r1[i];
                        }
                    }
                    else{
                        for(i=tmp1;i<n[m];i++){
                            r1[i] = r2[i];
                        }
                    }
                }
                loss2 = _loss(n, M, p, q, lambda2, r1, betahat);
                tmp = (loss1-loss2)/epsilon;
                lambda1[k+1]    = tmp < tmp2 ? tmp : tmp2;
            }
        }
        // printf("%lf\n",loss2);
        for(j=0;j<P;j++){
            betahat1_j = betahat1 + j*max_iter;
            tmp2 = betahat1_j[k];
            if(tmp2!=0.0){
                betahat1_j[k+1] = tmp2;
            }
        }
        betahat1[count*max_iter+k+1] = betahat[count];
        for(l=0;l<number;l++){
            tmp1 = active[l];
            if((tmp1%p)==0){
                m = tmp1/p;
                num[m]--;
            }
        }
        bic = 0.0;
        for(m=0;m<M;m++){
            if(m==0){tmp1=0;}
            else{tmp1 = n[m-1];}
            tmp3 = n[m];
            tmp5 = num[m];
            nm = tmp3 - tmp1;
            bic = bic + 2*nm*_loss1(tmp1, tmp3, q, r2) + tmp5*log(nm) + 2.0*lchoose(p-1, tmp5);
        }
        // printf("%lf\n",bic);
        // printf("%d\n",num[M-1]);
        for(l=0;l<number;l++){
            tmp1 = active[l];
            if((tmp1%p)==0){
                m = tmp1/p;
                num[m]++;
            }
        }
        vbic[k+1] = bic;
        if(bic<temp){
            temp = bic;
            location = k;
        }
        if(lambda1[k+1] < xi){break;}
    }
    for(j=0;j<P;j++){
        betahat1_j = betahat1 + j*max_iter;
        betahat[j] = betahat1_j[location+1];
        betahat[j] /= sdx[j];
    }

    free(r1);
    free(derivative);
    free(r2);
    free(betahat1);
    free(num);

    return(k);
}

SEXP _ESTIMATE(SEXP Y, SEXP X, SEXP N, SEXP M, SEXP P, SEXP Q, SEXP LAMBDA2, SEXP EPSILON, SEXP DELTA, SEXP XI, SEXP MAX_ITER){
    int p = INTEGER(P)[0],i, m=INTEGER(M)[0], dimension = m*p;
    int max_iter = INTEGER(MAX_ITER)[0], iter, *active;
    double *sdx;

    // active = (int*)malloc(sizeof(int) *p);
    // meanx = (double*)malloc(sizeof(double) *p);
    sdx = (double*)malloc(sizeof(double) *dimension);
    active = (int*)malloc(sizeof(int) *max_iter);

    SEXP BETAHAT, LAMBDA, BIC, Iter;
    SEXP Result, R_names;
    char *names[3] = { "iter", "bic", "lambda"};
    // printf("%s\n",names[0]);

    // PROTECT(ACTIVE = allocVector(INTSXP, dimension));
    PROTECT(BETAHAT = allocVector(REALSXP, dimension));
    PROTECT(LAMBDA = allocVector(REALSXP, max_iter));
    PROTECT(BIC = allocVector(REALSXP, max_iter));
    PROTECT(Iter = allocVector(INTSXP, 1));
    // PROTECT(Location = allocVector(REALSXP, 1));
    PROTECT(Result = allocVector(VECSXP, 3));
    PROTECT(R_names   = allocVector(STRSXP,  3));
    
    for(i = 0; i < 3; ++i) SET_STRING_ELT(R_names, i,  mkChar(names[i]));
    // printf("%d\n",m);
    
    iter = estimate(REAL(Y), REAL(X), sdx, active, INTEGER(N), INTEGER(M)[0], p, REAL(Q)[0], REAL(LAMBDA2)[0], REAL(EPSILON)[0], REAL(DELTA)[0], REAL(XI)[0], max_iter, REAL(BIC), REAL(LAMBDA), REAL(BETAHAT));
    // printf("%d\n",iter);
    INTEGER(Iter)[0] = iter;
    // SET_VECTOR_ELT(Result, 0, LAMBDA);
    // SET_VECTOR_ELT(Result, 0, BETAHAT);
    SET_VECTOR_ELT(Result, 0, Iter);
    SET_VECTOR_ELT(Result, 1, BIC);
    SET_VECTOR_ELT(Result, 2, LAMBDA);
    setAttrib(Result, R_NamesSymbol, R_names);

    free(active);
    // free(meanx);
    free(sdx);

    UNPROTECT(6);
    return Result;
}

SEXP ESTIMATE_1(SEXP Y, SEXP X, SEXP N, SEXP M, SEXP P, SEXP Q, SEXP LAMBDA2, SEXP EPSILON, SEXP DELTA, SEXP XI, SEXP MAX_ITER){
    int p = INTEGER(P)[0],i, m=INTEGER(M)[0], dimension = m*p;
    int max_iter = INTEGER(MAX_ITER)[0], iter, *active;
    double *sdx;

    // active = (int*)malloc(sizeof(int) *p);
    // meanx = (double*)malloc(sizeof(double) *p);
    sdx = (double*)malloc(sizeof(double) *dimension);
    active = (int*)malloc(sizeof(int) *max_iter);

    SEXP BETAHAT, LAMBDA, BIC, Iter;
    SEXP Result, R_names;
    char *names[4] = {"betahat", "iter", "bic", "lambda"};
    // printf("%s\n",names[0]);

    // PROTECT(ACTIVE = allocVector(INTSXP, dimension));
    PROTECT(BETAHAT = allocVector(REALSXP, dimension));
    PROTECT(LAMBDA = allocVector(REALSXP, max_iter));
    PROTECT(BIC = allocVector(REALSXP, max_iter));
    PROTECT(Iter = allocVector(INTSXP, 1));
    // PROTECT(Location = allocVector(REALSXP, 1));
    PROTECT(Result = allocVector(VECSXP, 4));
    PROTECT(R_names   = allocVector(STRSXP,  4));
    
    for(i = 0; i < 4; ++i) SET_STRING_ELT(R_names, i,  mkChar(names[i]));
    
    iter = estimate1(REAL(Y), REAL(X), sdx, active, INTEGER(N), INTEGER(M)[0], p, REAL(Q)[0], REAL(LAMBDA2)[0], REAL(EPSILON)[0], REAL(DELTA)[0], REAL(XI)[0], max_iter, REAL(BIC), REAL(LAMBDA), REAL(BETAHAT));
    // printf("%d\n",iter);
    INTEGER(Iter)[0] = iter;
    // SET_VECTOR_ELT(Result, 0, LAMBDA);
    SET_VECTOR_ELT(Result, 0, BETAHAT);
    SET_VECTOR_ELT(Result, 1, Iter);
    SET_VECTOR_ELT(Result, 2, BIC);
    SET_VECTOR_ELT(Result, 3, LAMBDA);
    setAttrib(Result, R_NamesSymbol, R_names);

    free(active);
    // free(meanx);
    free(sdx);

    UNPROTECT(6);
    return Result;
}
