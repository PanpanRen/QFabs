#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
// #include "testhb.h"

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

void sortN(int *ind, double *x, int n, int dd){
    int i, j, MaxInd, d;
    double tmp;

    d = (dd==n?dd-1:dd);
    for(i=0; i<d; i++){
        tmp = x[0]; MaxInd = ind[0];
        for(j=1; j<n-i; j++)
        {
            if(x[j]<tmp){
                x[j-1]      = x[j];
                x[j]        = tmp;
                ind[j-1]    = ind[j];
                ind[j]      = MaxInd;
            }
            else{
                tmp     = x[j];
                MaxInd  = ind[j];
            }
        }
    }
}

double calculate_bic(double loss, int *param, double gamma){
    int n       = param[0];
    int p       = param[1];
    int number  = param[2];

    return 2.0*n*(loss) + number*log(n) + 2.0*gamma*lchoose(p, number);
}

void preprocess(double *x, int n, int p, double *sdx){
	int i, j;
	double cor, var, *x_j, tmp1;

	for(j=0; j<p; j++){
		cor = 0.0;
		var = 0.0;
		x_j = x + j*n;
		for(i=0; i<n; i++){
			tmp1 = x_j[i];
			cor += tmp1;
			var += tmp1*tmp1;
		}
		tmp1    = cor/n;
		sdx[j]  = sqrt((var-n*tmp1*tmp1)/(n-1));
	}
}

double _loss(int n, double tau, double *r){
	int i;
	double tmp, loss1 = 0.0;

	for(i=0; i<n; i++){
		tmp     = r[i];
		loss1  += tmp*(tau-(tmp<0.0));
	}
	loss1 /= n;

	return loss1;
}

double _Intercept(double *r, int n, double tau, double delta){
    int i;
    double tmp, tmp1 = 0.0, tmp2 = 0.0, tmp3;
    for(i=0; i<n; i++){
        tmp = r[i];
        tmp3 = delta+fabs(tmp);
        tmp1 += tmp/tmp3;
        tmp2 += 1.0/tmp3;
    }
    return (tmp1-(1.0-2.0*tau)*n)/tmp2;
}

void der(double *x, double *sdx, int n, int p, double tau, double *r, double delta, double *derivative, double *d, double *residual){
    int i, j;
    double tmp = 0.0, tmp1, *x_j;

    tmp1 = 1.0-2.0*tau;
	for(i=0; i<n; i++){
		tmp         = r[i];
		residual[i] = tmp1 - tmp/(delta + fabs(tmp));
	}

    for(j=0; j<p; j++){
        tmp = 0.0;
        x_j = x + j*n;
        for(i=0; i<n; i++){
            tmp += residual[i]*x_j[i];
        }
        tmp             *= 0.5/(n*sdx[j]);
        derivative[j]    = tmp;
        d[j]             = fabs(tmp);
    }
}

void der1(double *x, double *sdx, int n, double tau, double *r, double delta, double *derivative, double *d1, double *residual, int max_s, int *index_max_s){
    int i, j, m;
    double tmp = 0.0, tmp1, *x_j;

    tmp1 = 1.0-2.0*tau;
	for(i=0; i<n; i++){
		tmp         = r[i];
		residual[i] = tmp1 - tmp/(delta + fabs(tmp));
	}

    for(m=0; m<max_s; m++){
        tmp = 0.0;
        j   = index_max_s[m];
        x_j = x + j*n;
        for(i=0; i<n; i++){
            tmp += residual[i]*x_j[i];
        }
        tmp            *= 0.5/(n*sdx[j]);
        derivative[m]   = tmp;
        d1[m]           = fabs(tmp);
    }
}

void Initial(double *x, double *sdx, int *param, double tau, double delta, double epsilon, double gamma, double *loss, double *lambda, 
            int *direction, int *active, double *bic, int *df, double *derivative, double *d, double *beta, double *r, double *residual, int max_s, int *index_max_s)
{
    int i, j, m, n = *param, p = param[1], ind[p];
    double tmp, tmp1, value, *x_j;

    der(x, sdx, n, p, tau, r, delta, derivative, d, residual);
    for(j=0; j<p; j++){
        ind[j] = j;
    }

    sortN(ind, d, p, max_s);
    for(m=0; m<max_s; m++){
        index_max_s[m] = ind[p-1-m];
    }

    m       = ind[p-1];
    value   = epsilon;
    if(derivative[m] > 0.0){value *= -1.0;}
    *beta = value;

    tmp     = _loss(n, tau, r);
    x_j     = x + m*n;
    tmp1    = value/sdx[m];
    for(i=0; i<n; i++){
		r[i] -= x_j[i]*tmp1;
	}
    value       = _loss(n, tau, r);
    *loss       = value;
    *lambda     = (tmp - value)/epsilon;
    *direction  = 1;
    param[2]    = 1; // param[0] = n; param[1] = p; param[2] = active
    *active     = m;
    *bic        = calculate_bic(value, param, gamma);
    *df         = 1;
}

void usi(int *s_i, int *s_j, int *tt_b, int *tt_a, int *act, int ns, int iter)
{
    int i;
    s_i += *tt_a;
    s_j += *tt_a;
    for (i = 0; i < ns; ++i)
    {
        s_i[i] = act[i];
        s_j[i] = iter;
    }
    *tt_b  = *tt_a;
    *tt_a += ns;
}

int Backward(double *x, double *sdx, int *param, double tau, double delta, double epsilon, double gamma, double xi, double *loss, double *lambda, 
            int *direction, int *active, double *bic, int *df, double *derivative, double *d, double *derivative1, double *beta, double *betanew, double *r, double *residual, int max_s, int *index_max_s, int k)
{
    int ind1, i, j, l = 0, m, count = 0, count1, n = *param, p = param[1], number = param[2], c = -1;
    double value, *x_j, tmp, tmp1 = 0.0, d1[max_s], temp1;

    value = 1e14;
    if((k+1)%20 != 0){
        der1(x, sdx, n, tau, r, delta, derivative1, d1, residual, max_s, index_max_s);

        for(l=0; l<number; l++){
            tmp         = beta[l];
            betanew[l]  = tmp;
            ind1        = active[l];
            for(m=0; m<max_s; m++){
                if(index_max_s[m] == ind1){
                    tmp1 = derivative1[m];
                    break;
                }
            }
            if(tmp > 0.0){tmp1 *= -1.0;}
            if(tmp1 < value){value = tmp1; count = ind1; c = l;}
        }

        tmp     = *d1;
        ind1    = 0;
        for(m=1; m<max_s; m++){
            value = d1[m];
            if(value > tmp){
                tmp     = value;
                ind1    = m;
            }
        }
        count1                  = index_max_s[ind1];
        index_max_s[ind1]       = index_max_s[max_s-1];
        index_max_s[max_s-1]    = count1;
        derivative1[max_s-1]    = derivative1[ind1];
    }
    else{
        der1(x, sdx, n, tau, r, delta, derivative1, d1, residual, max_s, index_max_s);
        sortN(index_max_s, d1, max_s, max_s);
        der(x, sdx, n, p, tau, r, delta, derivative, d, residual);
        for(j=0;j<p;j++){
            tmp = d[j];
            if(tmp>d1[0]){
                c = 0;
                for(m=0;m<max_s;m++){
                    if(j==index_max_s[m]){c++;break;}
                }
                if(c==0){
                    if(d1[max_s-1]<=tmp){
                        l = max_s-1;
                        for(i=0;i<max_s-1;i++){
                            d1[i]           = d1[i+1];
                            index_max_s[i]  = index_max_s[i+1];
                        }
                    }
                    else{
                        for(m=1;m<max_s;m++){
                            if(d1[m]>tmp){
                                l = m;
                                for(i=0;i<m;i++){
                                    d1[i]           = d1[i+1];
                                    index_max_s[i]  = index_max_s[i+1];
                                }
                                break;
                            }
                        }
                    }
                    index_max_s[l]  = j;
                    d1[l]           = tmp;
                }
            }
        }
        count1                  = index_max_s[max_s-1];
        derivative1[max_s-1]    = derivative[count1];

        for(l=0; l<number; l++){
            tmp         = beta[l];
            betanew[l]  = tmp;
            ind1        = active[l];
            tmp1        = derivative[ind1];
            if(tmp > 0.0){tmp1 *= -1.0;}
            if(tmp1 < value){value = tmp1; count = ind1; c = l;}
        }
    }

    value   = epsilon;
    tmp1    = beta[c];
    if (tmp1 > 0) value *= -1.0;
    betanew[c] += value;
    temp1       = value/sdx[count];
    x_j         = x + count*n;
    for(ind1=0; ind1<n; ind1++){
        residual[ind1] = r[ind1] - x_j[ind1]*temp1;
    }
    tmp = _loss(n, tau, residual);

    if (tmp - *loss - (*lambda)*epsilon < -1.0*xi){
        loss[1]         = tmp;
        direction[1]    = -1;
        df[1]           = *df;
        lambda[1]       = *lambda;
        for (ind1 = 0; ind1 < n; ++ind1)
            r[ind1] = residual[ind1];
        if(fabs(betanew[c]) < xi){
            for(l=c; l<number-1; l++){ 
                active[l]   = active[l+1];
                betanew[l]  = betanew[l+1];
            }
            param[2]--;
            df[1]--;
        }
        bic[1] = calculate_bic(tmp, param, gamma);
        return 0;
    }

    betanew[c] = tmp1;
    return 1;
}

void Forward(double *x, double *sdx, int *param, double tau, double epsilon, double gamma, double *loss, double *lambda, 
            int *direction, int *active, double *bic, int *df, double *derivative1, double *beta, double *r, int max_s, int *index_max_s)
{
    int i, l, m, n = *param, number = param[2], count1 = 0;
    double tmp, value, *x_j;

    count1  = index_max_s[max_s-1];
    value   = epsilon;
    if(derivative1[max_s-1] > 0.0){value *= -1.0;}

    m = 0;
    for(l=0; l<number; l++){
        if(count1 == active[l]){
            m++;
            beta[l] += value;
            df[1]    = *df;
            break;}
    }
    if(m==0){
        active[number]  = count1;
        beta[number]    = value;
        param[2]       += 1;
        df[1]           = *df + 1;
    }
    x_j = x + count1*n;
    tmp = value/sdx[count1];
    for(i=0; i<n; i++){
        r[i] -= x_j[i]*tmp;
    }
    value           = _loss(n, tau, r);
    loss[1]         = value;
    bic[1]          = calculate_bic(value, param, gamma);
    tmp             = (*loss - value)/epsilon;
    lambda[1]       = tmp < *lambda ? tmp : *lambda;
    direction[1]    = 1;
}

int quantile_EST(double *x, double *sdx, int *param, double tau, double delta, double epsilon, double gamma, double xi, int *max_iter, double *loss, double *lambda, int *direction, 
                int *active, double *bic, int *df, double *derivative, double *d, double *derivative1, int *sparse_i, int *sparse_j, double lam_m, int max_s, double *beta, double *r, double *residual, int *index_max_s, double *intercept, int update)
{
    int i, k, l, n = param[0], p = param[1];
    double tmp;

    preprocess(x, n, p, sdx);
    
    Initial(x, sdx, param, tau, delta, epsilon, gamma, loss, lambda, direction, active, bic, df, derivative, d, beta, r, residual, max_s, index_max_s);
    int tt_act_b = 0;
    int tt_act_a = 0;
    usi(sparse_i, sparse_j, &tt_act_b, &tt_act_a, active, param[2], 0);

    for(k=0; k<*max_iter-1; k++){
        l = Backward(x, sdx, param, tau, delta, epsilon, gamma, xi, loss+k, lambda+k, 
                    direction+k, active, bic+k, df+k, derivative, d, derivative1, beta+tt_act_b, beta+tt_act_a, r, residual, max_s, index_max_s, k);
        if((k+1)%update==0){
            tmp = _Intercept(r, n, tau, delta);
            for(i=0; i<n; i++){
                r[i] -= tmp;
            }
            *intercept += tmp;
            loss[k] = _loss(n, tau, r);
            l = 1;
        }
        if(l){
            Forward(x, sdx, param, tau, epsilon, gamma, loss+k, lambda+k, direction+k, active, bic+k, df+k, derivative1, beta+tt_act_a, r, max_s, index_max_s);
        }
        usi(sparse_i, sparse_j, &tt_act_b, &tt_act_a, active, param[2], k+1);

        tmp = lambda[k+1];
        if ((tmp <= lambda[0] * lam_m) ) {
            *max_iter = k+2;
            if (tmp < 0) 
            {
                (*max_iter)--;
                tt_act_a -= param[2];
            }
            break;
        }
        if (param[2] > max_s) {
            Rprintf("Warning! Max nonzero number is larger than predetermined threshold. Program ended early.\n");
            *max_iter = k+2;
            break;
        }
        if (k == *max_iter-2) {
            Rprintf("Solution path unfinished, more iterations are needed.\n");
            break;
        }
    }

	return tt_act_a;
}

SEXP _QFABS(SEXP Y, SEXP X, SEXP PARA_INT_, SEXP PARA_DOUBLE_){
    double *sdx, *loss, *lambda, *bic, *derivative, *d, *derivative1, *beta, *residual, *r, intercept = 0.0, tmp;
    int *direction, *active, *df, *sparse_i, *sparse_j, *index_max_s, temp;
	int i, tt_act_a, para[3];

	int n     	= INTEGER(PARA_INT_)[0];
	int p     	= INTEGER(PARA_INT_)[1];
	int iter	= INTEGER(PARA_INT_)[2];
    int max_s   = INTEGER(PARA_INT_)[3];
    int update  = INTEGER(PARA_INT_)[4];
    para[0] = n;
    para[1] = p;
    para[2] = 0;

    double *para1 	= REAL(PARA_DOUBLE_);
    double tau      = para1[0];
	double eps    	= para1[1];
    double delta    = para1[2];
	double xi		= para1[3];
    double lam_m    = para1[4];
    double gamma    = para1[5];

    sdx         = (double*)malloc(sizeof(double)*p);
    loss        = (double*)malloc(sizeof(double)*iter);
    lambda      = (double*)malloc(sizeof(double)*iter);
    direction   =    (int*)malloc(sizeof(int)   *iter);
    // active    =    (int*)calloc(max_s+2, sizeof(int));
    active      = (int*)malloc(sizeof(int) *(max_s+2));
    bic         = (double*)malloc(sizeof(double)*iter);
    df          =    (int*)malloc(sizeof(int)   *iter);
    derivative  = (double*)malloc(sizeof(double)  *p);
    d           = (double*)malloc(sizeof(double)  *p);
    derivative1 = (double*)malloc(sizeof(double)  *max_s);
    sparse_i    =    (int*)calloc(iter*max_s, sizeof(int));
    sparse_j    =    (int*)calloc(iter*max_s, sizeof(int));
    beta        = (double*)malloc(sizeof(double)*iter*max_s);
    residual    = (double*)malloc(sizeof(double)  *n);
    index_max_s = (int*)malloc(sizeof(int) *(max_s));
    r           = (double*)malloc(sizeof(double)*n);

    for(i=0; i<n; i++){
        r[i] = REAL(Y)[i];
    }

    tt_act_a = quantile_EST(REAL(X), sdx, para, tau, delta, eps, gamma, xi, &iter, loss, lambda, direction, active, bic, df,
                            derivative, d, derivative1, sparse_i, sparse_j, lam_m, max_s, beta, r, residual, index_max_s, &intercept, update);
    tmp = _Intercept(r, n, tau, delta);
    intercept += tmp;

    SEXP Beta, Lambda, Direction, Loops, BIC, Loss, Df, Indexi, Indexj;
    SEXP Intercept, Result, R_names;
    char *names[10] = {"beta", "lambda", "direction", "iter", "bic", "loss", 
    "df", "index_i", "index_j", "intercept"};
    PROTECT(Beta      = allocVector(REALSXP, tt_act_a));
    PROTECT(Indexi    = allocVector(INTSXP,  tt_act_a));
    PROTECT(Indexj    = allocVector(INTSXP,  tt_act_a));
    PROTECT(Lambda    = allocVector(REALSXP, iter));
    PROTECT(BIC       = allocVector(REALSXP, iter));
    PROTECT(Loss      = allocVector(REALSXP, iter));
    PROTECT(Direction = allocVector(INTSXP,  iter));
    PROTECT(Df        = allocVector(INTSXP,  iter));
    PROTECT(Loops     = allocVector(INTSXP,  1));
    PROTECT(Intercept = allocVector(REALSXP, 1));
    PROTECT(Result    = allocVector(VECSXP,  10));
    PROTECT(R_names   = allocVector(STRSXP,  10));

    for(i = 0; i < 10; ++i) SET_STRING_ELT(R_names, i,  mkChar(names[i]));
    INTEGER(Loops)[0] = iter;
    REAL(Intercept)[0] = intercept;
    for (i = 0; i < tt_act_a; ++i) 
    {
        temp = sparse_i[i];
        INTEGER(Indexi)[i]  = temp;
        INTEGER(Indexj)[i]  = sparse_j[i];
        REAL(Beta)[i]       = beta[i]/sdx[temp];
    }
    for (i = 0; i < iter; ++i) 
    {
        REAL(BIC)[i]          = bic[i];
        REAL(Loss)[i]         = loss[i];
        REAL(Lambda)[i]       = lambda[i];
        INTEGER(Direction)[i] = direction[i];
        INTEGER(Df)[i]        = df[i];
    }
    
    SET_VECTOR_ELT(Result, 0, Beta);
    SET_VECTOR_ELT(Result, 1, Lambda);
    SET_VECTOR_ELT(Result, 2, Direction);
    SET_VECTOR_ELT(Result, 3, Loops); 
    SET_VECTOR_ELT(Result, 4, BIC);
    SET_VECTOR_ELT(Result, 5, Loss);  
    SET_VECTOR_ELT(Result, 6, Df);  
    SET_VECTOR_ELT(Result, 7, Indexi);  
    SET_VECTOR_ELT(Result, 8, Indexj); 
    SET_VECTOR_ELT(Result, 9, Intercept);    
    setAttrib(Result, R_NamesSymbol, R_names); 


    free(sdx);
    free(loss);
    free(lambda);
    free(direction);
    free(active);
    free(bic);
    free(df);
    free(derivative);
    free(d);
    free(derivative1);
    free(sparse_i);
    free(sparse_j);
    free(beta);
    free(residual);
    free(index_max_s);
    free(r);

	UNPROTECT(12);
	return Result;
}