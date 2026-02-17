#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "proto.h"


 //TOOLS//



int isSymmetric(int* ia, int* ja, double* a, int *n) {
    for (int i = 0; i < *n; i++) {
        for (int j = ia[i]; j < ia[i + 1]; j++) {
            int col = ja[j];
            int row = i;
            double value = a[j];

            // Find the corresponding element in the transpose
            int found = 0;
            for (int k = ia[col]; k < ia[col + 1]; k++) {
                if (ja[k] == row) {
                    if (fabs(a[k] - value) > 1e-10) {
                        return 0; // Not symmetric
                    }
                    found = 1;
                    break;
                }
            }
            if (!found) {
                return 0; // Corresponding element not found, not symmetric
            }
        }
    }
    return 1; // Matrix is symmetric
}

int addVect(int *n , double *v1, double *v2){
	for (int i = 0; i < *n; i++){
		v1[i] += v2[i];
	}
	return 0;
}

double computeBound(double x, double y){
    //compute the value of a point on the bound of the domain
	//with a function BOUND defined in the macros
    return BOUND(x,y);
}

double computeResNorm(int *n, int *ia, int *ja, double *a, 
						double *b, double *u,  double *r){

	double rn = 0;

	computeRes(n,ia,ja,a,u,b,r);
	for (int i = 0; i < *n; i++){
		rn += r[i] * r[i];
	}
	rn = sqrt(rn);
	return rn;//*100000000000000;
}

int computeRes(int *n, int *ia, int *ja, double *a, 
			   		double *u, double *b, double *r){
	
	//r = b -Au
	int i = 0;
	int jai = 0;
	while (i < *n){
		r[i] = b[i];
		int ite = ia[i + 1] - ia[i];
		int j = 0;
		while (j < ite){
			r[i] -= a[jai + j] * u[ja[jai + j]]; 		
			j += 1;
		}
		jai += ite;
		i += 1;
	}	
	//trop gourmant  de faire memoire multMatVectCsr(n, ia, ja, a, u, au);
	//soustVect(n, b, au, r);

	return 0;
}

int copy( int *n, double *toCopy, double *copy){
    for (int i = 0; i < *n; i++)
    {
        copy[i] = toCopy[i];
    }
    
    return 0;
}



int scalProd(int *n, double *v1, double *v2, double *res){

    *res = 0;
    int i = 0;
    while (i < *n){
        *res += v1[i] * v2[i];
        i += 1;
    }

    return 0;
}

int subVectProd(int *n, double *alpha, double *v2, double *v1){

    int i = 0;
    while (i < *n){
        v2[i] -= (*alpha) * v1[i];
        i += 1;
    }
    return 0;
}

int addVectProd(int *n, double *alpha, double *v2, double *v1){

    int i = 0;
    while (i < *n){
        v2[i] += (*alpha) * v1[i];
        i += 1;
    }
    return 0;
}

int dSum(int *n, double *beta, double *d, double *z){

    int i = 0;
    while (i < *n){
        d[i] = z[i] + (*beta) * d[i];
        i += 1;
    }
    return 0;
}

int multCsrVector(int *n, int *ia, int *ja, 
                 double *a, double *Ax, double *x){
	//do v2 = A*v1 
	
	int i = 0;
	int jai = 0;
	while (i < *n){
		int ite = ia[i + 1] - ia[i];
		int j = 0;
        Ax[i] = 0;
		while (j < ite){
			Ax[i] += a[jai + j] * x[ja[jai + j]]; 		
			j += 1;
		}
		jai += ite;
		i += 1;
	}

	return 0;
}

int preInitialization(double *u0){
    //at the begining of mg method u0 must be set at 0
    //at lower level (ci), initialisation at 0 will be done 
    // at the begining of EACH iteration
    //for the record we valid that only u need to be initialized
    //other will be set (r[i] = b[i] in computeres
    //and in gaussres d is x and x[i] = b[i]) 
	//arbitrary u0
	for (int i=0; i < globVal.vectStart[1]; i++){
		u0[i] = 0; 
	}
	
	return 0;
}

//opti while ? rec ?
int reInitialization(double *ul){
    //ci are all set to 0 (ci are ul except top level)
    for (int i = globVal.vectStart[1]; i < globVal.vectStart[globVal.levelMax+1]; i++){
		ul[i] = 0; 
	}
    return 0;
}