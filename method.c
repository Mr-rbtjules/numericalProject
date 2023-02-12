#include <stdlib.h>
#include <stdio.h>


//pener a compute B que 1 fois dans le multigrid
int forwardGS(int iter, int n, int *ia, int *ja, double *a, double *b, double *u){

	int *il;
	int *jl;
	double *l;

	//preconditionner B

	stationaryIter(iter, n, ia, ja, a, b, u, il,jl,l);

	

	return 0;
}

int computeL(int n, int nnz, int *ia, int *ja, double *a, int n, int *il, int *jl, double *l){

	int nnzl = ((nnz-n)/2) + n;
	il  = malloc((n + 1) * sizeof(int));
    jl  = malloc(nnzl * sizeof(int));
    l   = malloc(nnzl * sizeof(double));

	int i = 0;
	int jai = 0;
	il[0] = 0;
	while (i < n){
		
		int ite = ia[i + 1] - ia[i];
		int j = 0;
		while (j < ite){
			if (ja[jai] <= n){
				l[ite + j] = a[ite + j];
				jl[ite + j] = ja[ite + j];
			}
			else {
				il[n+1] = jai + 1; 
				j = ite; // break car arrivee juste avant milieu de ligne
			}
		}
		jai += ite;
		
		i++;
	}

	return 0;
}

int computeU(int n, int nnz, int *ia, int *ja,
			 double *a, int *iu, int *ju, double *u){

	//sizel = sizea/2 car symetrique
	int nnzu = ((nnzu-n)/2) + n; //retire elem diag divise par 2 car sym puis rajoute elem diag
	iu  = malloc((n + 1) * sizeof(int));
    ju  = malloc(nnzu * sizeof(int));
    u   = malloc(nnzu * sizeof(double));
	
	int i = 0;
	int jai = 0;
	iu[0] = 0;
	while (i < n){
		
		int ite = ia[i + 1] - ia[i];
		int j = 0;
		while (j < ite){
			if (ja[jai] <= n){
				u[ite + j] = a[ite + j];
				ju[ite + j] = ja[ite + j];
			}
			else {
				iu[n+1] = jai + 1; 
				j = ite; // break car arrivee juste avant milieu de ligne
			}
		}
		jai += ite;
		i++;
	}

	return 0;
}

int backwardGS(){

	return 0;
}

int stationaryIter(int iter, int n, int *ia, int *ja, double *a,
					 double *b, double *u, int *iB, int *jB, double *B){
	double *r;

	//initialisation
	*u = calloc(n * sizeof(double));
	computeRes(n, ia, ja, a, u, b, r);

	for (int m = 0; m < iter; m++){
		
	}
}

int computeRes(int n, int *ia, int *ja, double *a, double *u, double *b, double *r){

	double *au;
	multMatVectCsr(n, ia, ja, a, u, au);
	soustVect(n, b, au, r);

	return 0;
}

int multMatVectCsr(int n , int *ia, int *ja, double *a, double *u, double *b){
	if (b != NULL){
		printf("Error b de syst non null\n");
		return 1;
	}
	b = calloc(n ,sizeof(double));
	int i = 0;
	int jai = 0;
	while (i < n){
		int ite = ia[i + 1] - ia[i];
		int j = 0;
		while (j < ite){
			b[i] += a[jai + j] * u[ja[jai + j]]; 		
			j += 1;
		}
		jai += ite;
		i += 1;
	}

	return 0;
}

int addVect(int n , double *v1, double *v2, double *vout){
	vout = malloc(n * sizeof(double));
	for (int i = 0; i < n; i++){
		vout[i] = v1[i] + v2[i];
	}
	return 0;
}

int soustVect(int n , double *v1, double *v2, double *vout){
	vout = malloc(n * sizeof(double));
	for (int i = 0; i < n; i++){
		vout[i] = v1[i] - v2[i];
	}
	return 0;
}

int gaussResD(int n , int *il, int *jl, double *l, double *u, double *b){) {

	int i = 0;
	if (u == NULL){
		u = malloc(sizeof(double) * n);
	}
	while (i < n){

		int start_ind_jl = il[i];
		int end_ind_jl = il[i+1] - 1;
		u[i] = b[i]; //copie b sur u
		while (start_ind_jl < end_ind_jl){

			u[i] -= (l[start_ind_jl] 
							* u[jl[start_ind_jl]]);
			start_ind_jl += 1;
		}
		//car start == end et donc fin de ligne => elem diag
		u[i] /= l[end_ind_jl];
		i += 1;
	}
	return 0;
}


int invL(int n , int *il, int *jl, double *l, double *u, double *b){) {

		int i = 0;
		if (u == NULL){
			u = malloc(sizeof(double) * n);
		}
		while (i < n){

			int start_ind_jl = il[i];
			int end_ind_jl = il[i+1] - 1;
			u[i] = b[i]; //copie b sur u
			while (start_ind_jl < end_ind_jl){

				u[i] -= (l[start_ind_jl] 
								* u[jl[start_ind_jl]]);
				start_ind_jl += 1;
			}
			//car start == end et donc fin de ligne => elem diag
			u[i] /= l[end_ind_jl];
			i += 1;
		}
}
	return 0;
}

