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

int stationaryIter(int iter, int n, int *ia, int *ja,
					 double *a,double *b, double *u, int forward){
	
	if (iter != 0){
		iter -= 1;
		stationaryIter(iter, n, ia, ja, a,
					 b, u, iB, jB, B);

	}
	else{
		//init qqchose
		//initialisation u0
		u = calloc(n * sizeof(double));
		r = malloc(n *sizeof(double));
		//qqchose
	}
	//rm
	computeRes(n, ia, ja, a, u, b, r);
	if (forward){
		gaussResL()
	}
	//dm == solve B*d = r

	//rm+1

	//um+1
	

	for (int m = 0; m < iter; m++){
		
	}
}

int plot(double *r)

int computeRes(int n, int *ia, int *ja, double *a, double *u, double *b, double *r){
	//r = b -Au
	if (b != NULL){
		printf("Error b de syst non null\n");
		return 1;
	}
	
	int i = 0;
	int jai = 0;
	while (i < n){
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

	
	
	//trop gourmant memoire multMatVectCsr(n, ia, ja, a, u, au);
	//soustVect(n, b, au, r);

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




//pas besoin de crer L, juste prendre ds A trou de balle
int gaussResL(int n , int *il, int *jl, double *l, double *x, double *b){) {
	//resoud Lx = b trouve x
	int i = 0;
	if (x == NULL){
		x = malloc(sizeof(double) * n);
	}
	while (i < n){

		int start_ind_jl = il[i];
		int end_ind_jl = il[i+1] - 1;
		x[i] = b[i]; //copie b sur u
		while (start_ind_jl < end_ind_jl && jl[start_ind_jl] < i){

			x[i] -= (l[start_ind_jl] 
							* x[jl[start_ind_jl]]);
			start_ind_jl += 1;
		}
		//car start == end et donc fin de ligne => elem diag
		x[i] /= l[end_ind_jl];
		i += 1;
	}
	return 0;
}

int gaussResU(int n , int *iu, int *ju, double *u, double *x, double *b){) {
	//resoud Lx = b trouve x
	int i = 0;
	if (x == NULL){
		x = malloc(sizeof(double) * n);
	}
	while (i < n){

		int start_ind_ju = iu[i];
		int end_ind_ju = iu[i+1] - 1;
		x[i] = b[i]; //copie b sur u
		while (start_ind_ju < end_ind_ju){

			x[i] -= (l[start_ind_ju] 
							* x[ju[start_ind_ju]]);
			start_ind_ju += 1;
		}
		//car start == end et donc fin de ligne => elem diag
		x[i] /= u[end_ind_ju];
		i += 1;
	}
	return 0;
}


