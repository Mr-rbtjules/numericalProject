#include <stdlib.h>
#include <stdio.h>


//pener a compute B que 1 fois dans le multigrid
int forwardGS(int iter, int n, int *ia, int *ja, double *a, double *b, double **u, double **r, double **d){

	
	stationaryIter(iter, n, ia, ja, a, b, &u, &r, &d, 1);

	return 0;
}

int backwardGS(int iter, int n, int *ia, int *ja, double *a, double *b, double **u, double **r, double **d){

	
	stationaryIter(iter, n, ia, ja, a, b, &u, &r, &d, 0);

	return 0;
}



int stationaryIter( int iter, int n, int *ia, int *ja,
					double *a,double *b, double **u, 
					double **r, double **d, int forward){
	
	if (iter != 1){
		
		stationaryIter(iter-1, n, ia, ja, a,
					 b, &u, &r, &d, forward);

	}
	else{
		//pour etre sur d'utiliser l'approx
		if (*u == NULL && *r == NULL && *d == NULL){
			//initialisation u0
			*u = calloc(n * sizeof(double));
			*r = malloc(n *sizeof(double));
			*d = malloc(n *sizeof(double));
		}
		else {
			//cas ou existe deja une appprox pour u0
			//pas besoin d'allouer de nouveau de la memoire
			//u r d existe deja
		}
		
	}
	//rm
	computeRes(n, ia, ja, a, *u, *b, *r);
	
	// dm == solve B*d = r (B construit a partir de A direct 
	// ds la methode)
	if (forward){
		gaussResL(n , ia, ja, a, *d, *r);
	}
	else{
		gaussResR(n , ia, ja, a, *d, *r);
	}
	
	//um+1 (ajout de la correction)
	addVect(n , *u, *d)
	
	return 0;
}

int plot(double *r);

int computeRes(int n, int *ia, int *ja, double *a, double *u, double *b, double *r){
	
	//r = b -Au
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







//pas besoin de crer L, juste prendre ds A trou de balle
int gaussResL(int n , int *il, int *jl, double *l, double *x, double *b){) {
	//resoud Lx = b trouve x
	int i = 0;
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
	int i = n-1;
	
	while (i >= 0){

		int start_ind_ju = iu[i+1]-1;
		int end_ind_ju = iu[i];
		x[i] = b[i]; //copie b sur u
		while (start_ind_ju > end_ind_ju && ju[start_ind_ju] > i){

			x[i] -= (l[start_ind_ju] 
							* x[ju[start_ind_ju]]);
			start_ind_ju -= 1;
		}
		//car start == end et donc debut de ligne => elem diag
		x[i] /= u[end_ind_ju];
		i -= 1;
	}
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

int addVect(int n , double *v1, double *v2){
	for (int i = 0; i < n; i++){
		v1[i] += v2[i];
	}
	return 0;
}


