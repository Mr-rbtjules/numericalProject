#include <stdlib.h>
#include <stdio.h>
#include <math.h>


//ial liste de liste d'elem deja compute 
int tg_rec(int level, int m, int *nl, int mu1, int mu2, int **ial, 
					int **jal, double **al, double *b, double **ul, double **rl, double **dl){

	
	if (level != 0){
		forwardGS(mu1, nl[level] , ial[level], jal[level], al[level], b[level], u[level], r[level], d[level]);
		
		//restrict r
	//!!! pblme check b calculer que 1 er level et apres que des residus ? verif ça aussi dans l'algo
		computeRes(n, ial[level], jal[level], al[level], ul[level], bl[level], &rl[level]);
		restrictR(m, level, rl[level],rl[level-1], m, nl[level]); //modify malloc init
		//au dessus s'effecctue du haut du v vers le bas
		tg_rec(level-1, nl, mu1, mu2, ial, jal, al, rl[level], ul, rl);//recursivité , ici important on repart avec A c = r
		// on voit que ul est composé de u en 0 puis que des c puis en
		// remontant on va applique les corrections aux corrections, r est composé du 
		//'vrai residu en 0 pus des residu de residu etc
		//code en dessous s'effectue du bas vers le haut du v
		prolongR(level,m, ul);//modifer pour que multi + direct ajoute la correction quasi r a faire juste enlver le malloc
		backwardGS(mu1, nl[level] , ial[level], jal[level], al[level], b[level], u[level], r[level], d[level]);
		
		
	}
	else {
		//solve coarse pblm
		solve_umfpack(n, ial[level], jal[level], al[level], rl[level], ul[level]);
	}
	return 0;
}

int computeSize(int m, int level, int m, int **size_ial,int **size_jal,
				int **size_al, int **size_rl, int **size_dl, int **size_ul){

	int nx = m-2;

	*size_ial = malloc(level * sizeof(int));
	*size_jal = malloc(level * sizeof(int));
	*size_al = malloc(level * sizeof(int));
	*size_ul = malloc(level * sizeof(int));
	*size_dl = malloc(level * sizeof(int));

	for (int i = 1; level){
		
		int x0c = x0 / exp(2, i); // va arrondir au point grille coarse a droite 
		int x1c = ((x1+1)/exp(2, i)) - 1; // permet si x1 pair on retire 1 
		int y0c = (y0/exp(2, i)); // arrondu coarse au dessus (permet de pas ajouter des points dans le trou)
		int y1c = ((y1+1)/exp(2, i)) - 1;
		int pc = y1c - y0c + 1;
		int qc = x1c - x0c + 1;

		int nxc = nx/exp(2,i); // nb de points coars sur un ligne pas bord
		int nc = nxc * nxc - (pc * qc);
		int nnzc = 5 * nxc * nxc - 4 * nxc ; 
		//nb de points concernés dans le trou:(compliqué a comprendre sans shema) (marche que si aumoins 3 points sur la largeur)
		int trousc = (5 * (pc-2) * (qc-2) + 4 * 2 * (pc-2) + 4 * 2 * (qc-2) 
					+ 3 * 4 * 1 + 1 * 2 * pc + 1 * 2 * qc) + 2 * pc + 2 * qc;
		nnzc -= trousc;

		size_ial[0][i] = nc + 1;
		size_jal[0][i] = nnzc;
		size_al[0][i] = nnzc;
		size_dl[0][i] = nc;
		size_ral[0][i] = nc;
		size_ual[0][i] = nc;
	
	}
	return 0;
}

int mg_method(int iter, int level, int m){
	//initit memory and pointers
	//compute all the coarse matrix and nl
	int **ial;
	int **jal;
	double **al;
	double **rl;
	double **ul;
	double **dl;
	
	double *b;
	int *size_ial;
	int *size_jal;
	int *size_al;
	int *size_rl;
	int *size_dl;
	int *size_ul;

	computeSize(&size_ial, &size_jal, &size_al, &size_rl, &size_dl, &size_ul);

	b = malloc(size_ul[0]* sizeof(double));
	for (int i = 0; i < level; i++){
		ial[i] = malloc(size_ial[i] * sizeof(int));
		jal[i] = malloc(size_jal[i] * sizeof(int));
		al[i] = malloc(size_al[i] * sizeof(double));
		rl[i] = malloc(size_rl[i] * sizeof(double));
		dl[i] = malloc(size_dl[i] * sizeof(double));
		ul[i] = malloc(size_ul[i] * sizeof(double));
		//compute a and all ac
		probMg(m, i, ial[i], jal[i], al[i], b);
	}
	
	
	for (int i = 0; i < iter; i++){
		tg_rec(level, m, size_ul, mu1, mu2, ial, 
					jal, al, b, ul, rl, dl);
		//print(res et u) + save
	}
	//

	for (int j = 0; i < level; i++){
		free(ial[i]);
		free(jal[i]);
		free(al[i])
		free(rl[i]);
		free(dl[i]);
		free(ul[i]);
	}
	

}

//pener a compute B que 1 fois dans le multigrid
int forwardGS(int iter, int n, int *ia, int *ja, double *a, double *b, double **u, double **r, double **d){

	//resoud Au = b avec iter iteration
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
	computeRes(n, ia, ja, a, *u, b, *r);
	
	// dm == solve B*d = r (B construit a partir de A direct 
	// ds la methode)
	if (forward){
		gaussResL(n , ia, ja, a, *d, *r);
	}
	else{
		gaussResU(n , ia, ja, a, *d, *r);
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


