#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "proto.h"


#define RELAX 1 //
#define MODE 1
#define EXPLICIT 1

/*
        level 0
\      /1
 \    /2
  \  /3
   \/4 = level max = nb de fois qu'on restr
*/
int mg_method(int iter, int levelMax, int m){
	//initit memory and pointers
	//compute all the coarse matrix and nl
	int **ial; //liste donc chaque elem pointe vers 1 matrice d'un certain level (ici matrice == liste)
	int **jal;
	double **al;
	double **rl;
	double **ul;
	double **dl;
	
	double **bl; //attention aux autres b prc pt on les utilise pour autre
	
	int *nl; 

	int mu1 = 3;
	int mu2 = 3;

	allocGrids(m, levelMax, &nl, &ial, &jal, &al, &bl, &dl, &rl, &ul);

	//precomputation of all Ac and bc

	for (int l = 0; l <= levelMax; l++){
		
		probMg(m, l, &(nl[l]), ial[l], jal[l], al[l], bl[l]);

	}
	

	//start iterations of the multigrid cycle (tg_rec)
	int startLevelTg = 0; 

	initialization(startLevelTg, nl, ial, jal, al, bl, dl, rl, ul);

	for (int i = 0; i < iter; i++){
		//initialisation ici ?
		
		firstStep(startLevelTg, m, mu1, nl,
				  ial, jal, al, bl, ul, rl, dl);
		
		tg_rec( startLevelTg+1, levelMax, m, mu1, mu2, 
				nl, ial, jal, al, rl, ul, bl, dl);
		lastStep(startLevelTg, m, mu2, nl,
				  ial, jal, al, bl, ul, rl, dl);
		//print(res et u) + save
		
		
		
	}
	plot_res(rl[startLevelTg], m, startLevelTg);
	free(nl);
	for (int j = 0; j <= levelMax; j++){
		free(ial[j]);
		free(jal[j]);
		free(al[j]);
		free(rl[j]);
		free(dl[j]);
		free(ul[j]);
	}
	
	return 0;
}
int initialization(int level, int *nl, int **ial, int **jal, double **al,
						 double **bl, double **dl, double **rl, double **ul){

	//arbitrary u0
	for (int i=0; i< nl[level]; i++){
		ul[level][i] = 0; 
	}
	printf("\n initial res : %lf\n", computeResNorm(nl[level],
												ial[level],
												jal[level],
												al[level],
												bl[level],
												ul[level],
												rl[level]));
	return 0;
}

int firstStep(int startLevelTg, int m, int mu1, int *nl,
				int **ial, int **jal, double **al, double **bl,
				  double **ul, double **rl, double **dl){
	printf("\n First Step\n");
	forwardGS(mu1, nl[startLevelTg] , ial[startLevelTg], 
			  jal[startLevelTg], al[startLevelTg], bl[startLevelTg],
			  ul[startLevelTg], rl[startLevelTg], dl[startLevelTg]);//ici on utilise b pour stocker le residu de Ac=r et stock c dans u
	
	printf("\n after fgsres : %lf\n", computeResNorm(nl[startLevelTg],
												ial[startLevelTg],
												jal[startLevelTg],
												al[startLevelTg],
												bl[startLevelTg],
												ul[startLevelTg],
												rl[startLevelTg]));
		//restrict r : r-Ac stocker dans b !
	computeRes(nl[startLevelTg], ial[startLevelTg],
			   jal[startLevelTg], al[startLevelTg], 
			   ul[startLevelTg], bl[startLevelTg], rl[startLevelTg]);
		//restrict residu stocker dans b
	restrictR(startLevelTg, rl[startLevelTg],
			  rl[startLevelTg+1], m, &(nl[startLevelTg+1]));
	return 0;
}




int tg_rec(int level, int levelMax, int m, int mu1,
			int mu2, int *nl, int **ial, int **jal,
		   double **al, double **bl, double **ul, double **rl, double **dl){
		
		//des le level 1 on stock pas u2 mais c2 puis 'deviennent des u lorsqu'on ajoute la correction en remontant
		//pblm tt en haut smoothing Au=b mais en dessous sm Ac = r pareil computeres
	if (level < levelMax){

		//mu1 smoothing iterations
		forwardGS(mu1, nl[level] , ial[level], jal[level],
					al[level], bl[level], ul[level], rl[level], dl[level]);//ici on utilise b pour stocker le residu de Ac=r et stock c dans u
		
		//restrict r : r-Ac stocker dans b !
		computeRes(nl[level], ial[level], jal[level],
					al[level], ul[level], bl[level], rl[level]);
		//restrict residu stocker dans b
		restrictR(level, rl[level],rl[level+1], m, &(nl[level+1]));
		
		//au dessus s'effecctue du haut du v vers le bas
		tg_rec(level+1, levelMax, m, mu1, mu2,
		       nl, ial, jal, al, bl, ul, rl, dl);//recursivité , ici important on repart avec A c = r
		// on voit que ul est composé de u en 0 puis que des c puis en
		// remontant on va applique les corrections aux corrections, r est composé du 
		//'vrai residu en 0 pus des residu de residu etc
		//code en dessous s'effectue du bas vers le haut du v
		//prolongR(level,ul[level-1], ul[level], m, nl );// direct ajoute la correction quasi r a faire juste enlver le malloc
		
		addProlCorrection(level+1, ul[level], 
		      			  ul[level+1], m, &(nl[level]));//
		backwardGS(mu2, nl[level] , ial[level], jal[level], 
				   al[level], bl[level], ul[level], rl[level], dl[level]);
		
	}
	else {
		//solve coarse pblm
		printf("Solve at coarst level\n");
		solveAtCoarseLevel(MODE, nl[level], ial[level], jal[level], al[level],
						   		bl[level], ul[level], rl[level], dl[level]);
	}
	return 0;
}
int lastStep(int startLevelTg, int m, int mu2, 
			 int *nl, int **ial, int **jal, double **al, 
			 double **bl, double **ul, double **rl, double **dl){

	addProlCorrection(startLevelTg+1, ul[startLevelTg], 
					  ul[startLevelTg+1], m, &(nl[startLevelTg]));
	printf("\n after solve and prol : %lf\n", computeResNorm(nl[startLevelTg],
															 ial[startLevelTg],
															 jal[startLevelTg],
															 al[startLevelTg],
															 bl[startLevelTg],
															 ul[startLevelTg],
															 rl[startLevelTg]));
	backwardGS(mu2, nl[startLevelTg] , ial[startLevelTg], 
			   jal[startLevelTg], al[startLevelTg], bl[startLevelTg], 
			   ul[startLevelTg], rl[startLevelTg], dl[startLevelTg]);
	
	printf("\n res : %lf\n", computeResNorm(nl[startLevelTg],
											ial[startLevelTg],
											jal[startLevelTg],
											al[startLevelTg],
											bl[startLevelTg],
											ul[startLevelTg],
											rl[startLevelTg]));
	return 0;	
}




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
	//trop gourmant  de faire memoire multMatVectCsr(n, ia, ja, a, u, au);
	//soustVect(n, b, au, r);

	return 0;
}
double computeNorm(int n, double *v){

	double vn = 0;

	for (int i = 0; i < n; i++){
		vn += v[i] * v[i];
	}
	return sqrt(vn);
}

double computeResNorm(int n, int *ia, int *ja, double *a, double *b, double *u,  double *r){

	double rn = 0;

	computeRes(n,ia,ja,a,u,b,r);
	for (int i = 0; i < n; i++){
		rn += r[i] * r[i];
	}
	rn = sqrt(rn);
	return rn;
}

int addVect(int n , double *v1, double *v2){
	for (int i = 0; i < n; i++){
		v1[i] += v2[i];
	}
	return 0;
}



//pas besoin de crer L, juste prendre ds A trou de balle
int gaussResL(int n , int *il, int *jl, double *l, double *x, double *b){
	//resoud Lx = b trouve x
	int i = 0;
	while (i < n){

		int start_ind_jl = il[i];
		
		int end_ind_jl = il[i+1];
		
		x[i] = b[i]; //copie b sur u
		while (start_ind_jl < end_ind_jl && jl[start_ind_jl] < i){

			x[i] -= (l[start_ind_jl] 
							* x[jl[start_ind_jl]]);
			start_ind_jl += 1;
		}
		//car start == end et donc fin de ligne => elem diag
		x[i] /= l[start_ind_jl];
		i += 1;
	}
	return 0;
}

int gaussResU(int n , int *iu, int *ju, double *u, double *x, double *b){ 
	//resoud Ux = b trouve x
	int i = n-1;
	
	while (i >= 0){

		int start_ind_ju = iu[i+1]-1;
		int end_ind_ju = iu[i];
		x[i] = b[i]; //copie b sur u
		while (start_ind_ju > end_ind_ju && ju[start_ind_ju] > i){

			x[i] -= (u[start_ind_ju] 
							* x[ju[start_ind_ju]]);
			start_ind_ju -= 1;
		}
		//car start == end et donc debut de ligne => elem diag
		x[i] /= u[start_ind_ju];
		i -= 1;
	}
	return 0;
}

int gaussResD(int n , int *ia, int *ja, double *a, double *x, double *b){
	//resoud Lx = b trouve x
	int i = 0;
	while (i < n){

		int start_ind_j = ia[i];
		
		int end_ind_j = ia[i+1];
		
		x[i] = b[i]; //copie b sur u
		while (start_ind_j < end_ind_j && ja[start_ind_j] < i){
			start_ind_j += 1;
		}
		if (ja[start_ind_j] == i){
			x[i] /= a[start_ind_j];
		}
		else {
			printf("\n no diag elem\n");
			exit(0);
		}		
		i += 1;
	}
	return 0;
}



int backwardGS(int iter, int n, int *ia, int *ja, double *a,
				double *b, double *u, double *r, double *d){

	//printf("Backward GS\n");
	stationaryIter(iter, n, ia, ja, a, b, u, r, d, -1);

	return 0;
}
//penser a compute B que 1 fois dans le multigrid
int forwardGS(int iter, int n, int *ia, int *ja, double *a,
			    double *b, double *u, double *r, double *d){
	//printf("Forward GS\n");
	//resoud Au = b avec iter iteration
	stationaryIter(iter, n, ia, ja, a, b, u, r, d, 1);
	return 0;
}


int jacobiIter(int iter, double tol, int n, int *ia, int *ja, double *a,
			    double *b, double *u, double *r, double *d){
								
	stationaryIter(iter, n, ia, ja, a, b, u, r, d, 0);
	return 0;
}


//ial liste de liste d'elem deja compute 


int stationaryIter(int iter, int n, int *ia, int *ja,
                   double *a, double *b, double *u, double *r, double *d, int forward){
//forward = 0 => backwardGs
	
	if (iter != 1){
		stationaryIter(iter-1, n, ia, ja, a,
					 b, u, r, d, forward);
	}
	else{
		//initialiser verifier ? on initialise pas car considere que deja fait
	}
	
	//rm
	computeRes(n, ia, ja, a, u, b, r);
	
	// dm == solve B*d = r (B construit a partir de A direct 
	// ds la methode)
	if (forward == 1){
		
		
		gaussResL(n , ia, ja, a, d, r);
		
	}
	else if (forward == -1){
		gaussResU(n , ia, ja, a, d, r);
	}
	else if (forward == 0){
		gaussResD(n , ia, ja, a, d, r);	
	}
	relax(n, d);
	//um+1 (ajout de la correction)
	addVect(n , u, d);

	return 0;
}





int solveAtCoarseLevel(int mode, int n, int *ia, int *ja, double *a, double *b, double *u, double *r, double *d){

	if (mode == 0){
		solve_umfpack(n, ia, ja, a, b, u);
	}
	else if (mode == 1) {
		symGS(1, 0, n, ia,ja,a,b,u,r,d);
	}
	else {

	}

}

int symGS(int iter, double tol, int n, int *ia, int *ja, double *a,
			    double *b, double *u, double *r, double *d){
	//printf("\n Sym gs initial res=  %lf\n", computeResNorm(n,ia,ja,a,b,u,r));
								
	if (iter){
		int i = 0;
		while (i < iter){
			
			forwardGS( 1, n, ia, ja, a, b, u, r, d);
      		backwardGS( 1, n, ia, ja, a, b, u, r, d);
			i++;
			if (tol){//tol > computeNorm(n,r)){
				i == iter;
			}
		}
	}
	else{
		while (computeNorm(n,r) > tol){
			forwardGS( 1, n, ia, ja, a, b, u, r, d);
      		backwardGS( 1, n, ia, ja, a, b, u, r, d);
		}
	}
	printf("\n Sym gs end res=  %lf\n", computeResNorm(n,ia,ja,a,b,u,r));

	return 0;
}


int allocGrids(int m, int levelMax, int **nl, int ***ial,
               int ***jal, double ***al, double ***bl,
			   double ***dl, double ***rl, double ***ul){ /* *** var liste de liste et on acces la memoire donc 1 en plus */

	*ial = malloc((levelMax+1) * sizeof(int*));
	*jal = malloc((levelMax+1) * sizeof(int*));
	*al = malloc((levelMax+1) * sizeof(double*));
	*rl = malloc((levelMax+1) * sizeof(double*));
	*dl = malloc((levelMax+1) * sizeof(double*));
	*ul = malloc((levelMax+1) * sizeof(double*));
	*bl = malloc((levelMax+1) * sizeof(double*));

	*nl = malloc((levelMax+1) * sizeof(int));

	for (int l = 0; l <= levelMax; l++){	
		
		allocLevel(m, l, *nl, ial, jal, al, bl, dl, rl, ul);
	}


	return 0;
}
int allocLevel(int m, int level, int *nl, int ***ial,
                     int ***jal, double ***al, double ***bl,
                     double ***dl, double ***rl, double ***ul){

    double hl, invh2l;
    int x0l,x1l,y0l,y1l, nxl, nnzl;
    computeParamLevel(m, level, &hl,&invh2l,&y0l,&y1l,&x0l,&x1l,&nxl, &(nl[level]), &nnzl);
    
    (*ial)[level] = malloc((nl[level] + 1) * sizeof(int));
    (*jal)[level] = malloc(nnzl * sizeof(int));
    (*al)[level] = malloc(nnzl * sizeof(double));
    (*bl)[level] = malloc(nl[level] * sizeof(double));
    (*dl)[level] = malloc(nl[level] * sizeof(double));
    (*rl)[level] = malloc(nl[level] * sizeof(double));
    (*ul)[level] = malloc(nl[level] * sizeof(double));
    
    if (*bl == NULL || *ial == NULL || *jal == NULL || 
        *al == NULL ||*dl == NULL || *rl == NULL || *ul == NULL){
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n\n");
        return 1;
    }

    printf("\n Alloc level %d : hl = %lf nl ",level, hl);
    printf(" = %d nnzl = %d nxl = %d \n", nl[level], nnzl, nxl);
    printf("x0l = %d x1l = %d y0l = %d y1l = %d \n", x0l, x1l, y0l, y1l);


    return 0;
}

int relax(int n, double *d){

	for (int i=0; i < n; i++){
		d[i]*=RELAX;
	}
}
void printU(int n, double *u){
	printf("\n[");
	for (int i=0;i<n-1;i++){
		printf("%.16g; ", u[i]);
	}
	printf("%.16g]\n ", u[n-1]);
}

void printA(int n, int *ia, int *ja, double *a){

	printf("\n[");
	for (int i=0; i < n; i++){
		int s = ia[i];
		int end = ia[i+1]-1;
		int j = 0;
		while (j <n){
			if (j != ja[s]){
				printf(" 0");
			}
			else if (j == ja[s]){
				printf(" %.5g", a[s]);
				if (s < end){
					s++;
				}
			}
			j++;
		} 
		printf(";");
	}
	printf("]\n");

}