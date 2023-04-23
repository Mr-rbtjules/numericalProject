#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "proto.h"


#define COORD_X0 1.0 //col bord gauche
#define COORD_X1 2.5 // col bord droit
#define COORD_Y0 1.5 //ligne bord bas
#define COORD_Y1 2.0 //ligne bord haut


#define RELAX 1 //
#define MODE 1
#define EXPLICIT 1// add time launch when 0
#define ITERCORE 2

/*
        level 0
\      /1
 \    /2
  \  /3
   \/4 = level max = nb de fois qu'on restr
*/


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


//                                 METHOD                                     //

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



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

	int mu1 = 4;
	int mu2 = 4;

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
	//printf("\n Start smoothing \n");
	//symGS(5,0, nl[0],ial[0], jal[0], al[0], bl[0], ul[0], rl[0], dl[0]);
	//plot_res(rl[startLevelTg], m, startLevelTg);
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
	if (EXPLICIT){
		printf("\n initial res : %lf\n", computeResNorm(nl[level],
														ial[level],
														jal[level],
														al[level],
														bl[level],
														ul[level],
														rl[level]));
	}
	return 0;
}

int firstStep(int startLevelTg, int m, int mu1, int *nl,
				int **ial, int **jal, double **al, double **bl,
				  double **ul, double **rl, double **dl){

	if (EXPLICIT){ printf("\n First Step\n"); }

	forwardGS(mu1, nl[startLevelTg] , ial[startLevelTg], 
			  jal[startLevelTg], al[startLevelTg], bl[startLevelTg],
			  ul[startLevelTg], rl[startLevelTg], dl[startLevelTg]);//ici on utilise b pour stocker le residu de Ac=r et stock c dans u
	
	if (EXPLICIT){
		printf("\n after fgsres : %lf\n", computeResNorm(nl[startLevelTg],
														ial[startLevelTg],
														jal[startLevelTg],
														al[startLevelTg],
														bl[startLevelTg],
														ul[startLevelTg],
														rl[startLevelTg]));
	}
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
		if (EXPLICIT) {printf("Solve at coarst level\n");}
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
	if (EXPLICIT){
		printf("\n after solve and prol res: %lf\n", computeResNorm(nl[startLevelTg],
																ial[startLevelTg],
																jal[startLevelTg],
																al[startLevelTg],
																bl[startLevelTg],
																ul[startLevelTg],
																rl[startLevelTg]));
	}
	backwardGS(mu2, nl[startLevelTg] , ial[startLevelTg], 
			   jal[startLevelTg], al[startLevelTg], bl[startLevelTg], 
			   ul[startLevelTg], rl[startLevelTg], dl[startLevelTg]);
	
	if (EXPLICIT){
		printf("\n after backward res : %lf\n", computeResNorm(nl[startLevelTg],
															ial[startLevelTg],
															jal[startLevelTg],
															al[startLevelTg],
															bl[startLevelTg],
															ul[startLevelTg],
															rl[startLevelTg]));
	}
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


int stationaryIter(int iter, int n, int *ia, int *ja,
                   double *a, double *b, double *u, 
				   double *r, double *d, int forward){
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
		if (EXPLICIT){
			printf("using umfpack\n");	
		}
		solve_umfpack(n, ia, ja, a, b, u);
	}
	else if (mode == 1) {
		
		if (EXPLICIT){
			printf("using symGS\n");	
		}
		symGS(ITERCORE, 0, n, ia,ja,a,b,u,r,d);
	}
	else {
		if (EXPLICIT){
			printf("using jacobi\n");	
		}
		jacobiIter(1,0, n,ia,ja,a,b,u,r,d);
		
	}
	return 0;
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
	if (EXPLICIT){
		printf("\n Sym gs end res=  %lf\n", computeResNorm(n,ia,ja,a,b,u,r));
	}

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
    computeParamLevel(m, level, &hl,&invh2l,&y0l,&y1l,
					  &x0l,&x1l,&nxl, &(nl[level]), &nnzl);
    
    (*ial)[level] = malloc((nl[level] + 1) * sizeof(int));
    (*jal)[level] = malloc(nnzl * sizeof(int));
    (*al)[level] = malloc(nnzl * sizeof(double));
    (*bl)[level] = malloc(nl[level] * sizeof(double));
    (*dl)[level] = malloc(nl[level] * sizeof(double));
    (*rl)[level] = malloc(nl[level] * sizeof(double));
    (*ul)[level] = malloc(nl[level] * sizeof(double));
    
    if (*bl == NULL || *ial == NULL || *jal == NULL || 
        *al == NULL ||*dl == NULL || *rl == NULL || *ul == NULL){
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n");
        return 1;
    }
	if (EXPLICIT){
		printf("\n Alloc level %d : hl = %lf nl ",level, hl);
    	printf(" = %d nnzl = %d nxl = %d \n", nl[level], nnzl, nxl);
    	printf("x0l = %d x1l = %d y0l = %d y1l = %d \n", x0l, x1l, y0l, y1l);
	}
    return 0;
}

int relax(int n, double *d){

	for (int i=0; i < n; i++){
		d[i]*=RELAX;
	}
}



/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


//                                 PROB                                     //

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

int probMg(int m, int level, int *nl, int *ial, int *jal, double *al, double *bl){
    //ici level == celui dont on veut calculer A et b
    //mémoire deja allouée
    double hl, invh2l;
    int x0l,x1l,y0l,y1l, nxl, nnzl;

    computeParamLevel(m, level, &hl,&invh2l,&y0l,
                      &y1l,&x0l,&x1l,&nxl, nl, &nnzl);
    if (EXPLICIT){

    	printf("\n ProbMg level %d : hl = %lf ", level, hl);
    	printf("nl = %d nnzl = %d nxl = %d \n", *nl, nnzl, nxl);
    	printf("x0l = %d x1l = %d y0l = %d y1l = %d \n", x0l, x1l, y0l, y1l);
	}
	
    int nnzl_save = nnzl;
    nnzl = 0;
    
    //passage ligne suiv(plaque complete)
    int ind = 0;
    for (int iyl = 0; iyl < nxl; iyl++) { //iy ix indice sur grille hors bords  mais position (ix+1)*h
        for (int ixl = 0; ixl < nxl; ixl++) {
            //exclu interieur et bord du trou
            if(! in_hole(ixl,iyl,y0l,y1l,x0l,x1l)){
                //marquer le début de la ligne suivante dans le tableau 'ia'
                ial[ind] = nnzl;
                
                bl[ind] = 0;
               
                if (check_sud(ixl,iyl,y0l,y1l,x0l,x1l,nxl)){
                    al[nnzl] = -invh2l;
                    jal[nnzl] = indice(ixl,iyl-1,y0l,y1l,x0l,x1l, nxl);//ind - nx + skip_sud;
                    
                    //-nx car on regarde delui d'en bas(shema)
                    //+skip_sud car comme on a passe des points(trous)
                    // nx ramene trop loin en arrière
                    nnzl++;
                }
                else{
                    bl[ind] += computeBound((ixl+ 1)*hl, (iyl + 1 -1)*hl) * invh2l; 
                }

                //replissage de la ligne : voisin ouest 
                //si pas a droite d'un bord
               
                if (check_west(ixl,iyl,y0l,y1l,x0l,x1l,nxl)){
                    al[nnzl] = -invh2l;
                    jal[nnzl] = ind - 1;
                    nnzl++;
                }
                else{
                    bl[ind] += computeBound((ixl + 1 - 1)*hl, (iyl + 1)*hl) * invh2l;
                }
                

                // replissage de la ligne : élém. diagonal
                al[nnzl] = 4.0*invh2l;
                jal[nnzl] = ind;
                
                nnzl++;
                
                // replissage de la ligne : voisin est
                //si pas a gauche d'un bord
                
                if ( check_est(ixl,iyl,y0l,y1l,x0l,x1l,nxl) ){
                    al[nnzl] = -invh2l;
                    jal[nnzl] = ind + 1;
                    nnzl++;
                }
                else{
                    bl[ind] += computeBound((ixl + 1 +1)*hl, (iyl + 1)*hl) * invh2l;
                }

                // replissage de la ligne : voisin nord
                //si pas en dessous d'un bord
                
                if ( check_nord(ixl,iyl,y0l,y1l,x0l,x1l,nxl) ){
                        al[nnzl] = -invh2l;
                        jal[nnzl] = indice(ixl,iyl+1,y0l,y1l,x0l,x1l, nxl);
                        nnzl++;
                }
                else{
                    bl[ind] += computeBound((ixl + 1)*hl, (iyl + 1 +1)*hl) * invh2l;
                }
                // numéro de l'équation
                ind += 1;
            }
        }
    }

     if (*nl != ind){
        printf(" err nl %d ind %d\n", *nl, ind);
    }
    else if (nnzl != nnzl_save){
        printf(" err nnzl %d nnzl_save %d\n", nnzl, nnzl_save);
    }
    else {
        ial[ind] = nnzl;
    }

	return 0;
}	

void computeParamTop(int m, double *h, double *invh2, int *y0, 
                     int *y1, int *x0, int *x1, int *nx, int *n){

    *h = 3.0/(double)(m-1);
    *invh2 = 1.0/((*h)*(*h));

    computeHole(y0,y1,x0,x1, m);
    
    *nx = m-2;

    int p = (*y1) - (*y0) + 1;
    int q = (*x1) - (*x0) + 1;
    *n = ((*nx) * (*nx)) - (p * q);
}

void computeParamLevel(int m, int level, double *hl, double *invh2l, 
                        int *y0l, int *y1l, int *x0l, int *x1l,
                        int *nxl, int *nl, int *nnzl){
    //!!!! faux !
	/*int sizeRc;
	int mc = (m+1)/2;
	double hc = 3.0/(mc-1);
	int nxc = mc-2;
	int x1c = ((int)(COORD_X1 * (mc-1)) /3)  -1; 
    int x0c = (((int)(COORD_X0 * (mc-1)) + ((3 - ((int)(COORD_X0*(mc-1))%3))%3))/3)  -1;
    int y1c = ((int)(COORD_Y1 * (mc-1)) /3)  -1;
    int y0c = (((int)(COORD_Y0 * (mc-1)) + ((3 - ((int)(COORD_Y0*(mc-1))%3))%3))/3)  -1;
    */
	//peut pas calculer de nouveau m ou nx pour position du trou
    double h, invh2;
    int x0,x1,y0,y1, nx, n;

    computeParamTop(m, &h,&invh2,&y0,&y1,&x0,&x1,&nx, &n);

    *hl = h*pow(2, level);
	*invh2l = 3.0/((*hl)*(*hl)); 

    *x0l = x0 / pow(2, level); // va arrondir au point grille coarse a droite 
    *x1l = ((x1+1)/pow(2, level)) - 1; // permet si x1 pair on retire 1 
    *y0l = (y0/pow(2, level)); // arrondu coarse au dessus (permet de pas ajouter des points dans le trou)
    *y1l = ((y1+1)/pow(2, level)) - 1;
    int pl = *y1l - *y0l + 1;
    int ql = *x1l - *x0l + 1;
    

    *nxl = nx/pow(2,level); // nb de points coars sur un ligne pas bord
    *nl = ((*nxl) * (*nxl)) - (pl * ql);
    *nnzl =( 5 * (*nxl) * (*nxl)) - (4 * (*nxl)) ; 
    //nb de points concernés dans le trou:(compliqué a comprendre sans shema) (marche que si aumoins 3 points sur la largeur)
    int trousl = (5 * (pl-2) * (ql-2) + 4 * 2 * (pl-2) + 4 * 2 * (ql-2) 
                + 3 * 4 * 1 + 1 * 2 * pl + 1 * 2 * ql) + 2 * pl + 2 * ql;
    *nnzl -= trousl;
}

void computeHole(int *y0, int *y1, int *x0, int *x1, int m){

    *x0 = 0;
    *x1 = 0;
    *y0 = 0;
    *y1 = 0;
    //x0
    while ((*x0+1)*3 < COORD_X0*(m-1)){
        *x0 += 1;
    }
    
    //x1
    while ((*x1+1)*3 < COORD_X1*(m-1)){
        *x1 += 1;
    }
    if ((*x1+1)*3 != COORD_X1*(m-1)){
        *x1 -= 1;
    }
    //y0
    while ((*y0+1)*3 < COORD_Y0*(m-1)){
        *y0 += 1;
    }
    //y1
    while ((*y1+1)*3 < COORD_Y1*(m-1)){
        *y1 += 1;
    }
    if ((*y1+1)*3 != COORD_Y1*(m-1)){
        *y1 -= 1;
    }
}

double computeBound(double x, double y){
    
    return exp(sqrt(x*x + y*y));
}

int in_hole(int ix, int iy, int y0, int y1, int x0, int x1){
	
    if (ix >= x0 && ix <= x1 && iy >= y0 && iy <= y1){
		return 1;
	}
	else{
		return 0;
	}
}

int on_bound(int px, int py, int m){
	
    if (py == 0 || py == (m-1) || px == 0 || px == (m-1)){
		return 1;
	}
	else{
		return 0;
	}
}

//h1 et h2 car soit on check de grille a grille ou alors de grille a grille coarse
int check_nord(int ix, int iy, int y0, int y1, int x0, int x1, int nx){
    
	if (iy + 1 < nx && ! in_hole(ix, iy +1, y0,y1,x0,x1)){
		return 1;
	}
	else{
		return 0;
	}
}

int check_sud(int ix, int iy, int y0, int y1, int x0, int x1, int nx){
    
	if (iy > 0 && ! in_hole(ix, iy -1, y0,y1,x0,x1)){
		return 1;
	}
	else{
		return 0;
	}
}
int check_west(int ix, int iy, int y0, int y1, int x0, int x1, int nx){
    
	if (ix > 0 && ! in_hole(ix - 1, iy, y0,y1,x0,x1)){
		return 1;
	}
	else{
		return 0;
	}
}

int check_est(int ix, int iy, int y0, int y1, int x0, int x1, int nx){
    
	if (ix + 1 < nx && ! in_hole(ix+1, iy, y0,y1,x0,x1)){
		return 1;
	}
	else{
		return 0;
	}
}
int check_nw(int ixp, int iyp, int y0p, int y1p, int x0p, int x1p, int nxp){
    
	int ixnw = ixp - 1;//ind point nw sur la grille prolong 
    int iynw = iyp + 1;
    if (ixnw >= 0 && iynw < nxp && ! in_hole(ixnw, iynw, y0p,y1p,x0p,x1p)){
		return 1;
	}
	else{
		return 0;
	}
}

int check_ne(int ixp, int iyp, int y0p, int y1p, int x0p, int x1p, int nxp){
    
	int ixne =  ixp + 1;//ind point sw sur la grille prolong 
    int iyne = iyp + 1;
    if (ixne < nxp && iyne < nxp && ! in_hole(ixne, iyne, y0p,y1p,x0p,x1p)){
		return 1;
	}
	else{
		return 0;
	}
}

int check_sw(int ixp, int iyp, int y0p, int y1p, int x0p, int x1p, int nxp){
    //return 1/true lrsque point de la grille prolongé fait partie des variable
	int ixsw =  ixp - 1;//ind point sw sur la grille prolong 
    int iysw = iyp - 1;
    if (ixsw >= 0 && iysw >=0 && ! in_hole(ixsw, iysw, y0p,y1p,x0p,x1p)){
		return 1;
	}
	else{
		return 0;
	}
}

int check_se(int ixp, int iyp, int y0p, int y1p, int x0p, int x1p, int nxp){
    
	int ixse =  ixp+1;//ind point nw sur la grille prolong 
    int iyse = iyp - 1;
    if (ixse < nxp && iyse >= 0 && ! in_hole(ixse, iyse, y0p,y1p,x0p,x1p)){
		return 1;
	}
	else{
		return 0;
	}
}

int indice(int ix,int iy, int y0, int y1, int x0, int x1, int nx){ //ix iy -> indice dans matrice u (csr)
	
	int ind;
	int p = y1 - y0 + 1;
	int q = x1- x0 + 1;
	if (iy < y0 || (iy == y0 && ix < x0) ){//pas de retard
		ind = ix + (iy * nx);
	}
	else if ((ix > x1) && (iy >= y0 && iy < y1 )){ // retard en augmentation
		ind = ix + (iy * nx) - q*(iy - y0 +1) ;
	}
	else if ((ix < x0) && (iy >= y0 && iy <= y1 )){ // retard en augmentation
		ind = ix + (iy * nx) - q*(iy -1 - y0 +1) ; // -1 car retard ne change pas qd reste a gauche du trou et qu'on passe a la ligne sup
	}
	else {
		ind = ix + (iy * nx) - q*p; // retard constant
	}

	return ind;
}



/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


//                                 GRID                                     //

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


#define SCALE_FACT 0.5


/*fonctionnnement des level
-celui tt en haut avec le h le plus petit = level 0 puis coarse 1 coarse coarse 2 etc
- qd passe level en arg c'est tj le level auquel on est pas celui vers lequel on va
*/

int restrictR(int level, double *rp, double *rc, int m, int *nc){
    //level = 0 = ou on est ->1
    // level -> level +1
    

    //level -1 (prolonge donc monte dans la pyramyde)
    double hp, invh2p;
    int x0p,x1p,y0p,y1p, nxp, np, nnzp;
    computeParamLevel(m, level, &hp,&invh2p,&y0p,&y1p,&x0p,&x1p,&nxp, &np, &nnzp);
    
    //level ou on va =level+1
    double hc, invh2c;
    int x0c,x1c,y0c,y1c, nxc, nnzc;
    computeParamLevel(m, level+1, &hc,&invh2c,&y0c,&y1c,&x0c,&x1c,&nxc, nc, &nnzc);
    
    if (EXPLICIT){
		printf("\n Restrict -level %d:  hp = %lf ", level, hp);
		printf("np = %d nnzp = %d nxp = %d \n", np, nnzp, nxp);
		printf(" x0p = %d x1p = %d y0p = %d y1p = %d \n", x0p, x1p, y0p, y1p);
		printf("\n         -level %d hc = %lf ", level+1, hc);
		printf("nc = %d nnzc = %d nxc = %d \n", *nc, nnzc, nxc);
		printf("x0c = %d x1c = %d y0c = %d y1c = %d \n", x0c, x1c, y0c, y1c);
	}
    if (rc == NULL || rp == NULL){
        printf("r no memory \n");
        return 1;
    }
    int np_save = np;
	np = 0;
    int nc_save = *nc;
	*nc = 0;

	for (int iyp = 0; iyp < nxp; iyp++){
        //passage colonne suiv
        for (int ixp = 0; ixp < nxp; ixp++){      

            //exclu interieur et bord du trou
            if( ! in_hole(ixp,iyp,y0p,y1p,x0p,x1p)){
                
                //garde que ligne ind impair  & col impair 
                if (iyp % 2 == 1 && ixp % 2 == 1){  
                    rc[*nc] = 0;
                    //marquer le début de la ligne suivante dans le tableau 'ia'            
                    //replissage de la ligne : voisin sud //ui-1
                    // + verification si pas au dessus d'un bord
                    
                    if (check_sud(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){
                        int ind = indice(ixp,iyp-1,y0p,y1p,x0p,x1p, nxp);
                        rc[*nc] += 0.25 * rp[ind] * SCALE_FACT;
                    }    

                    else{
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT; 
                        //rc[*nc] += 0.25 * computeBound((ixp+1)*hp,(iyp+1 -1)*hp) * SCALE_FACT;//si fct qui calcule chaque fois :juste coord sinon
                    }

                    //replissage de la ligne : voisin ouest 
                    //si pas a droite d'un bord
                    if (check_west(ixp,iyp,y0p,y1p,x0p,x1p,nxp)){
                        int ind = np -1;
                        rc[*nc] += 0.25 * rp[ind] * SCALE_FACT;
                    }
                    else{
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT;
                        //rc[*nc] += 0.25 * computeBound((ixp+1-1)*hp,(iyp+1)*hp) * SCALE_FACT;
                    }

                    int ind = np;
                    rc[*nc] += rp[ind] * SCALE_FACT;
                   
                    // replissage de la ligne : élém. diagonal
                    

                    // replissage de la ligne : voisin est
                    //si pas a gauche d'un bord
                    if (check_est(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){
                        int ind = np +1;
                        rc[*nc] += 0.25 * rp[ind] * SCALE_FACT;
                         
                    }
                    else{
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT;
                        //rc[*nc] += 0.25 * computeBound((ixp+1+1)*hp ,(iyp+1)*hp) * SCALE_FACT;
                        
                    }

                    // replissage de la ligne : voisin nord
                    //si pas en dessous d'un bord
                    if ( check_nord(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){

                        rc[*nc] += 0.25 * rp[indice(ixp,iyp+1,y0p,y1p,x0p,x1p, nxp)] * SCALE_FACT;
                    }
                    else{
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT;
                        //rc[*nc] += 0.25 * computeBound((ixp+1)*hp, (iyp+1+1)*hp) * SCALE_FACT;
                    }
                    // numéro de l'équation
                    *nc += 1;
                }
                np += 1;
            }
        }
    }
	if (np_save != np){
		printf("nr != n\n");
		return 1;
	}
	if (nc_save != *nc){
		printf("nc_save %d != nc %d\n", nc_save, *nc );
		return 1;
	}
	return 0;
}

int addProlCorrection(int level, double *up, double *uc, int m, int *np){

    //level 1 = on ou on est -> on va vers 0
    // level -> level -1

    //level level
    double hc, invh2c;
    int x0c,x1c,y0c,y1c, nxc, nc, nnzc;

    computeParamLevel(m, level, &hc,&invh2c,&y0c,
                      &y1c,&x0c,&x1c,&nxc, &nc, &nnzc);


    //level -1 car on remonte dans la pyr inv
    double hp, invh2p;
    int x0p,x1p,y0p,y1p, nxp, nnzp;

    computeParamLevel(m, level-1, &hp,&invh2p,&y0p,
                      &y1p,&x0p,&x1p,&nxp, np, &nnzp);

	if (EXPLICIT){
		printf("\n Prolong -level %d:  hp = %lf ", level-1, hp);
		printf("np = %d nnzp = %d nxp = %d \n", *np, nnzp, nxp);
		printf(" x0p = %d x1p = %d y0p = %d y1p = %d \n", x0p, x1p, y0p, y1p);
		printf("\n         -level %d hc = %lf ", level, hc);
		printf("nc = %d nnzc = %d nxc = %d \n", nc, nnzc, nxc);
		printf("x0c = %d x1c = %d y0c = %d y1c = %d \n", x0c, x1c, y0c, y1c);
	}
    
    if (up == NULL || uc == NULL){
        printf("no memory  prol\n");
        return 1;
    }
	
    int nc_save = nc;
	nc = 0;
    int np_save = *np;
	*np = 0;

	int ind;
	for (int iyp = 0; iyp < nxp; iyp++){
        //passage colonne suiv
        for (int ixp = 0; ixp < nxp; ixp++){      

            //exclu interieur et bord du trou
            if(! in_hole(ixp,iyp,y0p,y1p,x0p,x1p)){
                
                //up[*np] = 0; ici on add la prolongation au vecteur de niveau au dessus

				//impair impair -> 1/4 somme des 4 autour
                
                if (iyp % 2 == 0 && ixp % 2 == 0){ 
                    //check si bord est bien la ou pense etre pour uc (pas faire une recherche alors que bord)
                    //somme coin gauche bas 
                    //check si point prol pas ds le trou
                    
                    if (check_sw(ixp,iyp,y0p,y1p,x0p,x1p,nxp)){ // question est ce que prolong peut etre domaine et coarse ds un board ? non ca depend que de p
                        
						ind = indice((ixp/2) - 1, (iyp/2)-1, y0c,y1c,x0c,x1c, nxc);//juste /2 -1 car point prol au milieu des 4 tjrs pair
                        up[*np] += 0.25 * uc[ind]; 
                    }
                    /*else{
                        double bound = computeBound((ixp+1-1)*hp,(iyp+1-1)*hp);
                        
                        up[*np] += 0.25 * bound;
                       
                    }*/
                    //coin droit bas
                    if (check_se(ixp,iyp,y0p,y1p,x0p,x1p,nxp)){ //cond droit
                        ind = indice((ixp/2) ,(iyp/2)-1, y0c,y1c,x0c,x1c, nxc);
                        up[*np] += 0.25 * uc[ind];
                        
                    }
                   /* else{
                        up[*np] += 0.25 * computeBound((ixp+1+1)*hp,(iyp+1-1)*hp);
                    }*/
                    //coin haut gauche
                    if (check_nw(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){ //cond gauche haut
                        ind = indice((ixp/2 - 1),(iyp/2), y0c,y1c,x0c,x1c, nxc);
                        up[*np] += 0.25 * uc[ind];
                    }
                    /*else{
                        double bound = computeBound((ixp+1-1)*hp,(iyp + 1+1)*hp);
                        up[*np] += 0.25 * bound;
                        
                    }*/
                    //coin droit haut
                    
                    if (check_ne(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){ //coin droit
                        ind = indice((ixp/2),(iyp/2), y0c,y1c,x0c,x1c, nxc);
                        up[*np] += 0.25 * uc[ind];

                    }
                    /*else{
                        up[*np] += 0.25 * computeBound((ixp+1+1)*hp,(iyp+1+1)*hp);
                    }*/
                }

				//impair pair => somme haut + bas
				else if (iyp % 2 == 0 && ixp % 2 == 1){
                    
                    //somme bas
                    if (check_sud(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){
                       
                        ind = indice((ixp/2),(iyp/2)-1, y0c,y1c,x0c,x1c, nxc);
                        up[*np] += 0.5 * uc[ind];
                    }
                    /*else{
                        up[*np] += 0.5 * computeBound((ixp+1)*hp,(iyp+1-1)*hp);
                    }*/
                    //somme haut
                    if(check_nord(ixp,iyp,y0p,y1p,x0p,x1p,nxp)){
                        ind = indice((ixp/2),(iyp/2), y0c,y1c,x0c,x1c, nxc);
                        up[*np] += 0.5 * uc[ind]; //pas de nx/2+1 car point uc[nc] deja ligne du haut
                    }
                    /*else{
                        up[*np] += 0.5 * computeBound((ixp+1)*hp,(iyp+1+1)*hp);
                    }*/
                }
                //pair impair 1/2 somme gauche droite
                else if (iyp % 2 == 1 && ixp % 2 == 0){

                    //somme gauche
                    //si pas a droite d'un bord
                    if (check_west(ixp,iyp,y0p,y1p,x0p,x1p,nxp)){
                        up[*np] += 0.5 * uc[nc - 1];
                    }
                    /*else{
                        up[*np] += 0.5 * computeBound((ixp+1-1)*hp,(iyp+1)*hp);
                    }*/

                    //somme droit
                    //si pas a gauche d'un bord
                    if (check_est(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){
                        up[*np] += 0.5 * uc[nc]; //nc car deja +1 car type point precedent + 1
                    }
                    /*else{
                        up[*np] += 0.5 * computeBound((ixp+1+1)*hp,(iyp+1)*hp);
                    }*/  
                }
                //ligne impair  & col impair => elem identique
                else if (iyp % 2 == 1 && ixp % 2 == 1){  //pour 1 meme ligne alterne entre type 1 et 2, commence par 1 finis par 1

                    up[*np] += uc[nc];
                    nc += 1;       //=> 2 choses, pr point type 2 droite = nc gauche == nc-1
                    // et aussi que pour type 3 et 4 nc rpz point au dessus tt a gauche
                }
                *np += 1;
            }
        }
    }
    
    if (*np != np_save){
        printf(" err np %d npcheck %d\n", *np, np_save);
    }
    if (nc != nc_save){
        printf(" err nc %d nccheck %d\n", nc, nc_save);
    } 
	return 0;
}



/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


//                                 TEST                                     //

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


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