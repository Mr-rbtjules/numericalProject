#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "proto.h"

#define BOUND1(x, y) (exp(sqrt((x)*(x) + (y)*(y))))
#define BOUND2(x, y) (sin(sqrt((x)*(x) + (y)*(y))))
#define BOUND0(x, y) (0)

//program parameters
//type of boundary conditions
#define BOUND BOUND2
//type of domain$
#define DOMAIN1 "[0, 6] × [0, 10] [2, 4] × [2, 5]"
#define DOMAIN2 "[0, 6] × [0, 10] [0, 2] × [3, 4]"
#define DOMAIN3 "[0, 6] × [0, 10] [4, 6] × [3, 6]"
#define DOMAIN4 "[0, 6] × [0, 10] [2, 3] × [5, 8]"
#define DOMAIN5 "[0, 6] × [0, 10] [3, 4] × [0, 4]"
#define DOMAIN6 "[0, 6] × [0, 10] [0, 4] × [0, 3]"
#define DOMAIN7 "[0, 6] × [0, 10] [2, 4] × [5, 10]"
#define DOMAIN8 "[0, 6] × [0, 10] [0, 4] × [3, 10]"
#define DOMAIN DOMAIN1
//discretisation 
#define M 13
#define LEVELMAX 1 
#define MU1 1
#define MU2 1
#define MODE 1


/*
        level 0
\      /1
 \    /2
  \  /3
   \/4 = level max = nb de fois qu'on restr
*/
//print parameters
#define EXPLICIT 1
#define LOAD 1
#define CHRONO 1

//start global variable (m and domain)
globVal_s globVal = { NULL, NULL, NULL, NULL}; 



/*
Details : accorder computeresnorm mgmethod et tt
pour avoir mm ordre ia ja u r, etc
*/




/// RESOL ///


int mg_method(int iter){
	
	
	//initit memory and pointers
	//compute all the coarse matrix and nl
    int *nl = NULL;

    //liste : chaque elem pointe vers 1 matrice 
    //d'un certain level (ici matrice == liste)
	int *ial = NULL; 
	int *jal = NULL;
	double *al = NULL;                                                         
	double *rl = NULL;
	double *ul = NULL;
	double *dl = NULL;
	
    //attention aux autres b prc pt on les utilise pour autre
	double *b = NULL;  

	//we pass the adress for each level
	allocGrids(&nl, &ial, &jal, &al, &b, &dl, &rl, &ul);
    
	//precomputation of all Ac and bc

	for (int l = 0; l <= LEVELMAX; l++){
		//adresses of each vector is simply the pointer 
        //shifted
		probMg(
            globVal.m[l], 
            nl + l, 
            ial + globVal.vectStart[l] + l,
            jal + globVal.matStart[l], 
            al + globVal.matStart[l],
            b
        );
	}
	
	int startLevelTg = 0; 

	//u0 = 0
	initialization(nl, ial, jal, al, b, rl, ul);


    
	
    if (LOAD){
		printf("\n Load percentage : \n");
	}
	//start iterations of the multigrid cycle (tg_rec)
	double t1 = mytimer();
    if( mg_iter(
            iter, LEVELMAX, globVal.m[0], 
            MU1, MU2, nl, ial ,jal,al,rl,ul,dl,b)){
        return 1;
    }
    double t2 = mytimer();
    if (CHRONO){
        printf("\nTemps de solution multigrid method (CPU): %5.1f sec\n",t2-t1);   
    }


    printf("\n final res : %lf\n", computeResNorm(
            *nl, ial, jal, al, b, ul, rl
        ));
    
        

    plot_res(rl + globVal.vectStart[startLevelTg], startLevelTg);
    free(nl);
    free(ial);
    free(jal);
    free(al);
    free(b);
    free(dl);
    free(ul);
    free(rl);
	return 0;
}

int mg_iter(int iter, int levelMax, int m, int mu1, int mu2, int *nl, int *ial,
             int *jal, double *al, double *rl, double *ul, double *dl, double *b){
	
	int startLevelTg = 0;
	for (int i = 0; i < iter; i++){
		//initialisation ici ?
		if (LOAD){
			printf("   %.3f   ", (double)(i+1)/iter);
        	fflush(stdout);
			printf("\r");
		}
        
		tg_rec( startLevelTg+1, levelMax, m, mu1, mu2, 
				nl, ial, jal, al, rl, ul, b, dl);
		
		//print(res et u) + save
		
	}
	return 0;
}

int tg_rec(int level, int levelMax, int m, int mu1,
			int mu2, int *nl, int **ial, int **jal,
		   double **al, double **bl, double **ul, double **rl, double **dl){
		
		//des le level 1 on stock pas u2 mais c2 puis 'deviennent des u lorsqu'on ajoute la correction en remontant
		//pblm tt en haut smoothing Au=b mais en dessous sm Ac = r pareil computeres
	
	if (level < levelMax){

        //mettre dans une fonction car si level = 0 applique syst Au=b et r = b-Au
        //sinon Ac = r 
        //deux possibilités soit :
        //1: on stock dans les b inf le res r-Ac et dans r les 'b' inf qui sont enft r dans le cours
        //2: on stock dans les b d'abord b puis
		/*//mu1 smoothing iterations
		forwardGS(mu1, nl[level] , ial[level], jal[level],
					al[level], bl[level], ul[level], rl[level], dl[level]);//ici on utilise b pour stocker le residu de Ac=r et stock c dans u
		
		//restrict r : r-Ac stocker dans b ! etape en trop non necessaire&
		computeRes(nl[level], ial[level], jal[level],
					al[level], ul[level], bl[level], rl[level]);
		//restrict residu stocker dans b
		restrictR(level, rl[level],rl[level+1], m, &(nl[level+1]));*/
		
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


/// PROB ///

int probMg(int m, int *n, 
		   int *ia, int *ja, double *a, double *b){
    
    //fait pour m/2
    //mémoire deja allouée
    double h, invh2;
    int x0,x1,y0,y1, nx, ny, nnz;

    computeParamLevel(m, &h,&invh2,&x0,&x1, &y0, &y1, &nx, &ny, n, &nnz);
    if (EXPLICIT){

    	printf("\n ProbMg h = %lf ", h);
    	printf("n = %d nnz = %d nx = %d ny= %d \n", *n, nnz, nx, ny);
    	printf("x0 = %d x1 = %d y0 = %d y1 = %d \n", x0, x1, y0, y1);
	}
	
    int nnz_save = nnz;
    nnz = 0;
    
    //passage ligne suiv(plaque complete)
    int ind = 0;
    //iy ix indice sur grille hors bords  mais position (ix+1)*h
    for (int iy = 0; iy < ny; iy++) { 
        for (int ix = 0; ix < nx; ix++) {
            //exclu interieur et bord du trou
            if(! in_hole(ix,iy,y0,y1,x0,x1)){
                //marquer le début de la ligne suivante dans le tableau 'ia'
                ia[ind] = nnz;
                
                b[ind] = 0;
               
                if (check_sud(ix,iy,y0,y1,x0,x1,ny)){
                    a[nnz] = -invh2;
                    ja[nnz] = indice(ix,iy-1,y0,y1,x0,x1, nx);
                    nnz++;
                }
                else{//besoin de b que pour top level
					if (m == globVal.m[0]){
						double x = (ix+ 1)*h;
						double y = (iy + 1 -1)*h;
						double bound = computeBound(x, y);
                    	b[ind] += bound * invh2; 
					}
                }

                //replissage de la ligne : voisin ouest 
                //si pas a droite d'un bord
                if (check_west(ix,iy,y0,y1,x0,x1,nx)){
                    a[nnz] = -invh2;
                    ja[nnz] = ind - 1;
                    nnz++;
                }
                else{
					if (m == globVal.m[0]){
						double x = (ix + 1 - 1)*h;
						double y = (iy + 1)*h;
						double bound = computeBound( x, y);
                    	b[ind] += bound * invh2;
					}
                }

                // replissage de la ligne : élém. diagonal
                a[nnz] = 4.0*invh2;
                ja[nnz] = ind;
                
                nnz++;
                
                // replissage de la ligne : voisin est
                //si pas a gauche d'un bord
                
                if ( check_est(ix,iy,y0,y1,x0,x1,nx)){
                    a[nnz] = -invh2;
                    ja[nnz] = ind + 1;
                    nnz++;
                }
                else{
					if (m == globVal.m[0]){
						double x = (ix + 1 +1)*h;
						double y = (iy + 1)*h;
						double bound = computeBound(x, y);	
						b[ind] += bound * invh2;
					}
                }

                // replissage de la ligne : voisin nord
                //si pas en dessous d'un bord
                if ( check_nord(ix,iy,y0,y1,x0,x1,ny) ){
                        a[nnz] = -invh2;
                        ja[nnz] = indice(ix,iy+1,y0,y1,x0,x1, nx);
                        nnz++;
                }
                else{
					if (m == globVal.m[0]){
						double x = (ix + 1)*h;
						double y = (iy + 1 +1)*h;
						double bound = computeBound(x, y);
						b[ind] += bound * invh2;
					}
                }
                // numéro de l'équation
                ind += 1;
            }
        }
    }

     if (*n != ind){
        printf(" err nl %d ind %d\n", *n, ind);
    }
    else if (nnz != nnz_save){
        printf(" err nnz %d nnz_save %d\n", nnz, nnz_save);
    }
    else {
        ia[ind] = nnz;
    }

	return 0;
}

int initialization(int *n0, int *ia0, int *ja0, double *a0,
						 double *b, double *r0, double *u0){

	//arbitrary u0
	for (int i=0; i< *n0; i++){
		u0[i] = 0; 
	}
	if (EXPLICIT){
		printf("\n initial res : %lf\n", computeResNorm(
            *n0, ia0, ja0, a0, b, u0, r0
        ));
	}
	return 0;
}

int allocGrids(int **nl, int **ial,
               int **jal, double **al, double **b,
			   double **dl, double **rl, double **ul){ 
    *nl = (int*)malloc((LEVELMAX +1)*sizeof(int));
    // pour taille de ia : c'est le nb de n total
    // (dernier de la list vectStart mais pour chaque level on a 1 eleme
    // de plus que n, nb de niveaux = levelmAX +1
    *ial = (int*)malloc(
        (globVal.vectStart[LEVELMAX+1] + (LEVELMAX+1)) * sizeof(int)
    ); 
    //levelmax = indice du debut du dernier niveau donc +1 p
    // our avoir le nb tot d'elem matriciel
	*jal = (int*)malloc(globVal.matStart[LEVELMAX+1] * sizeof(int));
	*al = (double*)malloc(globVal.matStart[LEVELMAX+1] * sizeof(double));
    //elem vectoriel
	*rl = (double*)malloc(globVal.vectStart[LEVELMAX+1] * sizeof(double));
	*dl = (double*)malloc(globVal.vectStart[LEVELMAX+1] * sizeof(double));
	*ul = (double*)malloc(globVal.vectStart[LEVELMAX+1] * sizeof(double));
    //on a besoin que de b au debut voir multigrid algo description
    //2eme elem de vectstart = n pour le top level
	*b = (double*)calloc(globVal.vectStart[1], sizeof(double));


    if (*nl == NULL || *b == NULL || *ial == NULL || *jal == NULL || 
        *al == NULL ||*dl == NULL || *rl == NULL || *ul == NULL){
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n");
        return 1;
    }
    if (EXPLICIT){
        printf("\nMemory allocated for level 0 to level %d\n", LEVELMAX);
    }
	
	return 0;
}

void computeParamLevel(int m, double *h, double *invh2, int *x0,
						int *x1, int *y0, int *y1, int *nx,
						int *ny, int *n, int *nnz){
    
    int perUnit = (m-1)/(globVal.domain[1] - globVal.domain[0]);
    
    *x0 = ((globVal.domain[4] - globVal.domain[0]) * perUnit) - 1;
    *x1 = ((globVal.domain[5] - globVal.domain[0]) * perUnit) - 1;
    *y0 = ((globVal.domain[6] - globVal.domain[2]) * perUnit) - 1;
    *y1 = ((globVal.domain[7] - globVal.domain[2]) * perUnit) - 1;
    *h = (double)1/perUnit;
    *invh2 = 1.0/((*h)*(*h));
    *nx = m-2; //par def
    *ny = ((globVal.domain[3] - globVal.domain[2])*perUnit) -1;
    getNnz(*nx, perUnit, *x0, *x1, *y0,* y1, n, nnz);
}

void getNnz(int nx, int perUnit, int x0, int x1, 
					int y0, int y1, int *n, int *nnz){

    int p = y1 - y0 + 1;
    int q = x1 - x0 + 1;
    
    //1 chaque point donne 5 element
    int dx = nx; //nb de points a compter
    int dy = ((globVal.domain[3] - globVal.domain[2]) * perUnit) - 1;
    
    //generalise pour trou sur bord exterieur
    //on rajoute ce qu'on a retiré en trop
    
    if ((x0<=0 && y0<=0) || (x0<=0 && y1>=dy) 
        || (x1>= dx && y0<=0) ||(x1>=dx && y1>=dy)) { // dans 1 coin
    //les q point de la largeur partage 1 point 
    //avec les p points de la longeur
        
        *nnz = 5*dx*dy - 5*(p-1) * (q-1);
        *n = (dx * dy) - ((p-1) * (q-1));
        
        //2 on retire les points qui depassent sur les bord exterieurs
        *nnz -= 2*dx + 2*dy;
    }
    else if (x0 <= 0 || x1 >= dx){
        *nnz = 5*dx*dy;
        *n = (dx * dy) - (p * q);
        //2 on retire les points qui depassent sur les bord exterieurs
        *nnz -= 2*dx + 2*dy;
        //3 pour ceux qui compte mais qui depasse sur le bord du trou
        *nnz -= 2*p;
        //4 
        *nnz -= 2*q;
        //5 les 5 genere par les points du trou
        *nnz -= 5*p*q;
        //revient sur 3 car deja été retiré au 2 
        //c'est comme si le trou devenait un bord ext
        *nnz += p;
        //revient de nouveau sur le 3 car pas besoin de retirer a gauche
        //car pas de point existant a gauche
        *nnz += p;
        //haut et bas à gauche on a 2 points qui font partie du bord
        //exterieur on revient sur le 4 
        *nnz += 2; 
        //revient sur le 5, car on retire des points
        //qui n'ont jamais ete compté dans 1
        *n += p; //car retiré a n p points qui ne font pas partie du domaine
        *nnz += 5*p;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
    }
    else if (y0 <= 0 || y1 >= dy){
        
        *nnz = 5*dx*dy;
        *n = (dx * dy) - (p * q);
        //2 on retire les points qui depassent sur les bord exterieurs
        *nnz -= 2*dx + 2*dy;
        //3 pour ceux qui compte mais qui depasse sur le bord du trou
        *nnz -= 2*p;
        //4 
        *nnz -= 2*q;
        //5 les 5 genere par les points du trou
        *nnz -= 5*p*q;
        //revient sur 3 car deja été retiré au 2 
        //c'est comme si le trou devenait un bord ext
        *nnz += q;
        //revient de nouveau sur le 4 car pas besoin de retirer a gauche
        //car pas de point existant a gauche
        *nnz += q;
        //haut et bas à gauche on a 2 points qui font partie du bord
        //exterieur on revient sur le 4 
        *nnz += 2; 
        //revient sur le 5, car on retire des points
        //qui n'ont jamais ete compté dans 1
        *nnz += 5*q; 
        *n += q;
    }
    else {

        *nnz = 5*dx*dy;
        *n = (dx * dy) - (p * q);
        //2 on retire les points qui depassent sur les bord exterieurs
        *nnz -= 2*dx + 2*dy;
        //3 pour ceux qui compte mais qui depasse sur le bord du trou
        *nnz -= 2*p;
        //4 
        *nnz -= 2*q;
        //5 les 5 genere par les points du trou
        *nnz -= 5*p*q;
    }
}

int on_bound(int px, int py, int mx, int my){
	
    if (py == 0 || py == (my-1) || px == 0 || px == (mx-1)){
		return 1;
	}
	else{
		return 0;
	}
}

int in_hole(int ix, int iy, int y0, int y1, int x0, int x1){
	//1 == True
    if (ix >= x0 && ix <= x1 && iy >= y0 && iy <= y1){
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

int indice(int ix,int iy, int y0, int y1, int x0, int x1, int nx){ 
    //ix iy -> indice dans matrice u (csr)
	
	int ind;
	int p = y1 - y0 + 1;
	int q = x1 - x0 + 1;
    //pour les trous speciaux
    if (x1 == nx){
        q -= 1;
    }
    if (y0 == -1){
        p -= 1;
        y0 = 0; //enleve probleme de retard au debut
    }
    if (x0 == -1){
        q -= 1;
        x0 = 0; //enleve probleme de retard au debut
    }
    
    if (iy < y0){//pas de retard
		ind = ix + (iy * nx);
	}
    //au niveau du debut du trou a gauche
	else if (iy == y0 && ix < x0 ){//pas de retard nn plus
		ind = ix + (iy * nx);
	}
	else if ((ix > x1) && (iy >= y0 && iy < y1 )){ // retard en augmentation
		
        ind = ix + (iy * nx) - q*(iy - y0 +1) ;
	}
	else if ((ix < x0) && (iy >= y0 && iy <= y1 )){ // retard en augmentation
    // -1 car retard ne change pas qd reste a gauche du
    // trou et qu'on passe a la ligne sup
		ind = ix + (iy * nx) - q*(iy -1 - y0 +1) ; 
	}
	else {
		ind = ix + (iy * nx) - q*p; // retard constant
	}
	return ind;
}








 //TOOLS//

double computeBound(double x, double y){
    //compute the value of a point on the bound of the domain
	//with a function BOUND defined in the macros
    return BOUND(x,y);
}

//pt chnger n en *n a voir
double computeResNorm(int n, int *ia, int *ja, double *a, 
						double *b, double *u,  double *r){

	double rn = 0;

	computeRes(n,ia,ja,a,u,b,r);
	for (int i = 0; i < n; i++){
		rn += r[i] * r[i];
	}
	rn = sqrt(rn);
	return rn*10000000000000000;
}

int computeRes(int n, int *ia, int *ja, double *a, 
			   		double *u, double *b, double *r){
	
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

void printVect(void *vect, int size, int type) {
    int i;
    printf("\n[");
    if (type == 0) {  // Assuming type 0 represents int
        int *intVect = (int *)vect;
        for (i = 0; i < size; i++) {
            printf("%d, ", intVect[i]);
        }
    } 
    else if (type == 1) {  // Assuming type 1 represents double
        double *doubleVect = (double *)vect;
        for (i = 0; i < size; i++) {
            printf("%lf, ", doubleVect[i]);
        }
    printf("]\n");
    }
}

void plot_res(double *r, int level){

    //cree un fichier texte (.dat pour gnuplot) contenant 2 colonnes
    //avec les coordonnées x et y et une colonnes avec les éléments de u
    // (bord et trous compris avec u = 0) et lance le script gnuplot pour
    // plot le resultat (nmp = numero de mode propre)

    double h, invh2;
    int x0,x1,y0,y1, nx, ny, nnz;
    int n;
    computeParamLevel(
        globVal.m[level], &h,&invh2,&x0,&x1, &y0, &y1, &nx, &ny, &n, &nnz
    );
    
    int mx = nx+2;
    int my = ny+2;

    //creation du ficher
    FILE* pointFile = NULL;                        
    const char* file_name = "coord_stat.dat";
    // ouverture en mode ecriture
    pointFile = fopen(file_name,"w+");           
    
    if (pointFile != NULL){

        //Calcul des constantes
        
        int ind = 0;
        //indice 00 different que pour prob, ici c'est le coin du domaine
        for (int py = 0; py < my; py++){ //passage ligne suivante
            for (int px = 0; px < mx; px++){      //passage colonne suivante
                
                //si sur bord ou trou
                if ( on_bound(px,py,mx,my) || in_hole(px-1,py-1,y0,y1,x0,x1)){
                    
                    //fprintf(pointFile, "%.16g %.16g 0\n", px*hl, py*hl);
                    fprintf(pointFile, "%.16g %.16g NaN\n", px*h, py*h);
                }
                else{
                    fprintf(
                        pointFile,
                        "%.16g %.16g %.16g\n",
                        (px*h),
                        (py*h),
                        r[ind]
                    );
                    ind += 1;
                }
            }
            fprintf(pointFile, "\n");
        }
        //fermeture du ficher
        fclose(pointFile);
    }

    //initialise fichier gnuplot
    //plot le graphique

    FILE *gnuplot = popen("gnuplot -persistent","w");
    fprintf(gnuplot, "set title 'm = %d level = %d'\n", mx, level);
    fprintf(gnuplot, "load 'gnuScript_stat.gnu'\n");
    fclose(gnuplot);

    //supprime coord_stat.dat
    int _ = system("rm coord_stat.dat");
}




/// GLOBAL VARIABLES ///

void initGlobVal(){
    int count;
    extractDomain(DOMAIN, &count, &(globVal.domain));
    
    int m = correctM(globVal.domain, M);
    
    initMlevels(m, &(globVal.m));
    
    initIndex(&(globVal.vectStart), &(globVal.matStart));
}

void extractDomain(const char* str_domain, int* count, int **domain) {
	//array with the int describing the domain
    *domain = (int*)malloc(16 * sizeof(int));
    *count = 0;
    //on all the string
    while (*str_domain) {
		//check if we start on a number
        if (*str_domain >= '0' && *str_domain <= '9') {
			//change all the digits present until another char into 1 number
			//(1 digit or more)
            (*domain)[*count] = atoi(str_domain);
            (*count)++;

            while (*str_domain >= '0' && *str_domain <= '9') {
                str_domain++;
            }
        } else {
            str_domain++;
        }
    }
}

int correctM(int *domain, int m){
    //correct m in order to obtain discretisation 
    //matching with the surface and for multigrid =>
    //domain ex : if "[0, 6] × [0, 10]  [2, 4] × [2, 5]" -> domain = [0,6,0,10,2,4,2,5]
    // L = smallest number of cases possible
    int L = domain[1] - domain[0];
    while (m > L+1){
        L *= 2;
    }
    return L + 1;
}

int initMlevels(int m, int **mLevels){
    
    *mLevels = (int *)malloc((LEVELMAX+1)*sizeof(int));
    (*mLevels)[0] = m;
    for (int l = 1; l <= LEVELMAX; l+=1){
        (*mLevels)[l] = m/(1 << l) + 1;
    }
    
    return 0;
}

int initIndex( int **vectStart, int **matStart){
    
    //contient les debuts de chaque vecteur (pour tt les n) et la fin 
    //du vecteur global donc nb de levels = levelmax +1 , et encore +1 
    //pour la fin de vecteur un peu comme dans ia
	*vectStart = (int*)malloc((LEVELMAX+1 + 1) * sizeof(int));
    *matStart = (int*)malloc((LEVELMAX+1 + 1) * sizeof(int));
    
    (*vectStart)[0] = 0;
    (*matStart)[0] = 0;
    int nnzTot = 0;
    int nTot = 0;
    for (int l = 0; l <= LEVELMAX; l++){
        double h, invh2;
        int x0,x1,y0,y1, nx, ny, nnz;
        int n;	
		computeParamLevel(globVal.m[l], &h, &invh2, &x0,
						&x1, &y0, &y1, &nx,
						&ny, &n, &nnz);
        
        nnzTot += nnz;
        nTot += n;
        (*vectStart)[l+1] = nTot;
        (*matStart)[l+1] = nnzTot;
	}
    return 0;
}

void freeGlobVal(){
    free(globVal.domain);
    free(globVal.vectStart);
    free(globVal.matStart);
    free(globVal.m);
}




//NOT USED //
int allocProb(int m, int *n, int **ia, int **ja, 
     		  double **a, double **b, double **u, double **r){

    double hl, invh2l;
    int x0l,x1l,y0l,y1l, nxl, nyl, nnzl;
    computeParamLevel(m, &hl,&invh2l,&y0l,&y1l,
					  &x0l,&x1l,&nxl, &nyl, n, &nnzl);
    
    (*ia) = malloc((*n + 1) * sizeof(int));
    (*ja) = malloc(nnzl * sizeof(int));
    (*a) = malloc(nnzl * sizeof(double));
    (*b) = malloc(*n * sizeof(double));
    (*u) = malloc(*n * sizeof(double));
    (*r) = malloc(*n * sizeof(double));
    
    if (*b == NULL || *ia == NULL || *ja == NULL || 
        *a == NULL || *r == NULL || *u == NULL){
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n");
        return 1;
    }
	if (EXPLICIT){
		printf("\n Alloc level %d : hl = %lf nl ",0, hl);
    	printf(" = %d nnzl = %d nxl = %d nyl = %d\n", *n, nnzl, nxl, nyl);
    	printf("x0l = %d x1l = %d y0l = %d y1l = %d \n", x0l, x1l, y0l, y1l);
	}
    return 0;
}