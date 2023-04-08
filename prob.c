/* les modifications sont marquées avec "<--" */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "proto.h"


#define COORD_X0 1.0 //col bord gauche
#define COORD_X1 2.5 // col bord droit
#define COORD_Y0 1.5 //ligne bord bas
#define COORD_Y1 2.0 //ligne bord haut

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

double computeResNorm(int n, int *ia, int *ja, double *a, double *u, double *b, double *r){

	double rn = 0;

	computeRes(n,ia,ja,a,u,b,r);
	for (int i = 0; i < n; i++){
		rn += r[i] * r[i];
	}
	rn = sqrt(rn);
	return rn;
}


int prob(int m, int *n, int **ia, int **ja, double **a, double **b)
/*
   But
   ===
   Génerer le système linéaire n x n 
                          
                             Au = b                                   

   qui correspond à la disrétisation sur une grille cartesienne 
   regulière m x m de l'équation de Poisson à deux dimensions
              
            d    d        d    d
         - == ( == u ) - == ( == u )  = 0     sur [0,3] x [0,3]
           dx   dx       dy   dy

  avec les conditions aux limites de Dirichlet
         
        u = exp(sqrt(x**2 + y**2)) R \ ([0,3]**2 \ [1,2.5] x [1.5,2])
  
  La numérotation des inconnues est lexicographique, la direction x étént 
  parcourue avant celle de y. La matrice A est retournée dans le format 
  CRS qui est défini via trois tableaux : 'ia', 'ja' et 'a'.

  Arguments
  =========
  m  (input)  - nombre de points par direction dans la grille 
  n  (output) - pointeur vers le nombre d'inconus dans le système
  ia (output) - pointeur vers le tableau 'ia' de la matrice A
  ja (output) - pointeur vers le tableau 'ja' de la matrice A
  a  (output) - pointeur vers le tableau 'a' de la matrice A
  b  (output) - pointeur vers le tableau 'b'

*/
{
    
    /*Le pblm est ici malheureusement u de probmg depend de m alors que devrait pas ?
    pq bonne forme mais pas bonne intensité ?s*/
    
    
    
    int  ix, iy, ind = 0;
    
    double h = 3.0/(double)(m-1);
    double invh2 = 1.0/(h*h);

    // coordonnees du trou sur la grille discrete
    
    /*int x1 = ((int)(COORD_X1 * (m-1)) /3)  -1; 
    int x0 = (((int)(COORD_X0 * (m-1)) + ((3 - ((int)(COORD_X0*(m-1))%3))%3))/3)  -1;
    int y1 = ((int)(COORD_Y1 * (m-1)) /3)  -1;
    int y0 = (((int)(COORD_Y0 * (m-1)) + ((3 - ((int)(COORD_Y0*(m-1))%3))%3))/3)  -1;
    */
    int x0,x1,y0,y1;
    computeHole(&x0,&x1,&y0,&y1, m);
    
    /*
    
    *m-1 car points sur la grille sont a des fraction de m-1

    /3 entier : 1 pierre 2 coup on arrondi au numerateur multiple de 3 le proche (floor
    pour les x1 y1 et pour les x0 y0  + (3 - (m-1) %3)%3 permet ceil au multiple de 3 exemple :
    6-> 6 , 7-> 9 ect)
     de notre point
     + on divise par 3 donc on a le direct le numero du point (0 -> m)
     -1 car on considere que les points pas sur le bords (nx) qui commence par l'indice 0 donc jusque nx-1
    */
    
    int nx = m-2;
    //nb de points sur la largeur du trous

    int p = y1 - y0 + 1;//m%2 + m/6; //== nb de points entre y1 et y0
    // "" sur la longueur
    int q = x1- x0 + 1;//m%2 + m/2;
    //plaque hors bords et trous
    *n = nx * nx - (p * q);
    //nombre d'element non nul pour membrane de base
    int nnz = 5 * nx * nx - 4 * nx ; 
    //nb de points concernés dans le trou:(compliqué a comprendre sans shema) (marche que si aumoins 3 points sur la largeur)
    int trous = (5 * (p-2) * (q-2) + 4 * 2 * (p-2) + 4 * 2 * (q-2) 
                + 3 * 4 * 1 + 1 * 2 * p + 1 * 2 * q) + 2 * p + 2 * q;
    nnz -= trous;
    
    
    /* allocation des tableaux */

    *ia  = malloc((*n + 1) * sizeof(int));
    *ja  = malloc(nnz * sizeof(int));
    *a   = malloc(nnz * sizeof(double));
    *b   = malloc(*n * sizeof(double));

    /* allocation réussite? */

    if (*ia == NULL || *ja == NULL || *a == NULL || *b == NULL ) {
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n\n");
        return 1;
    }

    /* partie principale : replissage de la matrice */

    int nnz_save = nnz;
    
    //passage ligne suiv(plaque complete)
    ind = 0;
    nnz = 0;
    for (iy = 0; iy < nx; iy++) { //iy ix indice sur grille hors bords  mais position (ix+1)*h
        for (ix = 0; ix < nx; ix++) {
            //exclu interieur et bord du trou
            if(! in_hole(ix,iy,y0,y1,x0,x1)){
                //marquer le début de la ligne suivante dans le tableau 'ia'
                (*ia)[ind] = nnz;
                (*b)[ind] = 0.0;
                //printf("ind : %d\n",ind);
               // printf("nnz : %d\n", nnz);
                //replissage de la ligne : voisin sud //ui-1
                // + verification si pas au dessus d'un bord
                
                if (check_sud(ix,iy,y0,y1,x0,x1,nx)){
                    (*a)[nnz] = -invh2;
                    (*ja)[nnz] = indice(ix,iy-1,y0,y1,x0,x1, nx);//ind - nx + skip_sud;
                    
                    //-nx car on regarde delui d'en bas(shema)
                    //+skip_sud car comme on a passe des points(trous)
                    // nx ramene trop loin en arrière
                    nnz++;
                }
                else{
                    (*b)[ind] += computeBound((ix+ 1)*h, (iy + 1 -1)*h); 
                }

                //replissage de la ligne : voisin ouest 
                //si pas a droite d'un bord
               
                if (check_west(ix,iy,y0,y1,x0,x1,nx)){
                    (*a)[nnz] = -invh2;
                    (*ja)[nnz] = ind - 1;
                    nnz++;
                }
                else{
                    (*b)[ind] += computeBound((ix+1 -1)*h, (iy + 1)*h);
                }

                // replissage de la ligne : élém. diagonal
                (*a)[nnz] = 4.0*invh2;
                (*ja)[nnz] = ind;
                
                nnz++;
                
                // replissage de la ligne : voisin est
                //si pas a gauche d'un bord
                
                if ( check_est(ix,iy,y0,y1,x0,x1,nx) ){
                    (*a)[nnz] = -invh2;
                    (*ja)[nnz] = ind + 1;
                    nnz++;
                }
                else{
                    (*b)[ind] += computeBound((ix+ 1 +1)*h, (iy+ 1)*h);
                }

                // replissage de la ligne : voisin nord
                //si pas en dessous d'un bord
                
                if ( check_nord(ix,iy,y0,y1,x0,x1,nx) ){
                        (*a)[nnz] = -invh2;
                        (*ja)[nnz] = indice(ix,iy+1,y0,y1,x0,x1, nx);
                        nnz++;
                }
                else{
                    (*b)[ind] += computeBound((ix+ 1)*h, (iy+ 1 +1)*h);
                }
                // numéro de l'équation
                ind += 1;
                
            }
        }
    }
    printf("%d \n", ind);
    if (nnz == nnz_save){
        (*ia)[ind] = nnz;
    }
    else{
        printf("Error nnz = %d!= nnz = %d\n", nnz, nnz_save);
        return 1;
    }
     for (int i=0; i<*n; i++){
            printf(" %lf", (*b)[i]);
        }
    


    /* retour de fonction habituel */
    return 0;
}
int ancienprob(int m, int *n, int **ia, int **ja, double **a, double **b)
/*
   But
   ===
   Génère la matrice n x n qui correspond à la disrétisation sur une grille 
   cartesienne regulière m x m de l'operateur de Laplace à deux dimensions
              
            d    d        d    d
         - == ( == u ) - == ( == u )        sur [0,3] x [0,3] \ [ 1, 5/2]×[ 3/2, 2]
           dx   dx       dy   dy

  avec la fonction u qui satisfait les conditions aux limites de Dirichlet
         
         u = 0  sur (0,y), (3,y), (x,0) et (x,3), pour tt x | 0<=x<=3 et tt y | 0<=y<=3
  
  La numérotation des inconnues est lexicographique, la direction x étant 
  parcourue avant celle de y. La matrice est retournée dans le format CRS
  qui est défini par le scalaire 'n' et les trois tableaux 'ia, 'ja' et 'a'.

  Arguments
  =========
  m (input)   - nombre de points par direction dans la grille 
  n  (output) - pointeur vers le nombre d'inconus dans le système
  ia (output) - pointeur vers le tableau 'ia' de la matrice A
  ja (output) - pointeur vers le tableau 'ja' de la matrice A
  a  (output) - pointeur vers le tableau 'a' de la matrice A
  nnz (output) - pointeur vers le nombre d'elements non nuls de A

*/
{

    int ix, iy, nx, ind = 0;
    double h;
    double invh2;

    nx = m - 2; /* noeuds de Dirichlet ne sont pas pris en compte */ 
    
    h = 3.0/(double)(m-1);
    invh2 = (double)((m-1)*(m-1)) / 9.0; /* h^-2 pour L=3 */

    int p = 1 + m/6;           //nb de points sur la largeur du trous
    int q = 1 + m/2;           // "" sur la longueur

    *n  = nx * nx - (p * q); 
    int nnz = 5 * nx * nx - 4 * nx ; //nombre d'element non nul pour membrane sans trous
    int nnz_save = nnz;
    //avec trous:
    int trous = (5*(p-2)*(q-2) + 4*2*(p-2) + 4*2*(q-2) +3*4*1 + 1*2*p + 1*2*q) +2*p +2*q;

    nnz -= trous;

    /* allocation des tableaux */

    *ia  = malloc((*n + 1) * sizeof(int));
    *ja  = malloc((nnz) * sizeof(int));
    *a   = malloc((nnz) * sizeof(double));
    *b = calloc(*n, sizeof(double));
    /* allocation réussite? */

    if (*ia == NULL || *ja == NULL || *a == NULL ) {
        printf("\n ERREUR : pas assez de mémoire pour générer la matrice\n\n");
        return 1;
    }

    /** partie principale : replissage de la matrice **/

    int skip_sud = 0;
    int skip_nord = 0;
    int ecart = 0;

    /*Calcul des constantes*/
    int y0 = (int)ceil(COORD_Y0/h) - 1;    //ceil arrondis pour ne pas avoir par exemple (int)6.99 = 6
    int y1 = (int)ceil(COORD_Y1/h) - 1;
    int x0 = (int)ceil(COORD_X0/h) - 1;
    int x1 = (int)ceil(COORD_X1/h) - 1;


    /*calcul du nb de point sur le bord le plus long*/

    for (ix = 0; ix < nx; ix++){
        if ( ix >= x0 && ix<= x1 ){ //nombre de points sur un bord interne(du trou)
                ecart += 1;
        }
    }

    (nnz) = 0;
    for (iy = 0; iy < nx; iy++){          //passage ligne suiv
        for (ix = 0; ix < nx; ix++){      //passage colonne suiv

            //initialisation du retard a ajouter à cause du bord (depend d'ou on se trouve sur la plaque)
            if ( iy >= y0-1 && ix >= x1){
                skip_nord = ecart;
            }
            if ( iy >= y1 && ix >= x1){
                skip_nord = 0;
            }
            if ( iy >= y0 && ix >= x1){
                skip_sud = ecart;
            }
            if ( iy > y1 && ix >= x1){
                skip_sud = 0;
            }
            
            
            if(iy > y1  || iy < y0  || ix < x0 || ix > x1){   //exclu interieur et bord du trou
            
                /* marquer le début de la ligne suivante dans le tableau 'ia' */
                (*ia)[ind] = (nnz);// val in = 0
            
                /* replissage de la ligne : voisin sud */ //ui-1
                if (iy > 0 && ( iy != y1+1  || ix < x0 || ix > x1) ){        //si pas au dessus d'un bord
                    (*a)[(nnz)] = -invh2;
                    (*ja)[(nnz)] = ind - nx + skip_sud;
                    /*
                    -nx car on regarde delui d'en bas(shema)
                    +skip_sud car comme on a passe des points(trous) nx ramene trop loin en arrière 
                    */
                    (nnz)++;
                }
                else{
                    (*b)[ind] += computeBound(ix*h, (iy-1)*h); 
                }

                /* replissage de la ligne : voisin ouest */
                if (ix > 0 && ( ix != x1 + 1 || iy > y1 || iy < y0 )){     //si pas a droite d'un bord
                    (*a)[(nnz)] = -invh2;
                    (*ja)[(nnz)] = ind - 1;
                    (nnz)++;
                }
                else{
                    (*b)[ind] += computeBound((ix -1)*h, iy*h);
                }

                /* replissage de la ligne : élém. diagonal */
                (*a)[(nnz)] = 4.0*invh2;
                (*ja)[(nnz)] = ind;
                if (nnz == nnz_save-1){
                    printf("last nnz\n");
                }
                (nnz)++;

                /* replissage de la ligne : voisin est */
                if (ix < nx - 1 && ( iy > y1 || iy < y0 || ix != x0 - 1 ) ){ //si pas a gauche d'un bord
                    (*a)[(nnz)] = -invh2;
                    (*ja)[(nnz)] = ind + 1;
                    (nnz)++;
                    
                }
                else{
                    (*b)[ind] += computeBound((ix+1)*h, iy*h);
                }

                /* replissage de la ligne : voisin nord */
                if (iy < nx - 1 && ( iy < y0 - 1 || iy > y1 || ix < x0 || ix > x1 ) ){  //si pas en dessous d'un bord
                        (*a)[(nnz)] = -invh2;
                        (*ja)[(nnz)] = ind + nx - skip_nord;
                        (nnz)++;
                }
                else{
                    (*b)[ind] += computeBound(ix*h, (iy+1)*h);
                }
                
                
                /* numéro de l'équation */
                ind += 1;
            }
        }
    }
    /* dernier élément du tableau 'ia' */
    if (nnz == nnz_save){
        (*ia)[ind] = nnz;
    }
    else{
        printf("Error nnz != nnz\n");
        return 1;
    }
    
    printf("a0 %lf\n", (*a)[nnz-1]);
    printf("ja0 %d\n", ja[0][nnz-1]);
    printf("ia0 %d\n", ia[0][*n]);
    printf("b0 %lf\n", b[0][*n-1]);
    /* retour de fonction habituel */
    return 0;
}
