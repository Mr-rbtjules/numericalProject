#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "proto.h"



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
                    	b[ind] += bound; 
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
                    	b[ind] += bound;
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
						b[ind] += bound;
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
						b[ind] += bound;
					}
                }
                //scaling
                if (b[ind] != 0){
                    b[ind] *= invh2;
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



void computeParamLevel(int m, double *h, double *invh2, int *x0,
						int *x1, int *y0, int *y1, int *nx,
						int *ny, int *n, int *nnz){
    
    int perUnit = (m-1)/(globVal.domain[1] - globVal.domain[0]);
    
    *x0 = ((globVal.domain[4] - globVal.domain[0]) * perUnit) - 1;
    *x1 = ((globVal.domain[5] - globVal.domain[0]) * perUnit) - 1;
    *y0 = ((globVal.domain[6] - globVal.domain[2]) * perUnit) - 1;
    *y1 = ((globVal.domain[7] - globVal.domain[2]) * perUnit) - 1;
    *h = (double)1/perUnit;
    *invh2 = (double)perUnit*perUnit;//1.0/((*h)*(*h));
    *nx = m-2; //par def
    *ny = ((globVal.domain[3] - globVal.domain[2])*perUnit) -1;
    getNnz(*nx, perUnit, *x0, *x1, *y0,* y1, n, nnz);
    printf("perunit %d h : %e invh2 : %e\n", perUnit, *h, *invh2);
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
int check_nord(int ix, int iy, int y0, int y1, int x0, int x1, int ny){
    
	if (iy + 1 < ny && ! in_hole(ix, iy +1, y0,y1,x0,x1)){
		return 1;
	}
	else{
		return 0;
	}
}

int check_sud(int ix, int iy, int y0, int y1, int x0, int x1, int ny){
    
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

int check_nw(int ixp, int iyp, int y0p, int y1p, 
             int x0p, int x1p, int nxp, int nyp){
    
	int ixnw = ixp - 1;//ind point nw sur la grille prolong 
    int iynw = iyp + 1;
    //pas besoin de check si sur bord droit car pas possible
    if (ixnw >= 0 && iynw < nyp && ! in_hole(ixnw, iynw, y0p,y1p,x0p,x1p)){
		return 1;
	}
	else{
		return 0;
	}
}

int check_ne(int ixp, int iyp, int y0p, int y1p, 
             int x0p, int x1p, int nxp, int nyp){
    
	int ixne =  ixp + 1;//ind point sw sur la grille prolong 
    int iyne = iyp + 1;
    if (ixne < nxp && iyne < nyp && ! in_hole(ixne, iyne, y0p,y1p,x0p,x1p)){
		return 1;
	}
	else{
		return 0;
	}
}

int check_sw(int ixp, int iyp, int y0p, int y1p, 
             int x0p, int x1p, int nxp, int nyp){
    //return 1/true lorsque point de la grille prolongé fait partie des variable
	int ixsw =  ixp - 1;//ind point sw sur la grille prolong 
    int iysw = iyp - 1;
    if (ixsw >= 0 && iysw >=0 && ! in_hole(ixsw, iysw, y0p,y1p,x0p,x1p)){
		return 1;
	}
	else{
		return 0;
	}
}

int check_se(int ixp, int iyp, int y0p, int y1p, 
             int x0p, int x1p, int nxp, int nyp){
    
	int ixse =  ixp+1;//ind point nw sur la grille prolong 
    int iyse = iyp - 1;
    if (ixse < nxp && iyse >= 0 && ! in_hole(ixse, iyse, y0p,y1p,x0p,x1p)){
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

