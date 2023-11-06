#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "proto.h"

#define BOUND1(x, y) (exp(sqrt((x)*(x) + (y)*(y))))
#define BOUND2(x, y) (sin(sqrt((x)*(x) + (y)*(y))))


//program parameters
//type of boundary conditions
#define BOUND BOUND1
//type of domain
#define DOMAIN "[0, 6] × [0, 10]  [2, 4] × [2, 5]"
//discretisation 
#define M 7
//print parameters
#define EXPLICIT 1

//start global variable (m and domain)
globVal_s globVal = { 0, NULL}; 

/// PROBLEM DISCRETISATION ///
// a tester
void getIndexHole(int *domain, int mc, int *x0, int *x1, int *y0, int *y1){
    //on part du principe que le m est correct
    //perUnit = nombre de sous case dans 1 case de geo de base
    int perUnit = (mc-1)/(domain[1] - domain[0]);
    //(combien de case ds la geo de base) * nb de sous case par case ds la geo de base
    // - 1 pour index 0 
    *x0 = ((domain[4] - domain[0]) * perUnit) - 1;
    *x1 = ((domain[5] - domain[0]) * perUnit) - 1;
    *y0 = ((domain[6] - domain[2]) * perUnit) - 1;
    *y1 = ((domain[7] - domain[2]) * perUnit) - 1;
}

void computeParamLevel(int mc, double *hl, double *invh2l, int *x0l,
						int *x1l, int *y0l, int *y1l, int *nxl,
						int *nl, int *nnzl){
    
    int perUnit = (mc-1)/(globVal.domain[1] - globVal.domain[0]);
    *x0l = ((globVal.domain[4] - globVal.domain[0]) * perUnit) - 1;
    *x1l = ((globVal.domain[5] - globVal.domain[0]) * perUnit) - 1;
    *y0l = ((globVal.domain[6] - globVal.domain[2]) * perUnit) - 1;
    *y1l = ((globVal.domain[7] - globVal.domain[2]) * perUnit) - 1;
    *hl = (double)1/perUnit;
    *invh2l = 1.0/((*hl)*(*hl));
    *nxl = mc-2;
    getNnz(*nxl, perUnit, *x0l, *x1l, *y0l,* y1l, nl, nnzl);
    if (EXPLICIT){
        printf("\n ProbMg hl = %lf ", *hl);
        printf("nl = %d nnzl = %d nxl = %d \n", *nl, *nnzl, *nxl);
        printf("x0l = %d x1l = %d y0l = %d y1l = %d \n", *x0l, *x1l, *y0l, *y1l);
    }
}

void getNnz(int nx, int perUnit, int x0, int x1, int y0, int y1, int *n, int *nnz){

    
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






int probMg(int mc, int *nl, 
		   int *ial, int *jal, double *al, double *bl){
    
    //fait pour m/2
    //mémoire deja allouée
    double hl, invh2l;
    int x0l,x1l,y0l,y1l, nxl, nnzl;

    
    computeParamLevel(mc, &hl,&invh2l,&x0l,&x1l, &y0l, &y1l, &nxl, nl, &nnzl);
    if (EXPLICIT){

    	printf("\n ProbMg hl = %lf ", hl);
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
					if (mc == globVal.m){
						double x = (ixl+ 1)*hl;
						double y = (iyl + 1 -1)*hl;
						double bound = computeBound(x, y);
                    	bl[ind] += bound * invh2l; 
					}
                }

                //replissage de la ligne : voisin ouest 
                //si pas a droite d'un bord
               
                if (check_west(ixl,iyl,y0l,y1l,x0l,x1l,nxl)){
                    al[nnzl] = -invh2l;
                    jal[nnzl] = ind - 1;
                    nnzl++;
                }
                else{
					if (mc == globVal.m){
						double x = (ixl + 1 - 1)*hl;
						double y = (iyl + 1)*hl;
						double bound = computeBound( x, y);
                    	bl[ind] += bound * invh2l;
					}
					
                }
                

                // replissage de la ligne : élém. diagonal
                al[nnzl] = 4.0*invh2l;
                jal[nnzl] = ind;
                
                nnzl++;
                
                // replissage de la ligne : voisin est
                //si pas a gauche d'un bord
                
                if ( check_est(ixl,iyl,y0l,y1l,x0l,x1l,nxl)){
                    al[nnzl] = -invh2l;
                    jal[nnzl] = ind + 1;
                    nnzl++;
                }
                else{
					if (mc == globVal.m){
						double x = (ixl + 1 +1)*hl;
						double y = (iyl + 1)*hl;
						double bound = computeBound(x, y);	
						bl[ind] += bound * invh2l;
					}
                }

                // replissage de la ligne : voisin nord
                //si pas en dessous d'un bord
                
                if ( check_nord(ixl,iyl,y0l,y1l,x0l,x1l,nxl) ){
                        al[nnzl] = -invh2l;
                        jal[nnzl] = indice(ixl,iyl+1,y0l,y1l,x0l,x1l, nxl);
                        nnzl++;
                }
                else{
					if (mc == globVal.m){
						double x = (ixl + 1)*hl;
						double y = (iyl + 1 +1)*hl;
						double bound = computeBound(x, y);
						bl[ind] += bound * invh2l;
					}

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


int in_hole(int ix, int iy, int y0, int y1, int x0, int x1){
	
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





/// GLOBAL VARIABLES ///
void initGlobVal(){
    int count;
    extractDomain(DOMAIN, &count, &(globVal.domain));
    globVal.m = correctM(globVal.domain, M);
}

void freeGlobVal(){
    free(globVal.domain);
}

double computeBound(double x, double y){
    //compute the value of a point on the bound of the domain
	//with a function BOUND defined in the macros
    return BOUND(x,y);
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

//a tester
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


