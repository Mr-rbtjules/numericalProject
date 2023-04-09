#include <stdlib.h>
#include <stdio.h>
#include "proto.h"
#include <math.h>


#define SCALE_FACT 0.5




/*fonctionnnement des level
-celui tt en haut avec le h le plus petit = level 0 puis coarse 1 coarse coarse 2 etc
- qd passe level en arg c'est tj le level auquel on est pas celui vers lequel on va
*/

int restrictR(int level, double *rp, double **rc, int m, int *nc){
    //level = 0 = ou on est ->1
    // lvel -> level +1
    

    //level -1 (prolonge donc monte dans la pyramyde)
    double hp, invh2p;
    int x0p,x1p,y0p,y1p, nxp, np, nnzp;
    computeParamLevel(m, level, &hp,&invh2p,&y0p,&y1p,&x0p,&x1p,&nxp, &np, &nnzp);
    
    printf("\n restrict level : hp = %lf np = %d nnzp = %d nxp = %d \n", hp, np, nnzp, nxp);
    printf("x0p = %d x1p = %d y0p = %d y1p = %d \n", x0p, x1p, y0p, y1p);

    
   

    //level ou on va =level+1
    double hc, invh2c;
    int x0c,x1c,y0c,y1c, nxc, nnzc;
    computeParamLevel(m, level+1, &hc,&invh2c,&y0c,&y1c,&x0c,&x1c,&nxc, nc, &nnzc);
    
    printf("\n restrict level+1 : hc = %lf nc = %d nnzc = %d nxc = %d \n", hc, *nc, nnzc, nxc);
    printf("x0c = %d x1c = %d y0c = %d y1c = %d \n", x0c, x1c, y0c, y1c);


    if (*rc == NULL){
        *rc = malloc(*nc *sizeof(double));
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
                    (*rc)[*nc] = 0;
                    //marquer le début de la ligne suivante dans le tableau 'ia'            
                    //replissage de la ligne : voisin sud //ui-1
                    // + verification si pas au dessus d'un bord
                    
                    if (check_sud(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){
                        int ind = indice(ixp,iyp-1,y0p,y1p,x0p,x1p, nxp);
                        (*rc)[*nc] += 0.25 * rp[ind] * SCALE_FACT;
                    }    

                    else{
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT; 
                        (*rc)[*nc] += 0.25 * computeBound((ixp+1)*hp,(iyp+1 -1)*hp) * SCALE_FACT;//si fct qui calcule chaque fois :juste coord sinon
                    }

                    //replissage de la ligne : voisin ouest 
                    //si pas a droite d'un bord
                    if (check_west(ixp,iyp,y0p,y1p,x0p,x1p,nxp)){
                        int ind = np -1;
                        (*rc)[*nc] += 0.25 * rp[ind] * SCALE_FACT;
                    }
                    else{
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT;
                        (*rc)[*nc] += 0.25 * computeBound((ixp+1-1)*hp,(iyp+1)*hp) * SCALE_FACT;
                    }

                    int ind = np;
                    (*rc)[*nc] += rp[ind] * SCALE_FACT;
                   
                    // replissage de la ligne : élém. diagonal
                    

                    // replissage de la ligne : voisin est
                    //si pas a gauche d'un bord
                    if (check_est(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){
                        int ind = np +1;
                        (*rc)[*nc] += 0.25 * rp[ind] * SCALE_FACT;
                         
                    }
                    else{
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT;
                        (*rc)[*nc] += 0.25 * computeBound((ixp+1+1)*hp ,(iyp+1)*hp) * SCALE_FACT;
                        
                    }

                    // replissage de la ligne : voisin nord
                    //si pas en dessous d'un bord
                    if ( check_nord(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){

                        (*rc)[*nc] += 0.25 * rp[indice(ixp,iyp+1,y0p,y1p,x0p,x1p, nxp)] * SCALE_FACT;
                    }
                    else{
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT;
                        (*rc)[*nc] += 0.25 * computeBound((ixp+1)*hp, (iyp+1+1)*hp) * SCALE_FACT;
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


int prolongR(int level, double **up, double *uc, int m, int *np){ //ici m de u pas de uc !
    //level 1 = on ou on est -> on va vers 0
    // level -> level -1
    
    //level level
    double hc, invh2c;
    int x0c,x1c,y0c,y1c, nxc, nc, nnzc;

    computeParamLevel(m, level, &hc,&invh2c,&y0c,&y1c,&x0c,&x1c,&nxc, &nc, &nnzc);


    //level -1 car on remonte dans la pyr inv
    double hp, invh2p;
    int x0p,x1p,y0p,y1p, nxp, nnzp;

    computeParamLevel(m, level-1, &hp,&invh2p,&y0p,&y1p,&x0p,&x1p,&nxp, np, &nnzp);

    printf("\n prolong : hp = %lf np = %d nnzp = %d nxp = %d \n", hp, *np, nnzp, nxp);
    printf(" x0p = %d x1p = %d y0p = %d y1p = %d \n", x0p, x1p, y0p, y1p);
    printf("\nand hc = %lf nc = %d nnzc = %d nxc = %d \n", hc, nc, nnzc, nxc);
    printf("x0c = %d x1c = %d y0c = %d y1c = %d \n", x0c, x1c, y0c, y1c);


    if (*up == NULL){
        *up = malloc(*np * sizeof(double));
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
                
                (*up)[*np] = 0;

				//impair impair -> 1/4 somme des 4 autour
                
                if (iyp % 2 == 0 && ixp % 2 == 0){ 
                    //check si bord est bien la ou pense etre pour uc (pas faire une recherche alors que bord)
                    //somme coin gauche bas 
                    //check si point prol pas ds le trou
                    
                    if (check_sw(ixp,iyp,y0p,y1p,x0p,x1p,nxp)){ // question est ce que prolong peut etre domaine et coarse ds un board ? non ca depend que de p
                        
						ind = indice((ixp/2) - 1, (iyp/2)-1, y0c,y1c,x0c,x1c, nxc);//juste /2 -1 car point prol au milieu des 4 tjrs pair
                        (*up)[*np] += 0.25 * uc[ind]; 
                    }
                    else{
                        double bound = computeBound((ixp+1-1)*hp,(iyp+1-1)*hp);
                        
                        (*up)[*np] += 0.25 * bound;
                       
                    }
                    //coin droit bas
                    if (check_se(ixp,iyp,y0p,y1p,x0p,x1p,nxp)){ //cond droit
                        ind = indice((ixp/2) ,(iyp/2)-1, y0c,y1c,x0c,x1c, nxc);
                        (*up)[*np] += 0.25 * uc[ind];
                        
                    }
                    else{
                        (*up)[*np] += 0.25 * computeBound((ixp+1+1)*hp,(iyp+1-1)*hp);
                    }
                    //coin haut gauche
                    if (check_nw(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){ //cond gauche haut
                        ind = indice((ixp/2 - 1),(iyp/2), y0c,y1c,x0c,x1c, nxc);
                        (*up)[*np] += 0.25 * uc[ind];
                    }
                    else{
                        double bound = computeBound((ixp+1-1)*hp,(iyp + 1+1)*hp);
                        (*up)[*np] += 0.25 * bound;
                        
                    }
                    //coin droit haut
                    
                    if (check_ne(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){ //coin droit
                        ind = indice((ixp/2),(iyp/2), y0c,y1c,x0c,x1c, nxc);
                        (*up)[*np] += 0.25 * uc[ind];

                    }
                    else{
                        (*up)[*np] += 0.25 * computeBound((ixp+1+1)*hp,(iyp+1+1)*hp);
                    }
                }

				//impair pair => somme haut + bas
				else if (iyp % 2 == 0 && ixp % 2 == 1){
                    
                    //somme bas
                    if (check_sud(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){
                       
                        ind = indice((ixp/2),(iyp/2)-1, y0c,y1c,x0c,x1c, nxc);
                        (*up)[*np] += 0.5 * uc[ind];
                    }
                    else{
                        (*up)[*np] += 0.5 * computeBound((ixp+1)*hp,(iyp+1-1)*hp);
                    }
                    //somme haut
                    if(check_nord(ixp,iyp,y0p,y1p,x0p,x1p,nxp)){
                        ind = indice((ixp/2),(iyp/2), y0c,y1c,x0c,x1c, nxc);
                        (*up)[*np] += 0.5 * uc[ind]; //pas de nx/2+1 car point uc[nc] deja ligne du haut
                    }
                    else{
                        (*up)[*np] += 0.5 * computeBound((ixp+1)*hp,(iyp+1+1)*hp);
                    }
                }
                //pair impair 1/2 somme gauche droite
                else if (iyp % 2 == 1 && ixp % 2 == 0){

                    //somme gauche
                    //si pas a droite d'un bord
                    if (check_west(ixp,iyp,y0p,y1p,x0p,x1p,nxp)){
                        (*up)[*np] += 0.5 * uc[nc - 1];
                    }
                    else{
                        (*up)[*np] += 0.5 * computeBound((ixp+1-1)*hp,(iyp+1)*hp);
                    }

                    //somme droit
                    //si pas a gauche d'un bord
                    if (check_est(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){
                        (*up)[*np] += 0.5 * uc[nc]; //nc car deja +1 car type point precedent + 1
                    }
                    else{
                        (*up)[*np] += 0.5 * computeBound((ixp+1+1)*hp,(iyp+1)*hp);
                    }  
                }
                //ligne impair  & col impair => elem identique
                else if (iyp % 2 == 1 && ixp % 2 == 1){  //pour 1 meme ligne alterne entre type 1 et 2, commence par 1 finis par 1

                    (*up)[*np] = uc[nc];
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


int probMg(int m, int level, int *nl, int *ial, int *jal, double *al, double *bl){
    //ici level == celui dont on veut calculer A et b
    //mémoire deja allouée

    
    double hl, invh2l;
    int x0l,x1l,y0l,y1l, nxl, nnzl;

    computeParamLevel(m, level, &hl,&invh2l,&y0l,&y1l,&x0l,&x1l,&nxl, nl, &nnzl);
    
    printf("\n ProbMg : hl = %lf nl = %d nnzl = %d nxl = %d \n", hl, *nl, nnzl, nxl);
    printf("x0l = %d x1l = %d y0l = %d y1l = %d \n", x0l, x1l, y0l, y1l);

	
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


int allocGridLevel(int m, int level, int *nl, int **ial,
                     int **jal, double **al, double **bl){

    double hl, invh2l;
    int x0l,x1l,y0l,y1l, nxl, nnzl;
    computeParamLevel(m, level, &hl,&invh2l,&y0l,&y1l,&x0l,&x1l,&nxl, nl, &nnzl);
    
    *ial = malloc((*nl + 1) * sizeof(int));
    *jal = malloc(nnzl * sizeof(int));
    *al = malloc(nnzl * sizeof(double));
    *bl = malloc(*nl * sizeof(double));
    
    if (*bl == NULL || *ial == NULL || *jal == NULL || *al == NULL){
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n\n");
        return 1;
    }

    printf("\n Allow : hl = %lf nl = %d nnzl = %d nxl = %d \n", hl, *nl, nnzl, nxl);
    printf("x0l = %d x1l = %d y0l = %d y1l = %d \n", x0l, x1l, y0l, y1l);


    return 0;
}






/*
int returnAc(int ix, int iy, int m, int level, double h, int *ia, int *ja, double *a){

    //return Ac on the base of A

	//level 1 : 2grid method
	int multi = level -1;
	//valeurs de u
	double h = 3.0/(double)(m-1);
    double invh2 = 1.0/(h*h);
    int x1 = ((int)(COORD_X1 * (m-1)) /3)  -1; 
    int x0 = (((int)(COORD_X0 * (m-1)) + ((3 - ((int)(COORD_X0*(m-1))%3))%3))/3)  -1;
    int y1 = ((int)(COORD_Y1 * (m-1)) /3)  -1;
    int y0 = (((int)(COORD_Y0 * (m-1)) + ((3 - ((int)(COORD_Y0*(m-1))%3))%3))/3)  -1;
    
    int nx = m-2;

    int p = y1 - y0 + 1;
    int q = x1 - x0 + 1;
    int n = nx * nx - (p * q);
   
	double hc = 2*h;
	double invh2c = 1.0/(hc*hc);  // ou 1/(4* h**2) = 0.25 * invh2

	int x0c = (x0 / 2*multi); // va arrondir au point grille coarse a droite 
	int x1c = ((x1+1)/2*multi) - 1; // permet si x1 pair on retire 1 
	int y0c = (y0/2*multi); // arrondu coarse au dessus (permet de pas ajouter des points dans le trou)
    int y1c = ((y1+1)/2*multi) - 1;
	int pc = y1c - y0c + 1;
    int qc = x1c - x0c + 1;

	int nxc = nx/2; // nb de points coars sur un ligne pas bord
    int nc = nxc * nxc - (pc * qc);
	int nnzc = 5 * nxc * nxc - 4 * nxc ; 
    //nb de points concernés dans le trou:(compliqué a comprendre sans shema) (marche que si aumoins 3 points sur la largeur)
    int trousc = (5 * (pc-2) * (qc-2) + 4 * 2 * (pc-2) + 4 * 2 * (qc-2) 
                + 3 * 4 * 1 + 1 * 2 * pc + 1 * 2 * qc) + 2 * pc + 2 * qc;
    nnzc -= trousc;

	int nnz_save = nnz;
    
    //passage ligne suiv(plaque complete)
    ind = 0;
    nnzc = 0;
    for (iyc = 0; iyc < nxc; iyc++) { //iy ix indice sur grille hors bords  mais position (ix+1)*h
        for (ixc = 0; ixc < nxc; ixc++) {
            //exclu interieur et bord du trou
            if(! in_hole(ixc,iyc,y0c,y1c,x0c,x1c)){
                //marquer le début de la ligne suivante dans le tableau 'ia'
                (*iac)[ind] = nnzc;
               
                if (check_sud(ixc,iyc,y0c,y1c,x0c,x1c,nxc)){
                    (*ac)[nnzc] = -invh2c;
                    (*jac)[nnzc] = indice(ixc,iyc-1,y0c,y1c,x0c,x1c, nxc);//ind - nx + skip_sud;
                    
                    //-nx car on regarde delui d'en bas(shema)
                    //+skip_sud car comme on a passe des points(trous)
                    // nx ramene trop loin en arrière
                    nnzc++;
                }

                //replissage de la ligne : voisin ouest 
                //si pas a droite d'un bord
               
                if (check_west(ixc,iyc,y0c,y1c,x0c,x1c,nxc)){
                    (*ac)[nnzc] = -invh2c;
                    (*jac)[nnzc] = ind - 1;
                    nnzc++;
                }
                

                // replissage de la ligne : élém. diagonal
                (*ac)[nnzc] = 4.0*invh2c;
                (*jac)[nnzc] = ind;
                
                nnzc++;
                
                // replissage de la ligne : voisin est
                //si pas a gauche d'un bord
                
                if ( check_est(ixc,iyc,y0c,y1c,x0c,x1c,nxc) ){
                    (*ac)[nnzc] = -invh2c;
                    (*jac)[nnzc] = ind + 1;
                    nnzc++;
                }

                // replissage de la ligne : voisin nord
                //si pas en dessous d'un bord
                
                if ( check_nord(ixc,iyc,y0c,y1c,x0c,x1c,nxc) ){
                        (*ac)[nnzc] = -invh2c;
                        (*jac)[nnzc] = indice(ixc,iyc+1,y0c,y1c,x0c,x1c, nxc);
                        nnzc++;
                }
                // numéro de l'équation
                ind += 1;
                
            }
        }
    }

	return 0;
}	
*/
