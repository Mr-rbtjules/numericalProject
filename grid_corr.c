#include <stdlib.h>
#include "proto.h"


#define SCALE_FACT 0.5
#define COORD_X0 1.0 //col bord gauche
#define COORD_X1 2.5 // col bord droit
#define COORD_Y0 1.5 //ligne bord bas
#define COORD_Y1 2.0 //ligne bord haut



/*fonctionnnement des level
celui tt en haut avec le h le plus petit = level 0 puis coarse 1 coarse coarse 2 etc*/

int restrictR(double **r, double **rc, int m, int *n){
//ajouter mode multi

	int ix, iy;

	//valeurs de r
	double h = 3.0/(double)(m-1);
    double invh2 = 1.0/(h*h);
    int x1 = ((int)(COORD_X1 * (m-1)) /3)  -1; 
    int x0 = (((int)(COORD_X0 * (m-1)) + ((3 - ((int)(COORD_X0*(m-1))%3))%3))/3)  -1;
    int y1 = ((int)(COORD_Y1 * (m-1)) /3)  -1;
    int y0 = (((int)(COORD_Y0 * (m-1)) + ((3 - ((int)(COORD_Y0*(m-1))%3))%3))/3)  -1;
    
    int nx = m-2;

    int p = y1 - y0 + 1;
    int q = x1 - x0 + 1;
    *n = nx * nx - (p * q);
    
	//valeurs de rc
	int x0c = (x0 / 2); // va arrondir au point grille coarse a droite 
	int x1c = ((x1+1)/2) - 1; // permet si x1 par on retire 1 et pair on ajoute 1
	int y0c = (y0/2); // arrondu coarse au dessus (permet de pas ajouter des points dans le trou)
    int y1c = ((y1+1)/2) - 1;
	int pc = y1c - y0c + 1;
    int qc = x1c - x0c + 1;

	int nxc = nx/2; // nb de points coars sur un ligne pas bord
    int nc = nxc * nxc - (pc * qc);

    if (*rc == NULL){
        *rc = malloc(nc * sizeof(double));
    }
	

	int skip_sud = 0;
    int skip_nord = 0;

	int nr = 0;
	int nrc = 0;


	for (iy = 0; iy < nx; iy++){
        //passage colonne suiv
        for (ix = 0; ix < nx; ix++){      

            //exclu interieur et bord du trou
            if( ! in_hole(ix,iy,y0,y1,x0,x1)){
                
                //garde que ligne ind impair  & col impair 
                if (iy % 2 == 1 && ix % 2 == 1){  

                    //marquer le début de la ligne suivante dans le tableau 'ia'            
                    //replissage de la ligne : voisin sud //ui-1
                    // + verification si pas au dessus d'un bord
                    if (check_sud(ix,iy,y0,y1,x0,x1,nx) ){
                        
                        (*rc)[nrc] += 0.25 * (*r)[indice(ix,iy-1,y0,y1,x0,x1, nx)] * SCALE_FACT;
                    }
                    else{
                        (*rc)[nrc] += 0.25 * computeBound((ix+1)*h,(iy+1 -1)*h) * SCALE_FACT;//si fct qui calcule chaque fois :juste coord sinon
                    }

                    //replissage de la ligne : voisin ouest 
                    //si pas a droite d'un bord
                    if (check_west(ix,iy,y0,y1,x0,x1,nx)){
                        (*rc)[nrc] += 0.25 * (*r)[nr - 1] * SCALE_FACT;
                    }
                    else{
                        (*rc)[nrc] += 0.25 * computeBound((ix+1-1)*h,(iy+1)*h) * SCALE_FACT;
                    }


                    (*rc)[nrc] += (*r)[nr] * SCALE_FACT;
                    // replissage de la ligne : élém. diagonal
                    

                    // replissage de la ligne : voisin est
                    //si pas a gauche d'un bord
                    if (check_est(ix,iy,y0,y1,x0,x1,nx) ){
                        (*rc)[nrc] += 0.25 * (*r)[nr + 1] * SCALE_FACT;
                    }
                    else{
                        (*rc)[nrc] += 0.25 * computeBound((ix+1+1)*h ,(iy+1)*h) * SCALE_FACT;
                    }

                    // replissage de la ligne : voisin nord
                    //si pas en dessous d'un bord
                    if ( check_nord(ix,iy,y0,y1,x0,x1,nx) ){

                        (*rc)[nrc] += 0.25 * (*r)[indice(ix,iy+1,y0,y1,x0,x1, nx)] * SCALE_FACT;
                    }
                    else{
                        (*rc)[nrc] += 0.25 * computeBound((ix+1)*h, (iy+1+1)*h) * SCALE_FACT;
                    }
                    // numéro de l'équation
                    nrc += 1;
                
                }
                nr += 1;
            }
        }
    }

	if (nr != *n){
		printf("nr != n\n");
		return 1;
	}
	if (nrc != nc){
		printf("nrc != nc\n");
		return 1;
	}


	return 0;
}


int prolongR(double **u, double **uc, int m){ //ici m de u pas de uc !
//ajouter mode multi

	int ixc, iyc;

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

	int x0c = (x0 / 2); // va arrondir au point grille coarse a droite 
	int x1c = ((x1+1)/2) - 1; // permet si x1 par on retire 1 et pair on ajoute 1
	int y0c = (y0/2); // arrondu coarse au dessus (permet de pas ajouter des points dans le trou)
    int y1c = ((y1+1)/2) - 1;
	int pc = y1c - y0c + 1;
    int qc = x1c - x0c + 1;

	int nxc = nx/2; // nb de points coars sur un ligne pas bord
    int nc = nxc * nxc - (pc * qc);

	*u = malloc(n * sizeof(double));

	
	int nu = 0;
	int nuc = 0;
	
	for (iy = 0; iy < nx; iy++){
        //passage colonne suiv
        for (ix = 0; ix < nx; ix++){      

            //exclu interieur et bord du trou
            if(! in_hole(ix,iy,y0,y1,x0,x1)){
                
				//impair impair -> 1/4 somme des 4 autour
                if (iy % 2 == 0 && ix % 2 == 0){
                    //check si bord est bien la ou pense etre pour uc (pas faire une recherche alors que bord)
                    //somme coin gauche bas 
                    //check si pas ds le trou
                    if (check_sw(ix,iy,y0,y1,x0,x1,nx)){
                        
						int ind = indice((ix/2)-1,(iy/2-1), y0c,y1c,x0c,x1c, nxc);//indicie point sw coarse
                        (*u)[n] += 0.25 * (*uc)[ind]; 
                    }
                    else{
                        (*u)[n] += 0.25 * computeBound((ix+1-1)*h,(iy+1-1)*h);
                    }
                    //coin droit bas
                    if (check_se(ix,iy,y0,y1,x0,x1,nx)){ //cond droit
                        int ind = indice((ix/2),(iy/2-1), y0c,y1c,x0c,x1c, nxc);
                        (*u)[n] += 0.25 * (*uc)[ind];
                    }
                    else{
                        (*u)[n] += 0.25 * computeBound((ix+1+1)*h,(iy+1-1)*h);
                    }
                    //coin haut gauche
                    if (check_nw(ix,iy,y0,y1,x0,x1,nx) ){ //cond gauche haut
                        int ind = indice((ix/2 - 1),(iy/2), y0c,y1c,x0c,x1c, nxc);
                        (*u)[n] += 0.25 * (*uc)[ind];
                    }
                    else{
                        (*u)[n] += 0.25 * computeBound((ix+1-1)*h,(iy + 1+1)*h);
                    }
                    //coin droit haut
                    
                    if (check_ne(ix,iy,y0,y1,x0,x1,nx) ){ //cond droit
                        int ind = indice((ix/2),(iy/2), y0c,y1c,x0c,x1c, nxc);
                        (*u)[n] += 0.25 * (*uc)[ind];
                    }
                    else{
                        (*u)[n] += 0.25 * computeBound((ix+1+1)*h,(iy+1+1)*h);
                    }
                }

				//impair pair => somme haut + bas
				else if (iy % 2 == 0 && ix % 2 == 1){
                    
                    //somme bas
                    if (check_sud(ix,iy,y0,y1,x0,x1,nx) ){
                       
                        int ind = indice((ix/2),(iy/2)-1, y0c,y1c,x0c,x1c, nxc);
                        (*u)[n] += 0.5 * (*uc)[ind];
                    }
                    else{
                        (*u)[n] += 0.5 * computeBound((ix+1)*h,(iy+1-1)*h);
                    }
                    //somme haut
                    if(check_nord(ix,iy,y0,y1,x0,x1,nx)){
                        int ind = indice((ix/2),(iy/2), y0c,y1c,x0c,x1c, nxc);
                        (*u)[n] += 0.5 * (*uc)[ind]; //pas de nx/2+1 car point uc[nc] deja ligne du haut
                    }
                    else{
                        (*u)[n] += 0.5 * computeBound((ix+1)*h,(iy+1+1)*h);
                    }
                }
                //pair impair 1/2 somme gauche droite
                else if (iy % 2 == 1 && ix % 2 == 0){

                    //somme gauche
                    //si pas a droite d'un bord
                    if (check_west(ix,iy,y0,y1,x0,x1,nx)){
                        (*u)[n] += 0.5 * (*uc)[nc - 1];
                    }
                    else{
                        (*u)[n] += 0.5 * computeBound((ix+1-1)*h,(iy+1)*h);
                    }

                    //somme droit
                    //si pas a gauche d'un bord
                    if (check_est(ix,iy,y0,y1,x0,x1,nx) ){
                        (*u)[n] += 0.5 * (*uc)[nc]; //nc car deja +1 car type point precedent + 1
                    }
                    else{
                        (*u)[n] += 0.5 * computeBound((ix+1+1)*h,(iy+1)*h);
                    }  
                }
                //ligne impair  & col impair => elem identique
                else if (iy % 2 == 1 && ix % 2 == 1){  //pour 1 meme ligne alterne entre type 1 et 2, commence par 1 finis par 1

                    (*u)[n] = (*uc)[nc];
                    nc += 1;       //=> 2 choses, pr point type 2 droite = nc gauche == nc-1
                    // et aussi que pour type 3 et 4 nc rpz point au dessus tt a gauche
                }
                
                n += 1;
            }
        }
    }

	return 0;
}


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

int probMg(int m, int level, int *iac, int *jac, double *ac, double *b){


	//level 1 : 2grid method
	int ixc, iyc;
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

	double hc = h*exp(2, level);
	double invh2c = 3.0/(hc*hc); 

    int x0c = x0 / exp(2, level); // va arrondir au point grille coarse a droite 
    int x1c = ((x1+1)/exp(2, level)) - 1; // permet si x1 pair on retire 1 
    int y0c = (y0/exp(2, level)); // arrondu coarse au dessus (permet de pas ajouter des points dans le trou)
    int y1c = ((y1+1)/exp(2, i)) - 1;
    int pc = y1c - y0c + 1;
    int qc = x1c - x0c + 1;

    int nxc = nx/exp(2,level); // nb de points coars sur un ligne pas bord
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
                iac[ind] = nnzc;
               
                if (check_sud(ixc,iyc,y0c,y1c,x0c,x1c,nxc)){
                    ac[nnzc] = -invh2c;
                    jac[nnzc] = indice(ixc,iyc-1,y0c,y1c,x0c,x1c, nxc);//ind - nx + skip_sud;
                    
                    //-nx car on regarde delui d'en bas(shema)
                    //+skip_sud car comme on a passe des points(trous)
                    // nx ramene trop loin en arrière
                    nnzc++;
                }
                else{
                    if (level ==0){
                        b[ind] += computeBound((ix+ 1)*h, (iy + 1 -1)*h); 
                    }
                }

                //replissage de la ligne : voisin ouest 
                //si pas a droite d'un bord
               
                if (check_west(ixc,iyc,y0c,y1c,x0c,x1c,nxc)){
                    ac[nnzc] = -invh2c;
                    jac[nnzc] = ind - 1;
                    nnzc++;
                }
                else{
                    if (level ==0){
                        b[ind] += computeBound((ix + 1 - 1)*h, (iy + 1)*h);
                    }
                }
                

                // replissage de la ligne : élém. diagonal
                ac[nnzc] = 4.0*invh2c;
                jac[nnzc] = ind;
                
                nnzc++;
                
                // replissage de la ligne : voisin est
                //si pas a gauche d'un bord
                
                if ( check_est(ixc,iyc,y0c,y1c,x0c,x1c,nxc) ){
                    ac[nnzc] = -invh2c;
                    jac[nnzc] = ind + 1;
                    nnzc++;
                }
                else{
                    if (level ==0){
                        b[ind] += computeBound((ix + 1 +1)*h, (iy + 1)*h);
                    }
                }

                // replissage de la ligne : voisin nord
                //si pas en dessous d'un bord
                
                if ( check_nord(ixc,iyc,y0c,y1c,x0c,x1c,nxc) ){
                        ac[nnzc] = -invh2c;
                        jac[nnzc] = indice(ixc,iyc+1,y0c,y1c,x0c,x1c, nxc);
                        nnzc++;
                }
                else{
                    if (level ==0){
                        b[ind] += computeBound((ix + 1)*h, (iy + 1 +1)*h);
                    }
                }
                // numéro de l'équation
                ind += 1;
                
            }
        }
    }

	return 0;
}	


