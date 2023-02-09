#include <stdlib.h>
#include "proto.h"


#define SCALE_FACT 0.5
#define COORD_X0 1.0 //col bord gauche
#define COORD_X1 2.5 // col bord droit
#define COORD_Y0 1.5 //ligne bord bas
#define COORD_Y1 2.0 //ligne bord haut



int restrictR(double **r, double **rc, int m, int *n){


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
	int sizeRc;
	int mc = (m+1)/2;
	int nxc = mc-2;
	int x1c = ((int)(COORD_X1 * (mc-1)) /3)  -1; 
    int x0c = (((int)(COORD_X0 * (mc-1)) + ((3 - ((int)(COORD_X0*(mc-1))%3))%3))/3)  -1;
    int y1c = ((int)(COORD_Y1 * (mc-1)) /3)  -1;
    int y0c = (((int)(COORD_Y0 * (mc-1)) + ((3 - ((int)(COORD_Y0*(mc-1))%3))%3))/3)  -1;
    

    int pc = y1c - y0c + 1;
    int qc = x1c - x0c + 1;
    int nc = nxc * nxc - (pc * qc);

	*r = malloc(nc * sizeof(double));

	int skip_sud = 0;
    int skip_nord = 0;

	int nr = 0;
	int nrc = 0;


	for (iy = 0; iy < nx; iy++){
        //passage colonne suiv
        for (ix = 0; ix < nx; ix++){      

            //initialisation du retard a ajouter à cause du bord 
            //(depend d'ou on se trouve sur la plaque)
            if ( iy >= y0-1 && ix >= x1){
                skip_nord = q;
            }
            if ( iy >= y1 && ix >= x1){
                skip_nord = 0;
            }
            if ( iy >= y0 && ix >= x1){
                skip_sud = q;
            }
            if ( iy > y1 && ix >= x1){
                skip_sud = 0;
            }
            //exclu interieur et bord du trou
            if( ! in_hole((ix+1)*h, (iy+1)*h)){
                
                //garde que ligne ind impair  & col impair 
                if (iy % 2 == 1 && ix % 2 == 1){  

                    //marquer le début de la ligne suivante dans le tableau 'ia'            
                    //replissage de la ligne : voisin sud //ui-1
                    // + verification si pas au dessus d'un bord
                    if (check_sud(ix,iy,h,h) ){
                        
                        (*rc)[nrc] += 0.25 * (*r)[nr - nx + skip_sud] * SCALE_FACT;
                        //-nx car on regarde delui d'en bas(shema)
                        //+skip_sud car comme on a passe des points(trous)
                        // nx ramene trop loin en arrière
                    }
                    else{
                        (*rc)[nrc] += 0.25 * computeBound((ix+1)*h,(iy+1 -1)*h) * SCALE_FACT;//si fct qui calcule chaque fois :juste coord sinon
                    }

                    //replissage de la ligne : voisin ouest 
                    //si pas a droite d'un bord
                    if (check_west(ix,iy, h,h)){
                        (*rc)[nrc] += 0.25 * (*r)[nr - 1] * SCALE_FACT;
                    }
                    else{
                        (*rc)[nrc] += 0.25 * computeBound((ix+1-1)*h,(iy+1)*h) * SCALE_FACT;
                    }


                    (*rc)[nrc] += (*r)[nr] * SCALE_FACT;
                    // replissage de la ligne : élém. diagonal
                    

                    // replissage de la ligne : voisin est
                    //si pas a gauche d'un bord
                    if (check_est(ix,iy,h,h) ){
                        (*rc)[nrc] += 0.25 * (*r)[nr + 1] * SCALE_FACT;
                    }
                    else{
                        (*rc)[nrc] += 0.25 * computeBound((ix+1+1)*h ,(iy+1)*h) * SCALE_FACT;
                    }

                    // replissage de la ligne : voisin nord
                    //si pas en dessous d'un bord
                    if (check_nord(ix,iy,h,h) ){

                        (*rc)[nrc] += 0.25 * (*r)[nr + nx - skip_nord] * SCALE_FACT;
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


int prolongR(double **u, double **uc, int m){

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
    *n = nx * nx - (p * q);
    
	//valeurs de uc
	int sizeRc;
	int mc = (m+1)/2;
	double hc = 3.0/(mc-1);
	int nxc = mc-2;
	int x1c = ((int)(COORD_X1 * (mc-1)) /3)  -1; 
    int x0c = (((int)(COORD_X0 * (mc-1)) + ((3 - ((int)(COORD_X0*(mc-1))%3))%3))/3)  -1;
    int y1c = ((int)(COORD_Y1 * (mc-1)) /3)  -1;
    int y0c = (((int)(COORD_Y0 * (mc-1)) + ((3 - ((int)(COORD_Y0*(mc-1))%3))%3))/3)  -1;
    

    int pc = y1c - y0c + 1;
    int qc = x1c - x0c + 1;
    int nc = nxc * nxc - (pc * qc);

	*r = malloc(nc * sizeof(double));

	int skip_sud = 0;
    int skip_nord = 0;

	int nu = 0;
	int nuc = 0;

	for(ixc = 0; ixc < nxc; ixc++){
		for (iyc = 0; iyc < nxc; iyc++){
			if ( iyc >= y0c-1 && ixc >= x1c){//ecart reste jusque cond suivante appliquée, quand ? qd sorti des alentours du trous (condition suiv s'applique en compensation)
                skip_nord = qc;
            }
            if ( iyc >= y1c && ixc >= x1c){
                skip_nord = 0;
            }
            if ( iyc >= y0c && ixc >= x1c){
                skip_sud = qc;
            }
            if ( iyc > y1c && ixc >= x1c){
                skip_sud = 0;
            }
			if(! in_hole((ixc+1)*hc, (iyc+1)*hc)){

				if (check_sw(ixc,iyc, hc, h)){//(u)[transfert de coord uc a u pout le point origin(!! ) puis soustrait ce qui faut pour arriver au point target]
					(*u)[] = 0.25 * (*uc)[nuc];
				}
				else{

				}
				nuc += 1;
			}
		}
	}

	return 0;
}


	


/*
int prolongR(double **u, double **uc, int m){

	int ix, iy;

	//valeurs de u
	double h = 3.0/(double)(m-1);
    double invh2 = 1.0/(h*h);
    int x1 = ((int)(COORD_X1 * (m-1)) /3)  -1; 
    int x0 = ((int)(COORD_X0 * (m-1)) /3)  -1;
    int y1 = ((int)(COORD_Y1 * (m-1)) /3)  -1;
    int y0 = ((int)(COORD_Y0 * (m-1)) /3)  -1;
  
    int nx = m-2;

    int p = y1 - y0 + 1;
    int q = x1 - x0 + 1;
    int n = nx * nx - (p * q);
    
	*u = malloc(n * sizeof(double));

	//valeurs de uc -> les plus determinantes mais ici pas suffisante
	int sizeRc;
	int mc = (m+1)/2;
	int nxc = mc-2;
	int x1c = ((int)(COORD_X1 * (mc-1)) /3)  -1; 
    int x0c = ((int)(COORD_X0 * (mc-1)) /3)  -1;
    int y1c = ((int)(COORD_Y1 * (mc-1)) /3)  -1;
    int y0c = ((int)(COORD_Y0 * (mc-1)) /3)  -1;
  

    int pc = y1c - y0c + 1;
    int qc = x1c - x0c + 1;
    int nc = nxc * nxc - (pc * qc);

	

	int skip_sud = 0;
    int skip_nord = 0;

	int nr = n;
	int nrc = nc;

	int skip_sud = 0;
    int skip_nord = 0;
    int skip_coin_so = 0;
    int skip_coin_se = 0;
    int skip_coin_no = 0;
    int skip_coin_ne = 0;


	int n = 0;
    int nc = 0;
    //passage ligne suiv(plaque complete)
    //pour chaque ligne:
    for (iy = 0; iy < nx; iy++){
        //passage colonne suiv
        for (ix = 0; ix < nx; ix++){      

            //initialisation du retard a ajouter à cause du bord 
            //(depend d'ou on se trouve sur la plaque)
            
            if ( iy >= y0-1 && ix >= x1){
                skip_nord = q;
            }
            if ( iy >= y1 && ix >= x1){ //pblm : pense que y a un skip car point ds zone 
                skip_nord = 0;
            }
            if ( iy >= y0 && ix >= x1){
                skip_sud = q;
            }
            if ( iy > y1 && ix >= x1){
                skip_sud = 0;
            }
            
            //exclu interieur et bord du trou
            if(iy > y1  || iy < y0  || ix < x0 || ix > x1){
                
				//impair impair -> 1/4 somme des 4 autour
                if (iy % 2 == 0 && ix % 2 == 0){
                    //check si bord est bien la ou pense etre pour uc (pas faire une recherche alors que bord)
                    //somme coin gauche bas 
                    //check si pas ds le trou
                    if (iy -1 >= 0  && (iy-1 > y1  || ix -1 < x0 || x1 < ix -1) //condition bas
                            && ix -1 >= 0 && (ix - 1 > x1 || iy -1 > y1 || iy -1 < y0 )){//condition gauche
                        //doit chercher dans la liste donc check si ecart
                        if ( ix -1 < x0 && iy-1 >= y0 && iy -1 <= y1){
                            skip_coin_so = q;
                        }
                        else{
                            skip_coin_so = 0;
                        }
                        (*u)[n] += 0.25 * (*uc)[nc + ix/2 - (nx/2 + 1) + skip_coin_so]; 
                    }
                    else{
                        (*u)[n] += 0.25 * computeBound((ix-1+1)*h,(iy-1+1)*h);
                    }
                    //coin droit bas
                    if (iy - 1 >= 0  && (iy-1 > y1  || ix +1 < x0 || ix +1 > x1) //cond bas
                            && ix + 1 < nx && ( ix + 1 < x0 || iy-1 > y1 || iy-1  < y0)){ //cond droit
                        if ( ix +1 < x0 && iy-1 >= y0 && iy -1 <= y1){
                            skip_coin_se = q;
                        }
                        else{
                            skip_coin_se = 0;
                        }
                        (*u)[n] += 0.25 * (*uc)[nc + 1 + ix/2 - (nx/2 + 1) + skip_coin_se];
                    }
                    else{
                        (*u)[n] += 0.25 * computeBound((ix+1+1)*h,(iy-1+1)*h);
                    }
                    //coin haut gauche
                    if ( iy+1 < nx && ( iy+1 < y0 || ix -1 < x0 || ix -1 > x1) //cond haut gauche
                            && ix -1 >= 0 && (ix - 1 > x1 || iy +1 > y1 || iy +1 < y0 )){ //cond gauche haut
                        if (ix -1 > x1 && iy + 1 >= y0 && iy +1 <= y1){ //qd point no entre ds zone
                            skip_coin_no = q;
                        }
                        else{//qd point centrale en dehors zone
                            skip_coin_no = 0;
                        }
                        (*u)[n] += 0.25 * (*uc)[nc + ix/2 - skip_coin_no];
                    }
                    else{
                        (*u)[n] += 0.25 * computeBound((ix-1+1)*h,(iy + 1+1)*h);
                    }
                    //coin droit haut
                    
                    if ( iy+1 < nx && ( iy+1 < y0 || ix + 1 < x0 || ix + 1 > x1) //cond haut
                            && ix + 1 < nx && ( ix + 1 < x0 || iy + 1 > y1 || iy + 1  < y0)){ //cond droit
                        if (ix + 1 > x1 && iy + 1 >= y0 && iy +1 <= y1){
                            skip_coin_ne = q;
                        }
                        else{
                            skip_coin_ne = 0;
                        }
                        (*u)[n] += 0.25 * (*uc)[nc + 1 + ix/2 - skip_coin_ne];
                    }
                    else{
                        (*u)[n] += 0.25 * computeBound((ix+1+1)*h,(iy+1+1)*h);
                    }
                }

				//impair pair => somme haut + bas
				else if (iy % 2 == 0 && ix % 2 == 1){
                    
                    //somme bas
                    if (iy-1 >= 0 && ( iy-1 > y1  || ix < x0 || ix > x1) ){
                        if ( ix < x0 && iy-1 >= y0 && iy -1 <= y1){ //que dans zone gauche, quand point sud niveau y0 et jusqu' a ce que le point sud depasse y1
                            skip_s = q;
                        }
                        else {
                            skip_s = 0; //cas ou on est dans la zone de droite car passe sur le trou dans 1 sens puis dans l'autre donc s'annule

                        }
                        //+ix/2 car tt les 2 points il faut realigner le point de uc qui de base est ligne au dessus tt a gauche
                        (*u)[n] += 0.5 * (*uc)[nc + ix/2 - (nx/2 + 1) + skip_sud];//-nx/2+1 pour aller ligne en dessous
                    }
                    else{
                        (*u)[n] += 0.5 * computeBound((ix+1)*h,(iy-1+1)*h);
                    }
                    //somme haut
                    if(iy+1 < nx && ( iy+1 < y0 || ix < x0 || ix > x1 )){
                        if ( ix > x1 && iy + 1 >= y0 && iy +1 <= y1){//que dans zone de droite , a partir de point nord sur y0, et que 
                            skip_n = q;
                        }
                        else{
                            skip_n = 0;
                        }
                        (*u)[n] += 0.5 * (*uc)[nc + ix/2 - skip_n]; //pas de nx/2+1 car point uc[nc] deja ligne du haut
                    }
                    else{
                        (*u)[n] += 0.5 * computeBound((ix+1)*h,(iy+1+1)*h);
                    }
                }

                
                //pair impair 1/2 somme gauche droite
                else if (iy % 2 == 1 && ix % 2 == 0){

                    //somme gauche
                    //si pas a droite d'un bord
                    if (ix -1 >= 0 && ( ix - 1 > x1 || iy > y1 || iy < y0 )){
                        
                        (*u)[n] += 0.5 * (*uc)[nc - 1];
                    }
                    else{
                        (*u)[n] += 0.5 * computeBound((ix-1+1)*h,(iy+1)*h);
                    }

                    //somme droit
                    //si pas a gauche d'un bord
                    if (ix + 1 < nx && ( ix + 1 < x0 || iy > y1 || iy < y0) ){
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
*/