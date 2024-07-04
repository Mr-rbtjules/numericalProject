#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "proto.h"

 //V par defaut


/*
cycle V
ce qui fonctionne pour 1500 : 7 4 4 1 20

6145

300
  
cycle w 10 iter 5 4 4 10
cycle v 13 iter 5 4 4 20 un peu plus rapide 1.3

cycle f  10 iter 5 4 4 20 1.3s
*/



/*
        level 0
\      /1
 \    /2
  \  /3
   \/4 = level max = nb de fois qu'on restr

0               0
 \             / 
  1           1
   \         /
    2       2
     \     /
      3   3
       \ /
        4 
*/
//print parameters

//start global variable (m and domain)
globVal_s globVal = { 0, NULL, NULL, NULL, NULL}; 



/*simillar(level+1, fig, step);t tt
pour avoir mm ordre ia ja u r, etc
*/




/// RESOL ///


int mg_method(int m, int iter, int levelMax, int mu, int mode, int cycle, double *res){
	
	initGlobVal(m, levelMax);
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
	
    //'b' conceptuel, cad b au niveau 0 mais ri au niveau inf
	double *bl = NULL;  

	//we pass the adress for each level
	allocGrids(levelMax, &nl, &ial, &jal, &al, &bl, &dl, &rl, &ul);
    
	//precomputation of all Ac and bc

	for (int l = 0; l <= levelMax; l++){
		//adresses of each vector is simply the pointer shifted
        // bl not shifted because only the firs
		probMg(
            globVal.m[l], 
            nl + l, 
            ial + globVal.vectStart[l] + l,
            jal + globVal.matStart[l],
            al + globVal.matStart[l],
            bl
        );
	}
	
    
	int startLevelTg = 0; 

	//u0 = 0
	preInitialization(ul);
    if (LOAD){
		printf("\n Load percentage : \n");
	}
	//start iterations of the multigrid cycle (tg_rec)
	double t1 = mytimer();

    
    if( mg_iter(
            iter, levelMax, globVal.m[0], 
            mu, mu, cycle, nl, ial ,jal,al,rl,ul,dl,bl, res)){
        return 1;
    }
    double t2 = mytimer();
    if (CHRONO){
        printf("\nTemps de solution multigrid method (CPU): %5.1f sec\n",t2-t1);   
    }

    

    //top level so no shift
    double resi = computeResNorm(
            nl, ial, jal, al, bl, ul, rl
        );
    printf("\n final res : %.16f\n", resi);

    if (PLOT){
        plotCycle(levelMax, cycle);
        plot_res(rl + globVal.vectStart[startLevelTg], startLevelTg);
    }
    

    free(nl), free(ial), free(jal), free(al);
    free(bl), free(dl), free(ul) , free(rl);
    freeGlobVal();

    return 0;
}


int mg_iter(int iter, int levelMax, int m, int mu1, int mu2, int cycle, int *nl, int *ial,
             int *jal, double *al, double *rl, double *ul, double *dl, double *bl, double *res){
	double *nres = malloc(*nl * sizeof(double));
	int startLevelTg = 0;
	for (int i = 0; i < iter; i++){
		//initialisation ici ?
		if (LOAD && res != NULL){
            
			printf("   %.3f   ", (double)(i+1)/iter);
        	fflush(stdout);
			printf("\r");
		}   
        //ci = 0 pr tt i != top level
        //ligne tres opti !
        reInitialization(ul);
        //solve_umfpack(*nl, ial, jal, al, bl, ul);
		tg_rec( startLevelTg, m, mu1, mu2, cycle, 
				nl, ial, jal, al, bl, ul, rl, dl);
		if (res != NULL){
            
            res[i] = computeResNorm(
              nl, ial, jal, al, bl, ul, nres
            );
            
        }
		//print(res et u) + save
	}
    free(nres);
	return 0;
}

//on considere qu'on donne quoi a tg rec ?
// -> pointer sur la liste glob et on s'occupe des details inside
int tg_rec(int level, int m, int mu1, int mu2, 
		   int cycle, int *nl, int *ial, int *jal,
		   double *al, double *bl, double *ul, double *rl, double *dl){
		
		//des le level 1 on stock pas u2 mais c2 puis 'deviennent des u lorsqu'on ajoute la correction en remontant
		//pblm tt en haut smoothing Au=b mais en dessous sm Ac = r pareil computeres
	
	if (level < globVal.levelMax){

        //mettre dans une fonction car si level = 0 applique syst Au=b et r = b-Au
        //sinon Ac = r 
        //deux possibilités soit :
        //1: on stock dans les b inf le res r-Ac et dans r les 'b' inf qui sont enft r dans le cours
        //2: on stock dans les b d'abord b puis
        //mu1 smoothing iterations
        //ici on utilise b pour stocker le residu de Ac=r et stock c dans u
        //ici bien le system Au=b dont on s'occupe

        if (EXPLICIT){
            printf("\n\n  LEVEL : %d\n\n", level);
        }

        //getdown

        getDown( level, m, mu1, mu2, 
				nl, ial, jal, al, bl, ul, rl, dl);
		
		//au dessus s'effecctue du haut du v vers le bas
		tg_rec(level+1, m, mu1, mu2, cycle,
		       nl, ial, jal, al, bl, ul, rl, dl);
        
        if (cycle == 1){
            tg_rec(level+1, m, mu1, mu2, cycle,
		       nl, ial, jal, al, bl, ul, rl, dl);
        }
        
        //recursivité , ici important on repart avec Au = b avec u = c 
        //et b = r = R(res syst level-1)
        //
		// on voit que ul est composé de u en 0 puis que des c puis en
		// remontant on va applique les corrections aux corrections,
		//code en dessous s'effectue du bas vers le haut du v
		
		
        //attention revoir b qd prolong
        //et la boucle for ny !!
            //pblm dans indice check
        getUp( level + 1, m, mu1, mu2, 
				nl, ial, jal, al, bl, ul, rl, dl);
        if (EXPLICIT){
            printf("\n\n  LEVEL : %d\n\n", level);
        }

        if (level > 0){
            if (cycle == -1){
                FCycleLoop(level, m, mu1, mu2,
		       nl, ial, jal, al, bl, ul, rl, dl);
            }
        }

        
	}
	else {
        
		//solve coarse pblm
		if (EXPLICIT) {printf("\n\n LEVEL : %d, solve at bottom\n\n", level);}
		
        solveAtCoarseLevel(
            MODE, 
            nl + level, 
            ial + globVal.vectStart[level] + level, //car pour chaque level ia possede n+1 elem
            jal + globVal.matStart[level], //matstart contient la somme des nnz
            al + globVal.matStart[level],
            bl + globVal.vectStart[level],
            ul + globVal.vectStart[level],
            rl + globVal.vectStart[level],
            dl + globVal.vectStart[level]
        );
        
	}
	return 0;
}

void FCycleLoop(int level, int m, int mu1,
			int mu2, int *nl, int *ial, int *jal,
		   double *al, double *bl, double *ul, double *rl, double *dl){

    int initLevel = level;
    while (level < globVal.levelMax){
        getDown( level, m, mu1, mu2, 
				nl, ial, jal, al, bl, ul, rl, dl);
        level += 1;
        if (EXPLICIT){
            printf("\n\n  LEVEL : %d\n\n", level);
        }
    }

    solveAtCoarseLevel(
        MODE, 
        nl + level, 
        ial + globVal.vectStart[level] + level, //car pour chaque level ia possede n+1 elem
        jal + globVal.matStart[level], //matstart contient la somme des nnz
        al + globVal.matStart[level],
        bl + globVal.vectStart[level],
        ul + globVal.vectStart[level],
        rl + globVal.vectStart[level],
        dl + globVal.vectStart[level]
    );

    while (level > initLevel){
        getUp( level, m, mu1, mu2, 
				nl, ial, jal, al, bl, ul, rl, dl);
        level -= 1;
        if (EXPLICIT){
            printf("\n\n  LEVEL : %d\n\n", level);
        }
    } 
}

//on considere qu'on donne quoi a tg rec ?
// -> pointer sur la liste glob et on s'occupe des details inside




int getDown(int level, int m, int mu1,
			int mu2, int *nl, int *ial, int *jal,
		   double *al, double *bl, double *ul, double *rl, double *dl){

    forwardGS(
        mu1, 
        nl + level, 
        ial + globVal.vectStart[level] + level, //car pour chaque level ia possede n+1 elem
        jal + globVal.matStart[level], //matstart contient la somme des nnz
        al + globVal.matStart[level],
        bl + globVal.vectStart[level],
        ul + globVal.vectStart[level],
        rl + globVal.vectStart[level],
        dl + globVal.vectStart[level]
    );

    computeRes(
        nl + level, 
        ial + globVal.vectStart[level] + level, //car pour chaque level ia possede n+1 elem
        jal + globVal.matStart[level], //matstart contient la somme des nnz
        al + globVal.matStart[level],
        ul + globVal.vectStart[level],
        bl + globVal.vectStart[level],
        rl + globVal.vectStart[level]
    );

    //ici ligne cruciale avant de restrict residu et on met dans b, 
    //cad rc = b(level en dessous)
    restrictR(
        level,
        rl + globVal.vectStart[level],
        bl + globVal.vectStart[level+1]
    );

    return 0;
}

int getUp(int level, int m, int mu1,
			int mu2, int *nl, int *ial, int *jal,
		   double *al, double *bl, double *ul, double *rl, double *dl){

    addProlCorrection(
        level,
        ul + globVal.vectStart[level-1], //ajoute la corr ici
        ul + globVal.vectStart[level]
    );
        
    backwardGS(
        mu2, 
        nl + level-1, 
        ial + globVal.vectStart[level-1] + level-1, //car pour chaque level ia possede n+1 elem
        jal + globVal.matStart[level-1], //matstart contient la somme des nnz
        al + globVal.matStart[level-1],
        bl + globVal.vectStart[level-1],
        ul + globVal.vectStart[level-1],
        rl + globVal.vectStart[level-1],
        dl + globVal.vectStart[level-1]
    );
    return 0;
}

int forwardGS(int iter, int *n, int *ia, int *ja, double *a,
			    double *b, double *u, double *r, double *d){
    //se passe sur un mm level donc pointer change pas
     
    if (EXPLICIT){
        printf("Forward GS, start stationnary iteration recursively with B = L(A)\n");
    }
	
	//resoud Au = b avec iter iteration
    //forward == 1 => B = L(A)
    //deja initialiser à 0
    stationaryIter(iter, n, ia, ja, a, b, u, r, d, 1);
	
    return 0;
}
int backwardGS(int iter, int *n, int *ia, int *ja, double *a,
				double *b, double *u, double *r, double *d){
    
	if (EXPLICIT){
        printf("Backward GS, start stationnary iteration recursively with B = U(A)\n");
    }
    //pas besoin d'init utilise le c precedent
	stationaryIter(iter, n, ia, ja, a, b, u, r, d, -1);

	return 0;
}

int jacobiIter(int iter, int *n, int *ia, int *ja, double *a,
			    double *b, double *u, double *r, double *d){
	if (EXPLICIT){
        printf("Jacobi , start stationnary iteration recursively with B = D(A)\n");
    }						
	stationaryIter(iter, n, ia, ja, a, b, u, r, d, 0);
	return 0;
}


#define SCALE_FACT 0.5

/*fonctionnnement des level
-celui tt en haut avec le h le plus petit = level 0 puis coarse 1 coarse coarse 2 etc
- qd passe level en arg c'est tj le level auquel on est pas celui vers lequel on va
*/

//on ne passe pas juste rl, car rp rc pointe reellement dans des listes diff
//rc = bl (level +1)
int restrictR(int level, double *rp, double *rc){
    //level = 0 = ou on est ->1
    // level -> level +1
    int mp = globVal.m[level];
    int mc = globVal.m[level+1];
    double hp, invh2p;
    int x0p,x1p,y0p,y1p, nxp, nyp, np, nnzp;
    computeParamLevel(mp, &hp, &invh2p, &x0p, &x1p, 
                      &y0p, &y1p, &nxp, &nyp, &np, &nnzp);
    
    //level ou on va =level+1
    double hc, invh2c;
    int x0c,x1c,y0c,y1c, nxc, nyc, nc, nnzc;
    computeParamLevel(mc, &hc, &invh2c, &x0c, &x1c, 
                      &y0c, &y1c, &nxc, &nyc, &nc, &nnzc);
    
    if (EXPLICIT){
		printf("\nRestrict -level %d: mp = %d hp = %lf ", level, mp, hp);
		printf("np = %d nnzp = %d nxp = %d nyp = %d ", np, nnzp, nxp, nyp);
		printf(" x0p = %d x1p = %d y0p = %d y1p = %d \n", x0p, x1p, y0p, y1p);
		printf("\n         -level %d mc = %d hc = %lf ",level+1, mc, hc);
		printf("nc = %d nnzc = %d nxc = %d nyc = %d ", nc, nnzc, nxc, nyc);
		printf("x0c = %d x1c = %d y0c = %d y1c = %d \n", x0c, x1c, y0c, y1c);
	}
    if (rc == NULL || rp == NULL){
        printf("r no memory \n");
        return 1;
    }
    int np_save = np;
	np = 0;
    int nc_save = nc;
	nc = 0;


    

    int test = 1;

	int ind;
	for (int iyp = 0; iyp < nyp; iyp++){
        //passage colonne suiv
        for (int ixp = 0; ixp < nxp; ixp++){      

            //exclu interieur et bord du trou
            if( ! in_hole(ixp,iyp,y0p,y1p,x0p,x1p)){
                
                //garde que ligne ind impair  & col impair 
                if (iyp % 2 == 1 && ixp % 2 == 1){  
                    double nb_voisins = 0;
                    rc[nc] = 0;
                    //marquer le début de la ligne suivante dans le tableau 'ia'            
                    //replissage de la ligne : voisin sud //ui-1
                    // + verification si pas au dessus d'un bord
                    
                    if (check_sud(ixp,iyp,y0p,y1p,x0p,x1p,nyp) ){
                        ind = indice(ixp,iyp-1,y0p,y1p,x0p,x1p, nxp);
                        rc[nc] += rp[ind];
                        nb_voisins += 1;
                    }    

                    else{
                        if (!test){
                            nb_voisins += 1;
                        }
                        //residus sur bord censé etre 0 ?
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT; 
                        //rc[*nc] += 0.25 * computeBound((ixp+1)*hp,(iyp+1 -1)*hp) * SCALE_FACT;//si fct qui calcule chaque fois :juste coord sinon
                    }

                    //replissage de la ligne : voisin ouest 
                    //si pas a droite d'un bord
                    if (check_west(ixp,iyp,y0p,y1p,x0p,x1p,nxp)){
                        ind = np -1;
                        rc[nc] += rp[ind];
                        nb_voisins += 1;
                    }
                    else{
                        if (!test){
                            nb_voisins += 1;
                        }
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT;
                        //rc[*nc] += 0.25 * computeBound((ixp+1-1)*hp,(iyp+1)*hp) * SCALE_FACT;
                    }

                    
                    
                    
                   
                    // replissage de la ligne : élém. diagonal
                    

                    // replissage de la ligne : voisin est
                    //si pas a gauche d'un bord
                    if (check_est(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){
                        ind = np +1;
                        rc[nc] += rp[ind];
                        nb_voisins += 1;
                    }
                    else{
                        if (!test){
                            nb_voisins += 1;
                        }
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT;
                        //rc[*nc] += 0.25 * computeBound((ixp+1+1)*hp ,(iyp+1)*hp) * SCALE_FACT;
                        
                    }

                    // replissage de la ligne : voisin nord
                    //si pas en dessous d'un bord
                    if ( check_nord(ixp,iyp,y0p,y1p,x0p,x1p,nyp) ){
						ind = indice(ixp,iyp+1,y0p,y1p,x0p,x1p, nxp);
                        rc[nc] += rp[ind];
                        nb_voisins += 1;
                    }

                    else{
                        if (!test){
                            nb_voisins += 1;
                        }
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT;
                        //rc[*nc] += 0.25 * computeBound((ixp+1)*hp, (iyp+1+1)*hp) * SCALE_FACT;
                    }
                    rc[nc] /= nb_voisins;
                    rc[nc] += rp[np];
                    rc[nc] *= SCALE_FACT;

                    // numéro de l'équation
                    nc += 1;
                }
                np += 1;
            }
        }
    }
	if (np_save != np){
		printf("nr != n\n");
		return 1;
	}
	if (nc_save != nc){
		printf("nc_save %d != nc %d\n", nc_save, nc );
		return 1;
	}
	return 0;
}

int addProlCorrection(int level, double *up, double *uc){

    //level 1 = on ou on est -> on va vers 0
    // level -> level -1

    int mp = globVal.m[level-1];
    int mc = globVal.m[level];
    

    double hp, invh2p;
    int x0p,x1p,y0p,y1p, nxp, nyp, np, nnzp;
    computeParamLevel(mp, &hp, &invh2p, &x0p, &x1p, 
                      &y0p, &y1p, &nxp, &nyp, &np, &nnzp);
    
    //level ou on va =level+1
    double hc, invh2c;
    int x0c,x1c,y0c,y1c, nxc, nyc, nc, nnzc;
    computeParamLevel(mc, &hc, &invh2c, &x0c, &x1c, 
                      &y0c, &y1c, &nxc, &nyc, &nc, &nnzc);

    if (EXPLICIT){
		printf("\nProlong -level %d: mp = %d hp = %lf ", level-1, mp, hp);
		printf("np = %d nnzp = %d nxp = %d nyp = %d ", np, nnzp, nxp, nyp);
		printf(" x0p = %d x1p = %d y0p = %d y1p = %d \n", x0p, x1p, y0p, y1p);
		printf("\n         -level %d mc = %d hc = %lf ", level, mc, hc);
		printf("nc = %d nnzc = %d nxc = %d nyc = %d ", nc, nnzc, nxc, nyc);
		printf("x0c = %d x1c = %d y0c = %d y1c = %d \n", x0c, x1c, y0c, y1c);
	}

    
    
    if (up == NULL || uc == NULL){
        printf("no memory  prol\n");
        return 1;
    }
	
    int nc_save = nc;
	nc = 0;
    int np_save = np;
	np = 0;

	int ind;

    int test = 0;
    double nb;
    double corr;

	for (int iyp = 0; iyp < nyp; iyp++){
        //passage colonne suiv
        for (int ixp = 0; ixp < nxp; ixp++){      

            //exclu interieur et bord du trou
            if(! in_hole(ixp,iyp,y0p,y1p,x0p,x1p)){
                
                //up[*np] = 0; ici on add la prolongation au vecteur de niveau au dessus

				//impair impair -> 1/4 somme des 4 autour
                
                if (iyp % 2 == 0 && ixp % 2 == 0){ 
                    nb = 0;
                    corr = 0;

                    //check si bord est bien la ou pense etre pour uc (pas faire une recherche alors que bord)
                    //somme coin gauche bas 
                    //check si point prol pas ds le trou
                    
                    if (check_sw(ixp,iyp,y0p,y1p,x0p,x1p,nxp,nyp)){ // question est ce que prolong peut etre domaine et coarse ds un board ? non ca depend que de p
                        
						ind = indice((ixp/2)-1,(iyp/2)-1, 
									 y0c,y1c,x0c,x1c, nxc);//juste /2 -1 car point prol au milieu des 4 tjrs pair
                        corr += uc[ind]; 
                        nb +=1;
                    }
                    //on met qqchose ici ?
                    /*else{
                        double bound = computeBound((ixp+1-1)*hp,(iyp+1-1)*hp);
                        
                        up[*np] += 0.25 * bound;
                       
                    }*/
                    //coin droit bas
                    if (check_se(ixp,iyp,y0p,y1p,x0p,x1p,nxp,nyp)){ //cond droit
                        ind = indice((ixp/2) ,(iyp/2) - 1, 
									 y0c,y1c,x0c,x1c, nxc);
                        corr += uc[ind];
                        nb +=1;
                    }
                   /* else{
                        up[*np] += 0.25 * computeBound((ixp+1+1)*hp,(iyp+1-1)*hp);
                    }*/
                    //coin haut gauche
                    if (check_nw(ixp,iyp,y0p,y1p,x0p,x1p,nxp, nyp) ){ //cond gauche haut
                        ind = indice((ixp/2 - 1),(iyp/2), 
									 y0c,y1c,x0c,x1c, nxc);
                        corr += uc[ind];
                        nb +=1;
                    }
                    /*else{
                        double bound = computeBound((ixp+1-1)*hp,(iyp + 1+1)*hp);
                        up[*np] += 0.25 * bound;
                        
                    }*/
                    //coin droit haut
                    
                    if (check_ne(ixp,iyp,y0p,y1p,x0p,x1p,nxp, nyp) ){ //coin droit
                        ind = indice((ixp/2),(iyp/2), y0c,y1c,x0c,x1c, nxc);
                        corr += uc[ind];
                        nb +=1;
                    }
                    up[np] += corr/nb;
                    /*else{
                        up[*np] += 0.25 * computeBound((ixp+1+1)*hp,(iyp+1+1)*hp);
                    }*/
                }

				//impair pair => somme haut + bas
				else if (iyp % 2 == 0 && ixp % 2 == 1){
                    nb = 0;
                    corr = 0;
                    //somme bas
                    if (check_sud(ixp,iyp,y0p,y1p,x0p,x1p,nyp) ){
                       
                        ind = indice((ixp/2),(iyp/2)-1,
									 y0c,y1c,x0c,x1c, nxc);
                        corr += uc[ind];
                        nb +=1;
                    }
                    /*else{
                        up[*np] += 0.5 * computeBound((ixp+1)*hp,(iyp+1-1)*hp);
                    }*/
                    //somme haut
                    if(check_nord(ixp,iyp,y0p,y1p,x0p,x1p,nyp)){
                        ind = indice((ixp/2),(iyp/2), y0c,y1c,x0c,x1c, nxc);
                        corr += uc[ind]; //pas de nx/2+1 car point uc[nc] deja ligne du haut
                        nb +=1 ;
                    }

                    up[np] += corr/nb;
                    /*else{
                        up[*np] += 0.5 * computeBound((ixp+1)*hp,(iyp+1+1)*hp);
                    }*/
                }
                //pair impair 1/2 somme gauche droite
                else if (iyp % 2 == 1 && ixp % 2 == 0){
                    nb = 0;
                    corr = 0;
                    //somme gauche
                    //si pas a droite d'un bord
                    if (check_west(ixp,iyp,y0p,y1p,x0p,x1p,nxp)){
                        corr += uc[nc - 1];
                        nb +=1;
                    }
                    /*else{
                        up[*np] += 0.5 * computeBound((ixp+1-1)*hp,(iyp+1)*hp);
                    }*/

                    //somme droit
                    //si pas a gauche d'un bord
                    if (check_est(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){
                        corr += uc[nc]; //nc car deja +1 car type point precedent + 1
                        nb +=1;
                    }
                    up[np] += corr/nb;
                    /*else{
                        up[*np] += 0.5 * computeBound((ixp+1+1)*hp,(iyp+1)*hp);
                    }*/  
                }
                //ligne impair  & col impair => elem identique
                else if (iyp % 2 == 1 && ixp % 2 == 1){  //pour 1 meme ligne alterne entre type 1 et 2, commence par 1 finis par 1

                    up[np] += uc[nc];
                    nc += 1;       //=> 2 choses, pr point type 2 droite = nc gauche == nc-1
                    // et aussi que pour type 3 et 4 nc rpz point au dessus tt a gauche
                }
                np += 1;
            }
        }
    }
    
    if (np != np_save){
        printf(" err np %d npcheck %d\n", np, np_save);
    }
    if (nc != nc_save){
        printf(" err nc %d nccheck %d\n", nc, nc_save);
    } 
	return 0;
}


int stationaryIter(int iter, int *n, int *ia, int *ja,
                   double *a, double *b, double *u, 
				   double *r, double *d, int forward){
//attention considere que deja init
//forward = 0 => backwardGs
	
	if (iter != 1){
		stationaryIter(iter-1, n, ia, ja, a,
					 b, u, r, d, forward);
	}
	//rm (aussi initialisation de r0)
	computeRes(n, ia, ja, a, u, b, r);
	
	// dm == solve B*d = r (B construit a partir de A direct 
	// ds la methode)
	if (forward == 1){
        //ici A correspond a L(A)= B
        //on ne modifie que d (ici qu'on initialise d)
		gaussResL(n , ia, ja, a, d, r);
	}
	else if (forward == -1){
		gaussResU(n , ia, ja, a, d, r);
	}
	else if (forward == 0){
		gaussResD(n , ia, ja, a, d, r);	
	}


	//um+1 (ajout de la correction)
    //peut opti ici en modifiant les gauss de sorte a ajouter 
    //direct a u et pas apres
	addVect(n , u, d);

	return 0;
}


//pas besoin de crer L, juste prendre ds A 
int gaussResL(int *n , int *il, int *jl, double *l, double *x, double *b){
	//resoud Lx = b trouve x
    if (EXPLICIT){
        printf("  Gauss resolution for L(A)d=r\n");
    }
	int i = 0;
	while (i < *n){

		int start_ind_jl = il[i];
		
		int end_ind_jl = il[i+1];
		
		x[i] = b[i]; //copie b sur u
		while (start_ind_jl < end_ind_jl && jl[start_ind_jl] < i){

			x[i] -= (l[start_ind_jl] 
							* x[jl[start_ind_jl]]);
			start_ind_jl += 1;
		}
		//car start == end et donc fin de ligne => elem diag
        if (RELAX){
            double relax = 1;//relax function ?
            x[i] *= relax/l[start_ind_jl];
        }
        else{
            x[i] /= l[start_ind_jl];
        }
		
		i += 1;
	}
	return 0;
}

int gaussResU(int *n , int *iu, int *ju, double *u, double *x, double *b){ 
	//resoud Ux = b trouve x
    if (EXPLICIT){
        printf("  Gauss resolution for U(A)d=r\n");
    }
	int i = *n-1;
	
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

int gaussResD(int *n , int *ia, int *ja, double *a, double *x, double *b){
	//resoud Lx = b trouve x
    if (EXPLICIT){
        printf("  Gauss resolution for D(A)d=r\n");
    }
	int i = 0;
	while (i < *n){

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


int solveAtCoarseLevel(int mode, int *n, int *ia, int *ja, double *a, 
						double *b, double *u, double *r, double *d){

	if (mode == 0){
		
		solve_umfpack(*n, ia, ja, a, b, u);
		if (EXPLICIT){
			printf("\n umfpack end res(de Ac=r)=  %lf\n", computeResNorm(n,ia,ja,a,b,u,r));

		}
	}
	else if (mode == 1) {
		
        symGS(ITERCORE, 0, n, ia, ja, a, b, u, r, d);
		if (EXPLICIT){
			printf("\n symGS end res(de Ac=r)=  %lf\n", computeResNorm(n,ia,ja,a,b,u,r));	
		}
	}
	else {
		
		
        jacobiIter(ITERCORE, n, ia, ja, a, b, u, r, d);
		if (EXPLICIT){
			printf("\n jacobi end res(de Ac=r)=  %lf\n", computeResNorm(n,ia,ja,a,b,u,r));	
		}
        
	}
	return 0;
}

//verif si pas jabobi inside
int symGS(int iter, double tol, int *n, int *ia, int *ja, double *a,
			    double *b, double *u, double *r, double *d){
	//printf("\n Sym gs initial res=  %lf\n", computeResNorm(n,ia,ja,a,b,u,r));
								
	if (iter){
		int i = 0;
		while (i < iter){
			
			forwardGS( 1, n, ia, ja, a, b, u, r, d);
            //jacobiIter(1, n, ia, ja, a, b, u, r, d);
      		backwardGS( 1, n, ia, ja, a, b, u, r, d);
			i++;
			if (tol){//tol > computeNorm(n,r)){
				i = iter;
			}
		}
	}
	/*
    else{
		while (computeNorm(n,r) > tol){
			forwardGS( 1, n, ia, ja, a, b, u, r, d);
      		backwardGS( 1, n, ia, ja, a, b, u, r, d);
		}
	}
    */

	return 0;
}




int allocGrids(int levelMax, int **nl, int **ial,
               int **jal, double **al, double **bl,
			   double **dl, double **rl, double **ul){ 
    *nl = (int*)malloc((levelMax +1)*sizeof(int));
    // pour taille de ia : c'est le nb de n total
    // (dernier de la list vectStart mais pour chaque level on a 1 eleme
    // de plus que n, nb de niveaux = levelmAX +1
    *ial = (int*)malloc(
        (globVal.vectStart[levelMax+1] + (levelMax+1)) * sizeof(int)
    ); 
    //levelmax = indice du debut du dernier niveau donc +1 p
    // our avoir le nb tot d'elem matriciel
	*jal = (int*)malloc(globVal.matStart[levelMax+1] * sizeof(int));
	*al = (double*)malloc(globVal.matStart[levelMax+1] * sizeof(double));
    //elem vectoriel
	*rl = (double*)malloc(globVal.vectStart[levelMax+1] * sizeof(double));
	*dl = (double*)malloc(globVal.vectStart[levelMax+1] * sizeof(double));
	*ul = (double*)malloc(globVal.vectStart[levelMax+1] * sizeof(double));
    //on a besoin que de b au debut voir multigrid algo description
    //2eme elem de vectstart = n pour le top level
	*bl = (double*)malloc(globVal.vectStart[levelMax+1] * sizeof(double));


    if (*nl == NULL || *bl == NULL || *ial == NULL || *jal == NULL || 
        *al == NULL ||*dl == NULL || *rl == NULL || *ul == NULL){
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n");
        return 1;
    }
    if (EXPLICIT){
        printf("\nMemory allocated for level 0 to level %d\n", levelMax);
    }
	
	return 0;
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

void plotVCycle(int levelMax, int cycle){
    printf("\nTop Level\n");
    for (int i = 0; i <= 2*levelMax; i++){
       //espace(=i) nombre espace(1+(levmax*2)-1)*2 nombre
        if ( i % 2 == 0){
            int nombre = i/2;
            //espace:
            for (int j=0; j < i; j++){
                printf(" ");
            }
            printf("%d", nombre);
            if (i != 2*levelMax){
                int spaces = 1+((levelMax*2)-1)*2 -2*i; 
                for (int j=0; j < spaces; j++){
                    printf(" ");
                }
                printf("%d\n", nombre);
            }
        }
        // espace \ espace /
        else {
            //espace:
            for (int j=0; j < i; j++){
                printf(" ");
            }
            printf("\\");
            
            int spaces = 1+((levelMax*2)-1)*2 -2*i; 
            for (int j=0; j < spaces; j++){
                printf(" ");
            }
            printf("/\n");
            

        }
    }
    printf("\nBottom Level\n");
}


//lanczos iter genre 50 max c bon si n 1000 ou plus a tester
int computeExtremeEv(int *n, int lanczosIter, int *ia, int *ja, 
                    double *a, double *lowest, double *largest, int forward){
    //fairte list param necessaire : 
    //-> iter la
    //faire list memoire necessaire

    //tridiag de A avec lanczos
    double *alphal; //size ?
    double *betal; // size ?
    //add forward and Bx=Av gauss res
    lanczosAlgo(lanczosIter, n, ia, ja, a, &alphal, &betal); //->3 malloc v + 2 malloc al bet

    //seek for lowest e vector (succesive thomas resol for Tx=v v=x)
    //inverse iteration method
    double *evect; //init
    double mu; 
    computeMu(alphal, betal, n , &mu);
    int max_iter_Thom = 40; 
    double evTol = 1e-6; //?
    invIter(alphal, betal, evect, n, &mu, max_iter_Thom, evTol);
    //rayleigh quot for eigenvalue
    double *Ax; //init
    rayleighQuot(n, lowest, ia, ja, a, Ax, evect);

    //seek for largest e vector (T-muI)x=v mu = μ=max⁡(∣αj∣)+max⁡(∣βj∣)
    //shifted inverse iteration method
    //double *evect; //reinit
    mu = 0;
    invIter(alphal, betal, evect, n, &mu, max_iter_Thom, evTol);
    //rayleigh quot for eigenvalue
    rayleighQuot(n, largest, ia, ja, a, Ax, evect);

    free(evect);
    free(Ax);
    return 0;
}


//1000X1000 seul 50 iter suffisant !!mettre en param to see speed
int lanczosAlgo(int iter, int *n, int *ia, int *ja, 
                 double *a, double **alpha, double **beta){

    double tol = 1e-12;

    double *v_current = malloc(*n * sizeof(double));
    double *v_previous = calloc(*n, sizeof(double)); // Initializes to zero
    double *v_next = malloc(*n * sizeof(double));

    if (v_current == NULL || v_next == NULL || v_previous == NULL){
        return 1;
    }

    initRandomV(v_current, n);

    (*alpha) = (double *)malloc(iter * sizeof(double));
    (*beta) = (double *)malloc((iter+1) * sizeof(double));

    (*beta)[0] = 0; //b1=0
    for (int j = 1; j <= iter; j++){
        //vj+1=next = v[2] au debut
        //direct wj dans vj+1 = A*vj - betaj*vj-1
        multCsrSubvector(n, ia, ja, a, v_next, v_current, v_previous, &((*beta)[j-1]));
        //alphaj = <wj,vj> = <vj+1,vj>
        scalProd(n, v_next, v_current, &((*alpha)[j-1]));//j-1 car ind 0
        //wj=vj+1=vj+1 - alj*vj
        subVectProd(n, &(*alpha)[j-1], v_next, v_current);
        //bj+1 = ||wj==vj+1 ||
        computeVectNorm2(n,&((*beta)[j]), v_next);

        //tol ?
        if ( (*beta)[j] < tol){
            j = iter +1;
        }

        double *temp = v_previous;
        v_previous = v_current;
        v_current = v_next;
        v_next = temp;
        return 0;
    }


    free(v_current);
    free(v_next);
    free(v_previous);
    return 0;
}

//allocateoutside ?
int invIter(double *alpha, double *beta, double *v, 
            int *n, double *mu, int max_iter, double tol){

    //iterative thomas algo with mu=0
    double *y = malloc(*n * sizeof(double));
    if (y == NULL) {
        // Handle memory allocation failure
        return 1;
    }

    double norm;
    double error;
    for (int iter = 0; iter < max_iter; iter++) {
        // Solve Ty = v
        thomas_algorithm(alpha, beta, v, y, mu, n);

        // Normalize y and copy it back to v
        normalize(n, y);
        //replace v par y
        for (int i = 0; i < *n; i++) {
            v[i] = y[i];
        }

        // Check for convergence
        computeVectNorm2( n, &norm, v);
        error = norm - 1.0;
        if ((error < 0 && error > -tol) 
            || (error >= 0 && error < tol)) {
            iter = max_iter;
        }
    }

    // Copy the result to x
    

    free(y);
    return 0;
}



int thomas_algorithm(double *alpha, double *beta, 
                      double *v, double *x, double *mu, int *n){
    //solve(T-uI)x = v
    double *c_prime = malloc(*n * sizeof(double));
    if (c_prime == NULL) {
        // Handle memory allocation failure
    }

    c_prime[0] = beta[0] / (alpha[0] - *mu);
    v[0] = v[0] / (alpha[0] - *mu);

    // Forward elimination
    for (int i = 1; i < *n; i++) {
        double m = 1.0 / (alpha[i] - *mu - beta[i-1] * c_prime[i-1]);
        c_prime[i] = beta[i] * m;
        v[i] = (v[i] - beta[i-1] * v[i-1]) * m;
    }

    // Backward substitution
    x[*n-1] = v[*n-1];
    for (int i = *n-2; i >= 0; i--) {
        x[i] = v[i] - c_prime[i] * x[i+1];
    }

    free(c_prime);
    return 0;
}

//check randomness
int initRandomV(double *v, int *n){
    
    // init randomness
    srand((unsigned int)time(NULL));

    double norm = 0;
    for (int i = 0; i < *n; i++) {
        v[i] = (double)rand() / RAND_MAX * 2.0 - 1.0; // Random values bet -1 and 1
        norm += v[i] * v[i];
    }

    norm = sqrt(norm);

    // Normalize the vector
    for (int i = 0; i < *n; i++) {
        v[i] /= norm;
    }
    return 0;
}


int multCsrSubvector(int *n, int *ia, int *ja, double *a, 
                     double *v2, double *v1, double *v0, double *beta){
	//do v2 = A*v1 - beta*v0
	
	int i = 0;
	int jai = 0;
	while (i < *n){
		int ite = ia[i + 1] - ia[i];
		int j = 0;
        v2[i] = 0;
		while (j < ite){
			v2[i] += a[jai + j] * v1[ja[jai + j]]; 		
			j += 1;
		}
        //comme ca direct 
        v2[i] -= (*beta) * v0[i];

		jai += ite;
		i += 1;
	}

	return 0;
}



int rayleighQuot(int *n, double *res, int *ia, int *ja, 
                 double *a, double *Ax, double *x){

    multCsrVector(n, ia, ja, a, Ax, x);
    double xAx;
    double xx;
    scalProd( n, x, Ax, &xAx);
    scalProd( n , x, x, &xx);

    if (xx != 0){
        *res = xAx/xx;
    }
    else {
        return 1;
    }

    return 0;
}




int computeVectNorm2(int *n, double *norm, double *v){

    *norm = 0;
    for (int i = 0; i < *n; i++) {
        *norm += v[i] * v[i];
    }
    *norm = sqrt(*norm);
    return 0;
}

int normalize(int *n, double *v){

    double norm;
    computeVectNorm2(n, &norm, v);
    if (norm > 0){
        for(int i = 0; i < *n; i++){
            v[i] /= norm;
        }
    }
    else{
        return 1;
    }
    return 0;
}

int vectScalDivide(int *n, double *v, double *beta){

    int i = 0;
    while (i < *n){
        v[i] /= *beta;
        i += 1;
    }
    return 0;
}



int computeMu(double *alpha, double *beta, int *n, double *mu) {
    double max_alpha = fabs(alpha[0]);
    double max_beta = fabs(beta[0]);

    for (int i = 1; i < *n; i++) {
        if (fabs(alpha[i]) > max_alpha) {
            max_alpha = fabs(alpha[i]);      //!!!!al bet not same size
        }
        if (i < *n - 1 && fabs(beta[i]) > max_beta) {
            max_beta = fabs(beta[i]);
        }
    }

    *mu =  max_alpha + max_beta;
    return 0;
}

//inverse iteration -> smallest
//random v, u=0

//solve (T)x = v
// v =x and repeat

//shifted inverse iteration --> largest
// μ=max⁡(∣αj∣)+max⁡(∣βj∣)
////solve (T-u)x = v
//norm v
//v = x repeat
//rayleigh quotient


