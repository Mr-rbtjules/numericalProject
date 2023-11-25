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
//type of domain
//si marche pas check indice restrict prolong
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
#define M 6145
#define LEVELMAX  10
#define MU1 4
#define MU2 4
#define MODE 1 //0umpf 1 sym 2jacobi 
#define ITERCORE 20


/*
ce qui fonctionne pour 1500 : 7 4 4 1 20



*/



/*
        level 0
\      /1
 \    /2
  \  /3
   \/4 = level max = nb de fois qu'on restr
*/
//print parameters
#define EXPLICIT 0
#define LOAD 1
#define CHRONO 1
#define RELAX 0
#define PLOT 0

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
	
    //'b' conceptuel, cad b au niveau 0 mais ri au niveau inf
	double *bl = NULL;  

	//we pass the adress for each level
	allocGrids(&nl, &ial, &jal, &al, &bl, &dl, &rl, &ul);
    
	//precomputation of all Ac and bc

	for (int l = 0; l <= LEVELMAX; l++){
		//adresses of each vector is simply the pointer shifted
        // bl not shifted because only the first needed in prob
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
            iter, LEVELMAX, globVal.m[0], 
            MU1, MU2, nl, ial ,jal,al,rl,ul,dl,bl)){
        return 1;
    }
    double t2 = mytimer();
    if (CHRONO){
        printf("\nTemps de solution multigrid method (CPU): %5.1f sec\n",t2-t1);   
    }

    //top level so no shift
    printf("\n final res : %lf\n", computeResNorm(
            nl, ial, jal, al, bl, ul, rl
        ));

    if (PLOT){
        plot_res(rl + globVal.vectStart[startLevelTg], startLevelTg);
    }
    

    free(nl);
    free(ial);
    free(jal);
    free(al);
    free(bl);
    free(dl);
    free(ul);
    free(rl);

	return 0;
}

int mg_iter(int iter, int levelMax, int m, int mu1, int mu2, int *nl, int *ial,
             int *jal, double *al, double *rl, double *ul, double *dl, double *bl){
	
	int startLevelTg = 0;
	for (int i = 0; i < iter; i++){
		//initialisation ici ?
		if (LOAD){
			printf("   %.3f   ", (double)(i+1)/iter);
        	fflush(stdout);
			printf("\r");
		}
        //ci = 0 pr tt i != top level
        //ligne tres opti !
        reInitialization(ul);

		tg_rec( startLevelTg, m, mu1, mu2, 
				nl, ial, jal, al, bl, ul, rl, dl);
		
		//print(res et u) + save
	}
	return 0;
}

//on considere qu'on donne quoi a tg rec ?
// -> pointer sur la liste glob et on s'occupe des details inside
int tg_rec(int level, int m, int mu1,
			int mu2, int *nl, int *ial, int *jal,
		   double *al, double *bl, double *ul, double *rl, double *dl){
		
		//des le level 1 on stock pas u2 mais c2 puis 'deviennent des u lorsqu'on ajoute la correction en remontant
		//pblm tt en haut smoothing Au=b mais en dessous sm Ac = r pareil computeres
	
	if (level < LEVELMAX){

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
		
		//au dessus s'effecctue du haut du v vers le bas
		tg_rec(level+1, m, mu1, mu2,
		       nl, ial, jal, al, bl, ul, rl, dl);
        
        
        if (EXPLICIT){
            printf("\n\n  LEVEL : %d\n\n", level);
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
        addProlCorrection(
            level+1,
            ul + globVal.vectStart[level], //ajoute la corr ici
            ul + globVal.vectStart[level+1]
        );
        
        backwardGS(
            mu2, 
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

	int ind;
	for (int iyp = 0; iyp < nyp; iyp++){
        //passage colonne suiv
        for (int ixp = 0; ixp < nxp; ixp++){      

            //exclu interieur et bord du trou
            if( ! in_hole(ixp,iyp,y0p,y1p,x0p,x1p)){
                
                //garde que ligne ind impair  & col impair 
                if (iyp % 2 == 1 && ixp % 2 == 1){  
                    rc[nc] = 0;
                    //marquer le début de la ligne suivante dans le tableau 'ia'            
                    //replissage de la ligne : voisin sud //ui-1
                    // + verification si pas au dessus d'un bord
                    
                    if (check_sud(ixp,iyp,y0p,y1p,x0p,x1p,nyp) ){
                        ind = indice(ixp,iyp-1,y0p,y1p,x0p,x1p, nxp);
                        rc[nc] += 0.25 * rp[ind] * SCALE_FACT;
                    }    

                    else{
                        //residus sur bord censé etre 0 ?
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT; 
                        //rc[*nc] += 0.25 * computeBound((ixp+1)*hp,(iyp+1 -1)*hp) * SCALE_FACT;//si fct qui calcule chaque fois :juste coord sinon
                    }

                    //replissage de la ligne : voisin ouest 
                    //si pas a droite d'un bord
                    if (check_west(ixp,iyp,y0p,y1p,x0p,x1p,nxp)){
                        ind = np -1;
                        rc[nc] += 0.25 * rp[ind] * SCALE_FACT;
                    }
                    else{
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT;
                        //rc[*nc] += 0.25 * computeBound((ixp+1-1)*hp,(iyp+1)*hp) * SCALE_FACT;
                    }

                    ind = np;
                    rc[nc] += rp[ind] * SCALE_FACT;
                   
                    // replissage de la ligne : élém. diagonal
                    

                    // replissage de la ligne : voisin est
                    //si pas a gauche d'un bord
                    if (check_est(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){
                        ind = np +1;
                        rc[nc] += 0.25 * rp[ind] * SCALE_FACT;
                         
                    }
                    else{
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT;
                        //rc[*nc] += 0.25 * computeBound((ixp+1+1)*hp ,(iyp+1)*hp) * SCALE_FACT;
                        
                    }

                    // replissage de la ligne : voisin nord
                    //si pas en dessous d'un bord
                    if ( check_nord(ixp,iyp,y0p,y1p,x0p,x1p,nyp) ){
						ind = indice(ixp,iyp+1,y0p,y1p,x0p,x1p, nxp);
                        rc[nc] += 0.25 * rp[ind] * SCALE_FACT;
                    }
                    else{
                        //(*rc)[*nc] += 0.25 * rp[np] * SCALE_FACT;
                        //rc[*nc] += 0.25 * computeBound((ixp+1)*hp, (iyp+1+1)*hp) * SCALE_FACT;
                    }
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
	for (int iyp = 0; iyp < nyp; iyp++){
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
                    
                    if (check_sw(ixp,iyp,y0p,y1p,x0p,x1p,nxp,nyp)){ // question est ce que prolong peut etre domaine et coarse ds un board ? non ca depend que de p
                        
						ind = indice((ixp/2)-1,(iyp/2)-1, 
									 y0c,y1c,x0c,x1c, nxc);//juste /2 -1 car point prol au milieu des 4 tjrs pair
                        up[np] += 0.25 * uc[ind]; 
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
                        up[np] += 0.25 * uc[ind];
                        
                    }
                   /* else{
                        up[*np] += 0.25 * computeBound((ixp+1+1)*hp,(iyp+1-1)*hp);
                    }*/
                    //coin haut gauche
                    if (check_nw(ixp,iyp,y0p,y1p,x0p,x1p,nxp, nyp) ){ //cond gauche haut
                        ind = indice((ixp/2 - 1),(iyp/2), 
									 y0c,y1c,x0c,x1c, nxc);
                        up[np] += 0.25 * uc[ind];
                    }
                    /*else{
                        double bound = computeBound((ixp+1-1)*hp,(iyp + 1+1)*hp);
                        up[*np] += 0.25 * bound;
                        
                    }*/
                    //coin droit haut
                    
                    if (check_ne(ixp,iyp,y0p,y1p,x0p,x1p,nxp, nyp) ){ //coin droit
                        ind = indice((ixp/2),(iyp/2), y0c,y1c,x0c,x1c, nxc);
                        up[np] += 0.25 * uc[ind];

                    }
                    /*else{
                        up[*np] += 0.25 * computeBound((ixp+1+1)*hp,(iyp+1+1)*hp);
                    }*/
                }

				//impair pair => somme haut + bas
				else if (iyp % 2 == 0 && ixp % 2 == 1){
                    
                    //somme bas
                    if (check_sud(ixp,iyp,y0p,y1p,x0p,x1p,nyp) ){
                       
                        ind = indice((ixp/2),(iyp/2)-1,
									 y0c,y1c,x0c,x1c, nxc);
                        up[np] += 0.5 * uc[ind];
                    }
                    /*else{
                        up[*np] += 0.5 * computeBound((ixp+1)*hp,(iyp+1-1)*hp);
                    }*/
                    //somme haut
                    if(check_nord(ixp,iyp,y0p,y1p,x0p,x1p,nyp)){
                        ind = indice((ixp/2),(iyp/2), y0c,y1c,x0c,x1c, nxc);
                        up[np] += 0.5 * uc[ind]; //pas de nx/2+1 car point uc[nc] deja ligne du haut
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
                        up[np] += 0.5 * uc[nc - 1];
                    }
                    /*else{
                        up[*np] += 0.5 * computeBound((ixp+1-1)*hp,(iyp+1)*hp);
                    }*/

                    //somme droit
                    //si pas a gauche d'un bord
                    if (check_est(ixp,iyp,y0p,y1p,x0p,x1p,nxp) ){
                        up[np] += 0.5 * uc[nc]; //nc car deja +1 car type point precedent + 1
                    }
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
				i == iter;
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


int preInitialization(double *u0){
    //at the begining of mg method u0 must be set at 0
    //at lower level (ci), initialisation at 0 will be done 
    // at the begining of EACH iteration
    //for the record we valid that only u need to be initialized
    //other will be set (r[i] = b[i] in computeres
    //and in gaussres d is x and x[i] = b[i]) 
	//arbitrary u0
	for (int i=0; i < globVal.vectStart[1]; i++){
		u0[i] = 0; 
	}
	/*
        if (EXPLICIT){
		printf("\n initial res : %lf\n", computeResNorm(
            *n0, ia0, ja0, a0, b, u0, r0
        ));
	}
    */
	return 0;
}

//opti while ? rec ?
int reInitialization(double *ul){
    //ci are all set to 0 (ci are ul except top level)
    for (int i = globVal.vectStart[1]; i < globVal.vectStart[LEVELMAX+1]; i++){
		ul[i] = 0; 
	}
    return 0;
}

int allocGrids(int **nl, int **ial,
               int **jal, double **al, double **bl,
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
	*bl = (double*)malloc(globVal.vectStart[LEVELMAX+1] * sizeof(double));


    if (*nl == NULL || *bl == NULL || *ial == NULL || *jal == NULL || 
        *al == NULL ||*dl == NULL || *rl == NULL || *ul == NULL){
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n");
        return 1;
    }
    if (EXPLICIT){
        printf("\nMemory allocated for level 0 to level %d\n", LEVELMAX);
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








 //TOOLS//

int addVect(int *n , double *v1, double *v2){
	for (int i = 0; i < *n; i++){
		v1[i] += v2[i];
	}
	return 0;
}

double computeBound(double x, double y){
    //compute the value of a point on the bound of the domain
	//with a function BOUND defined in the macros
    return BOUND(x,y);
}

double computeResNorm(int *n, int *ia, int *ja, double *a, 
						double *b, double *u,  double *r){

	double rn = 0;

	computeRes(n,ia,ja,a,u,b,r);
	for (int i = 0; i < *n; i++){
		rn += r[i] * r[i];
	}
	rn = sqrt(rn);
	return rn*10000000000000000;
}

int computeRes(int *n, int *ia, int *ja, double *a, 
			   		double *u, double *b, double *r){
	
	//r = b -Au
	int i = 0;
	int jai = 0;
	while (i < *n){
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
    printf("\nmi = %d,", m);
    for (int l = 1; l <= LEVELMAX; l+=1){
        (*mLevels)[l] = m/(1 << l) + 1;
        printf(" %d,", (*mLevels)[l]);
    }
    printf("\n");
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