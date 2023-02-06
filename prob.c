/* les modifications sont marquées avec "<--" */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define COORD_X0 1.0 //col bord gauche
#define COORD_X1 2.5 // col bord droit
#define COORD_Y0 1.5 //ligne bord bas
#define COORD_Y1 2.0 //ligne bord haut



void initBaseVal(int m){
    while( (m-1)%6 != 0 && m > 7){
        m++;
    }
    int nx = m-2;
    //nb de points sur la largeur du trous
    int p = 1 + m/6; //== nb de points entre y1 et y0
    // "" sur la longueur
    int q = 1 + m/2;
    //plaque hors bords et trous
    int dim = nx * nx - (p * q);
    //nombre d'element non nul pour membrane de base
    int nnz = 5 * nx * nx - 4 * nx ; 
    //nb de points dans le trou:(compliqué a comprendre sans shema)
    int trous = (5 * (p-2) * (q-2) + 4 * 2 * (p-2) + 4 * 2 * (q-2) 
                + 3 * 4 * 1 + 1 * 2 * p + 1 * 2 * q) + 2 * p + 2 * q;
    nnz -= trous;

    double h = 3.0/(double)(m-1);


    int y0 = (int)ceil(COORD_Y0/h) - 1;
    int y1 = (int)ceil(COORD_Y1/h) - 1;
    int x0 = (int)ceil(COORD_X0/h) - 1;
    int x1 = (int)ceil(COORD_X1/h) - 1;
}

double computeBound(double x, double y){//mettra ix*h +-1
    return exp(sqrt(x*x + y*y));
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
    int  nnz, ix, iy, ind, nx;
    double xx, yy, h, invh2;

    
    int nx = m-2;
    //nb de points sur la largeur du trous
    int p = 1 + m/6; //== nb de points entre y1 et y0
    // "" sur la longueur
    int q = 1 + m/2;
    //plaque hors bords et trous
    int dim = nx * nx - (p * q);
    //nombre d'element non nul pour membrane de base
    int nnz = 5 * nx * nx - 4 * nx ; 
    //nb de points dans le trou:(compliqué a comprendre sans shema)
    int trous = (5 * (p-2) * (q-2) + 4 * 2 * (p-2) + 4 * 2 * (q-2) 
                + 3 * 4 * 1 + 1 * 2 * p + 1 * 2 * q) + 2 * p + 2 * q;
    nnz -= trous;

    double h = 3.0/(double)(m-1);


    int y0 = (int)ceil(COORD_Y0/h) - 1;
    int y1 = (int)ceil(COORD_Y1/h) - 1;
    int x0 = (int)ceil(COORD_X0/h) - 1;
    int x1 = (int)ceil(COORD_X1/h) - 1;
    
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


    int skip_sud = 0;
    int skip_nord = 0;
    int ecart = 0;

    //Calcul des constantes (indice des points du bord du trou)
    //ceil arrondis pour ne pas avoir par exemple (int)6.99 = 6 
   

    //calcul du nb de point sur le bord le plus long
    //quand on passe par odre lexicographique et qu'on
    //arrive au trou les points voisin ne sont plus voisins
    //dans la matrice et sont 'decalé' il y a un ecart
    for (ix = 0; ix < nx; ix++){
        if ( ix >= x0 && ix<= x1 ){ //nombre de points sur un bord interne(du trou)
                ecart += 1;
        }
    }
    int nnz_save = nnz;
    
    printf(" nnz %d  ecart %d nx %d\n", nnz, ecart, nx);
    //passage ligne suiv(plaque complete)

    nnz = 0;
    for (iy = 0; iy < nx; iy++) {
        for (ix = 0; ix < nx; ix++) {
            /* numéro de l'équation */
            //initialisation du retard a ajouter à cause du bord 
            //(depend d'ou on se trouve sur la plaque)
            if ( iy >= y0-1 && ix >= x1){//ecart reste jusque cond suivante appliquée, quand ? qd sorti des alentours du trous
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
            //exclu interieur et bord du trou
            if(iy > y1  || iy < y0  || ix < x0 || ix > x1){
                //marquer le début de la ligne suivante dans le tableau 'ia'
                (*ia)[ind] = nnz;
                (*b)[ind] = 0.0;
            
                //replissage de la ligne : voisin sud //ui-1
                // + verification si pas au dessus d'un bord
                if (iy > 0 && ( iy-1 != y1 || ix < x0 || ix > x1) ){
                    (*a)[nnz] = -invh2;
                    (*ja)[nnz] = ind - nx + skip_sud;
                    
                    //-nx car on regarde delui d'en bas(shema)
                    //+skip_sud car comme on a passe des points(trous)
                    // nx ramene trop loin en arrière
                    nnz++;
                }
                else{
                    (*b)[ind] += computeBound(ix*h, (iy-1)*h); 
                }

                //replissage de la ligne : voisin ouest 
                //si pas a droite d'un bord
                if (ix > 0 && ( ix-1 != x1 || iy > y1 || iy < y0 )){
                    (*a)[nnz] = -invh2;
                    (*ja)[nnz] = ind - 1;
                    nnz++;
                }
                else{
                    (*b)[ind] += computeBound((ix -1)*h, iy*h);

                }

                // replissage de la ligne : élém. diagonal
                (*a)[nnz] = 4.0*invh2;
                (*ja)[nnz] = ind;
                nnz++;

                // replissage de la ligne : voisin est
                //si pas a gauche d'un bord
                if (ix < nx - 1 && ( iy > y1 || iy < y0 || ix != x0 - 1 ) ){
                    (*a)[nnz] = -invh2;
                    (*ja)[nnz] = ind + 1;
                    nnz++;
                }
                else{
                    (*b)[ind] += computeBound((ix+1)*h, iy*h);
                }

                // replissage de la ligne : voisin nord
                //si pas en dessous d'un bord
                if (iy < nx - 1 && ( iy +1  != y0 || ix < x0 || ix > x1 ) ){
                        (*a)[nnz] = -invh2;
                        (*ja)[nnz] = ind + nx - skip_nord;
                        nnz++;
                }
                else{
                    (*b)[ind] += computeBound(ix*h, (iy+1)*h);
                }
                // numéro de l'équation
                ind += 1; 
            }
        }
    }

    if (nnz == nnz_save){
        (*ia)[ind] = nnz;
    }
    else{
        printf("Error nnz != nnz\n");
        return 1;
    }

    /* dernier élément du tableau 'ia' */
    (*ia)[ind + 1] = nnz;

    /* retour de fonction habituel */
    return 0;
}

