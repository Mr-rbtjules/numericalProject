/* les modifications sont marquées avec "<--" */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int prob(int m, int *n, int **ia, int **ja, double **a, double **b)
/*
   But
   ===
   Génerer le système linéaire n x n 
                          
                             Au = b                                   

   qui correspond à la disrétisation sur une grille cartesienne 
   regulière m x m de l'équation de Poisson à deux dimensions
              
            d    d        d    d
         - == ( == u ) - == ( == u )  = 0     sur [0,1] x [0,1]
           dx   dx       dy   dy

  avec les conditions aux limites de Dirichlet
         
         u = 0  sur (x,0)                 , avec 0 <=  y  <= 1 ,
         u = 1  sur (0,y), (1,y) et (x,1) , avec 0 <= x,y <= 1 . 
  
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

    nx = m - 2; /* noeuds de Dirichlet ne sont pas pris en compte */
    h = 1.0/(m-1); /* pas de discrétisation */
    invh2 = (m-1)*(m-1); /* h^-2 pour L=1 */
    *n  = nx * nx; /* nombre d'inconnues */
    nnz = 5 * nx * nx - 4 * nx; /* nombre d'éléments non nuls */

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

    nnz = 0;
    for (iy = 0; iy < nx; iy++) {
        for (ix = 0; ix < nx; ix++) {
            /* numéro de l'équation */
            ind = ix + nx * iy;

            /* marquer le début de la ligne suivante dans le tableau 'ia' */
            (*ia)[ind] = nnz;

            /* calculer le membre de droite */
            (*b)[ind] = 0.0; /* rho = 0.0 */
   
            /* remplissage de la ligne : voisin sud */
            if (iy > 0)  {
                (*a)[nnz] = -invh2; /* pour D=1 */
                (*ja)[nnz] = ind - nx;
                nnz++;
            }

            /* remplissage de la ligne : voisin ouest */
            if (ix > 0)  {
                (*a)[nnz] = -invh2; /* pour D=1 */
                (*ja)[nnz] = ind - 1;
                nnz++;
            } else
               (*b)[ind] += invh2; /* conditions de Dirichlet, bord ouest */

            /* remplissage de la ligne : élém. diagonal */
            (*a)[nnz] = 4.0*invh2; /* pour D=1 */
            (*ja)[nnz] = ind;
            nnz++;

            /* remplissage de la ligne : voisin est */
            if (ix < nx - 1) {
                (*a)[nnz] = -invh2; /* pour D=1 */
                (*ja)[nnz] = ind + 1;
                nnz++;
            } else
                (*b)[ind] += invh2; /* conditions de Dirichlet, bord est */

            /* remplissage de la ligne : voisin nord */
            if (iy < nx - 1) {
                (*a)[nnz] = -invh2; /* pour D=1 */
                (*ja)[nnz] = ind + nx;
                nnz++;
            } else
               (*b)[ind] += invh2; /* conditions de Dirichlet, bord nord */
        }
    }

    /* dernier élément du tableau 'ia' */
    (*ia)[ind + 1] = nnz;

    /* retour de fonction habituel */
    return 0;
}
