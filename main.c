#include <stdio.h>
#include <stdlib.h>

/* Declaratiopn des prototypes */

double mytimer();
double dnrm2_(int*, double*, int*); /* norme euclidienne dans BLAS */
prob(int m, int *n, int **ia, int **ja, double **a, double **b);
/* Fonction main */

int main(int argc, char *argv[])
{

  /* déclarer les variables */

  int m = 100;
  int n, *ia, *ja; 
  double *a, *b, *x;
  double t1, t2;

  /* générér le problème */


   while( (m-1)%6 != 0 && m > 7){
        m++;
   }
  if (prob(m, &n, &ia, &ja, &a, &b))
     return 1;
  printf("\nPROBLEM: ");
  printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n] );


   
  /* allouer la mémoire pour le vecteur de solution */

  x = malloc(n * sizeof(double));
  if ( x == NULL ) {
  	printf("\n ERREUR : pas de mémoire pour vecteur des solutions\n\n");
        return 1;
  }

  /* résoudre et mesurer le temps de solution */

  t1 = mytimer();
  if( solve_umfpack(n, ia, ja, a, b, x) )
     return 1;
  t2 = mytimer();
  printf("\nTemps de solution (CPU): %5.1f sec\n",t2-t1);

  /* libérér la mémoire */
  free(ia); free(ja); free(a); free(b); free(x);
  return 0;
}

