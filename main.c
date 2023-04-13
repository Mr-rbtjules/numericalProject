#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "proto.h"
/* Declaratiopn des prototypes */


/* Fonction main */

int main(int argc, char *argv[])
{


mg_method(1000, 2, 120);




	

   /*int m = 100;
   int levelMax = 1;
   int startLevel = 1;

   int *nl;
   int **ial = NULL;
   int **jal = NULL; 
   double **al = NULL;
   double **bl = NULL;
   double **ul = NULL;
   double **dl = NULL;
   double **rl = NULL;

   
   free(nl);
	for (int j = 0; j <= levelMax; j++){
		free(ial[j]);
		free(jal[j]);
		free(al[j]);
		free(rl[j]);
		free(dl[j]);
		free(ul[j]);
	}*/
	
 
   
/*
  int m = 26; //  >= 13 et impaire pour la restriction et conserver distance par rapport au bords
  int level = 1; //level 1 min 26 , 2 52,3 104, level max 6-7

  int n;
  int *ia = NULL;
  int *ja = NULL; 
  double *a = NULL;
  double *b = NULL;
  double *x = NULL;
  double t1, t2;

allocGridLevel(m, level, &n, &ia, &ja, &a, &b);

  if (probMg(m, level, &n, ia, ja, a, b))
     return 1;
  printf("\nPROBLEM: ");
  printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n] );
   
  x = malloc(n * sizeof(double));
  if ( x == NULL ) {
  	printf("\n ERREUR : pas de mémoire pour vecteur des solutions\n\n");
        return 1;
  }


  t1 = mytimer();
  if( solve_umfpack(n, ia, ja, a, b, x) )
     return 1;
  t2 = mytimer();
  printf("\nTemps de solution (CPU): %5.1f sec\n",t2-t1);


   
  free(ia); free(ja); free(a); free(b); free(x);
  */
  return 0;
}

