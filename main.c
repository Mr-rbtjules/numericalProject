#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "proto.h"
/* Declaratiopn des prototypes */


/* Fonction main */

int main(int argc, char *argv[])
{
/*
   double a[] = {4,3,6,1,9,7,1,5};
   int ia[] = {0,3,5,7,8};
   int ja[] = {0,2,3,2,3,0,1,1};
   int n = 4;

   printA(n,ia,ja,a);
*/


   //initit memory and pointers
	//compute all the coarse matrix and nl
	int **ial; //liste donc chaque elem pointe vers 1 matrice d'un certain level (ici matrice == liste)
	int **jal;
	double **al;
	double **rl;
	double **ul;
	double **dl;
	
	double **bl; //attention aux autres b prc pt on les utilise pour autre
	
	int *nl;
   int levelMax = 0;
   int m = 1000;

   int iter = 500;
   double tol = 1;

	allocGrids(m, levelMax, &nl, &ial, &jal, &al, &bl, &dl, &rl, &ul);

	//precomputation of all Ac and bc

	for (int l = 0; l <= levelMax; l++){
		
		probMg(m, l, nl, ial[l], jal[l], al[l], bl[l]);
		
	}
	

	//start iterations of the multigrid cycle (tg_rec)
	int startLevelTg = 0; 

	initialization(startLevelTg, nl, ial, jal, al, bl, dl, rl, ul);
   
   
	

	free(nl);
	for (int j = 0; j <= levelMax; j++){
		free(ial[j]);
		free(jal[j]);
		free(al[j]);
		free(rl[j]);
		free(dl[j]);
		free(ul[j]);
	}
	

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

