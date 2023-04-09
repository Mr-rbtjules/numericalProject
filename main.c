#include <stdio.h>
#include <stdlib.h>

#include "proto.h"
/* Declaratiopn des prototypes */


/* Fonction main */

int main(int argc, char *argv[])
{

   
   /*
   check si restrict -> prolong on revient a meme dim
   */

   int m1 = 100;
   int level1 = 1;

   int n1;
   int *ia1 = NULL;
   int *ja1 = NULL; 
   double *a1 = NULL;
   double *b1 = NULL;
   double *x1 = NULL;
   
   allocGridLevel(m1, level1, &n1, &ia1, &ja1, &a1, &b1);
   probMg(m1, level1, &n1, ia1, ja1, a1, b1);

   x1 = malloc(n1 * sizeof(double));
   
   solve_umfpack(n1, ia1, ja1, a1, b1, x1);

   
   
   double *r1 = malloc(n1*sizeof(double));
   double rn1 = computeResNorm(n1,ia1,ja1,a1,x1, b1,r1);
   printf(" \nnorme résidu r:%.16g\n", rn1);
   free(r1);

   /*printf("base :\n");
   for( int i = 0; i < n1; i++){
         printf(" %lf ", x1[i]);
   }*/
   plot_static(x1, m1, level1);

   double *uc2 = NULL;
   int nc2;
   restrictR(level1, x1, &uc2, m1, &nc2);
  
   /*printf("restricted :\n");
   for( int i = 0; i < nc2; i++){
         printf(" %lf ", uc2[i]);
   }*/
   plot_static(uc2, m1, level1+1);



   
   double *up = NULL;
   int np;
   prolongR(level1+1, &up, uc2, m1, &np);

   

   /*
   printf("prolonged :\n");
   for( int i = 0; i < np; i++){
         printf(" %lf ", up[i]);
   }
   printf("\n");*/
   plot_static(up, m1, level1);



   free(ia1);free(ja1); free(a1); free(b1); free(x1); free(uc2); free(up);


   
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


   int plotx = 1;
   if (plotx){
      for( int i = 0; i < n; i++){
         printf(" %lf ", x[i]);
      }
   }

   int n2;
  int *ia2 = NULL;
  int *ja2 = NULL; 
  double *a2 = NULL;
  double *b2 = NULL;
  double *x2 = NULL;


 


allocGridLevel(m, 0, &n2, &ia2, &ja2, &a2, &b2);




  if (probMg(m, 0, &n2, ia2, ja2, a2, b2))
     return 1;
  printf("\nPROBLEM to coarse: ");
  printf("m = %5d   n = %8d  nnz = %9d\n", m, n2, ia2[n2] );

   int nc;  
   double *uc;
   restrictR(1, x, &uc, m, &nc);
   //res
   double *r = malloc(n*sizeof(double));
   double rn = computeResNorm(n,ia,ja,a,x, b,r);
   printf(" \nnorme résidu:%.16g\n", rn);


   



//plot
   plot_static(uc, 13, 0);

free(uc);
   
  free(ia); free(ja); free(a); free(b); free(x);
  */
  return 0;
}

