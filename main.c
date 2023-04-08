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

   int m2 = 26;
   int level2 = 1;

   int n2;
   int *ia2 = NULL;
   int *ja2 = NULL; 
   double *a2 = NULL;
   double *b2 = NULL;
   double *x2 = NULL;
   double *uc2 = NULL;
   
   allocGridLevel(m2, level2, &n2, &ia2, &ja2, &a2, &b2);
   probMg(m2, level2, &n2, ia2, ja2, a2, b2);

   x2 = malloc(n2 * sizeof(double));
   
   solve_umfpack(n2, ia2, ja2, a2, b2, x2);
   double *r2 = malloc(n2*sizeof(double));
   double rn2 = computeResNorm(n2,ia2,ja2,a2,x2, b2,r2);
   printf(" \nnorme résidu r:%.16g\n", rn2);
   free(r2);

   /*printf("base :\n");
   for( int i = 0; i < n2; i++){
         printf(" %lf ", x2[i]);
   }*/
   //plot_static(x2, m2, level2);

   /*
   int nc2;
   restrictR(level2, x2, &uc2, m2, &nc2);
   */
   /*printf("restricted :\n");
   for( int i = 0; i < nc2; i++){
         printf(" %lf ", uc2[i]);
   }*/
   //plot_static(uc2, m2, level2+1);


   //plot
   
   double *up = NULL;
   int np;
   prolongR(level2, &up, x2, m2, &np);

   /*
   

   */

   
   printf("prolonged :\n");
   for( int i = 0; i < np; i++){
         printf(" %lf ", up[i]);
   }
   printf("\n");
   plot_static(up, m2, 0);


   int m1 = m2;
   int level1 = 0;

   int n1;
   int *ia1 = NULL;
   int *ja1 = NULL; 
   double *a1 = NULL;
   double *b1 = NULL;
   
   allocGridLevel(m1, level1, &n1, &ia1, &ja1, &a1, &b1);
   probMg(m1, level1, &n1, ia1, ja1, a1, b1);

   double *r3 = malloc(n1*sizeof(double));
   double rn3 = computeResNorm(n1,ia1,ja1,a1,up, b1,r3);
   printf(" \nnorme résidu r:%.16g\n", rn3);
   free(r3);



   free(up);

   free(uc2);
 free(ia2);
  free(ja2); 
  free(a2);
   free(b2); 
   free(x2);


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

