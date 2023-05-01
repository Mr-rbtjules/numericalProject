#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "proto.h"
/* Declaratiopn des prototypes */


/* Fonction main */

int main(int argc, char *argv[])
{
   int m = 150; //  >= 13 et impaire pour la restriction et conserver distance par rapport au bords
   mg_method(1000, 3, m);
   int level = 0; //level 1 min 26 , 2 52,3 104, level max 6-7

   int n;
   int *ia = NULL;
   int *ja = NULL; 
   double *a = NULL;
   double *b = NULL;
   double *x = NULL;
   double *r = NULL;
   double t1, t2;

   allocProb(m, &n, &ia, &ja, &a, &b, &x, &r);

   if (probMg(m, level, &n, ia, ja, a, b))
      return 1;
   printf("\nPROBLEM: ");
   printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n] );
  

   t1 = mytimer();
   if( solve_umfpack(n, ia, ja, a, b, x) )
      return 1;
   t2 = mytimer();
   printf("\nTemps de solution umpfpack(CPU): %5.1f sec\n",t2-t1);

   printf("\n res = %lf\n",computeResNorm(n, ia, ja, a, b, x, r));

   free(ia); free(ja); free(a); free(b); free(x);



   return 0;
}

