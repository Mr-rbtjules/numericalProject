#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "proto.h"

//must be the form of boundaries \ holes



int main() {
    //int* numbers = extract_numbers(TEST, &count);
    
    initGlobVal();  
    
    int *n;
    int *nnz;
    int *ia = NULL;
    int *ja = NULL; 
    double *a = NULL;
    double *b = NULL;
    double *x = NULL;
    double *r = NULL;
    double *d;
    
    double t1, t2;

    //allocProb(globVal.m, &n, &ia, &ja, &a, &b, &x, &r);
    int level = 1;
    allocGrids(globVal.m, 1, &n, &nnz, &ia, &ja, &a, &b, &d, &r, &x);
    
    if (probMg(globVal.m/2 +1, n+1, ia+n[0]+1, ja + nnz[0], a+ nnz[0], b + n[0]))
        return 1;
        //ia+n[0] =? nnz[0]
    printf("\nPROBLEM: ");
    printf("m = %5d   n = %8d  nnz = %9d\n", globVal.m/2 + 1, *(n+1), *(ia +n[0]+1 +n[1] ));


    t1 = mytimer();
    if( solve_umfpack(*(n+1), ia + n[0] + 1, ja + nnz[0], a+ nnz[0], b+ n[0], x+n[0]) )
        return 1;
    t2 = mytimer();
    printf("\nTemps de solution umpfpack(CPU): %5.1f sec\n",t2-t1);

    printf("\n res = %lf\n",computeResNorm(*(n+1), ia + n[0] + 1, ja + nnz[0],
     a + nnz[0], b + n[0], x + n[0], r + n [0]));

    free(n); free(nnz);free(ia); free(ja); free(a); free(b); free(x);
    free(d); free(r);

    freeGlobVal();
   return 0;
}