#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "proto.h"

//must be the form of boundaries \ holes



int main() {
    //int* numbers = extract_numbers(TEST, &count);
    
    initGlobVal();  
    int n;
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
    allocGrids(globVal.m, level, &ia, &ja, &a, &b, &d, &r, &x);
    //globVal.m, n, ia, ja, a, b
    if (probMg(globVal.m[level], &n, ia+globVal.vectStart[level] + level, ja + globVal.matStart[level], a + globVal.matStart[level], b))
        return 1;
        //ia+n[0] =? nnz[0]
    printf("\nPROBLEM: ");
    printf("m = %5d   n = %8d  nnz = %9d\n", globVal.m[level], n, *(ia + globVal.vectStart[level] + n +1));
    


    t1 = mytimer();
    //*(n+1), ia + n[0] + 1, ja + nnz[0], a+ nnz[0], b+ n[0], x+n[0]
    if( solve_umfpack(n, ia+globVal.vectStart[level] + level, ja + globVal.matStart[level], a + globVal.matStart[level], b, x + globVal.vectStart[level]), x)
        return 1;
    t2 = mytimer();
    printf("\nTemps de solution umpfpack(CPU): %5.1f sec\n",t2-t1);

    /*printf("\n res = %lf\n",computeResNorm(*(n+1), ia + n[0] + 1, ja + nnz[0],
     a + nnz[0], b + n[0], x + n[0], r + n [0]));
*/
    free(ia); free(ja); free(a); free(b); free(x);
    free(d); free(r);

    freeGlobVal();
   return 0;
}