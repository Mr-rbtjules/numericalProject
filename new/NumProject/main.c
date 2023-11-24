#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "proto.h"

//must be the form of boundaries \ holes



int main() {
    //int* numbers = extract_numbers(TEST, &count);
    
    initGlobVal();

    mg_method(1);
    
    freeGlobVal();

/*    initGlobVal();  
    //pour remplacer nnz on peut utiliser ia dernier elem ? 
    int *nl;
    int *ia = NULL;
    int *ja = NULL; 
    double *a = NULL;
    double *b = NULL;
    double *x = NULL;
    double *r = NULL;
    double *d;
    
    double t1, t2;

    //allocProb(globVal.m, &n, &ia, &ja, &a, &b, &x, &r);
    int level = 0;
    allocGrids(&nl, &ia, &ja, &a, &b, &d, &r, &x);
    //globVal.m, n, ia, ja, a, b
    //on decide de pointer la ou le level vecteur commence (translation du pointer)
    if (probMg(
            globVal.m[level], 
            (nl + level), 
            ia+globVal.vectStart[level] + level, 
            ja + globVal.matStart[level],
            a + globVal.matStart[level], 
            b
        )){
        return 1;
    }
    
    printf("\nPROBLEM: ");
    printf("m = %5d   n = %8d  nnz = %9d\n", globVal.m[level], *(nl+level), *(ia + globVal.vectStart[level] + *(nl+level) + level));
    


    t1 = mytimer();
    //*(n+1), ia + n[0] + 1, ja + nnz[0], a+ nnz[0], b+ n[0], x+n[0]
    if( solve_umfpack(
            *(nl + level),
            ia+globVal.vectStart[level] + level,
            ja + globVal.matStart[level], 
            a + globVal.matStart[level], 
            b, 
            x + globVal.vectStart[level]
        )){
        return 1;
    }
    t2 = mytimer();
    printf("\nTemps de solution umpfpack(CPU): %5.1f sec\n",t2-t1);

    printf("\n res = %lf\n",computeResNorm(
            *(nl+level), ia + globVal.vectStart[level]+ level, 
            ja + globVal.matStart[level], a + globVal.matStart[level], 
            b, x + globVal.vectStart[level], r +globVal.vectStart[level]
        )
    );

    free(nl); free(ia); free(ja); free(a); free(b); free(x);
    free(d); free(r);

    freeGlobVal();*/
   return 0;
}