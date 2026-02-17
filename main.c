#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "proto.h"
#include <float.h>



 
int main() {
    
    
    int ms[9] = {1537, 769, 385, 193, 97, 49, 25, 13, 7};
    int iter = 15;
    int m = ms[3];
    double *res = malloc(iter * sizeof(double));
      
    //conjugateGradientCSR();
    mg_method(m, iter, 5, 4, 0, 0, res); 
    
    //CGmethod(iter, m, 7, res);
    plotIter(res, iter, m);

    free(res);
    
    

    
/*

probleme potentiel : 
-double x = (ix+ 1)*h;
double y = (iy + 1 -1)*h;
double bound = computeBound(x, y);

-enleve -invh2 mais alors 
*/



   return 0;
}