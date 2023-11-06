#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "proto.h"

//must be the form of boundaries \ holes



int main() {
    //int* numbers = extract_numbers(TEST, &count);
    
    initGlobVal();  
    int nl;
    int *ial, *jal, *al, *bl; 
    double hl, invh2l;
    int x0l,x1l,y0l,y1l, nxl, nnzl;
    computeParamLevel(globVal.m, &hl,&invh2l,&x0l,&x1l, &y0l, &y1l, &nxl, &nl, &nnzl);
    //probMg(m, &nl, ial,jal, al, bl);
    
    freeGlobVal();
   return 0;
}