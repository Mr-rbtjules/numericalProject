#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "proto.h"

//must be the form of boundaries \ holes


int main() {
    //int* numbers = extract_numbers(TEST, &count);
    
    
    initGlobVal();

    int per = (globVal.m[0]-1)/(globVal.domain[1] - globVal.domain[0]);
    printf("h per unit %d\n", per);

    mg_method(13);
    //CGmethod(12);
    

    


    freeGlobVal();


   return 0;
}