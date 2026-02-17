#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "proto.h"


/// GLOBAL VARIABLES ///

void initGlobVal(int m , int levelMax){
    int count;
	globVal.levelMax = levelMax;
    extractDomain(DOMAIN, &count, &(globVal.domain));
    
    int corrm = correctM(globVal.domain, m);
    
    initMlevels(corrm, &(globVal.m), levelMax);
    
    initIndex(levelMax, &(globVal.vectStart), &(globVal.matStart));
}

void extractDomain(const char* str_domain, int* count, int **domain) {
	//array with the int describing the domain
    *domain = (int*)malloc(16 * sizeof(int));
    *count = 0;
    //on all the string
    while (*str_domain) {
		//check if we start on a number
        if (*str_domain >= '0' && *str_domain <= '9') {
			//change all the digits present until another char into 1 number
			//(1 digit or more)
            (*domain)[*count] = atoi(str_domain);
            (*count)++;

            while (*str_domain >= '0' && *str_domain <= '9') {
                str_domain++;
            }
        } else {
            str_domain++;
        }
    }
}

int correctM(int *domain, int m){
    //correct m in order to obtain discretisation 
    //matching with the surface and for multigrid =>
    //domain ex : if "[0, 6] × [0, 10]  [2, 4] × [2, 5]" -> domain = [0,6,0,10,2,4,2,5]
    // L = smallest number of cases possible
    int L = domain[1] - domain[0];
    while (m > L+1){
        L *= 2;
    }

    return L + 1;
}

int initMlevels(int m, int **mLevels, int levelMax){
    
    *mLevels = (int *)malloc((levelMax+1)*sizeof(int));
    (*mLevels)[0] = m;
    printf("\nmi = %d,", m);
    for (int l = 1; l <= levelMax; l+=1){
        (*mLevels)[l] = m/(1 << l) + 1;
        printf(" %d,", (*mLevels)[l]);
    }
    printf("\n");
    return 0;
}

int initIndex( int levelMax, int **vectStart, int **matStart){
    
    //contient les debuts de chaque vecteur (pour tt les n) et la fin 
    //du vecteur global donc nb de levels = levelmax +1 , et encore +1 
    //pour la fin de vecteur un peu comme dans ia
	*vectStart = (int*)malloc((levelMax+1 + 1) * sizeof(int));
    *matStart = (int*)malloc((levelMax+1 + 1) * sizeof(int));
    
    (*vectStart)[0] = 0;
    (*matStart)[0] = 0;
    int nnzTot = 0;
    int nTot = 0;
    for (int l = 0; l <= levelMax; l++){
        double h, invh2;
        int x0,x1,y0,y1, nx, ny, nnz;
        int n;	
		computeParamLevel(globVal.m[l], &h, &invh2, &x0,
						&x1, &y0, &y1, &nx,
						&ny, &n, &nnz);
        
        nnzTot += nnz;
        nTot += n;
        (*vectStart)[l+1] = nTot;
        (*matStart)[l+1] = nnzTot;
	}
    return 0;
}

void freeGlobVal(){
    free(globVal.domain);
    free(globVal.vectStart);
    free(globVal.matStart);
    free(globVal.m);
}
