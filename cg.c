#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "proto.h"


int CGmethod(int iter, int m, int levelMax, double *res){


    //INIT PRECONDITION PART
    initGlobVal(m,levelMax);
	//initit memory and pointers
	//compute all the coarse matrix and nl
    int *nl = NULL;

    //liste : chaque elem pointe vers 1 matrice 
    //d'un certain level (ici matrice == liste)
	int *ial = NULL; 
	int *jal = NULL;
	double *al = NULL;                                                         
	double *rl = NULL;
	double *ul = NULL;
	double *dl = NULL;
	
    //'b' conceptuel, cad b au niveau 0 mais ri au niveau inf
	double *bl = NULL;  
	//we pass the adress for each level
	allocGrids(levelMax, &nl, &ial, &jal, &al, &bl, &dl, &rl, &ul);
    
	//precomputation of all Ac and bc

	for (int l = 0; l <= levelMax; l++){
		//adresses of each vector is simply the pointer shifted
        // bl not shifted because only the firs
		probMg(
            globVal.m[l], 
            nl + l, 
            ial + globVal.vectStart[l] + l,
            jal + globVal.matStart[l], 
            al + globVal.matStart[l],
            bl
        );
	}
	
	int startLevelTg = 0; 

	//u0 = 0 
	preInitialization(ul);

    //INIT CONJUGUATE GRADIANT PART
    double *uCG = malloc(*nl * sizeof(double));
    double *rCG = malloc(*nl * sizeof(double));
    double *dCG = malloc(*nl * sizeof(double)); 
    double *AdCG = malloc(*nl * sizeof(double));
    double *bCG = malloc(*nl * sizeof(double));
    double save_old;
    double save_new; 
    double alpha;
    double beta;
    double dAd;

    
    

    int iterPrecond = 2;
    int mu = 4;
    int mode = 0;
    int cycle = 0;

    
    //uCG0=0
    copy(nl, bl, bCG); // comme ça on conserve b ds qq chose de fixe
    preInitialization(uCG);
    computeRes(nl, ial, jal, al, uCG, bCG, rCG);
    //u0
    copy(nl, rCG, bl);
    mg_iter( iterPrecond, levelMax, globVal.m[0], 
            mu, mu, cycle, nl, ial ,jal,al, rl,ul,dl,bl, NULL); // dans ul on a z
    //saveold - r0z0
    //scalProd(nl, rCGl, ul, &save_old);
    //d0 = u0
    
    copy(nl, ul, dCG);

    scalProd(nl, rCG, ul, &save_old);

    if (LOAD){
		printf("\n Load percentage : \n");
	}
	//start iterations of CG method
	double t1 = mytimer();
    
    double *nres = malloc(*nl * sizeof(double));

    for (int i = 0; i < iter; i++){
        if (LOAD){
            
			printf("   %.3f   ", (double)(i+1)/iter);
        	fflush(stdout);
			printf("\r");
		}
        //Ad = A*d
        
        multCsrVector( nl, ial, jal, al, AdCG, dCG); 

        //save new = num alpha and denom beta+1
        //rappel rcGl est le bl de mgmethod precond
        //save new = <r , B-1 r>
        

        //alpham = savenew/ <d,Ad>
        scalProd(nl, dCG, AdCG, &dAd);
        alpha = save_old/dAd;


        //u = u+ alpha*d
        addVectProd(nl, &alpha, uCG, dCG);
        
        //r = r - alpha*d
        subVectProd(nl, &alpha, rCG, AdCG);

        copy(nl, rCG, bl);
        preInitialization(ul);
        mg_iter( iterPrecond, levelMax, globVal.m[0], 
            mu, mu,0, nl, ial ,jal,al, rl,ul,dl,bl, NULL); 
        scalProd(nl, rCG, ul, &save_new);
        //check conv

        //betam+1 save new/saveold
        beta = save_new/save_old; 

        //ul+1 = B-1rm+1
        //part de quel ul+1 ? tester les 2
        //preInitialization(ul);
        
        
        //dm+1 = ul+1 + beta*dm
        dSum( nl, &beta, dCG, ul);


        save_old = save_new;

        //peut utiliser rl pour mettre ce qu'on veut vu que reset 
        //a chaque debut de precond
        
        res[i] = computeResNorm(
            nl, ial, jal, al, bCG, uCG, nres
        );
        
    }
    free(nres);
    
    
    double t2 = mytimer();
    if (CHRONO){
        printf("\nTemps de solution CG method (CPU): %5.1f sec\n",t2-t1);   
    }


    

    //top level so no shift
    
    printf("\n final res : %.16f\n", computeResNorm(
            nl, ial, jal, al, bCG, uCG, rCG
    ));
    

    if (PLOT){
        plotCycle(levelMax, cycle);
        plot_res(rCG + globVal.vectStart[startLevelTg], startLevelTg);
    }
    

    free(nl);
    free(ial);
    free(jal);
    free(al);
    free(rCG);
    free(dl);
    free(ul);
    free(rl);

    free(uCG);
    free(bCG);
    free(dCG);
    free(AdCG);

    freeGlobVal();

	return 0;
}













// CSR matrix-vector multiplication
void csrMatVecMult(int* ia, int* ja, double* a, double* x, double* result, int N) {
    for (int i = 0; i < N; i++) {
        result[i] = 0;
        for (int j = ia[i]; j < ia[i + 1]; j++) {
            result[i] += a[j] * x[ja[j]];
        }
    }
}

double dotProduct(double* vec1, double* vec2, int N) {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += vec1[i] * vec2[i];
    }
    return sum;
}

// Function to add two vectors: result = vec1 + alpha * vec2
void vecAdd(double* vec1, double* vec2, double* result, double alpha, int N) {
    for (int i = 0; i < N; i++) {
        result[i] = vec1[i] + alpha * vec2[i];
    }
}

// Function to subtract two vectors: result = vec1 - alpha * vec2
void vecSub(double* vec1, double* vec2, double* result, double alpha, int N) {
    for (int i = 0; i < N; i++) {
        result[i] = vec1[i] - alpha * vec2[i];
    }
}

// Conjugate Gradient Method adapted for CSR format
void conjugateGradientCSR(int iter, int* ia, int* ja, double* a, double* b, double* x, int N, double *res) {
    double *r, *p, *Ap;
    r = (double*) malloc(N * sizeof(double));
    p = (double*) malloc(N * sizeof(double));
    Ap = (double*) malloc(N * sizeof(double));

    // Initial guess x = 0
    for (int i = 0; i < N; i++) {
        x[i] = 0;
    }

    // r = b - Ax
    csrMatVecMult(ia, ja, a, x, Ap, N);
    vecSub(b, Ap, r, 1, N);

    // p = r
    for (int i = 0; i < N; i++) {
        p[i] = r[i];
    }

    double rs_old = dotProduct(r, r, N);
    for (int i = 0; i < iter; i++){
        csrMatVecMult(ia, ja, a, p, Ap, N);
        double alpha = rs_old / dotProduct(p, Ap, N);
        vecAdd(x, p, x, alpha, N);
        vecAdd(r, Ap, r, -alpha, N);
        double rs_new = dotProduct(r, r, N);
        if (sqrt(rs_new) < 1e-10) {
			printf(" i = %d\n", N);
            break;
        }
        vecAdd(r, p, p, rs_new / rs_old, N);
        rs_old = rs_new;
    }

    free(r);
    free(p);
    free(Ap);
}