#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "proto.h"

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
void conjugateGradientCSR(int* ia, int* ja, double* a, double* b, double* x, int N) {
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
    for (int i = 0; i < N; i++) {
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