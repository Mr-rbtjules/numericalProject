#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 4 // Size of the matrix A
#define MAX_ITER 1000 // Maximum number of iterations for Lanczos

void matvec(double A[N][N], double *x, double *result) {
    for (int i = 0; i < N; i++) {
        result[i] = 0.0;
        for (int j = 0; j < N; j++) {
            result[i] += A[i][j] * x[j];
        }
    }
}

void lanczos(double A[N][N], int max_iter, double *lmin, double *lmax) {
    double *v = (double *)malloc(N * sizeof(double));
    double *w = (double *)malloc(N * sizeof(double));
    double alpha, beta, old_beta;

    for (int i = 0; i < N; i++) {
        v[i] = rand() % 10; // Random starting vector
    }

    *lmin = 1e10;
    *lmax = -1e10;

    beta = 0.0;

    for (int j = 0; j < max_iter; j++) {
        matvec(A, v, w);

        alpha = 0.0;
        for (int i = 0; i < N; i++) {
            alpha += v[i] * w[i];
        }

        for (int i = 0; i < N; i++) {
            w[i] = w[i] - alpha * v[i];
            if (j > 0) {
                w[i] = w[i] - beta * v[i];
            }
        }

        old_beta = beta;
        beta = sqrt(0.0);
        for (int i = 0; i < N; i++) {
            beta += w[i] * w[i];
        }
        beta = sqrt(beta);

        // Update approximations of the eigenvalues
        if (alpha + beta > *lmax) *lmax = alpha + beta;
        if (alpha - beta < *lmin) *lmin = alpha - beta;

        for (int i = 0; i < N; i++) {
            v[i] = w[i] / beta;
        }
    }

    free(v);
    free(w);
}

int main() {
    double A[N][N] = {
        {4, 1, 0, 0},
        {1, 3, 1, 0},
        {0, 1, 2, 1},
        {0, 0, 1, 4}
    };

    double lmin, lmax;
    lanczos(A, MAX_ITER, &lmin, &lmax);

    printf("Min eigenvalue approximation: %f\n", lmin);
    printf("Max eigenvalue approximation: %f\n", lmax);

    return 0;
}
