#ifndef PROTO
#define PROTO


double mytimer();
double dnrm2_(int*, double*, int*); /* norme euclidienne dans BLAS */

int solve_umfpack(int n, int *ia, int *ja, double *a, 
                  double *b, double *x);

//prob.c
int prob(int m, int *n, int **ia, int **ja, double **a, double **b);
double computeBound(double x, double y);
int in_hole(int ix, int iy, int y0, int y1, int x0, int x1);
int check_nord(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_sud(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_west(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_est(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_nw(int ixc, int iyc, int y0, int y1, int x0, int x1, int nx);
int check_ne(int ixc, int iyc, int y0, int y1, int x0, int x1, int nx);
int check_sw(int ixc, int iyc, int y0, int y1, int x0, int x1, int nx);
int check_se(int ixc, int iyc, int y0, int y1, int x0, int x1, int nx);
int indice(int ix,int iy, int m, int y0, int y1, int x1, int x0);


//grid_corr.c
int restrictR(double **r, double **rc, int m, int *n);

#endif