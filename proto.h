#ifndef PROTO
#define PROTO


double mytimer();
double dnrm2_(int*, double*, int*); /* norme euclidienne dans BLAS */

int solve_umfpack(int n, int *ia, int *ja, double *a, 
                  double *b, double *x);

//prob.c
void computeParamTop(int m, double *h, double *invh2, int *y0,
                     int *y1, int *x0, int *x1, int *nx, int *n);
void computeParamLevel(int m, int level, double *hl, double *invh2l, 
                        int *y0l, int *y1l, int *x0l, int *x1l,
                        int *nxl, int *nl,int *nnzl);

int prob(int m, int *n, int **ia, int **ja, double **a, double **b);
double computeBound(double x, double y);
int in_hole(int ix, int iy, int y0, int y1, int x0, int x1);
int check_nord(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_sud(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_west(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_est(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_nw(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_ne(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_sw(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_se(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int indice(int ix,int iy, int y0, int y1, int x0, int x1, int nx);
int on_bound(int ix, int iy, int m);
void computeHole(int *y0, int *y1, int *x0, int *x1, int m);
//grid_corr.c
int restrictR(int level, double *rp, double **rc, int m, int *nc);
int prolongR(int level, double **up, double *uc, int m, int *np);
int probMg(int m, int level, int *nl, int *ial, int *jal, double *al, double *b);
int allocGridLevel(int m, int level, int *nl, int **ial,
                     int **jal, double **al, double **bl);
//plot.c
void plot_static(double *x, int m, int level);

//method.c
int mg_method(int iter, int levelMax, int m);
int tg_rec(int level, int levelMax, int m, int mu1, int mu2, 
		   int *nl, int mu1, int mu2, int **ial, int **jal,
		   double **al, double *b, double **ul, double **rl, double **dl);
int allocGrids(int m, int levelMax, int **nl, int ***ial,
               int ***jal, double ***al, double ***bl,
			   double ***dl, double ***rl, double ***ul);
int allocLevel(int m, int level, int *nl, int ***ial,
                     int ***jal, double ***al, double ***bl,
                     double ***dl, double ***rl, double ***ul);

int stationaryIter(int iter, int n, int *ia, int *ja, double *a,
					 double *b, double *u, int forward);
int gaussResL(int n , int *il, int *jl, double *l, double *x, double *b);

//ds prob a changer
int computeRes(int n, int *ia, int *ja, double *a, double *u, double *b, double *r);
double computeResNorm(int n, int *ia, int *ja, double *a, double *u, double *b, double *r);


#endif