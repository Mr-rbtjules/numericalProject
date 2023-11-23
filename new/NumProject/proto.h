#ifndef PROTO
#define PROTO



typedef struct globVal_s globVal_s;
struct globVal_s{

	int *m;
	int *domain;
	int *vectStart;
	int *matStart;

};
extern globVal_s globVal;



// PROB //
int probMg(int m, int *n, 
		   int *ia, int *ja, double *a, double *b);
int mg_method(int iter);
int initialization(int *n0, int *ia0, int *ja0, double *a0,
						 double *b, double *r0, double *u0);
int allocProb(int m, int *n, int **ia, int **ja, 
     		  double **a, double **b, double **u, double **r);
int allocGrids(int **nl, int **ial,
               int **jal, double **al, double **bl,
			   double **dl, double **rl, double **ul);
void computeParamLevel(int m, double *h, double *invh2, int *x0,
						int *x1, int *y0, int *y1, int *nx,
						int *ny, int *n, int *nnz);
    
void getNnz(int nx, int perUnit, int x0, int x1, int y0, int y1, int *n, int *nnz);
double computeBound(double x, double y);
int on_bound(int px, int py, int mx, int my);
int in_hole(int ix, int iy, int y0, int y1, int x0, int x1);
int check_nord(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_sud(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_west(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_est(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int indice(int ix,int iy, int y0, int y1, int x0, int x1, int nx);

// Global variable //

void extractDomain(const char* str_domain, int* count, int **domain);
int correctM(int *domain, int m);
void initGlobVal();
int initMlevels(int m, int **mLevels);
void freeGlobVal();

double mytimer();
int solve_umfpack(int n, int *ia, int *ja, double *a, 
                  double *b, double *x);
double computeResNorm(int n, int *ia, int *ja, double *a,
                             double *b,double *u, double *r);
int computeRes(int n, int *ia, int *ja, double *a,
                    double *u, double *b, double *r);
void printVect(void *vect, int size, int type);

//plot//
void plot_res(double *r, int level);

#endif