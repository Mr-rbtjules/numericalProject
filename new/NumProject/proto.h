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



//changer type de retour

// RESOL//

int mg_method(int iter);
int mg_iter(int iter, int levelMax, int m, int mu1, int mu2, int *nl, int *ial,
             int *jal, double *al, double *rl, double *ul, double *dl, double *b);

// PROB //

int probMg(int m, int *n, 
		   int *ia, int *ja, double *a, double *b);
int initialization(int *n0, int *ia0, int *ja0, double *a0,
						 double *b, double *r0, double *u0);
int allocGrids(int **nl, int **ial,
               int **jal, double **al, double **bl,
			   double **dl, double **rl, double **ul);
void computeParamLevel(int m, double *h, double *invh2, int *x0,
						int *x1, int *y0, int *y1, int *nx,
						int *ny, int *n, int *nnz);
void getNnz(int nx, int perUnit, int x0, int x1, 
					int y0, int y1, int *n, int *nnz);
int on_bound(int px, int py, int mx, int my);
int in_hole(int ix, int iy, int y0, int y1, int x0, int x1);
int check_nord(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_sud(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_west(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_est(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int indice(int ix,int iy, int y0, int y1, int x0, int x1, int nx);


// TOOLS //
double computeBound(double x, double y);
double computeResNorm(int n, int *ia, int *ja, double *a,
                             double *b,double *u, double *r);
int computeRes(int n, int *ia, int *ja, double *a,
                    double *u, double *b, double *r);
void printVect(void *vect, int size, int type);
void plot_res(double *r, int level);

// Global variable //
void initGlobVal();
void extractDomain(const char* str_domain, int* count, int **domain);
int correctM(int *domain, int m);
int initMlevels(int m, int **mLevels);
int initIndex( int **vectStart, int **matStart);
void freeGlobVal();


// OTHER //
double mytimer();
int solve_umfpack(int n, int *ia, int *ja, double *a, 
                  double *b, double *x);


//not used
int allocProb(int m, int *n, int **ia, int **ja, 
     		  double **a, double **b, double **u, double **r);


#endif