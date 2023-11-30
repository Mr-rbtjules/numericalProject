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
             int *jal, double *al, double *rl, double *ul, double *dl, double *bl);

int tg_rec(int level, int m, int mu1,
			int mu2, int *nl, int *ial, int *jal,
		   double *al, double *bl, double *ul, double *rl, double *dl);
void FCycleLoop(int level, int m, int mu1,
			int mu2, int *nl, int *ial, int *jal,
		   double *al, double *bl, double *ul, double *rl, double *dl);
int getDown(int level, int m, int mu1,
			int mu2, int *nl, int *ial, int *jal,
		   double *al, double *bl, double *ul, double *rl, double *dl);
int getUp(int level, int m, int mu1,
			int mu2, int *nl, int *ial, int *jal,
		   double *al, double *bl, double *ul, double *rl, double *dl);
int forwardGS(int iter, int *n, int *ia, int *ja, double *a,
			    double *b, double *u, double *r, double *d);
int backwardGS(int iter, int *n, int *ia, int *ja, double *a,
				double *b, double *u, double *r, double *d);
int jacobiIter(int iter, int *n, int *ia, int *ja, double *a,
			    double *b, double *u, double *r, double *d);
int stationaryIter(int iter, int *n, int *ia, int *ja,
                   double *a, double *b, double *u, 
				   double *r, double *d, int forward);
int gaussResL(int *n , int *il, int *jl, double *l, double *x, double *b);
int gaussResU(int *n , int *iu, int *ju, double *u, double *x, double *b);
int gaussResD(int *n , int *ia, int *ja, double *a, double *x, double *b);
int solveAtCoarseLevel(int mode, int *n, int *ia, int *ja, double *a, 
						double *b, double *u, double *r, double *d);
int symGS(int iter, double tol, int *n, int *ia, int *ja, double *a,
			    double *b, double *u, double *r, double *d);
int restrictR(int level, double *rp, double *rc);
int addProlCorrection(int level, double *up, double *uc);
int preInitialization(double *u0);
int reInitialization(double *ul);
int allocGrids(int **nl, int **ial,
               int **jal, double **al, double **bl,
			   double **dl, double **rl, double **ul);
// PROB //

int probMg(int m, int *n, 
		   int *ia, int *ja, double *a, double *b);

void computeParamLevel(int m, double *h, double *invh2, int *x0,
						int *x1, int *y0, int *y1, int *nx,
						int *ny, int *n, int *nnz);
void getNnz(int nx, int perUnit, int x0, int x1, 
					int y0, int y1, int *n, int *nnz);
int on_bound(int px, int py, int mx, int my);
int in_hole(int ix, int iy, int y0, int y1, int x0, int x1);
int check_nord(int ix, int iy, int y0, int y1, int x0, int x1, int ny);
int check_sud(int ix, int iy, int y0, int y1, int x0, int x1, int ny);
int check_west(int ix, int iy, int y0, int y1, int x0, int x1, int nx);
int check_est(int ix, int iy, int y0, int y1, int x0, int x1, int nx);

int check_nw(int ixp, int iyp, int y0p, int y1p, 
             int x0p, int x1p, int nxp, int nyp);
int check_ne(int ixp, int iyp, int y0p, int y1p, 
             int x0p, int x1p, int nxp, int nyp);
int check_sw(int ixp, int iyp, int y0p, int y1p, 
             int x0p, int x1p, int nxp, int nyp);
int check_se(int ixp, int iyp, int y0p, int y1p, 
             int x0p, int x1p, int nxp, int nyp);

int indice(int ix,int iy, int y0, int y1, int x0, int x1, int nx);


// TOOLS //
int addVect(int *n , double *v1, double *v2);


double computeBound(double x, double y);
double computeResNorm(int *n, int *ia, int *ja, double *a,
                             double *b,double *u, double *r);
int computeRes(int *n, int *ia, int *ja, double *a,
                    double *u, double *b, double *r);
void printVect(void *vect, int size, int type);
void plot_res(double *r, int level);
void simillar(int level, char **fig, int *step);
void plotCycle();

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

//extreme ev 
int computeExtremeEv(int *n, int lanczosIter, int *ia, int *ja, 
                    double *a, double *lowest, double *largest);
//1000X1000 seul 50 iter suffisant !!mettre en param to see speed
int lanczosAlgo(int iter, int *n, int *ia, int *ja, 
                 double *a, double **alpha, double **beta);
int invIter(double *alpha, double *beta, double *v, 
            int *n, double *mu, int max_iter, double tol)
void thomas_algorithm(double *alpha, double *beta, 
					  double *v, double *x, double mu, int *n);
//check randomness
int initRandomV(double *v, int *n);
int multCsrSubvector(int *n, int *ia, int *ja, 
                 double *a, double *v2, double *v1, double *v0, double *beta);
int multCsrVector(int *n, int *ia, int *ja, 
                 double *a, double *Ax, double *x);
int rayleighQuot(int *n, double *res, int *ia, int *ja, 
                 double *a, double *Ax, double *x);
int scalProd(int *n, double *v1, double *v2, double *res);
int computeVectNorm2(int *n, double *norm, double *v);
int normalize(int *n, double *v);
int vectScalDivide(int *n, double *v, double *beta);
int subVectProd(int *n, double *alpha, double *v2, double *v1);
#endif