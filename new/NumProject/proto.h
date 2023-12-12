#ifndef PROTO
#define PROTO

#define BOUND1(x, y) (exp(sqrt((x)*(x) + (y)*(y))))
#define BOUND2(x, y) (sin(sqrt(pow(x,2) + pow(y,2))))
#define BOUND0(x, y) (0)

//program parameters
//type of boundary conditions
#define BOUND BOUND2
//type of domain
//si marche pas check indice restrict prolong
#define DOMAIN1 "[0, 6] × [0, 10] [2, 4] × [2, 5]"
#define DOMAIN2 "[0, 6] × [0, 10] [0, 2] × [3, 4]"
#define DOMAIN3 "[0, 6] × [0, 10] [4, 6] × [3, 6]"
#define DOMAIN4 "[0, 6] × [0, 10] [2, 3] × [5, 8]"
#define DOMAIN5 "[0, 6] × [0, 10] [3, 4] × [0, 4]"
#define DOMAIN6 "[0, 6] × [0, 10] [0, 4] × [0, 3]"
#define DOMAIN7 "[0, 6] × [0, 10] [2, 4] × [5, 10]"
#define DOMAIN8 "[0, 6] × [0, 10] [0, 4] × [3, 10]"
#define DOMAIN DOMAIN2
//discretisation 


#define MODE 1 //0umpf 1 sym 2jacobi 
#define ITERCORE 20
#define CYCLEW 0
#define CYCLEF 0

#define EXPLICIT 0
#define LOAD 1
#define CHRONO 1
#define RELAX 0
#define PLOT 0


typedef struct globVal_s globVal_s;
struct globVal_s{

	int levelMax;
	int *m;
	int *domain;
	int *vectStart;
	int *matStart;

};
extern globVal_s globVal;



//changer type de retour


int CGmethod(int iter, int m, int levelMax, double *res);
void conjugateGradientCSR(int iter, int* ia, int* ja, double* a, double* b, double* x, int N, double *res);
// RESOL//

int mg_method(int m, int iter, int levelMax, int mu, int mode, int cycle, double *res);
int mg_iter(int iter, int levelMax, int m, int mu1, int mu2, int cycle, int *nl, int *ial,
             int *jal, double *al, double *rl, double *ul, double *dl, double *bl, double *res);

int tg_rec(int level, int m, int mu1,
			int mu2, int cycle, int *nl, int *ial, int *jal,
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

int allocGrids(int levelMax, int **nl, int **ial,
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
int isSymmetric(int* ia, int* ja, double* a, int *n);
int addVect(int *n , double *v1, double *v2);


double computeBound(double x, double y);
double computeResNorm(int *n, int *ia, int *ja, double *a,
                             double *b,double *u, double *r);
int computeRes(int *n, int *ia, int *ja, double *a,
                    double *u, double *b, double *r);
int copy( int *n, double *toCopy, double *copy);
int scalProd(int *n, double *v1, double *v2, double *res);
int subVectProd(int *n, double *alpha, double *v2, double *v1);
int addVectProd(int *n, double *alpha, double *v2, double *v1);
int dSum(int *n, double *beta, double *d, double *z);
int multCsrVector(int *n, int *ia, int *ja, 
                 double *a, double *Ax, double *x);
int preInitialization(double *u0);
int reInitialization(double *ul);

//PLOT//
void plotCycle(int levelMax, int cycle);
void simillar(int levelMax, int cycle, int level, char **fig, int *step);
void FPart(int levelMax, int level, char **fig, int *step);
void plotUp(int level, char **fig, int *step);
void plotDown(int level, char **fig, int *step);
void plot_res(double *r, int level);
void printVect(void *vect, int size, int type);
void printMatrix(int* ia, int* ja, double* a, int *n);
int plotIter(double *res, int iter, int m);

// Global variable //
void initGlobVal(int m, int levelMax);
void extractDomain(const char* str_domain, int* count, int **domain);
int correctM(int *domain, int m);
int initMlevels(int m, int **mLevels, int levelMax);
int initIndex( int levelMax, int **vectStart, int **matStart);
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
                    double *a, double *lowest, double *largest, int forward);
//1000X1000 seul 50 iter suffisant !!mettre en param to see speed
int lanczosAlgo(int iter, int *n, int *ia, int *ja, 
                 double *a, double **alpha, double **beta);
int invIter(double *alpha, double *beta, double *v, 
            int *n, double *mu, int max_iter, double tol);

int thomas_algorithm(double *alpha, double *beta, 
					  double *v, double *x, double *mu, int *n);
//check randomness
int initRandomV(double *v, int *n);
int multCsrSubvector(int *n, int *ia, int *ja, 
                 double *a, double *v2, double *v1, double *v0, double *beta);

int rayleighQuot(int *n, double *res, int *ia, int *ja, 
                 double *a, double *Ax, double *x);


int computeVectNorm2(int *n, double *norm, double *v);
int normalize(int *n, double *v);
int vectScalDivide(int *n, double *v, double *beta);

int computeMu(double *alpha, double *beta, int *n, double *mu);



void csrMatVecMult(int* ia, int* ja, double* a, double* x, double* result, int N);
double dotProduct(double* vec1, double* vec2, int N);
void vecAdd(double* vec1, double* vec2, double* result, double alpha, int N);
void vecSub(double* vec1, double* vec2, double* result, double alpha, int N);
#endif