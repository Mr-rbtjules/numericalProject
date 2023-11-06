#ifndef PROTO
#define PROTO



typedef struct globVal_s globVal_s;
struct globVal_s{

	int m;
	int *domain;

};
extern globVal_s globVal;



// PROB //
int probMg(int mc, int *nl, 
		   int *ial, int *jal, double *al, double *bl);

void computeParamLevel(int mc, double *hl, double *invh2l, int *x0l,
						int *x1l, int *y0l, int *y1l, int *nxl,
						int *nl, int *nnzl);
void getNnz(int nx, int perUnit, int x0, int x1, int y0, int y1, int *n, int *nnz);
double computeBound(double x, double y);
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
void freeGlobVal();


#endif