#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "proto.h"


void plotCycle(int levelMax, int cycle){
    int numberOfLines = (levelMax)*2 +1; 
    char **paragraph = malloc(numberOfLines * sizeof(char*)); 
    int sizemax = 200;
    for (int i = 0; i < numberOfLines; i++) {
        paragraph[i] = malloc(sizemax * sizeof(char)); // Allocate memory for each line (100 characters)
        for (int j = 0; j < sizemax; j++){
            paragraph[i][j] = ' ';
        }
    }
    
    int step = 0;
    simillar(levelMax, cycle, 0, paragraph, &step);

    // Print the paragraph
    printf("\n             MULTIGRID CYCLE \n");
    printf("Top Level\n");
    for (int i = 0; i < numberOfLines; i++) {
        printf("%s\n", paragraph[i]);
    }
    printf("Bottom Level\n");

    // Free the allocated memory
    for (int i = 0; i < numberOfLines; i++) {
        free(paragraph[i]);
    }
    free(paragraph);
}

void simillar(int levelMax, int cycle, int level, char **fig, int *step){

    if (level < levelMax){
        plotDown(level, fig, step);
        /*for (int i = 0; i < (LEVELMAX)*2 +1; i++) {
            printf("%s\n", fig[i]);
        }  */ 
        
        simillar(levelMax, cycle, level+1, fig, step);
        if (cycle == 1){
            simillar(levelMax, cycle,level+1, fig, step);
        }
        
        plotUp(level+1, fig,step); 
        if (level > 0){
            if (cycle == -1){
                FPart(levelMax,level,fig,step); 
            }
        }
    }
    else{            
        if (fig[level*2][*step - 1] != level + '0'){
            fig[level*2][*step] = level + '0';
            *step +=1;
        }
    }

}

void FPart(int levelMax, int level, char **fig, int *step){
    int initLevel = level;
    printf("level %d\n", level);
    while (level < levelMax){
        level += 1;
        fig[level*2-1][*step] = '\\';
        *step +=1;
        if (fig[level*2][*step - 1] != level + '0'){
            fig[level*2][*step] = level + '0';
            *step +=1;
        }
        
    }
    while (level > initLevel){
        plotUp(level, fig,step);
        level -= 1;
    } 
}

void plotUp(int level, char **fig, int *step){

        //getup 
        // print / step+1
        //print level step+1
        fig[(level-1)*2 + 1][*step] = '/';
        *step +=1;
        fig[(level-1)*2][*step] = (level-1) + '0';
        *step +=1;
}


void plotDown(int level, char **fig, int *step){
    if (fig[level*2][*step - 1] != level + '0'){
        fig[level*2][*step] = level + '0';
        *step +=1;
    }
    
    fig[level*2+1][*step] = '\\';
    *step +=1;
}

void plot_res(double *r, int level){

    //cree un fichier texte (.dat pour gnuplot) contenant 2 colonnes
    //avec les coordonnées x et y et une colonnes avec les éléments de u
    // (bord et trous compris avec u = 0) et lance le script gnuplot pour
    // plot le resultat (nmp = numero de mode propre)

    double h, invh2;
    int x0,x1,y0,y1, nx, ny, nnz;
    int n;
    computeParamLevel(
        globVal.m[level], &h,&invh2,&x0,&x1, &y0, &y1, &nx, &ny, &n, &nnz
    );
    
    int mx = nx+2;
    int my = ny+2;

    //creation du ficher
    FILE* pointFile = NULL;                        
    const char* file_name = "coord_stat.dat";
    // ouverture en mode ecriture
    pointFile = fopen(file_name,"w+");           
    
    if (pointFile != NULL){

        //Calcul des constantes
        
        int ind = 0;
        //indice 00 different que pour prob, ici c'est le coin du domaine
        for (int py = 0; py < my; py++){ //passage ligne suivante
            for (int px = 0; px < mx; px++){      //passage colonne suivante
                
                //si sur bord ou trou
                if ( on_bound(px,py,mx,my) || in_hole(px-1,py-1,y0,y1,x0,x1)){
                    
                    //fprintf(pointFile, "%.16g %.16g 0\n", px*hl, py*hl);
                    fprintf(pointFile, "%.16g %.16g NaN\n", px*h, py*h);
                }
                else{
                    fprintf(
                        pointFile,
                        "%.16g %.16g %.16g\n",
                        (px*h),
                        (py*h),
                        r[ind]
                    );
                    ind += 1;
                }
            }
            fprintf(pointFile, "\n");
        }
        //fermeture du ficher
        fclose(pointFile);
    }

    //initialise fichier gnuplot
    //plot le graphique

    FILE *gnuplot = popen("gnuplot -persistent","w");
    fprintf(gnuplot, "set title 'm = %d level = %d'\n", mx, level);
    fprintf(gnuplot, "load 'gnuScript_stat.gnu'\n");
    fclose(gnuplot);

    //supprime coord_stat.dat
    int _ = system("rm coord_stat.dat");
}

void printVect(void *vect, int size, int type) {
    int i;
    printf("\n[");
    if (type == 0) {  // Assuming type 0 represents int
        int *intVect = (int *)vect;
        for (i = 0; i < size; i++) {
            printf("%d, ", intVect[i]);
        }
    } 
    else if (type == 1) {  // Assuming type 1 represents double
        double *doubleVect = (double *)vect;
        for (i = 0; i < size; i++) {
            printf("%lf, ", doubleVect[i]);
        }
    printf("]\n");
    }
}

void printMatrix(int* ia, int* ja, double* a, int *n) {
    for (int i = 0; i <= *n; i++) {
        int rowStart = ia[i];
        int rowEnd = ia[i + 1];
        int jaIndex = rowStart;

        for (int j = 0; j <= *n; j++) {
            if (jaIndex < rowEnd && ja[jaIndex] == j) {
                printf("%.2f\t", a[jaIndex]);
                jaIndex++;
            } else {
                printf("0.00\t");
            }
        }
        printf("\n");
    }
}

int plotIter(double *res, int iter, int m){

    // Step 1: Write the data to a text file
    FILE *file = fopen("dataIter.txt", "w");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }
    for (int i = 0; i < iter; i++) {
        fprintf(file, "%d %.16f\n", i + 1, res[i]);
    }
    fclose(file);

   char title[50]; // Ensure this is large enough to hold the entire title
sprintf(title, "set title 'm = %d';", m);

// Step 2: Use gnuplot to plot the data
char command[1024]; // Ensure this is large enough to hold the entire gnuplot command
sprintf(command,
       "gnuplot -p -e \""
       "set terminal png size 800,600; "
       "set output 'iter.png'; "
       "set xlabel 'Iterations'; "
       "set ylabel 'Res'; "
       "%s "                   // This inserts the title with the value of m
       "set grid; "
       "set logscale y; "
       "unset key; "
       "plot 'dataIter.txt' using 1:2 with linespoints; "
       "\"", title);

// Execute the gnuplot command
system(command);
    printf("Plot generated as 'iter.png'\n");

    return 0;
}


