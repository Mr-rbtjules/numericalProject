#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "proto.h"

#define COORD_X0 1.0
#define COORD_X1 2.5
#define COORD_Y0 1.5
#define COORD_Y1 2.0

#define TAILLE_MAX 200


int add_plot(double *x, int n){
    //faudra juste savoir qu'on alterne 
}

void plot_static(double *x, int m, int level){

    /*cree un fichier texte (.dat pour gnuplot) contenant 2 colonnes
    avec les coordonnées x et y et une colonnes avec les éléments de u (bord et trous compris avec u = 0)
    et lance le script gnuplot pour plot le resultat (nmp = numero de mode propre)

    Arguments
    =========
    x (input) - pointeur vers le vecteur a plot
    m (input)     - nombre de points par directions dans la grille
    */
    
    double h = 3.0/(double)(m-1);
    double invh2 = 1.0/(h*h);

    // coordonnees du trou sur la grille discrete
    
    /*int x1 = ((int)(COORD_X1 * (m-1)) /3)  -1; 
    int x0 = (((int)(COORD_X0 * (m-1)) + ((3 - ((int)(COORD_X0*(m-1))%3))%3))/3)  -1;//4 alors que 3
    int y1 = ((int)(COORD_Y1 * (m-1)) /3)  -1;
    int y0 = (((int)(COORD_Y0 * (m-1)) + ((3 - ((int)(COORD_Y0*(m-1))%3))%3))/3)  -1; //6 alors qu'on veut 5
    */
    int x0,x1,y0,y1;
    computeHole(&x0,&x1,&y0,&y1, m); //y0 = 7 mais devrait 6 comme 16

    
    int nx = m-2;
    //nb de points sur la largeur du trous

    int p = y1 - y0 + 1;//m%2 + m/6; //== nb de points entre y1 et y0
    // "" sur la longueur
    int q = x1- x0 + 1;//m%2 + m/2;
    //plaque hors bords et trous
    int n = nx * nx - (p * q);
   

   /*level*/
    double hl = h*pow(2, level);
	double invh2l = 3.0/(hl*hl); 

    int x0l = x0 / pow(2, level); // va arrondir au point grille coarse a droite 
    int x1l = ((x1+1)/pow(2, level)) - 1; // permet si x1 pair on retire 1 
    int y0l = (y0/pow(2, level)); // arrondu coarse au dessus (permet de pas ajouter des points dans le trou)
    int y1l = ((y1+1)/pow(2, level)) - 1;
    int pl = y1l - y0l + 1;
    int ql = x1l - x0l + 1;

    int nxl = nx/pow(2,level); // nb de points coars sur un ligne pas bord
    int nl = nxl * nxl - (pl * ql);

    int ml = nxl+2;

    /*creation du ficher*/
    FILE* pointFile = NULL;                        //renvoi un pointeur pointant vers un type FILE
    const char* file_name = "coord_stat.dat";

    pointFile = fopen(file_name,"w+");           // ouverture en mode ecriture
    
    if (pointFile != NULL){

        /*Calcul des constantes*/
        
        int ind = 0;
        for (int py = 0; py < ml; py++){ //passage ligne suivante
            for (int px = 0; px < ml; px++){      //passage colonne suivante
                
                //si sur bord ou trou
                if ( on_bound(px,py,ml) ){
                    
                    fprintf(pointFile, "%.16g %.16g 0\n", (px*hl), py*hl);
                }
                else if ( in_hole(px-1,py-1,y0l,y1l,x0l,x1l)){
                    fprintf(pointFile, "%.16g %.16g 0\n", (px*hl), py*hl);
                }
                else{
                    fprintf(pointFile, "%.16g %.16g %.16g\n", (px*hl), (py*hl), x[ind]);
                    ind += 1;
                }
            }
            fprintf(pointFile, "\n");
        }
        //fermeture du ficher
        fclose(pointFile);
        
    }

    //plot le graphique

    FILE *gnuplot = popen("gnuplot gnuScript_stat.gnu -persistent","w");

    fclose(gnuplot);

    //supprime coord_stat.dat
    system("rm coord_stat.dat");

}
/*
void plot_dyn(double **u_l, int m, int *size_ul, double time)
{

    
    
    double h = 3.0 / (double)(m-1);
    
    int x0,x1,y0,y1;
    computeHole(&x0,&x1,&y0,&y1, m);
    FILE* pointFile = NULL;                        //renvoi un pointeur pointant vers un type FILE
    const char* file_name = "coord_dyn.dat";

    pointFile = fopen(file_name,"w+");           // ouverture en mode ecriture
    
    if (pointFile != NULL){

        printf("Creation du fichier .dat\n");

        int y0 = (int)ceil(COORD_Y0/h);
        int y1 = (int)ceil(COORD_Y1/h);
        int x0 = (int)ceil(COORD_X0/h);
        int x1 = (int)ceil(COORD_X1/h);

        int nx = m;
        int ind = 0;
        for (int iy = 0; iy < nx; iy++){ //passage ligne suivante
            for (int ix = 0; ix < nx; ix++){      //passage colonne suivante
                
                //2colonnes pour coord x y
                fprintf(pointFile, "%f %f ", (ix*h), (iy)*h);

                //ajoute temperature au cours du temps

                //si sur bord ou interieur au trou
                if (iy == 0 || iy == (nx-1) || ix == 0 || ix == (nx-1) || ( iy <= y1  && iy >= y0  && ix >= x0 && ix <= x1 ) ){
                    for (int i = 0; i < size_ut; i++){
                        fprintf(pointFile, "0.0000000 ");
                    }
                }
                else{
                    for (int i = 0; i < size_ut; i++){
                        fprintf(pointFile, "%f ", u_t[i][ind]);
                    }
                    ind += 1;
                }
                
                //passe a la ligne pour initialiser la coord suivante:
                fprintf(pointFile,"\n");
            }
            fprintf(pointFile,"\n"); //passe a la ligne a chaque fois que la coord y change pour être adapté a gnu
        }
        //fermeture du ficher
        fclose(pointFile);
    }
    


    printf("Creation du gif...\n");

    FILE* gnuplot = popen("gnuplot -persistent", "w");

    //deverse le contenu du script ligne par ligne dans gnuplot pour alleger le code
    FILE *script = fopen("gnuScript_dyn.gnu", "r");
    char ligne[TAILLE_MAX] = "";                      //prend max 50 caract par ligne

    while (fgets(ligne, TAILLE_MAX, script) != NULL){
        fprintf(gnuplot, "%s\n",ligne);
    }

    //initialise le nombre d'image 
    fprintf(gnuplot,"\nimax = %d\n", size_ut);
    
    //initialise pas et m(pour l'indiquer sur le graphique)
    fprintf(gnuplot, "pas = %f\nm = %d\n", dt, m);

    //lance la boucle
    fprintf(gnuplot, "load 'boucle_dyn.gnu'\n");

    //fermeture 
    fclose(gnuplot);
    
    //supprime coord_dyn.dat
    system("rm coord_dyn.dat");
}
*/


