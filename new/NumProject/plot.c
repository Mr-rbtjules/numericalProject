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



 /*ajouter fonctionnalité ou trou en couleur diff*/



    /*cree un fichier texte (.dat pour gnuplot) contenant 2 colonnes
    avec les coordonnées x et y et une colonnes avec les éléments de u (bord et trous compris avec u = 0)
    et lance le script gnuplot pour plot le resultat (nmp = numero de mode propre)

    Arguments
    =========
    x (input) - pointeur vers le vecteur a plot
    m (input)     - nombre de points par directions dans la grille
    */
    double hl, invh2l;
    int x0l,x1l,y0l,y1l, nxl, nl, nnzl;

    computeParamLevel(m, level, &hl,&invh2l,&y0l,&y1l,&x0l,&x1l,&nxl, &nl, &nnzl);
   
    
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
                if ( on_bound(px,py,ml) ||  in_hole(px-1,py-1,y0l,y1l,x0l,x1l) ){
                    

                    //fprintf(pointFile, "%.16g %.16g 0\n", px*hl, py*hl);
                    fprintf(pointFile, "%.16g %.16g %.16g\n", px*hl, py*hl, computeBound(px*hl, py*hl));
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

    //initialise fichier gnuplot


    //plot le graphique

    FILE *gnuplot = popen("gnuplot -persistent","w");
    fprintf(gnuplot, "set title 'm = %d level = %d'\n", m, level);
    fprintf(gnuplot, "load 'gnuScript_stat.gnu'\n");
    fclose(gnuplot);

    //supprime coord_stat.dat
    system("rm coord_stat.dat");

}


void plot_res(double *r, int m, int level){



 /*ajouter fonctionnalité ou trou en couleur diff*/



    /*cree un fichier texte (.dat pour gnuplot) contenant 2 colonnes
    avec les coordonnées x et y et une colonnes avec les éléments de u (bord et trous compris avec u = 0)
    et lance le script gnuplot pour plot le resultat (nmp = numero de mode propre)

    Arguments
    =========
    x (input) - pointeur vers le vecteur a plot
    m (input)     - nombre de points par directions dans la grille
    */
    double hl, invh2l;
    int x0l,x1l,y0l,y1l, nxl, nl, nnzl;

    computeParamLevel(m, level, &hl,&invh2l,&y0l,&y1l,&x0l,&x1l,&nxl, &nl, &nnzl);
   
    
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
                if ( on_bound(px,py,ml) ||  in_hole(px-1,py-1,y0l,y1l,x0l,x1l) ){
                    

                    //fprintf(pointFile, "%.16g %.16g 0\n", px*hl, py*hl);
                    fprintf(pointFile, "%.16g %.16g 0\n", px*hl, py*hl);
                }
                else{
                    fprintf(pointFile, "%.16g %.16g %.16g\n", (px*hl), (py*hl), r[ind]);
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
    fprintf(gnuplot, "set title 'm = %d level = %d'\n", m, level);
    fprintf(gnuplot, "load 'gnuScript_stat.gnu'\n");
    fclose(gnuplot);

    //supprime coord_stat.dat
    system("rm coord_stat.dat");

}