#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <proto.h>

#define COORD_X0 1.0
#define COORD_X1 2.5
#define COORD_Y0 1.5
#define COORD_Y1 2.0


void plot_static(double *x, int m){

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
    
    int x1 = ((int)(COORD_X1 * (m-1)) /3)  -1; 
    int x0 = (((int)(COORD_X0 * (m-1)) + ((3 - ((int)(COORD_X0*(m-1))%3))%3))/3)  -1;
    int y1 = ((int)(COORD_Y1 * (m-1)) /3)  -1;
    int y0 = (((int)(COORD_Y0 * (m-1)) + ((3 - ((int)(COORD_Y0*(m-1))%3))%3))/3)  -1;
    /*
    
    /*creation du ficher*/
    FILE* pointFile = NULL;                        //renvoi un pointeur pointant vers un type FILE
    const char* file_name = "coord_stat.dat";

    pointFile = fopen(file_name,"w+");           // ouverture en mode ecriture
    
    if (pointFile != NULL){

        /*Calcul des constantes*/
        
        int ind = 0;
        for (int iy = 0; iy < m; iy++){ //passage ligne suivante
            for (int ix = 0; ix < m; ix++){      //passage colonne suivante
                
                //si sur bord ou trou
                if (iy == 0 || iy == (m-1) || ix == 0 || ix == (m-1) || ( iy <= y1  && iy >= y0  && ix >= x0 && ix <= x1 )){
                    fprintf(pointFile, "%f %f 0\n", (ix*h), iy*h);
                }
                else{
                    fprintf(pointFile, "%f %f %f\n", (ix*h), (iy*h), x[ind]);
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

