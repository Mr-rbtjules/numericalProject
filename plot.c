#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <proto.h>

#define COORD_X0 1.0
#define COORD_X1 2.5
#define COORD_Y0 1.5
#define COORD_Y1 2.0


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
                if ( on_bound(ix,iy,m) || in_hole(ix,iy,y0,y1,x0,x1)){
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
void plot_dyn(double **u_t, int m, int size_ut, double dt)
{

    /*cree un fichier texte contenant 2 colonnes
    avec les coordonnées x et y et des colonnes sucessive donnant la temperature au cours du temps
    
    Arguments
    =========
    u_t (input) - pointeur vers les vecteurs u_ti
    m (input)   - nombre de points par directions dans la grille
    ite (input) - nombre d'elements de u_t
    dt (input)  - pas
    */
    
    double h = 3.0 / (double)(m-1);
    
    
    /*creation du ficher*/
    FILE* pointFile = NULL;                        //renvoi un pointeur pointant vers un type FILE
    const char* file_name = "coord_dyn.dat";

    pointFile = fopen(file_name,"w+");           // ouverture en mode ecriture
    
    if (pointFile != NULL){

        printf("Creation du fichier .dat\n");

        /*Calcul des constantes*/
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
    

    /*script gnuplot*/

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

