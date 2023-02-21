#include <stdlib.h>
#include <stdio.h>

int plot(int m, double *x) {
	/* Affiche le graphe via gnuplot */
	int ind = 0, l_window, r_window, down_door, up_door;
	double Tp = 20.0, Tf = 0.0,h;
	r_window = (m-1)*7/8;
    l_window = (m-1)*3/8;
    up_door = (m-1)*3/8;
    down_door = (m-1)/8; 
    h = 4.0/(m-1);
	FILE *f = NULL;
	f = fopen("data.txt","w"); 
	
	for (int iy = 0; iy < m; iy++) { // y va bien de 0 à m en terme d'indices, pour x c'est différent
		for (int ix = 0;ix < m; ix++){
			fprintf(f, "%f %f ",ix*h,iy*h);
			
			//Gèrer les cas où on doit afficher une autre valeur que celle venant de x;
			if (iy == 0 && ix >= l_window && ix <= r_window){
				fprintf(f,"%f",Tf);
			} else if (ix == 0 && iy >= down_door && iy <= up_door){
				fprintf(f,"%f",Tp);
			} else if (iy > (m-1)*5/8 && ix > (m-1)*5/8){
				fprintf(f,"NaN");
			} else {
				fprintf(f,"%f",x[ind]);
				ind++;
			}
			fprintf(f,"\n");
		}
		fprintf(f,"\n");
	}
	fclose(f);
	system("gnuplot -p gnuplot.gp\\  ");
	return 1;
}
	
