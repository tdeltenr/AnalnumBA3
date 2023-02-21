#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double residu(int n,int *ia, int *ja, double *a, double *b, double *x){
	/* Renvoie la norme du résidu*/
	double Ax, r, normr, normb;
	normr = 0.0;
	normb = 0.0;
	for(int i = 0; i < n; i++){
		Ax = 0.0;
		for(int nz_id = ia[i]; nz_id < ia[i+1]; nz_id++){
			Ax += a[nz_id]*x[ja[nz_id]];
		}
		r = Ax - b[i];
		normr += r*r;
		normb += b[i]*b[i];
	}
	return sqrt(normr/normb);
}

