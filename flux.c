#include <stdlib.h>
#include <stdio.h>
#include <math.h>


double mean(double *x, int n){
	/* Calcul la moyenne */
	double mean_value = 0.0;
	for(int k = 0; k < n ; k++){
		mean_value += x[k];
	}
	mean_value = mean_value/n;
	return mean_value;
}

double variance(double *x, double mean, int n){
	/* Calcul la variance via la formule la moins couteuse */
	double variance_value = 0.0;
	for (int k = 0; k<n; k++){
		variance_value += x[k]*x[k];
	}
	variance_value /= n;
	variance_value -= mean*mean; 
	variance_value = sqrt(variance_value); 
	return variance_value;
}


double flux(int m, int n, double *x, double rho_value){
	/*déterminer le flux passant par la fenêtre */
	
	double h = 4.0/(m-1);
	int ind_l_w = (m-1)*3/8-1; //Plus grand indice à la gauche de la porte
	int ind_r_w = (m-1)*3/8; //Plus petit indice à la droite de la porte
	
	/* Calcule la Puissance du radiateur (rho*S) et le flux sortant de la fenêtre via l'équation de fourier (flux = - lamda*grad(T)) */ 
	double window_flux,radiator_power,lambda, mean_x, variance_x;
	lambda = 0.026;
	window_flux = 0;
	
	// Puissance du radiateur (rho supposé constant) [rho*S*k]
	radiator_power = lambda*rho_value*2*0.3;
	printf("La puissance du radiateur : %f W/m \n",radiator_power);

	// Calcul du flux de chaleur le long d'un contour enserrant la fenêtre [flux = integrale(k*gradT*dL)]
	if (m<33){
		printf("Le flux n'est pas défini pour un m < 33\n");
	} 
	else { 
		double grad1,grad2;
	 	//flux normal en dx à gauche
		grad1 = 1/(2*h)*(x[ind_l_w -1] - 0);
		grad2 = 1/(2*h)*(x[ind_l_w + (m-1)/2 - 1 ] - x[ind_l_w + (m-1)/2 + 1]);
		window_flux += (grad2+grad1);
		
			
	 // flux normal en dy 
		grad1 = 1/(2*h)*(x[ind_l_w + (m-1)*3/2 +1] - x[ind_l_w]);
		grad2 = 1/(2*h)*(x[ind_l_w + (m-1)*3/2 + 2] - 0);
		window_flux += (grad2+grad1);
		
	 for(int k = (m-1)*7/8 ; k < (m-1)*11/8-1; k++){ 
		 grad1 = 1/(2*h)*(x[k+m] - 0);
		 grad2 = 1/(2*h)*(x[k+1+m] - 0);
		 window_flux += (grad2+grad1);
		
	}   
		grad1 = 1/(2*h)*(x[ind_r_w + 2*m -1] - 0);
		grad2 = 1/(2*h)*(x[ind_r_w + 2*m] - x[ind_r_w]);
		window_flux += (grad2+grad1);
	
	  // flux normal en dx à droite
		grad1 = 1/(2*h)*(x[ind_r_w +1] - 0);
		grad2 = 1/(2*h)*(x[ind_r_w + m + 1] - x[ind_r_w + m - 1]);
		window_flux += (grad2+grad1);
	
		
	window_flux *= lambda;
	window_flux *= h/2;
	printf("flux = %f W/m \n",window_flux);
 }
	mean_x = mean(x,n);
	variance_x = variance(x, mean_x,n);
		
	return variance_x/mean_x;
	
}

// Pour la configuration en face de la fenêtre on a un minimum de var/moy pour rho = 150 m = 19.3 e-t = 3.05
// Pour la seconde configuration à côté de la porte on trouve un minimum de var/moy pour rho = 7 m = 11.75 et e-t 4.82





/*
	//flux normal en dx à gauche
		grad1 = 1/(2*h)*(x[ind_l_w] - x[ind_l_w-2]);
		grad2 = 1/(2*h)*(x[ind_l_w + (m-1)/2] - x[ind_l_w -2 + (m-1)/2]);
		window_flux += (grad2+grad1);
		
		grad1 = grad2;
		grad2 = 1/(2*h)*(x[ind_l_w + (m-1)*3/2 +1 ] - x[ind_l_w -2 + (m-1)*3/2 +1]);
		window_flux += (grad2+grad1);
		
	 // flux normal en dy 
	 for(int k = (m-1)*7/8; k < (m-1)*11/8+1; k++){ 
		 grad1 = 1/(2*h)*(x[k+2*m] - x[k]);
		 grad2 = 1/(2*h)*(x[k+1+2*m] - x[k+1]);
		 window_flux += (grad2+grad1);
		
	}
	  // flux normal en dx à droite
	  for(int j =0; j <2; j++){
		  grad1 = 1/(2*h)*(x[ind_r_w +j*m] - x[ind_r_w+2+j*m]);
		  grad2 = 1/(2*h)*(x[ind_r_w +(j+1)*m] - x[ind_r_w+2+(j+1)*m]);
		  window_flux += (grad2+grad1);
		  
		}
*/



