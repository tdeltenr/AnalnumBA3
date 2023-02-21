#include <stdio.h>
#include <stdlib.h>

double rho(int ix, int iy, int m, double rho_value, int radiator_number){
	/* Calcule la valeur du second membre en (ix,iy) */
	double h,m_radiator; 
	h = 4.0/(m-1); 
	m_radiator = 0.15/h; //Calcul du nombre de point pour obtenir un radiateur de 30cm de "profondeur"
	// Radiateur proche de la fenêtre
	if(radiator_number == 0){
		if ( (ix >= (m-1)*3/8 && ix <= (m-1)*7/8) && (iy >= (m-1)/8 - m_radiator && iy <= (m-1)/8 + m_radiator)){
			return rho_value;
		}
	}

	// Radiateur proche de la porte
	else if(radiator_number == 1){
		if ((ix >= (m-1)/8 - m_radiator && ix <=  (m-1)/8 + m_radiator) && ((iy >= (m-1)*3/8 && iy <= (m-1)*7/8)) ){
			return rho_value;
		}
	}
	
	return 0.0;
}
