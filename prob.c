#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int up_skipped(int ix,int iy, int m){
	/* Calcule le nbre de points de différence avec un simple + m pour trouver l'indice nord dans prob.c*/
	if (iy == 0 && ix < (m-1)*3/8){
		if (m==9)
		 return (m-1)/2 +2; // points skippés par la fenêtre + la porte
		else
		return (m-1)/2 +1; // points skippés uniquement par la fenêtre la porte étant à iy =/= 1 
		
	 } else if (iy == 0 && ix > (m-1)*7/8 && m == 9){
		return 1;
	
	} else if (iy < (m-1)/8 -1 ){
		return 0;
		
	} else if (iy  < (m-1)*3/8 && iy >= (m-1)/8-1){  
		return 1;
	
	} else if (iy <= (m-1)*5/8){
		return 0;
	
	} else 
		return (m-1)*3/8;
}


int down_skipped(int ix,int iy, int m){
	/* Calcule le nbre de points de différence avec un simple - m pour trouver l'indice sud dans prob.c*/
	 if (iy == 1 && ix < (m-1)*3/8){
		 if (m==9)
		 return (m-1)/2 +2; // points skippés par la fenêtre + la porte
		else
		return (m-1)/2 +1; // points skippés uniquement par la fenêtre la porte étant à iy =/= 1 
	} else if (iy < ((m-1)/8 + 1)){
		return 0;
		
	} else if (iy >= ((m-1)/8+1) && iy <= ((m-1)*3/8)){
		return 1;
	
	} else if (iy <= ((m-1)*5/8 +1)){
		return 0; 
		
	} else 
		return (m-1)*3/8;
}


int prob(int m, int *n, int **ia, int **ja, double **a, double **b, double rho(int,int,int,double,int),double rho_value, int radiator_number)
/*
   But
   ===
   Générer le système linéaire n x n 
                          
                             Au = b                                   

   qui correspond à la discrétisation sur une grille cartésienne 
   régulière m x m de l'équation de Poisson à deux dimensions
              
            d    d        d    d
         - == ( == u ) - == ( == u )  = 0     sur [0,1] x [0,1]
           dx   dx       dy   dy

  avec les conditions aux limites de Dirichlet
         
         u = 20     sur (0,y)                 , avec 0.5 <=  y  <= 1.5
         du/dn = 1  sur (1,y), (x,0) et (x,1) , avec 0 <= x,y <= 1 .
         u = 0      sur (x,0)                 , avec 1.5 <= x <= 3.5

  La numérotation des inconnues est lexicographique, la direction x étant 
  parcourue avant celle de y. La matrice A est retournée dans le format 
  CRS qui est défini via trois tableaux : 'ia', 'ja' et 'a'.

  Arguments
  =========
  m  (input)  - nombre de points par direction dans la grille
                (les valeurs m inférieures à 2 ne sont pas valides) 
  n  (output) - pointeur vers le nombre d'inconnus dans le système
  ia (output) - pointeur vers le tableau 'ia' de la matrice A
  ja (output) - pointeur vers le tableau 'ja' de la matrice A
  a  (output) - pointeur vers le tableau 'a' de la matrice A
  b  (output) - pointeur vers le tableau 'b'

*/
{
    int  nnz, ix, iy, ind,next_ind,r_window,l_window,up_door,down_door;
    double invh2, Tp = 20.0, Tf = 0.0;

    if(((m-1)%8) != 0) { // Vérifier que le pas est une valeur valide !!! 
        printf("\n ERREUR : m = %d n'est pas une valeur valide\n\n",m);
        return 1;
    }
    
    invh2 = (m-1)*(m-1)/16; /* h^-2 pour L=4 */
    *n  = m*m - ((m-1)/2 + 1) - ((m-1)/4 +1) - (m-1)*(m-1)*9/64 ; /* nombre d'inconnues */
    nnz = 5*(*n) - 4*m - 4; /* nombre d'éléments non nuls 5n pour les pts intérieurs - 4m pour les murs - 4 pour les bords avec les fenêtres */
    
    /* allocation des tableaux */
    *ia  = (int*)malloc((*n + 1) * sizeof(int));
    *ja  = (int*)malloc(nnz * sizeof(int));
    *a   = (double*)malloc(nnz * sizeof(double));
    *b   = (double*)malloc(*n*sizeof(double));
    /* allocation réussite? */

    if (*ia == NULL || *ja == NULL || *a == NULL || *b == NULL ) {
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n\n");
        return 1;
    }

    /* partie principale : remplissage de la matrice */
	
    ind = 0; /* au cas où m<=1 */
    next_ind = 0;
    nnz = 0; 
    r_window = (m-1)*7/8;
    l_window = (m-1)*3/8;
    up_door = (m-1)*3/8;
    down_door = (m-1)*1/8; 
    int j = 0;
    for (iy =0 ; iy < m; iy++) { // y va bien de 0 à m en terme d'indices, pour x c'est différent
		for (ix =0;ix < m; ix++){
			
			// Supprimer les cas ou on ne doit pas tourner;
			if (iy == 0 && ix >= l_window && ix <= r_window){
				j++;
			} else if (ix == 0 && iy >= down_door && iy <= up_door){
				j++;
			} else if (iy > (m-1)*5/8 && ix > (m-1)*5/8){
				j++;
			} else {
				
				/* Trouver l'indice du point et l'incrémenter pour la suite */
				ind = next_ind;
				next_ind++;
				/* marquer le début de la ligne suivante dans le tableau 'ia' */
				(*ia)[ind] = nnz;
			
				/* calculer le membre de droite */   
			    (*b)[ind] = rho(ix, iy, m, rho_value, radiator_number);
   
				/* (v) remplissage de la ligne : voisin sud */
				if (iy == m-1 || (iy == (m-1)*5/8 && ix > (m-1)*5/8)) { /* condition de Neumann, bord nord    */         
					(*a)[nnz] = -2*invh2;
					(*ja)[nnz] = ind - m + down_skipped(ix,iy,m);
					nnz++;
				} else if (ix == 0 && iy == up_door +1){ /* Point au dessus de porte*/
					(*b)[ind] += Tp*invh2;
				}else if (iy == 1 && ix >= l_window && ix <= r_window) {/* Condition de Dirichlet de la fenêtre */
						(*b)[ind] += Tf*invh2;
	
				} else if (iy > 0) { /* reste du domaine */
					(*a)[nnz] = -invh2; 
					(*ja)[nnz] = ind - m + down_skipped(ix,iy,m);
					//printf("ind = %d\n",ind - m + down_skipped(ix,iy,m));
					nnz++;

				} 

				/* remplissage de la ligne : voisin ouest */
				if (ix == 1 && iy >= down_door && iy <= up_door){ /* points à droite de la porte  */
					(*b)[ind] += Tp*invh2;

				} else if (iy == 0 && ix == r_window + 1 ){ /* point à droite de la fenêtre*/
					(*b)[ind] += Tf*invh2;

				} else if ( (ix == m-1 && iy <= (m-1)*5/8) || (ix == (m-1)*5/8&& iy > (m-1)*5/8) )   { /* condition de Neumann, bord est */   		
					(*a)[nnz] = -2*invh2; 
					(*ja)[nnz] = ind - 1;
					nnz++;

				} else if (ix > 0) { /* milieu du domaine */
					(*a)[nnz] = -invh2; 
					(*ja)[nnz] = ind - 1;
					nnz++;

				}
			
     
				/* remplissage de la ligne : élément diagonal */
				(*a)[nnz] = 4.0*invh2;
				(*ja)[nnz] = ind;
				nnz++;
	

				/* remplissage de la ligne : voisin est */
				if (iy == 0 && ix == l_window-1){ /* point voisin de la fenêtre*/
					(*b)[ind] += Tf*invh2;
				} else if (ix == 0 && !( iy <= up_door && iy >= down_door)){ /* Condition de Neumann bord ouest */
					(*a)[nnz] = -2*invh2; 
					(*ja)[nnz] = ind + 1;
					nnz++;

				} else if ((iy <= (m-1)*5/8 && ix < m-1) || (iy > (m-1)*5/8 && ix <= ((m-1)*5/8)-1) )  { /* milieu du domaine */
					(*a)[nnz] = -invh2;
					(*ja)[nnz] = ind + 1;
					nnz++;

				} 

				/* remplissage de la ligne : voisin nord */
				if (ix == 0 && iy == down_door-1){
					if (m==9) {
					(*b)[ind] += 2*Tp*invh2; /* Point sous la fenêtre si m = 9 il coincide avec un terme de facteur 2 dans la dérivée*/
					}else 
					(*b)[ind] += Tp*invh2;
	
				} else if (iy == 0) { /* condition de Neumann, bord sud */         
					(*a)[nnz] = -2*invh2; 
					(*ja)[nnz] = ind + m - up_skipped(ix,iy,m);
					nnz++;
	
				} else if ((ix <= (m-1)*5/8 && iy < m-1) || (ix > (m-1)*5/8 && iy < (m-1)*5/8) ) { /* reste des points */
					(*a)[nnz] = -invh2;
					(*ja)[nnz] = ind + m - up_skipped(ix,iy,m);
					nnz++;
            }
         }
      }
    }	
    /* dernier élément du tableau 'ia' */
    (*ia)[ind + 1] = nnz;
    
    /* retour habituel de fonction */
    return 0;
}
