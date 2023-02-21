#include <stdio.h>
#include <stdlib.h>
#include "prob.h"
#include "umfpk.h"
#include "time.h"
#include "rho.h"
#include "find_norm.h"
#include "plot.h"
#include "flux.h"
#include "petsc_solve.h"
#include "petscksp.h"

/* Fonction main */

int main(int argc, char *argv[])
{
  /* déclarer les variables */
  char show_graph, solver_type;
  int m, m_user, radiator_number;
  int n, *ia, *ja; 
  double *a, *b, *x_direct, *x_iterative;
  double tc1, tc2, tw1, tw2, relative_error, rho_value; /* mis à jour le 13/10/22 */
  m = 161;
  
  // User Interface 
  printf("En cas de problème avec cette interface voir le readme\n");
  printf("Souhaitez vous afficher la solution y/n ?\n");
  scanf("%c",&show_graph);
  printf("Entrez le raffinement souhaité (automatiquement converti en un pas valide)\nPour une bonne préçision prendre n = 30\nPour afficher avec m >500 prendre n >= 63\n");
  scanf("%d",&m_user);
  m = 8*m_user + 1; 
  printf("Souhaitez vous le radiateur au niveau de la fenêtre (0) ou du mur (1) ?\n");
  scanf("%d",&radiator_number);
  if(radiator_number != 1 && radiator_number != 0){
	printf("La valeur entrée doit être : \n0 pour la fenêtre \n1 pour le mur \n");
	return 0;
}
  printf("Entrez la valeur de rho souhaité\nRho optimal pour la fenêtre : 150\nRho optimal pour la porte : 7\n");
  scanf("%lf",&rho_value);
  printf("Quel(s) solveur(s) utiliser ?\nItératif : i\nDirect : d\nLes deux : b\n");
  scanf("%s",&solver_type);
  
 
  
  /* générér le problème */
  if (prob(m, &n, &ia, &ja, &a, &b, rho,rho_value, radiator_number))
     return 1;
     	
  printf("\nPROBLEM: ");
  printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n] );

  /* allouer la mémoire pour le vecteur de solution */

  x_direct = malloc(n * sizeof(double));
  x_iterative =  malloc(n * sizeof(double));
  if ( x_direct == NULL || x_iterative == NULL ) {
  	printf("\n ERREUR : pas de mémoire pour les vecteurs des solutions\n\n");
        return 1;
  }
  if (solver_type == 'i' || solver_type == 'b'){
  // SOLVEUR ITERATIF
   /* résoudre et mesurer le temps de solution */
  
  tc1 = mytimer_cpu(); tw1 = mytimer_wall();
  if(petsc_solve(argc,argv,n, ia, ja, a, b,&x_iterative))
	return 1;
  tc2 = mytimer_cpu(); tw2 = mytimer_wall();
  printf("\nTemps de solution (CPU) [ITERATIF]: %5.1f sec",tc2-tc1); /* mis à jour le 13/10/22 */
  printf("\nTemps de solution (horloge) [ITERATIF]: %5.1f sec \n",tw2-tw1); /* mis à jour le 13/10/22 */
  
  // Obtenir l'erreur relative
  relative_error = residu(n, ia,ja,a,b,x_iterative); 
  printf("relative error [ITERATIF] = %.16e \n", relative_error); 
}
  
  if (solver_type == 'd' || solver_type == 'b'){
  // SOLVEUR DIRECT
  /* résoudre et mesurer le temps de solution */
  
  tc1 = mytimer_cpu(); tw1 = mytimer_wall(); 
  if( solve_umfpack(n, ia, ja, a, b, x_direct) )
     return 1;
  tc2 = mytimer_cpu(); tw2 = mytimer_wall(); 
  printf("\nTemps de solution (CPU) [DIRECT]: %5.1f sec",tc2-tc1); 
  printf("\nTemps de solution (horloge) [DIRECT]: %5.1f sec \n",tw2-tw1);
   
  // Obtenir l'erreur relative
  relative_error = residu(n, ia,ja,a,b,x_direct); 
  printf("relative error [DIRECT] = %.16e \n", relative_error); 
}
  
  //Plot le graphe
  if(show_graph == 'y'){
	if (solver_type == 'd' || solver_type == 'b'){
	plot(m,x_direct);
	}else {
	plot(m, x_iterative);
	}
}	

  // Calculer le flux 
  if (solver_type == 'd' || solver_type == 'b')
  flux(m,n,x_direct, rho_value);
  else 
  flux(m,n,x_iterative, rho_value);
  
  /* libérer la mémoire*/
  free(ia); free(ja); free(a); free(b); free(x_direct); free(x_iterative);
 
  return 0;
}

   
