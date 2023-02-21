/* prototype */
int prob(int m, int *n, int **ia, int **ja, double **a, double **b, double (*rho)(int,int,int,double,int), double rho_value, int radiator_number);

/* Utilisé pour calculer le voisin nord*/
int up_skipped(int ix,int iy, int m);

/* Utilisé pour calculer le voisin sud*/
int down_skipped(int ix,int iy, int m);
