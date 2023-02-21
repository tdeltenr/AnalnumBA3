static char help[] = "Solves a linear system\n\n";
/*T
   Concepts: KSP^solving a system of linear equations
   Processors: 1
T*/
/* DISCLAIMER : cette fonction est largement empruntée à la doc du solveur puis modifié pour s'adapter au problème. Le crédit leur revient donc ;) */


#include "petscksp.h"

int petsc_solve(int argc,char **args,int na, int *ia, int *ja, double *a, 
                  double *b, double **x_iterative)
/* Résoudre de manière itérative le problème pour trouver le vecteur x*/
{
  Vec            x,Vec_b;      /* approx solution, RHS, exact solution */
  Mat            A;            /* linear system matrix */
  KSP            ksp;          /* linear solver context */
  PC             pc;           /* preconditioner context */
  PetscReal      norm;         /* norm of solution error */
  PetscErrorCode ierr;
  PetscInt       tolerance = 1e-5, n = na, its;
  PetscScalar	 *x_array;
  PetscMPIInt    size;

ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

//Create vectors 
ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
ierr = VecSetFromOptions(x);CHKERRQ(ierr);
ierr = VecDuplicate(x,&Vec_b);CHKERRQ(ierr);
ierr = VecPlaceArray(Vec_b,b);

//Create the matrix using ia, ja, a
ierr = MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, na,na,ia,ja,a, &A);
  
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              Create the linear solver and set various options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

/*
   Set operators. Here the matrix that defines the linear system
   also serves as the matrix that defines the preconditioner.
*/
ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

/*
   Set linear solver defaults for this problem (optional).
   - By extracting the KSP and PC contexts from the KSP context,
     we can then directly call any KSP and PC routines to set
     various options.
   - The following four statements are optional; all of these
     parameters could alternatively be specified at runtime via
     KSPSetFromOptions();
*/
ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
ierr = KSPSetTolerances(ksp,tolerance,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

/*
  Set runtime options, e.g.,
      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
  These options will override those specified above as long as
  KSPSetFromOptions() is called _after_ any other customization
  routines.
*/
ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Solve the linear system
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
ierr = KSPSolve(ksp,Vec_b,x);CHKERRQ(ierr);

/*
   View solver info; we could instead use the option -ksp_view to
   print this info to the screen at the conclusion of KSPSolve().
*/
//ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

VecGetArray(x,&x_array);
ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
ierr = PetscPrintf(PETSC_COMM_WORLD,
       "La solution converge en %d itérations \n",its);CHKERRQ(ierr);


double value;
for(int k = 0; k<n; k++){
	value = x_array[k];
	(*x_iterative)[k] = value;
}

KSPDestroy(&ksp);
VecRestoreArray(x,&x_array);
VecDestroy(&Vec_b);
MatDestroy(&A);
VecDestroy(&x); 
PetscFinalize();
  return ierr;
}

