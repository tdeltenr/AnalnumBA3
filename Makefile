# librairies de SuiteSparse
L1 = SuiteSparse/UMFPACK/Lib/libumfpack.a
L2 = SuiteSparse/CHOLMOD/Lib/libcholmod.a 
L3 = SuiteSparse/AMD/Lib/libamd.a 
L4 = SuiteSparse/CAMD/Lib/libcamd.a  
L5 = SuiteSparse/COLAMD/Lib/libcolamd.a 
L6 = SuiteSparse/CCOLAMD/Lib/libccolamd.a 
L7 = SuiteSparse/metis-4.0/libmetis.a
L8 = SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a

# Préparation de PETSC 
PETSC_DIR = petsc

LIB = $(L1) $(L2) $(L3) $(L4) $(L5) $(L6) $(L7) $(L8) -lm -lblas -llapack 

 

COPT = -O3 -Wall 
#COPT = -g

default: main

clean: 
	rm *.o 
	rm main
	rm data.txt

main: main.c prob.o time.o umfpack.o rho.o find_norm.o plot.o flux.o petsc_solve.o
	cc $(COPT) $^ -o $@ $(LIB) $(LDLIBS) -Ipetsc/include  \
 -Ipetsc/arch-linux-c-debug/include -lpetsc

umfpack.o: umfpack.c
	cc $(COPT) -c $< -o $@ -ISuiteSparse/UMFPACK/Include \
  -ISuiteSparse/SuiteSparse_config  -ISuiteSparse/AMD/Include

petsc_solve.o : petsc_solve.c
	cc $(COPT) -c $< -o $@ -Ipetsc/include  \
 -Ipetsc/arch-linux-c-debug/include -lpetsc

%.o: %.c
	cc $(COPT) -c $< -o $@


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/test


