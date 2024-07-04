# librairies de SuiteSparse
PREF := /opt/homebrew/opt/suite-sparse
METIS_PREF := /opt/homebrew/opt/metis
OMP_PREF := /opt/homebrew/opt/libomp
L1 = $(PREF)/lib/libumfpack.a
L2 = $(PREF)/lib/libcholmod.a
L3 = $(PREF)/lib/libamd.a
L4 = $(PREF)/lib/libcamd.a  
L5 = $(PREF)/lib/libcolamd.a 
L6 = $(PREF)/lib/libccolamd.a 
L7 = -L$(METIS_PREF)/lib -lmetis
L8 = $(PREF)/lib/libsuitesparseconfig.a
L9 = $(OMP_PREF)/lib/libomp.a
LIB = $(L1) $(L2) $(L3) $(L4) $(L5) $(L6) $(L7) $(L8) $(L9) -lm -lblas -llapack

#L1 = SuiteSparse/UMFPACK/Lib/libumfpack.a
#L2 = SuiteSparse/CHOLMOD/Lib/libcholmod.a 
#L3 = SuiteSparse/AMD/Lib/libamd.a 
#L4 = SuiteSparse/CAMD/Lib/libcamd.a  
#L5 = SuiteSparse/COLAMD/Lib/libcolamd.a 
#L6 = SuiteSparse/CCOLAMD/Lib/libccolamd.a 
#L7 = SuiteSparse/metis-4.0/libmetis.a
#L8 = SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a
#LIB = $(L1) $(L2) $(L3) $(L4) $(L5) $(L6) $(L7) $(L8) -lm -lblas -llapack

INCLUDE_PATHS = -I$(PREF)/include/suitesparse -I$(PREF)/include -I$(METIS_PREFIX)/include 



COPT = -O3

default: main

clean: 
	rm *.o 
	rm main

main: main.c umfpack.o method.o time.o globVal.o plot.o tools.o prob.o cg.o
	cc $(COPT) $^ -o $@ $(LIB)

umfpack.o: umfpack.c
	cc $(COPT) -c $< -o $@ $(INCLUDE_PATHS)
#-ISuiteSparse/UMFPACK/Include -ISuiteSparse/SuiteSparse_config  -ISuiteSparse/AMD/Include

%.o: %.c proto.h
	cc $(COPT) -c $< -o $@ $(INCLUDE_PATHS)


