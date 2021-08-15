#
#        FILE: Makefile
# DESCRIPTION: Makefile for parallel PI prototype of ASiMoV-CCS.
#

FC = mpiifort
FFLAGS = -cpp -stand f18
MODFLAGS = -fsyntax-only

EXE = pi

MOD_SRC = parallel_mod.f90 parallel_mpi_env.f90 parallel_mpi_utils.f90 parallel_mpi_collectives.f90
SRC =   $(wildcard *.f90)
OBJ = $(SRC:%.f90=%.o)
MOD = $(MOD_SRC:%_mod.f90=%.mod)

all: $(MOD) $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(INC) $(LIB)

%.mod: %_mod.f90
	$(FC) $(FFLAGS) $(MODFLAGS) $< $(INC)

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $^ $(INC)

clean:
	rm -f $(EXE) *.o *.mod *.smod