#
#        FILE: Makefile
# DESCRIPTION: Makefile for PETSc ex3 prototype of ASiMoV-CCS, based on guide at
#              http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
#

.SUFFIXES:

FC = mpif90
FFLAGS = -cpp -std="f2018" -Wall -Wpedantic -Werror
MODFLAGS = -fsyntax-only
DEPFLAGS = -M -MD -MF

EXE = ex3

SRC = $(wildcard *.f90)
OBJ = $(SRC:%.f90=%.o)
DEP = $(SRC:%.f90=%.d)

BUILDOPTS = -DACCS_PETSC
INC = -I$(PETSC_DIR)/include
LIB = -L$(PETSC_DIR)/lib -lpetsc

all: $(EXE)

test: $(EXE)
	make -C ../tests all

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) $(BUILDOPTS) -o $@ $^ $(INC) $(LIB)

%.mod %.smod %._mod.o: %_mod.f90 %_mod.d
	$(FC) $(FFLAGS) $(BUILDOPTS) -c $< $(INC)

%.mod %._mod.o: %_mod.f90 %_mod.d
	$(FC) $(FFLAGS) $(BUILDOPTS) -c $< $(INC)

%.o: %.f90 %.d
	$(FC) $(FFLAGS) $(BUILDOPTS) -o $@ -c $< $(INC)

%.d: %.f90
	-@$(FC) $(FFLAGS) $(BUILDOPTS) $(DEPFLAGS) $@ $< $(INC)

clean:
	rm -f $(EXE) *.o *.mod *.smod
clean-full: clean
	rm -f *.d

$(DEP):

-include $(DEP)
