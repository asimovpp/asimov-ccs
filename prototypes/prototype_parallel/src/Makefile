
PY = python

ifeq ($(CMP),intel)
FC = mpiifort
FFLAGS = -fpp -O3 -stand f18 -fopenmp -warn all
CAFFLAGS = -coarray=single
MODFLAGS = -fsyntax-only
endif
ifeq ($(CMP),gnu)
FC = mpif90
FFLAGS = -cpp -O3 -std="f2018" -fopenmp -Wall
CAFFLAGS = -fcoarray=single
MODFLAGS = -fsyntax-only
DEPFLAGS = -M -MD -MF
endif
ifeq ($(CMP),cray)
FC = ftn
FFLAGS = -ef -en -eF -M969 -homp
CAFFLAGS = -hcaf
DEPFLAGS = -eP
endif


#FC = mpif90
#FFLAGS = -cpp -std="f2018" -fcoarray=single -fopenmp  #-Wall -Wpedantic -Werror
#MODFLAGS = -fsyntax-only

EXE = ccs_app
EXE_DEPS=ccs_app.deps
ALL_DEPS=all.deps

SRC = $(wildcard *.f90)
OBJ = $(SRC:%.f90=%.o)
DEP = $(SRC:%.f90=%.d)

FFLAGS += -DACCS_PETSC
#INC = -I$(PETSC_DIR)/include
#LIB = -L$(PETSC_DIR)/lib -lpetsc

all: $(OBJ) $(EXE)

test: $(EXE)
	make -C ../tests all

$(EXE_DEPS): config.json
	$(PY) generate_link_deps.py $<

$(EXE): $(EXE_DEPS)
	$(FC) $(FFLAGS) -o $@ $(filter-out $(EXE_DEPS),$^) $(INC) $(LIB)

%_caf.mod %_caf.o: %_caf.f90
	$(FC) $(FFLAGS) $(CAFFLAGS) -o $@ -c $< $(INC)

%.mod %.smod %_mod.o: %_mod.f90
	$(FC) $(FFLAGS) -c $< $(INC)

%.mod %_mod.o: %_mod.f90
	$(FC) $(FFLAGS) -c $< $(INC)

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(INC)

%_caf.o: %_caf.f90
	$(FC) $(FFLAGS) $(CAFFLAGS) -o $@ -c $< $(INC)

clean:
	rm -f $(EXE) *.o *.mod *.smod
clean-full: clean
	rm -f *.deps $(EXE_DEPS)

depend $(ALL_DEPS):
	makedepf90 *.f90 > $(ALL_DEPS)

include $(ALL_DEPS)

include $(EXE_DEPS)
