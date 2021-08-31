
PY = python
FC = mpif90
FFLAGS = -cpp -std="f2018" -fcoarray=single -fopenmp #-Wall -Wpedantic -Werror
MODFLAGS = -fsyntax-only
DEPFLAGS = -M -MD -MF

EXE = ccs_app
EXE_DEPS=ccs_app.deps

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

%_caf.mod %_caf.o: %_caf.f90 %_caf.d
	$(FC) $(FFLAGS) -fcoarray=single -o $@ -c $< $(INC)

%.mod %.smod %_mod.o: %_mod.f90 %_mod.d
	$(FC) $(FFLAGS) -c $< $(INC)

%.mod %_mod.o: %_mod.f90 %_mod.d
	$(FC) $(FFLAGS) -c $< $(INC)

%.o: %.f90 %.d
	$(FC) $(FFLAGS) -o $@ -c $< $(INC)

%_caf.o: %_caf.f90 %_caf.d
	$(FC) $(FFLAGS) -fcoarray=single -o $@ -c $< $(INC)

%.d: %.f90
	-@$(FC) $(FFLAGS) $(DEPFLAGS) $@ $< $(INC)

clean:
	rm -f $(EXE) *.o *.mod *.smod
clean-full: clean
	rm -f *.d $(EXE_DEPS)

$(DEP):

-include $(DEP)

-include $(EXE_DEPS)
