
PY = python3
FC = mpif90
FFLAGS = -cpp -std="f2018" -fopenmp #-Wall -Wpedantic -Werror
MODFLAGS = -fsyntax-only
CAFFLAGS = -fcoarray=single

EXE = ccs_app
EXE_DEPS=ccs_app.deps
ALL_DEPS=all.deps
TAG_DEPS=build_tags.deps

SRC = $(wildcard *.f90)
OBJ = $(SRC:%.f90=%.o)
DEP = $(SRC:%.f90=%.d)
include $(TAG_DEPS)

FFLAGS += -DACCS_PETSC
#INC = -I$(PETSC_DIR)/include
#LIB = -L$(PETSC_DIR)/lib -lpetsc

all: $(OBJ) $(EXE)

test: $(EXE)
	make -C ../tests all

$(EXE): $(EXE_DEPS)
	$(FC) $(FFLAGS) -o $@ $(filter-out $(EXE_DEPS),$^) $(INC) $(LIB)

%.mod %.smod %_mod.o: %_mod.f90
	$(FC) $(FFLAGS) -c $< $(INC)

%.mod %_mod.o: %_mod.f90
	$(FC) $(FFLAGS) -c $< $(INC)

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(INC)

$(CAF_OBJ): %.o: %.f90
	$(FC) $(FFLAGS) $(CAFFLAGS) -o $@ -c $< $(INC)

clean:
	rm -f $(EXE) *.o *.mod *.smod
clean-full: clean
	rm -f *.deps $(EXE_DEPS)

$(ALL_DEPS): $(SRC)
	makedepf90 $(SRC) > $(ALL_DEPS)

$(TAG_DEPS): $(SRC)
	$(PY) process_build_tags.py $(SRC) > $(TAG_DEPS)

$(EXE_DEPS): config.json $(ALL_DEPS)
	$(PY) generate_link_deps.py $^

include $(ALL_DEPS)

include $(EXE_DEPS)

