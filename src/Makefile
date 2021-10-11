PY = python3

ifeq ($(CMP),intel)
  FC = mpiifort
  FFLAGS = -cpp -fPIC -g -O3 -xHOST -warn all,noexternals -stand f18 -module obj
  CAFFLAGS = -coarray
endif
ifeq ($(CMP),gnu)
  FC = mpif90
  FFLAGS = -cpp -O3 -std="f2018" -fopenmp -Wall -Wpedantic -Werror -Wimplicit-interface -Wimplicit-procedure -Jobj
  CAFFLAGS = -fcoarray=single
endif
ifeq ($(CMP),cray)
  FC = ftn
  FFLAGS = -ef -en -eF -M969 -homp -Jobj
  CAFFLAGS = -hcaf
endif

# Only set this value if building a CAF binary
# otherwise keep unset
#ifeq ($(CMP),intel)
#CAFLINK= #-coarray
#endif
#ifeq ($(CMP),gnu)
#CAFLINK= #-fcoarray=single
#endif


EXE = ccs_app
TOOLS=build_tools/

DEP_PREFIX=obj/
EXE_DEPS=$(DEP_PREFIX)ccs_app.deps
ALL_DEPS=$(DEP_PREFIX)all.deps
TAG_DEPS=$(DEP_PREFIX)build_tags.deps
RULE_DEPS=$(DEP_PREFIX)rules.deps

OBJ_DIR=obj/

IGNORE = "create_boundary_conditions_bfs.f90 create_boundary_conditions_poisson.f90 discretisation_mod.f90 create_boundary_conditions_mod.f90 create_boundary_conditions_tgv.f90 spatial.f90"

#SRC = $(wildcard *.f90)
ALL_SRC = $(shell find . -type f -name '*.f90')
SRC = $(shell $(PY) $(TOOLS)filter-out.py $(IGNORE) "$(ALL_SRC)")
#OBJ = $(SRC:%.f90=%.o)
OBJ=$(addprefix $(OBJ_DIR), $(notdir $(SRC:.f90=.o)))

include $(TAG_DEPS)
include $(EXE_DEPS)

FFLAGS += -DACCS_PETSC
INC = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc

all: $(OBJ) $(EXE)

#test: $(EXE)
#	make -C ../tests all

$(EXE): $(EXE_DEPS)
	$(FC) $(FFLAGS) $(CAFLINK) -o $@ $(filter-out $(EXE_DEPS),$^) $(INC) $(LIB)

#%.o: %.f90
#§	$(FC) $(FFLAGS) -o $@ -c $< $(INC)

$(CAF_OBJ): %.o: %.f90
	$(FC) $(FFLAGS) $(CAFFLAGS) -o $@ -c $< $(INC)

clean:
	rm -f $(EXE) *.o *.mod *.smod
	rm -f $(OBJ_DIR)*.o $(OBJ_DIR)*.mod $(OBJ_DIR)*.smod
clean-full: clean
	rm -f $(OBJ_DIR)*.deps

$(ALL_DEPS): $(SRC)
	makedepf90 $(SRC) > $(DEP_PREFIX)tmp_all.deps
	$(PY) $(TOOLS)send_to_obj.py $(DEP_PREFIX)tmp_all.deps $(ALL_DEPS) 

$(RULE_DEPS): $(ALL_DEPS)
	$(PY) $(TOOLS)fix_makefile.py $(ALL_DEPS) $(RULE_DEPS) 

$(TAG_DEPS): $(SRC)
	$(PY) $(TOOLS)process_build_tags.py $(SRC) > $(TAG_DEPS)

$(EXE_DEPS): config.json $(ALL_DEPS)
	$(PY) $(TOOLS)generate_link_deps.py config.json $(ALL_DEPS) $(EXE_DEPS)

#include $(ALL_DEPS)

include $(RULE_DEPS)

docs: doxy 
doxy:
	doxygen doxygen.cfg
