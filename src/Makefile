PY = python3
CCS_DIR ?= ../
OBJ_DIR=$(CCS_DIR)src/obj/
BUILD ?= release

# these 2 variables are needed to parse "," and " " correctly
comma := ,
space := $(null) $(null)

ifndef CMP
$(warning "CMP is not set. Please define which compiler should be used.")
endif

ifeq ($(CMP),intel)
  FC = mpiifort -fc=ifort
  FFLAGS = -cpp -fPIC -g -O3 -stand f18 -qopenmp-simd
  ifeq ($(BUILD),debug)
    # Add debugging (i.e. expensive) flags
  endif
  ifeq ($(PROFILE),yes)
    FFLAGS += -qopt-report-per-object -qopt-report=5 -qopt-report-annotate=html -qopt-report-phase=ipo,loop,vec -qopt-report-annotate-position=both
  endif
  FFLAGS += -xHOST
  FFLAGS += -warn all,noexternals
  FFLAGS += -module $(OBJ_DIR)
  CAFFLAGS = -coarray
  MPIRUN = mpirun
endif
ifeq ($(CMP),intelx)
  FC = mpiifort -fc=ifx
  FFLAGS = -cpp -fPIC -g -O3 -stand f18
  ifeq ($(BUILD),debug)
    # Add debugging (i.e. expensive) flags
  endif
  FFLAGS += -xHOST
  FFLAGS += -warn all,noexternals
  FFLAGS += -module $(OBJ_DIR)
  CAFFLAGS = -coarray
  MPIRUN = mpirun
endif
ifeq ($(CMP),gnu)
  FC = mpif90
  FFLAGS = -cpp -O3 -std=f2018
  CAFFLAGS = -fcoarray=single
  ifeq ($(BUILD),debug)
    # Add debugging (i.e. expensive) flags
    FFLAGS += -fcheck=bounds
  endif
  FFLAGS += -fopenmp
  FFLAGS += -Wall -Wpedantic -Werror -Wimplicit-interface -Wimplicit-procedure
  FFLAGS += -J$(OBJ_DIR)
  MPIRUN = mpirun
endif
ifeq ($(CMP),cray)
  FC = ftn
  FFLAGS = -ef -en -eF -M969
  ifeq ($(BUILD),debug)
    # Add debugging (i.e. expensive) flags
  endif
  FFLAGS += -homp
  FFLAGS += -J$(OBJ_DIR)
  CAFFLAGS = -hcaf
  MPIRUN = srun
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
TOOLS=$(CCS_DIR)src/build_tools/

DEP_PREFIX=$(CCS_DIR)src/obj/
EXE_DEPS=$(DEP_PREFIX)ccs_app.deps
ALL_DEPS=$(DEP_PREFIX)all.deps
TAG_DEPS=$(DEP_PREFIX)build_tags.deps
RULE_DEPS=$(DEP_PREFIX)rules.deps

IGNORE = "create_boundary_conditions_bfs.f90 create_boundary_conditions_poisson.f90 discretisation_mod.f90 create_boundary_conditions_mod.f90 create_boundary_conditions_tgv.f90 spatial.f90"
CASES = poisson

ALL_SRC = $(shell find $(CCS_DIR)src/ -type f -name '*.f90')
SRC = $(shell $(PY) $(TOOLS)filter-out.py $(IGNORE) "$(ALL_SRC)")
OBJ=$(addprefix $(OBJ_DIR), $(notdir $(SRC:.f90=.o)))

include $(TAG_DEPS)
include $(EXE_DEPS)

FFLAGS += -DACCS_PETSC
INC = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc
INC += -I${FYAML}/build 
LIB += -Wl,-rpath,${FYAML}/build:${FYAML}/build/yaml-cpp -L${FYAML}/build -lfortran-yaml-cpp -L${FYAML}/build/yaml-cpp -lyaml-cpp

all: $(OBJ) $(EXE)

tests: all
	make -C ../testing tests

$(EXE): $(EXE_DEPS)
	$(FC) $(FFLAGS) $(CAFLINK) -o $@ $(filter-out $(EXE_DEPS),$^) $(INC) $(LIB)

#%.o: %.f90
#§	$(FC) $(FFLAGS) -o $@ -c $< $(INC)

$(CAF_OBJ): %.o: %.f90
	$(FC) $(FFLAGS) $(CAFFLAGS) -o $@ -c $< $(INC)

clean:
	rm -f $(EXE) *.o *.mod *.smod
	rm -f $(OBJ_DIR)*.o $(OBJ_DIR)*.mod $(OBJ_DIR)*.smod $(OBJ_DIR)*.html $(OBJ_DIR)*.optrpt 
clean-tests:
	make -C ../testing clean
clean-full: clean clean-tests clean-docs
	rm -f $(OBJ_DIR)*.deps

$(ALL_DEPS): $(SRC)
	makedepf90 -b $(DEP_PREFIX) $(SRC) > $(ALL_DEPS)

$(RULE_DEPS): $(ALL_DEPS)
	$(PY) $(TOOLS)fix_makefile.py $(ALL_DEPS) $(RULE_DEPS) 

$(TAG_DEPS): $(SRC)
	$(PY) $(TOOLS)process_build_tags.py $(SRC) > $(TAG_DEPS)

$(EXE_DEPS): config.yaml $(ALL_DEPS)
	$(PY) $(TOOLS)generate_link_deps.py config.yaml $(ALL_DEPS) $(EXE_DEPS)

#include $(ALL_DEPS)

include $(RULE_DEPS)

docs: doxy 
doxy:
	doxygen doxygen.cfg
docs-latex: doxy
	make -C latex
clean-docs:
	rm -rf html latex

#Needed to pass variables to children Makefiles, e.g. for the testing framework
export
