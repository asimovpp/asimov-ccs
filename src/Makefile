# these 2 variables are needed to parse "," and " " correctly
comma := ,
space := $(null) $(null)

# prevent the use of implicit suffix rules
.SUFFIXES:

<<<<<<< HEAD
ifeq ($(CMP),intel)
  FC = mpiifort -fc=ifort
  FFLAGS = -cpp -fPIC -g -O3 -stand f18
  ifeq ($(BUILD),debug)
    # Add debugging (i.e. expensive) flags
    -g -Og -traceback -check all -check bounds -check uninit -frapuv -debug all
  endif
  FFLAGS += -xHOST
  FFLAGS += -warn all,noexternals
  FFLAGS += -module $(OBJ_DIR)
  CAFFLAGS = -coarray
  MPIRUN = mpirun
=======
NEED_CMP = yes
ifneq (,$(filter $(MAKECMDGOALS),clean clean-tests clean-full clean-docs docs doxy docs-latex))
  NEED_CMP = no
>>>>>>> develop
endif

PY = python3
CCS_DIR ?= $(shell realpath $(PWD)/..)/
ARCH_DIR=$(CCS_DIR)/src/build_tools/archs/
OBJ_DIR=$(CCS_DIR)/src/obj/
BUILD ?= release

# this can be set to 'yes' in order to include proprietary code
CCS_PROPRIETARY ?= no

# path to proprietary repository
ifeq ($(CCS_PROPRIETARY),yes)
  ifndef CCS_PROPRIETARY_DIR
    $(error Error: Please set CCS_PROPRIETARY_DIR environment variable)
  endif
endif

ifeq ($(NEED_CMP),yes)
  ifndef CMP
    $(error "Error: CMP is not set. Please define which compiler should be used.")
  endif
  include $(ARCH_DIR)/Makefile.$(CMP)
endif


EXE = ccs_app
TOOLS=$(CCS_DIR)/src/build_tools/

DEP_PREFIX=$(OBJ_DIR)
EXE_DEPS=$(DEP_PREFIX)/ccs_app.deps
ALL_DEPS=$(DEP_PREFIX)/all.deps
TAG_DEPS=$(DEP_PREFIX)/build_tags.deps
RULE_DEPS=$(DEP_PREFIX)/rules.deps

IGNORE = "create_boundary_conditions_bfs.f90 create_boundary_conditions_poisson.f90 discretisation_mod.f90 create_boundary_conditions_mod.f90 create_boundary_conditions_tgv.f90 spatial.f90"
CASES = poisson

ALL_SRC = $(shell find $(CCS_DIR)/src/ -type f -name '*.f90')
ifeq ($(CCS_PROPRIETARY),yes)
  ALL_SRC += $(shell find $(CCS_PROPRIETARY_DIR)/src/ -type f -name '*.f90')
endif
SRC = $(shell $(PY) $(TOOLS)/filter-out.py $(IGNORE) "$(ALL_SRC)")
OBJ = $(addprefix $(OBJ_DIR), $(notdir $(SRC:.f90=.o)))

ifeq ($(NEED_CMP),yes)
  include $(TAG_DEPS)
  include $(EXE_DEPS)
endif

FFLAGS += -DACCS_PETSC
INC = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc
INC += -I${FYAML}/build 
LIB += -Wl,-rpath,${FYAML}/build:${FYAML}/build/yaml-cpp -L${FYAML}/build -lfortran-yaml-cpp -L${FYAML}/build/yaml-cpp -lyaml-cpp 
INC += $(shell $(ADIOS2)/bin/adios2-config --fortran-flags)
LIB += $(shell $(ADIOS2)/bin/adios2-config --fortran-libs)

# file containing Makefile commands particular to proprietary code
ifeq ($(CCS_PROPRIETARY),yes)
  CCS_PROPRIETARY_MAKEFILE=$(CCS_PROPRIETARY_DIR)/src/Makefile.inc
  include $(CCS_PROPRIETARY_MAKEFILE)
endif




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
	rm -f $(EXE) *.o *.mod *.smod *.deps
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.smod $(DEP_PREFIX)/*.deps $(OBJ_DIR)/*.html $(OBJ_DIR)/*.optrpt
clean-tests:
	make -C ../testing clean
clean-full: clean clean-tests clean-docs

$(ALL_DEPS): $(SRC)
	makedepf90 -b $(OBJ_DIR) $(SRC) > $(ALL_DEPS)

$(RULE_DEPS): $(ALL_DEPS)
	$(PY) $(TOOLS)/fix_makefile.py $(ALL_DEPS) $(RULE_DEPS) 

$(TAG_DEPS): $(SRC)
	$(PY) $(TOOLS)/process_build_tags.py $(SRC) > $(TAG_DEPS)

$(EXE_DEPS): config.yaml $(ALL_DEPS)
	$(PY) $(TOOLS)/generate_link_deps.py config.yaml $(ALL_DEPS) $(EXE_DEPS)

ifeq ($(NEED_CMP),yes)
  #include $(ALL_DEPS)
  include $(RULE_DEPS)
endif

docs: doxy 
doxy:
	doxygen doxygen.cfg
docs-latex: doxy
	make -C latex
clean-docs:
	rm -rf html latex

#Needed to pass variables to children Makefiles, e.g. for the testing framework
export
