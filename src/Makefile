# these 2 variables are needed to parse "," and " " correctly
comma := ,
space := $(null) $(null)

# prevent the use of implicit suffix rules
.SUFFIXES:

# NEED_CMP decides if full compilation is required or not
# e.g. if the build is a simple clean, makedepf90 does not need to be invoked
NEED_CMP = yes
ifneq (,$(filter $(MAKECMDGOALS),clean clean-tests clean-full clean-docs docs ford doxy docs-latex))
  NEED_CMP = no
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
SMOD_DEPS=$(DEP_PREFIX)/submodules.deps
TAG_DEPS=$(DEP_PREFIX)/build_tags.deps
RULE_DEPS=$(DEP_PREFIX)/rules.deps

# check makedepf90 version
MAKEDEPF90_SMODS=$(shell makedepf90 -h | grep -q '\-S PATH'; echo $$?)

IGNORE = " "

ALL_SRC = $(shell find $(CCS_DIR)/src/ -type f -name '*.f90')
ifeq ($(CCS_PROPRIETARY),yes)
  ALL_SRC += $(shell find $(CCS_PROPRIETARY_DIR)/src/ -type f -name '*.f90')
endif
SRC = $(shell $(PY) $(TOOLS)/filter_out.py $(IGNORE) "$(ALL_SRC)")
OBJ = $(addprefix $(OBJ_DIR), $(notdir $(SRC:.f90=.o)))

ifeq ($(NEED_CMP),yes)
  include $(TAG_DEPS)
  include $(EXE_DEPS)
endif

INC = -I${CCS_DIR}/src/include
FFLAGS += -DACCS_PETSC
INC += -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc
INC += -I${FYAML}/build 
LIB += -Wl,-rpath,${FYAML}/build:${FYAML}/build/yaml-cpp -L${FYAML}/build -lfortran-yaml-cpp -L${FYAML}/build/yaml-cpp -lyaml-cpp 

ifeq ($(NEED_CMP),yes)
  INC += $(shell $(ADIOS2)/bin/adios2-config --fortran-flags)
  LIB += $(shell $(ADIOS2)/bin/adios2-config --fortran-libs)
endif

# file containing Makefile commands particular to proprietary code
ifeq ($(CCS_PROPRIETARY),yes)
  CCS_PROPRIETARY_MAKEFILE=$(CCS_PROPRIETARY_DIR)/src/Makefile.inc
  include $(CCS_PROPRIETARY_MAKEFILE)
endif



all: obj app

app: $(EXE)

obj: $(OBJ)

tests: FFLAGS+=-DVERBOSE
tests: obj
	make -C $(CCS_DIR)/tests all

$(EXE): $(EXE_DEPS)
	$(FC) $(FFLAGS) $(CAFLINK) -o $@ $(filter-out $(EXE_DEPS),$^) $(INC) $(LIB)

#%.o: %.f90
#§	$(FC) $(FFLAGS) -o $@ -c $< $(INC)

$(CAF_OBJ): %.o: %.f90
	$(FC) $(FFLAGS) $(CAFFLAGS) -o $@ -c $< $(INC)

clean:
	rm -f $(EXE) *.o *.mod *.smod *.deps
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.smod $(DEP_PREFIX)/*.deps $(OBJ_DIR)/*.html $(OBJ_DIR)/*.optrpt $(OBJ_DIR)/*.lst opt_info.txt
clean-tests:
	make -C $(CCS_DIR)/tests clean
clean-full: clean clean-tests clean-docs

$(ALL_DEPS): $(SRC)
	if [ $(MAKEDEPF90_SMODS) = 0 ]; then makedepf90 -b $(OBJ_DIR) -S $(SMOD_DEPS) $(SRC) > $(ALL_DEPS); else makedepf90 -b $(OBJ_DIR) $(SRC) > $(ALL_DEPS); fi

$(RULE_DEPS): $(ALL_DEPS)
	$(PY) $(TOOLS)/fix_makefile.py $(ALL_DEPS) $(RULE_DEPS) 

$(TAG_DEPS): $(SRC)
	$(PY) $(TOOLS)/process_build_tags.py $(SRC) > $(TAG_DEPS)

$(EXE_DEPS): config.yaml $(ALL_DEPS)
	if [ $(MAKEDEPF90_SMODS) = 0 ]; then $(PY) $(TOOLS)/generate_link_deps.py config.yaml $(ALL_DEPS) $(EXE_DEPS) $(SMOD_DEPS); else $(PY) $(TOOLS)/generate_link_deps.py config.yaml $(ALL_DEPS) $(EXE_DEPS); fi

ifeq ($(NEED_CMP),yes)
  #include $(ALL_DEPS)
  include $(RULE_DEPS)
endif

docs: doxy 
doxy:
	doxygen .doxygen.cfg
ford:
	ford -d $(CCS_DIR)/src --exclude_dir doc .documentation_settings.ford
docs-latex: doxy
	make -C latex
clean-docs:
	rm -rf html latex

#Needed to pass variables to children Makefiles, e.g. for the testing framework
export
