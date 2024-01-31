# these 2 variables are needed to parse "," and " " correctly
comma := ,
space := $(null) $(null)

# print input and execute it as well
printdo = echo $(1); $1

# prevent the use of implicit suffix rules
.SUFFIXES:

# NEED_CMP decides if full compilation is required or not
# e.g. if the build is a simple clean, makedepf90 does not need to be invoked
NEED_CMP = yes
ifneq (,$(filter $(MAKECMDGOALS),clean clean-tests clean-full clean-docs docs ford dev_guide))
  NEED_CMP = no
endif

PY = python3
CCS_DIR ?= $(shell realpath $(PWD))
ARCH_DIR=$(CCS_DIR)/build_tools/archs
OBJ_DIR=$(CCS_DIR)/obj
BUILD ?= release

# this enables adding more dirs via `make SRC_DIRS="a b"` (must use absolute paths)
override SRC_DIRS += $(CCS_DIR)/src $(CCS_DIR)/build_tools

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
TOOLS=$(CCS_DIR)/build_tools

DEP_PREFIX=$(OBJ_DIR)
EXE_DEPS=$(DEP_PREFIX)/ccs_app.deps
ALL_DEPS=$(DEP_PREFIX)/all.deps
SMOD_DEPS=$(DEP_PREFIX)/submodules.deps
TAG_DEPS=$(DEP_PREFIX)/build_tags.deps
RULE_DEPS=$(DEP_PREFIX)/rules.deps

LIB_CCS=libccs.a

# check makedepf90 version
MAKEDEPF90_SMODS=$(shell makedepf90 -h | grep -q '\-S PATH'; echo $$?)

IGNORE = " "

ifeq ($(CCS_PROPRIETARY),yes)
  SRC_DIRS += $(CCS_PROPRIETARY_DIR)/src
endif

find_src_files = $(shell find $(dir) -type f -name '*.f90' -o -name '*.c')
ALL_SRC = $(foreach dir, $(SRC_DIRS), $(find_src_files))

SRC = $(shell $(PY) $(TOOLS)/filter_out.py $(IGNORE) "$(ALL_SRC)")
TMP_OBJ = $(addprefix $(OBJ_DIR)/, $(notdir $(SRC:.f90=.o)))
OBJ = $(TMP_OBJ:.c=.o)

ifeq ($(NEED_CMP),yes)
  include $(TAG_DEPS)
  include $(EXE_DEPS)
endif

INC = -I${CCS_DIR}/include
FFLAGS += -DACCS_PETSC
INC += -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc

INC += -I${RCMF90}/include
LIB += -L${RCMF90}/lib -lrcm

INC += -I${FYAMLC}/modules 
LIB += -Wl,-rpath,${FYAMLC}/lib -L${FYAMLC}/lib -lfortran-yaml-c

ifdef PARHIP
  INC += -I${PARHIP}/include
  LIB += -L${PARHIP}/lib -lparhip_interface -Wl,-rpath,${PARHIP}/lib
endif

ifdef PARMETIS
  INC += -I${PARMETIS}/include
  LIB += -L${PARMETIS}/lib -lGKlib -lmetis -lparmetis -Wl,-rpath,${PARMETIS}/lib
endif

ifeq ($(NEED_CMP),yes)
  INC += $(shell $(ADIOS2)/bin/adios2-config --fortran-flags)
  LIB += $(shell $(ADIOS2)/bin/adios2-config --fortran-libs)
endif

# file containing Makefile commands particular to proprietary code
ifeq ($(CCS_PROPRIETARY),yes)
  CCS_PROPRIETARY_MAKEFILE=$(CCS_PROPRIETARY_DIR)/src/Makefile.inc
  include $(CCS_PROPRIETARY_MAKEFILE)
endif

ifeq ($(BUILD),debug)
  FFLAGS += -DEXCLUDE_MISSING_INTERFACE
endif

all: obj lib app

lib: $(EXE_DEPS)
	rm -f $(LIB_CCS)
	ar rcs $(LIB_CCS) $(filter-out $(EXE_DEPS),$^)

app: $(EXE)

obj: $(OBJ)

all: obj app

$(EXE): $(EXE_DEPS)
	$(FC) $(FFLAGS) $(CAFLINK) -o $@ $(filter-out $(EXE_DEPS),$^) $(INC) $(LIB) 
	@echo -n "===> Built ccs_app with "
	@grep main $(CCS_DIR)/config.yaml

COMPILE_FORTRAN = $(call printdo, $(FC) $(FFLAGS) -o $@ -c $< $(INC))
COMPILE_C =       $(call printdo, $(CC) $(CFLAGS) -o $@ -c $< $(INC))
$(OBJ_DIR)/%.o: 
	@if [ $(suffix $<) = .f90 ]; then $(COMPILE_FORTRAN); elif [ $(suffix $<) = .c ]; then $(COMPILE_C); fi

$(CAF_OBJ): %.o: %.f90
	$(FC) $(FFLAGS) $(CAFFLAGS) -o $@ -c $< $(INC)

$(TAG_DEPS): $(SRC)
	$(PY) $(TOOLS)/process_build_tags.py $(OBJ_DIR) $(SRC) > $(TAG_DEPS)

GEN_DEPS         = $(call printdo, makedepf90 -b $(OBJ_DIR)                 $(SRC) > $(ALL_DEPS))
GEN_DEPS_W_SMODS = $(call printdo, makedepf90 -b $(OBJ_DIR) -S $(SMOD_DEPS) $(SRC) > $(ALL_DEPS))
$(ALL_DEPS): $(SRC)
	@if [ $(MAKEDEPF90_SMODS) = 0 ]; then $(GEN_DEPS_W_SMODS) ; else $(GEN_DEPS) ; fi

GEN_LINK_LINE         = $(call printdo, $(PY) $(TOOLS)/generate_link_deps.py config.yaml $(ALL_DEPS) $(EXE_DEPS)             )
GEN_LINK_LINE_W_SMODS = $(call printdo, $(PY) $(TOOLS)/generate_link_deps.py config.yaml $(ALL_DEPS) $(EXE_DEPS) $(SMOD_DEPS))
$(EXE_DEPS): config.yaml $(ALL_DEPS)
	@if [ $(MAKEDEPF90_SMODS) = 0 ]; then $(GEN_LINK_LINE_W_SMODS) ; else $(GEN_LINK_LINE) ; fi

.PHONY: tests
tests: FFLAGS+=-DVERBOSE
tests: obj
	make -C $(CCS_DIR)/tests all

ifeq ($(NEED_CMP),yes)
  include $(ALL_DEPS)
endif

.PHONY: docs dev_guide ford
docs: ford dev_guide
ford:
	ford .project_documentation_settings.md
dev_guide:
	make -C dev_guide all
clean-docs:
	rm -rf doc
	make -C dev_guide clean
  

clean:
	rm -f $(EXE) $(LIB_CCS) *.o *.mod *.smod *.deps
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.smod $(DEP_PREFIX)/*.deps $(OBJ_DIR)/*.html $(OBJ_DIR)/*.optrpt $(OBJ_DIR)/*.lst opt_info.txt
clean-tests:
	make -C $(CCS_DIR)/tests clean
clean-full: clean clean-tests clean-docs

#Needed to pass variables to children Makefiles, e.g. for the testing framework
export
