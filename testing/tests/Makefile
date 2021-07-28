#
#        FILE: Makefile
# DESCRIPTION: Builds/runs test cases for ASiMoV-CCS
#

# Test constants
TESTN = 100       # How many random inputs should a test generate?
TESTSCALE = 1.0e9 # How large to scale test inputs? This should cover any conceivable speed (past
                  # light-speed!) and pressures achieved by compressors with PR > 40 and 10^5 Pa
                  # input.

OBJDIR=obj
MODDIR=mod
INC = -I$(MODDIR)

OPT = -g

FC = gfortran
FFLAGS = -cpp -std=f2018 -Wall -Wpedantic -Werror $(OPT)
DECL = -DTESTSCALE=$(TESTSCALE) -DTESTN=$(TESTN)

# Dependencies
DEPDIR = ../../src
_DEPS = mod_constants.f90
DEPS = $(patsubst %, $(DEPDIR)/%, $(_DEPS))
_DEPMODS = $(_DEPS:mod_%.f90=%.mod)
DEPMODS = $(patsubst %, $(MODDIR)/%, $(_DEPMODS))
_DEPOBJ = $(_DEPS:%.f90=%.o)
DEPOBJ = $(patsubst %, $(OBJDIR)/%, $(_DEPOBJ))

# Test system modules
TESTDEPS = mod_accs_test_utils.f90
_TESTMODS = $(TESTDEPS:mod_%.f90=%.mod)
TESTMODS = $(patsubst %, $(MODDIR)/%, $(_TESTMODS))
_TESTOBJ = $(TESTDEPS:%.f90=%.o)
TESTOBJ = $(patsubst %, $(OBJDIR)/%, $(_TESTOBJ))

# Targets
all : deps testdeps

deps : depmods depobj
depmods : $(MODDIR) $(DEPMODS)
depobj : depmods $(OBJDIR) $(DEPOBJ)

testdeps : testmods testobj
testmods : $(MODDIR) $(TESTMODS)
testobj : testmods $(OBJDIR) $(TESTOBJ)

# Objects
$(OBJDIR) :
	mkdir $(OBJDIR)

$(OBJDIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(DECL) -c $^ -o $@ -I$(MODDIR)
$(OBJDIR)/%.o : $(DEPDIR)/%.f90
	$(FC) $(FFLAGS) $(DECL) -c $^ -o $@ -I$(MODDIR)

# Modules
$(MODDIR) :
	mkdir $(MODDIR)

$(MODDIR)/%.mod : %.mod
	mv $< $(MODDIR)/

%.mod : mod_%.f90
	$(FC) $(FFLAGS) $(DECL) -fsyntax-only $< $(INC)
%.mod : $(DEPDIR)/mod_%.f90
	$(FC) $(FFLAGS) $(DECL) -fsyntax-only $< $(INC)

clean :
	rm -r $(MODDIR) $(OBJDIR)
