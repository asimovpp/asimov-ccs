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

FC = gfortran
FFLAGS = -std=f2018 -Wall -Wpedantic -Werror

# Dependencies
DEPDIR = ../../src
_DEPS = mod_constants.f90
DEPS = $(patsubst %, $(DEPDIR)/%, $(_DEPS))
_DEPMODS = $(_DEPS:mod_%.f90=%.mod)
DEPMODS = $(patsubst %, $(MODDIR)/%, $(_DEPMODS))
_DEPOBJ = $(_DEPS:%.f90=%.o)
DEPOBJ = $(patsubst %, $(OBJDIR)/%, $(_DEPOBJ))

# Test system modules

# Tests

# Targets
all : deps

deps : depmods depobj
depmods : $(MODDIR) $(DEPMODS)
depobj : depmods $(OBJDIR) $(DEPOBJ)

# Objects
$(OBJDIR) :
	mkdir $(OBJDIR)

$(OBJDIR)/%.o : $(DEPDIR)/%.f90
	$(FC) $(FFLAGS) -c $^ -o $@ -I$(MODDIR)

# Modules
$(MODDIR) :
	mkdir $(MODDIR)

$(MODDIR)/%.mod : %.mod
	mv $< $(MODDIR)/

%.mod : $(DEPDIR)/mod_%.f90
	$(FC) $(FFLAGS) -fsyntax-only $<

clean :
	rm -r $(MODDIR) $(OBJDIR)
