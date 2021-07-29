#
#        FILE: Makefile
# DESCRIPTION: Builds/runs test cases for ASiMoV-CCS
#

# Test constants
TESTN = 100       # How many random inputs should a test generate?
TESTSCALE = 1.0e9 # How large to scale test inputs? This should cover any conceivable speed (past
                  # light-speed!) and pressures achieved by compressors with PR > 40 and 10^5 Pa
                  # input.

SRC = ../../../src
INC = ../mod
OBJ = $(wildcard ../obj/*.o)

OPT = -g

# pFUnit
include $(PFUNIT)/include/PFUNIT.mk

MK = make
FC = gfortran
FFLAGS = -cpp -std=f2018 -Wall -Wpedantic -Werror $(OPT)
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
DECL = -DTESTSCALE=$(TESTSCALE) -DTESTN=$(TESTN)

# Tests
discretisation_TESTS = spatial.pf
discretisation_OTHER_SRCS = $(SRC)/discretisation/spatial.f90
discretisation_OTHER_LIBRARIES = $(OBJ) # Yes, weird, I think interpret libraries==objects
$(eval $(call make_pfunit_test,discretisation))

# Targets
all : testsys discretisation run

run : testsys discretisation
	./discretisation

testsys :
	$(MK) -C ../

# Objects
%.o : %.F90
	$(FC) $(FFLAGS) $(DECL) -c $^ -o $@ -I$(INC)
%.o : %.f90
	$(FC) $(FFLAGS) $(DECL) -c $^ -o $@ -I$(INC)

# Modules

clean :
	rm -r *.F90 *.o *.mod *.inc
