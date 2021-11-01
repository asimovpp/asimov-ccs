#
# FILE: tests/linear_solvers/Makefile
# DESCRIPTION: Builds unit tests cases for ASiMoV-CCS linear_solvers/
#

CCS_DIR = ../../../src
OBJ_DIR = $(CCS_DIR)/obj

EXE = ccs_app
EXE_DEPS = $(OBJ_DIR)/ccs_app.deps

USEMPI = YES
include $(PFUNIT)/include/PFUNIT.mk

#
# Define test suites
#
TEST_SUITES =

TEST_SUITE = linear_solvers
$(TEST_SUITE)_TESTS = $(wildcard *.pf)
$(TEST_SUITE)_OTHER_LIBRARIES := $(addprefix $(CCS_DIR)/,$(filter-out obj/$(CASES).o $(EXE):,$(shell cat ${EXE_DEPS}))) -I$(OBJ_DIR) $(LIB)
$(eval $(call make_pfunit_test,$(TEST_SUITE)))
TEST_SUITES += $(TEST_SUITE)

#
# Should not need to edit below line
#
all : $(TEST_SUITES)

%.o : %.F90
	$(FC) $(FFLAGS) -o $@ -c $< $(INC) -I$(OBJ_DIR) $(LIB) $(PFUNIT_EXTRA_FFLAGS)

clean:
	rm *.o *.F90 *.inc *.mod
