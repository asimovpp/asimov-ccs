#
# FILE: tests/mesh/Makefile
# DESCRIPTION: Builds unit tests cases for ASiMoV-CCS mesh/
#

CCS_DIR = ../../../src
OBJ_DIR = $(CCS_DIR)/obj

EXE = ccs_app
EXE_DEPS = $(OBJ_DIR)/ccs_app.deps


#
# Define test suites
#
TEST_SUITE = mesh
TEST_SRC = test_mesh.f90
TEST_OBJ = $(TEST_SRC:%.f90=%.o)

#
# Should not need to edit below line
#
all : $(TEST_SUITE)
$(TEST_SUITE) : $(TEST_OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(addprefix $(CCS_DIR)/,$(filter-out obj/$(CASES).o $(EXE):,$(shell cat ${EXE_DEPS}))) $(CCS_OBJ) -I. -I$(OBJ_DIR) $(INC) $(LIB)

%.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $< -I$(OBJ_DIR) $(INC)

clean:
	rm *.o *.mod $(TEST_SUITE)
