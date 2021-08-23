#
# FILE: Makefile
# DESCRIPTION: Builds/runs test cases for ASiMoV-CCS PETSc/ex3 prototype
#

USEMPI=YES

FC = mpif90
FFLAGS = -cpp -std=f2018 -Wall -Wpedantic -Werror

SRC := $(wildcard ../src/*.f90)
SRC := $(filter-out ../src/ex3.f90, $(SRC)) # Exclude main program from sources
INC = ../src

INC += -I$(PETSC_DIR)/include

LPETSC = -L$(PETSC_DIR)/lib -lpetsc

include $(PFUNIT)/include/PFUNIT.mk
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)

test_vector_TESTS = test_vector.pf
test_vector_OTHER_SRCS = $(SRC)
test_vector_OTHER_LIBRARIES = $(LPETSC)
test_vector_EXTRA_INITIALIZE = accs_init
test_vector_EXTRA_FINIALIZE = accs_finalise
$(eval $(call make_pfunit_test,test_vector))

test_matrix_TESTS = test_matrix.pf
test_matrix_OTHER_SRCS = $(SRC)
test_matrix_OTHER_LIBRARIES = $(LPETSC)
test_matrix_EXTRA_INITIALIZE = accs_init
test_matrix_EXTRA_FINIALIZE = accs_finalise
$(eval $(call make_pfunit_test,test_matrix))

all: ex3 test_vector test_matrix

ex3:
	make -C ../src

run:
	mpirun ./test_vector
	wait
	mpirun ./text_matrix

%.o: %.F90
	$(FC) $(FFLAGS) -c $^ -o $@ -I$(INC)
%.o: %.f90
	$(FC) $(FFLAGS) -c $^ -o $@ -I$(INC)

clean:
	rm -f *.F90 *.o *.mod *.inc
