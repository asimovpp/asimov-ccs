#
# FILE: Makefile
# DESCRIPTION: Builds/runs test cases for ASiMoV-CCS PETSc/ex3 prototype
#

FC = gfortran
FFLAGS = -cpp

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
$(eval $(call make_pfunit_test,test_vector))

all: ex3 test_vector

ex3:
	make -C ../src

run:
	./test_vector

%.o: %.F90
	$(FC) $(FFLAGS) -c $^ -o $@ -I$(INC)
%.o: %.f90
	$(FC) $(FFLAGS) -c $^ -o $@ -I$(INC)

clean:
	rm *.F90 *.o *.mod *.inc
