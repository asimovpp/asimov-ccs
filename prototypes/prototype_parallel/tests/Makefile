#
# FILE: Makefile
# DESCRIPTION: Builds/runs test cases for parallel ASiMoV-CCS pi prototype
#

USEMPI=YES

FC = mpiifort
FFLAGS = -cpp -stand f18

SRC := $(wildcard ../src/*.f90)
SRC := $(filter-out ../src/pi.f90, $(SRC)) # Exclude main program from sources
INC = ../src

include $(PFUNIT)/include/PFUNIT.mk
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)

test_parallel_TESTS = test_parallel.pf
test_parallel_OTHER_SRCS = $(SRC)
$(eval $(call make_pfunit_test,test_parallel))

all: pi test_parallel

pi:
	make -C ../src

run:
	./test_parallel

%.o: %.F90
	$(FC) $(FFLAGS) -c $^ -o $@ -I$(INC)
%.o: %.f90
	$(FC) $(FFLAGS) -c $^ -o $@ -I$(INC)

clean:
	rm -f *.F90 *.o *.mod *.inc test_parallel
