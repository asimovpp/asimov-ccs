#FC=gfortran-8
FFLAGS=-std=f2018 -cpp

SRCDIR=../

.SUFFIXES:
.SUFFIXES: .f90 .o

SRC=$(wildcard $(SRCDIR)*.f90)
BASE=$(shell basename -a $(SRC))
OBJ=$(BASE:.f90=.o)

all: $(OBJ)

%.o: $(SRCDIR)%.f90
	$(FC) $(FFLAGS) -c $<

types.o: turbulence_types.o problem_former_types.o type_headers.o solver_types.o
problem_former_types.o: type_headers.o
turbulence_types.o: type_headers.o
main.o: stepper.o types.o
form_problem.o: problem_former.o
turbulence.o problem_former.o stepper.o: types.o
solver_amg.o solver_cgstab.o: solver.o
solver.o: types.o
stepper_2loop.o: stepper.o turbulence.o solver.o 
stepper_basic.o: stepper.o problem_former.o solver.o 
turbulence_ke.o turbulence_kw.o: turbulence.o

TEST_MAIN_FILES= \
				main.o \
				types.o \
				turbulence_types.o \
				type_headers.o \
				problem_former_types.o \
				problem_former.o \
				form_problem.o \
				solver.o \
				solver_types.o \
				solver_amg.o \
				stepper.o \
				stepper_basic.o \
				form_problem.o

exe:
	python3 link_executable.py config.json

test_main: $(TEST_MAIN_FILES)  
	$(FC) $(FFLAGS) -o $@ $^

clean:
	rm -f *.o *.mod *.smod test_main ccs_main

print:
	echo $(SRC)
	echo $(BASE)
	echo $(OBJ)
