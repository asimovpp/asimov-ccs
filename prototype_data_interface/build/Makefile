#FC=gfortran-8
FFLAGS=-std=f2018

SRCDIR=../

.SUFFIXES:
.SUFFIXES: .f90 .o

SRC=$(wildcard $(SRCDIR)*.f90)
BASE=$(shell basename -a $(SRC))
OBJ=$(BASE:.f90=.o)

all: $(OBJ)

%.o: $(SRCDIR)%.f90
	$(FC) $(FFLAGS) -c $<

main.o presolve.o flowsolve.o: types.o
main.o: presolve.o flowsolve.o
presolve_basic.o: presolve.o
flux_1order.o flux_2order.o particles_1.o particles_2.o solve_amg.o solve_cgstab.o turbulence_ke.o turbulence_kw.o stepper_2loop stepper_3loop: flowsolve.o presolve.o

TEST_MAIN_FILES= \
				main.o \
				types.o \
				presolve.o \
				flowsolve.o \
				presolve_basic.o \
				turbulence_ke.o \
				solve_amg.o \
				stepper_2loop.o

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
