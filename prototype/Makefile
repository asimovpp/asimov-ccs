#FC=gfortran-8
FFLAGS=

.SUFFIXES:
.SUFFIXES: .f90 .o

SRC=$(wildcard *.f90)
OBJ=$(SRC:.f90=.o)

all: $(OBJ)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

main.o: presolve.o postsolve.o flowsolve.o
postsolve_basic.o: postsolve.o
presolve_basic.o: presolve.o
flux_1order.o flux_2order.o particles_1.o particles_2.o solve_amg.o solve_cgstab.o turbulence_ke.o turbulence_kw.o: flowsolve.o

TEST_MAIN_FILES= \
				main.o \
				presolve.o \
				postsolve.o \
				flowsolve.o \
				presolve_basic.o \
				postsolve_basic.o \
				turbulence_ke.o \
				particles_1.o \
				flux_1order.o \
				solve_amg.o

test_main: $(TEST_MAIN_FILES)  
	$(FC) $(FFLAGS) -o $@ $^

clean:
	rm -f *.o *.mod *.smod test_main ccs_main
