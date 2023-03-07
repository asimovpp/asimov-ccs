FC = gfortran
KFC = gfortran
LKR = gfortran

all: kernel prog

kernel:
	$(KFC) -c poisson_discretisation_kernel.f90 -o poisson_discretisation_kernel.o

prog:
	$(FC) -c poisson_discretisation_bindc.f90 -o poisson_discretisation_bindc.o
	$(LKR) -o poisson_discretisation poisson_discretisation_bindc.o poisson_discretisation_kernel.o
