FC=gfortran-8
FFLAGS=

all: main1 main2

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

basic_sub1.o basic_sub2.o: basic.o
main.o: basic.o

main1: main.o basic.o basic_sub1.o
	$(FC) $(FFLAGS) -o $@ $^
main2: main.o basic.o basic_sub2.o
	$(FC) $(FFLAGS) -o $@ $^

clean:
	rm -f *.o *.mod *.smod main1 main2
