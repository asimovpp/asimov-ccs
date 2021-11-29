FYAML=/work/Projects/ASiMoV/fortran-yaml-cpp

rm tgv_intel *.o

ifort -I${FYAML}/build -c ../../io/read_yaml_mod.f90 -o read_yaml_mod.o
ifort -I${FYAML}/build -c ../../io/read_yaml_utils.f90 -o read_yaml_utils.o
ifort -I${FYAML}/build -c tgv.f90 -o tgv.o

ifort read_yaml_mod.o read_yaml_utils.o tgv.o -o tgv_intel -Wl,-rpath,${FYAML}/build:${FYAML}/build/yaml-cpp -L${FYAML}/build -lfortran-yaml-cpp -L${FYAML}/build/yaml-cpp -lyaml-cpp
