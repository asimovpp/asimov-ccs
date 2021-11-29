FYAML=/work/Projects/ASiMoV/fortran-yaml-cpp

rm tgv_intel *.o

ifort -cpp -I${FYAML}/build -c ../../ccs_core/kinds_mod.f90 -o kinds_mod.o
ifort -I${FYAML}/build -c ../../io/read_config_mod.f90 -o read_config_mod.o
ifort -I${FYAML}/build -c ../../io/read_config_yaml.f90 -o read_config_yaml.o
ifort -I${FYAML}/build -c tgv.f90 -o tgv.o

ifort kinds_mod.o read_config_mod.o read_config_yaml.o tgv.o -o tgv_intel -Wl,-rpath,${FYAML}/build:${FYAML}/build/yaml-cpp -L${FYAML}/build -lfortran-yaml-cpp -L${FYAML}/build/yaml-cpp -lyaml-cpp
