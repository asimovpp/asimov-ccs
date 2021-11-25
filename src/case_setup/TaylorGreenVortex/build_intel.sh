FYAML=/work/Projects/ASiMoV/fortran-yaml-cpp

rm test_yaml_intel example_intel *.o

ifort -I${FYAML}/build -c test_yaml.f90 -o test_yaml.o
ifort -I${FYAML}/build -c example.f90 -o example.o

ifort test_yaml.o -o test_yaml_intel -Wl,-rpath,${FYAML}/build:${FYAML}/build/yaml-cpp -L${FYAML}/build -lfortran-yaml-cpp -L${FYAML}/build/yaml-cpp -lyaml-cpp
ifort example.o -o example_intel -Wl,-rpath,${FYAML}/build:${FYAML}/build/yaml-cpp -L${FYAML}/build -lfortran-yaml-cpp -L${FYAML}/build/yaml-cpp -lyaml-cpp
