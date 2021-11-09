rm test_yaml_gnu *.o

gfortran -I/work/Projects/ASiMoV/fortran-yaml-cpp/build -c test_yaml.f90 -o test_yaml.o

gfortran test_yaml.o -o test_yaml -Wl,-rpath,/work/Projects/ASiMoV/fortran-yaml-cpp/build:/work/Projects/ASiMoV/fortran-yaml-cpp/build/yaml-cpp -L/work/Projects/ASiMoV/fortran-yaml-cpp/build -lfortran-yaml-cpp -L/work/Projects/ASiMoV/fortran-yaml-cpp/build/yaml-cpp -lyaml-cpp
