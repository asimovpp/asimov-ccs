rm test_yaml_intel *.o

ifort -I/work/Projects/ASiMoV/fortran-yaml-cpp/build -c test_yaml.f90 -o test_yaml.o
ifort -I/work/Projects/ASiMoV/fortran-yaml-cpp/build -c example.f90 -o example.o

ifort test_yaml.o -o test_yaml_intel -Wl,-rpath,/work/Projects/ASiMoV/fortran-yaml-cpp/build:/work/Projects/ASiMoV/fortran-yaml-cpp/build/yaml-cpp -L/work/Projects/ASiMoV/fortran-yaml-cpp/build -lfortran-yaml-cpp -L/work/Projects/ASiMoV/fortran-yaml-cpp/build/yaml-cpp -lyaml-cpp
ifort example.o -o example_intel -Wl,-rpath,/work/Projects/ASiMoV/fortran-yaml-cpp/build:/work/Projects/ASiMoV/fortran-yaml-cpp/build/yaml-cpp -L/work/Projects/ASiMoV/fortran-yaml-cpp/build -lfortran-yaml-cpp -L/work/Projects/ASiMoV/fortran-yaml-cpp/build/yaml-cpp -lyaml-cpp
