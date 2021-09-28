rm -R objs
rm -R bin

mkdir objs

cd ccs_core
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals -stand f18 -DACCS_PETSC -o kinds_mod.o -c kinds_mod.f90 -I/work/Projects/ASiMoV/petsc/include -I/work/Projects/ASiMoV/petsc/arch-linux-c-opt/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o constants_mod.o -c constants_mod.f90 -I/work/Projects/ASiMoV/petsc/include -I/work/Projects/ASiMoV/petsc/arch-linux-c-opt/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o types_mod.o -c types_mod.f90 -I/work/Projects/ASiMoV/petsc/include -I/work/Projects/ASiMoV/petsc/arch-linux-c-opt/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o mesh_utils_mod.o -c mesh_utils_mod.f90 -I/work/Projects/ASiMoV/petsc/include -I/work/Projects/ASiMoV/petsc/arch-linux-c-opt/include
mv *.o *.*mod ../objs

cd petsc
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o petsctypes_mod.o -c petsctypes_mod.f90 -I../../objs -I/work/Projects/ASiMoV/petsc/include -I/work/Projects/ASiMoV/petsc/arch-linux-c-opt/include
mv *.o *.*mod ../../objs

cd ../../linear_solvers
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o solver_mod.o -c solver_mod.f90 -I../objs -I/work/Projects/ASiMoV/petsc/include -I/work/Projects/ASiMoV/petsc/arch-linux-c-opt/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o vec_mod.o -c vec_mod.f90 -I../objs -I/work/Projects/ASiMoV/petsc/include -I/work/Projects/ASiMoV/petsc/arch-linux-c-opt/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o mat_mod.o -c mat_mod.f90 -I../objs -I/work/Projects/ASiMoV/petsc/include -I/work/Projects/ASiMoV/petsc/arch-linux-c-opt/include
mv *.o *.*mod ../objs

cd petsc
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o vec_petsc.o -c vec_petsc.f90 -I../../objs -I/work/Projects/ASiMoV/petsc/include -I/work/Projects/ASiMoV/petsc/arch-linux-c-opt/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o mat_petsc.o -c mat_petsc.f90 -I../../objs -I/work/Projects/ASiMoV/petsc/include -I/work/Projects/ASiMoV/petsc/arch-linux-c-opt/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o solver_petsc.o -c solver_petsc.f90 -I../../objs -I/work/Projects/ASiMoV/petsc/include -I/work/Projects/ASiMoV/petsc/arch-linux-c-opt/include
mv *.o *.*mod ../../objs

cd ../../ccs_core
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o utils_mod.o -c utils_mod.f90 -I../objs  -I/work/Projects/ASiMoV/petsc/include -I/work/Projects/ASiMoV/petsc/arch-linux-c-opt/include
mv *.o *.*mod ../objs

cd ../case_setup/Poisson
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o poisson.o -c poisson.f90 -I../../objs -I/work/Projects/ASiMoV/petsc/include -I/work/Projects/ASiMoV/petsc/arch-linux-c-opt/include
mv *.o *.*mod ../../objs

cd ../../
mkdir bin
cd bin
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o poisson ../objs/mat_mod.o ../objs/mat_petsc.o ../objs/solver_mod.o ../objs/solver_petsc.o ../objs/petsctypes_mod.o \
../objs/vec_mod.o ../objs/vec_petsc.o ../objs/poisson.o ../objs/types_mod.o ../objs/utils_mod.o ../objs/kinds_mod.o ../objs/mesh_utils_mod.o ../objs/constants_mod.o \
-I../objs -I/work/Projects/ASiMoV/petsc/include -I/work/Projects/ASiMoV/petsc/arch-linux-c-opt/include -L/work/Projects/ASiMoV/petsc/arch-linux-c-opt/lib -lpetsc
