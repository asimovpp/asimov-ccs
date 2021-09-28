PETSC_ARCH=arch-linux-c-opt-mpif08
PETSC_DIR=/work/Projects/ASiMoV/petsc

rm -R bin

mkdir objs

echo "Build ccs_core"
cd ccs_core
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals -stand f18 -DACCS_PETSC -o kinds_mod.o -c kinds_mod.f90 -I../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o constants_mod.o -c constants_mod.f90 -I../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o types_mod.o -c types_mod.f90 -I../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o mesh_utils_mod.o -c mesh_utils_mod.f90 -I../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mv *.o *.*mod ../objs

echo "Build ccs_core/petsc"
cd petsc
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o petsctypes_mod.o -c petsctypes_mod.f90 -I../../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mv *.o *.*mod ../../objs

echo "Build linear_solvers"
cd ../../linear_solvers
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o solver_mod.o -c solver_mod.f90 -I../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o vec_mod.o -c vec_mod.f90 -I../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o mat_mod.o -c mat_mod.f90 -I../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mv *.o *.*mod ../objs

echo "Build linear_solvers/petsc"
cd petsc
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o vec_petsc.o -c vec_petsc.f90 -I../../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o mat_petsc.o -c mat_petsc.f90 -I../../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o solver_petsc.o -c solver_petsc.f90 -I../../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mv *.o *.*mod ../../objs

echo "Build ccs_core utils"
cd ../../ccs_core
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o utils_mod.o -c utils_mod.f90 -I../objs  -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mv *.o *.*mod ../objs

echo "Build parallel"
cd ../parallel
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o parallel_types_mod.o -c parallel_types_mod.f90 -I../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o parallel_mod.o -c parallel_mod.f90 -I../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o parallel_types_mpi_mod.o -c parallel_types_mpi_mod.f90 -I../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o parallel_env_mpi.o -c parallel_env_mpi.f90 -I../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o parallel_errors_mpi.o -c parallel_errors_mpi.f90 -I../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o parallel_utils_mpi.o -c parallel_utils_mpi.f90 -I../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mv *.o *.*mod ../objs

echo "Build parallel/petsc"
cd petsc
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o parallel_env_mpi_petsc.o -c parallel_env_mpi_petsc.f90 -I../../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o parallel_errors_mpi_petsc.o -c parallel_errors_mpi_petsc.f90 -I../../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mv *.o *.*mod ../../objs

echo "Build Poisson"
cd ../../case_setup/Poisson
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o poisson.o -c poisson.f90 -I../../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
mv *.o ../../objs

cd ../../
mkdir bin
cd bin
echo "Link Poisson"
mpiifort -cpp -fPIC -g -O3 -xHOST -warn all,noexternals  -stand f18 -DACCS_PETSC -o poisson ../objs/mat_mod.o ../objs/mat_petsc.o ../objs/solver_mod.o ../objs/solver_petsc.o ../objs/petsctypes_mod.o \
../objs/vec_mod.o ../objs/vec_petsc.o ../objs/poisson.o ../objs/types_mod.o ../objs/utils_mod.o ../objs/kinds_mod.o ../objs/mesh_utils_mod.o ../objs/constants_mod.o \
../objs/parallel_types_mod.o ../objs/parallel_mod.o ../objs/parallel_types_mpi_mod.o ../objs/parallel_env_mpi_petsc.o ../objs/parallel_errors_mpi_petsc.o ../objs/parallel_utils_mpi.o \
-I../objs -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -L/work/Projects/ASiMoV/petsc/arch-linux-c-opt/lib -lpetsc
