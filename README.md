# ASiMoV-CCS

## Prerequisites

0 `MPI`
- `PETSc` - we are tracking the latest versions, however a bug in auto-generated Fortran interfaces means it is easiest to get started with the 3.13 series as found in the ARCHER2 modules.
- `makedepf90` - A version can be found on ARCHER2, otherwise get the patched version from Debian (if using Ubuntu/Debian this can be `apt-get install`'d)

## Building

With the prerequisites in place, ASiMoV-CCS can be built with
``
make CMP=<compiler>
``
where `<compiler>` is one of: `gnu`, `intel` or `cray` to build with GNU compilers (`mpif90`), Intel (`mpiifort`) or Cray (`ftn`) compilers, respectively. To override the compiler set `FC` on the command line, for instance to build with GNU compilers on a Cray machine
``
make CMP=gnu FC=ftn
``
