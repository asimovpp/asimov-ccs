
# How to run asimov-ccs
    make
    ./ccs_app

# Customising asimov-ccs
  Edit config.json
  
  Valid options:
  - main
    - pi_mpi
    - pi_caf
  -  parallel_env
    - mpi
    - caf
  -  parallel_utils
    - mpi
    - caf
  -  compute
    - mpi
    - mpi_omp
    - caf
  -  collectives
    - mpi
    - none
  -  parallel_types
    - mpi
    - none

  You cannot mix mpi and caf, but at the moment there will only the linker will throw an error if mixing does happen. But it does not catch all incorrect combinations!
  The only actually valid combinations are:
  - pi_mpi, mpi, mpi, mpi, mpi, mpi
  - pi_mpi, mpi, mpi, mpi_omp, mpi, mpi
  - pi_caf, caf, caf, caf, none, none
