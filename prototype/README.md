# How to run asimov-ccs
    export FC=${my_fortran_compiler}
    make all
    python3 link_executable.py config.json
    ./ccs_main

# Customising asimov-ccs
  Edit config.json
  
  Valid options:
  - turbulance
    - ke
    - kw 
  - particles
    - 1
    - 2
  - flux
    - 1st_order
    - 2nd_order
  - solver
    - amg
    - cgstab

