CC = mpicc
CFLAGS = -O3

FC = mpif90
FFLAGS = -cpp -std=f2018 -ffree-line-length-none
CAFFLAGS = -fcoarray=single
ifeq ($(VERBOSE),yes)
  FFLAGS += -DVERBOSE
endif
ifeq ($(BUILD),debug)
  # Add debugging (i.e. expensive) flags
  FFLAGS += -g -Og
  FFLAGS += -fcheck=bounds
  FFLAGS += -fbacktrace
  FFLAGS += -ffpe-trap=invalid,zero,overflow
  FFLAGS += -Wimplicit-interface -Wimplicit-procedure
else
  FFLAGS += -O3
endif
ifeq ($(PROFILE),yes)
  FFLAGS += -fopt-info-missed-optall=opt_info.txt
endif
FFLAGS += -fopenmp
FFLAGS += -Wall -Wpedantic -Werror 
FFLAGS += -J$(OBJ_DIR)
MPIRUN = mpirun

# Only set this value if building a CAF binary otherwise keep unset
#CAFLINK= #-fcoarray=single
