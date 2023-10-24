CC = mpicc
CFLAGS = -O3

FC = mpif90
FFLAGS = -cpp -std=f2018 -ffree-line-length-none -fimplicit-none -g
CAFFLAGS = -fcoarray=single
ifeq ($(VERBOSE),yes)
  FFLAGS += -DVERBOSE
endif
ifeq ($(BUILD),debug)
  # Add debugging (i.e. expensive) flags
  FFLAGS += -Og
  FFLAGS += -fcheck=bounds
  FFLAGS += -fbacktrace
  FFLAGS += -ffpe-trap=invalid,zero,overflow
  #FFLAGS += -Wimplicit-interface -Wimplicit-procedure
  FFLAGS += -Wall -Wpedantic -Wno-uninitialized -Werror # Wuninitialized is buggy
  FFLAGS += -DEXCLUDE_MISSING_INTERFACE
else
  FFLAGS += -O3
  # Warning, native flags don't allow cross compilation
  FFLAGS += -march=native -mtune=native
endif
ifeq ($(PROFILE),yes)
  FFLAGS += -fopt-info-missed-optall=opt_info.txt
endif

# Check for gfortran version >= 10
GFORTRAN_VER_GTE10 := $(shell echo `gfortran -dumpversion | cut -f1 -d.` \>= 10 | bc)
ifeq ($(GFORTRAN_VER_GTE10),1)
  FFLAGS += -fallow-argument-mismatch
endif

FFLAGS += -fopenmp
FFLAGS += -J$(OBJ_DIR)
MPIRUN = mpirun --oversubscribe

# Only set this value if building a CAF binary otherwise keep unset
#CAFLINK= #-fcoarray=single
