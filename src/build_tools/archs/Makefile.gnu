FC = mpif90
FFLAGS = -cpp -std=f2018
CAFFLAGS = -fcoarray=single
ifeq ($(BUILD),debug)
  # Add debugging (i.e. expensive) flags
  FFLAGS += -g -Og
  FFLAGS += -fcheck=bounds
  FFLAGS += -fbacktrace
else
  FFLAGS += -O3
endif
ifeq ($(PROFILE),yes)
  FFLAGS += -fopt-info-missed-optall=opt_info.txt
endif
FFLAGS += -fopenmp
FFLAGS += -Wall -Wpedantic -Werror -Wimplicit-interface -Wimplicit-procedure
FFLAGS += -J$(OBJ_DIR)
MPIRUN = mpirun

# Only set this value if building a CAF binary otherwise keep unset
#CAFLINK= #-fcoarray=single
