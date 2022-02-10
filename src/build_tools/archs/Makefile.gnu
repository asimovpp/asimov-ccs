FC = mpif90
FFLAGS = -cpp -O3 -std=f2018
CAFFLAGS = -fcoarray=single
ifeq ($(BUILD),debug)
  # Add debugging (i.e. expensive) flags
  FFLAGS += -fcheck=bounds
endif
FFLAGS += -fopenmp
FFLAGS += -Wall -Wpedantic -Werror -Wimplicit-interface -Wimplicit-procedure
FFLAGS += -J$(OBJ_DIR)
MPIRUN = mpirun

# Only set this value if building a CAF binary otherwise keep unset
#CAFLINK= #-fcoarray=single
