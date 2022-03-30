FC = gfortran
FFLAGS = -O3 -cpp -ffree-line-length-none -Wall \
	 -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface \
	 -Wunused-parameter -fwhole-file -fcheck=all -fbacktrace -mcmodel=large

# Print everything to .stdout file
FFLAGS += -DWRITESTDOUT

# Print debug messages
# FFLAGS += -DDEBUG

# Silence all messages (except warnings and errors)
# FFLAGS += -DSAY_LESS

# Disable system calls (call system(...))
# FFLAGS += -DNOSYSCALLS

# Do not allow too big indices (needed when the time step is big or domain has northern/eastern boundaries)
FFLAGS += -DSNAP_TO_BOUNDS

# Land interaction methods
# FFLAGS += -DPARTICLE_BOUNCE
FFLAGS += -DPARTICLE_REDIRECT
# FFLAGS += -DPARTICLE_BEACH_IMMEDIATELY

# Use OpenMP
FFLAGS += -DUSE_OMP -fopenmp

#--------------------------------------------------
MOD=./mod
OBJ=./obj
BIN=../bin

NETCDF_DIR=`nf-config --prefix`
NETCDF_INC=-I${NETCDF_DIR}/include
NETCDF_LIB=-L${NETCDF_DIR}/lib

FFLAGS += $(NETCDF_INC)

LIBS = ${NETCDF_LIB} -lnetcdff -lnetcdf

SRCS = precdefs.f90 errors.f90 \
       interp.f90 nc_manager.f90 field.f90 list.f90 \
	   domain.f90 \
	   datetime.f90 \
	   fieldset.f90 \
	   vars.f90 \
       particle.f90 output.f90 \
	   advection.f90 physics.f90 \
       init.f90 loop.f90 postprocessing.f90 \
	   main.f90 

OBJS = $(patsubst %.f90, $(OBJ)/%.o, $(SRCS))

#--------------------------------------------------
PROGRAM = $(BIN)/main

all: $(PROGRAM)

$(PROGRAM): $(OBJS) | $(BIN)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

$(OBJ)/%.o: %.f90 | $(OBJ)
	$(FC) $(FFLAGS) -c $< -o $@ $(LIBS)

$(BIN) $(OBJ):
	mkdir $@

clean:
	rm -rf ./obj
	rm *.mod

show:
	@ echo $(OBJS)
