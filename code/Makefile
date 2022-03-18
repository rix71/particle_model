FC = gfortran
FFLAGS = -O3 -cpp -ffree-line-length-none -Wall \
	 -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface \
	 -Wunused-parameter -fwhole-file -fcheck=all -fbacktrace -mcmodel=large

FFLAGS += -DWRITESTDOUT
# FFLAGS += -DDEBUG
# FFLAGS += -DNOSYSCALLS
FFLAGS += -DSNAP_TO_BOUNDS


MOD=./mod
OBJ=./obj
BIN=../bin

NETCDF_DIR=`nf-config --prefix`
NETCDF_INC=-I${NETCDF_DIR}/include
NETCDF_LIB=-L${NETCDF_DIR}/lib

FFLAGS += $(NETCDF_INC)
# FFLAGS += -Iinclude

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

# SRCS = Utils/precdefs.f90 \
#        Utils/errors.f90 \
#        Utils/interp.f90 \
# 	   Utils/nc_manager.f90 \
# 	   Utils/datetime.f90

# SRCS += Unknown/vars.f90 Unknown/physics.f90

# SRCS += Field/domain.f90 \
#         Field/fields.f90 

# SRCS += Particle/particle.f90

# SRCS += Unknown/output.f90 \
#         Core/init.f90 \
#         Core/loop_particle.f90 \
#         Unknown/mod_postprocessing.f90 \
#         Core/main.f90
        


OBJS = $(patsubst %.f90, $(OBJ)/%.o, $(SRCS))
# OBJS = $(patsubst %.f90, %.o, $(SRCS))
# OBJ_DIRS = $(OBJ) $(OBJ)/Utils $(OBJ)/Field $(OBJ)/Particle $(OBJ)/Core $(OBJ)/Unknown
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
	@ echo $(OBJ_DIRS)
