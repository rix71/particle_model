FC = gfortran
FFLAGS = -g -cpp -ffree-line-length-none -Wall -Wline-truncation \
	 -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface \
	 -Wunused-parameter -fwhole-file -fcheck=all -fbacktrace -mcmodel=large

FFLAGS += -DDEBUG

MOD=./mod
OBJ=./obj
BIN=../bin

NETCDF_DIR=`nf-config --prefix`
NETCDF_INC=-I${NETCDF_DIR}/include
NETCDF_LIB=-L${NETCDF_DIR}/lib

FFLAGS += $(NETCDF_INC)

LIBS = ${NETCDF_LIB} -lnetcdff -lnetcdf

SRCS = precdefs.f90 errors.f90 interp.f90 nc_manager.f90 modtime.f90 vars.f90 \
       physics.f90 domain.f90 fields.f90 particle.f90 output.f90 \
       init.f90 loop_particle.f90 main.f90 

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
