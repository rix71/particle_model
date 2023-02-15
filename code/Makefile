FC = gfortran
FFLAGS = -O3 -g -cpp -ffree-line-length-none -Wall \
	 -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface \
	 -Wunused-parameter -fwhole-file -fcheck=all -fbacktrace -mcmodel=large

#--------------------------------------------------
# COMPILER FLAGS
#--------------------------------------------------
# Profiling flags
# FFLAGS += -pg -no-pie

# Print everything to .stdout file
# FFLAGS += -DWRITESTDOUT

# Print debug messages
# FFLAGS += -DDEBUG

# Silence all messages (except warnings and errors)
# FFLAGS += -DSAY_LESS

# Disable system calls (call system(...))
# FFLAGS += -DNOSYSCALLS

# Do not allow too big indices (needed when the time step is big or domain has northern/eastern boundaries)
FFLAGS += -DSNAP_TO_BOUNDS

# Calculate the vertical displacement in 3D mode (enable to keep particles at initial depths) 
# FFLAGS += -DNO_ADVECT_VERTICAL

# Calculate the vertical displacement in 3D mode (enable to keep particles at initial depths)
# FFLAGS += -DNO_DIFFUSE_VERTICAL

# Smagorinsky parametrization method (full field much slower)
# FFLAGS += -DSMAGORINSKY_FULL_FIELD
FFLAGS += -DSMAGORINSKY_INTERP_UV

# Set particle to sealevel, if floating (disable to keep particles at initial depths (surface))
FFLAGS += -DPARTICLE_SNAP_SEALVL

# Land interaction methods
# FFLAGS += -DPARTICLE_BOUNCE
FFLAGS += -DPARTICLE_REDIRECT
# FFLAGS += -DPARTICLE_BEACH_IMMEDIATELY

# Ignore particles that are outside the domain
FFLAGS += -DIGNORE_BAD_PARTICLES

# Compile and run the postprocessing module
# FFLAGS += -DPOSTPROCESS

# Use OpenMP
FFLAGS += -DUSE_OMP -fopenmp

# Biofouling model (Kooi, Tsiaras or Simple)
BF_MODEL = Simple

#--------------------------------------------------
ifeq ($(BF_MODEL), Kooi)
	FFLAGS += -DBIOFOULING_KOOI
	BF_PATH = biofouling/kooi
	# BF_FILES = $(wildcard $(BF_PATH)/*.f90)
	BF_FILES = $(BF_PATH)/biofouling_kooi.f90 $(BF_PATH)/biofouling_kooi_funcs.f90
else ifeq ($(BF_MODEL), Tsiaras)
	FFLAGS += -DBIOFOULING_TSIARAS
	BF_PATH = biofouling/tsiaras
	BF_FILES = $(wildcard $(BF_PATH)/*.f90)
else ifeq ($(BF_MODEL), Simple)
	FFLAGS += -DBIOFOULING_SIMPLE
	BF_PATH = biofouling/simple
	# BF_FILES = $(wildcard $(BF_PATH)/*.f90)
	BF_FILES = $(BF_PATH)/biofouling_simple.f90 $(BF_PATH)/biofouling_simple_funcs.f90
else
	$(error "Unknown biofouling model")
endif

#--------------------------------------------------
# OBJ=./obj
BIN=../bin
BUILD=../build
MOD=$(BUILD)/mod

ifneq ($(wildcard $(MOD)),)
	FFLAGS += -I${MOD}
endif

FFLAGS += -I../include

NETCDF_DIR=`nf-config --prefix`
NETCDF_INC=-I${NETCDF_DIR}/include
NETCDF_LIB=-L${NETCDF_DIR}/lib

FFLAGS += $(NETCDF_INC)

LIBS = ${NETCDF_LIB} -lnetcdff -lnetcdf

SRCS = precdefs.f90 errors.f90 \
       utils.f90 interp.f90 nc_manager.f90 field.f90 list.f90 \
	   domain.f90 datetime.f90 \
	   fieldset.f90 \
	   vars.f90 \
       particle.f90 output.f90 \
	   advection.f90 physics.f90 \
	   vertical_motion.f90 \
       diffusion.f90 $(BF_FILES) \
	   init.f90 loop.f90 postprocessing.f90 \
	   main.f90 

OBJS = $(patsubst %.f90, $(BUILD)/%.o, $(SRCS))

#--------------------------------------------------
define fc-command
$(FC) $(FFLAGS) -c $< -o $@ $(LIBS)
endef
#--------------------------------------------------
PROGRAM = $(BIN)/main

all: $(PROGRAM)

$(PROGRAM): $(OBJS) | $(BIN) $(MOD)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)
	make remove_mods

$(BUILD)/%.o: %.f90 | $(BUILD)
	$(fc-command)

$(BUILD)/$(BF_PATH)/%.o: $(BF_PATH)/%.f90 | $(BUILD)/$(BF_PATH)
	$(fc-command)

$(BIN) $(BUILD) $(BUILD)/$(BF_PATH) $(MOD):
	mkdir -p $@

clean:
	rm -rf $(BUILD) $(MOD) 

remove_mods:
ifneq ($(wildcard ./*.mod),)	
	mv *.mod $(MOD)/.
endif
ifneq ($(wildcard ./*.smod),)	
	mv *.smod $(MOD)/.
endif

show:
	@ echo $(OBJS)
