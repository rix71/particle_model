#ifdef DEBUG
#define OUT_ID
#define OUT_VELOCITY
#define OUT_SETTLING_VELOCITY
#define OUT_DENSITY
#define OUT_DENSITY_PLASTIC
#define OUT_DELTA_RHO
#define OUT_RADIUS
#define OUT_RADIUS_PLASTIC
#define OUT_AGE
#define OUT_TRAJECTORY
#define OUT_TIME_ON_BEACH
#define OUT_H_BIOFILM
#define OUT_GROWTH_BIOFILM
#define OUT_STATE
#define OUT_KIN_VISCOSITY
#define OUT_FRICTION_VELOCITY
#else
#define OUT_ID
#define OUT_VELOCITY
#define OUT_SETTLING_VELOCITY
#define OUT_STATE
#endif


!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!(parallel)gathering directives to shorten code
#define GATHER_LOOP_RV(varname) do ipart = 1, nwrite; var(ipart,1) = particles(ipart)%varname; end do
#define GATHER_LOOP_IV(varname) do ipart = 1, nwrite; var_int(ipart,1) = particles(ipart)%varname; end do

#ifdef USE_OMP
#define GATHER_OMP_START_RV START_OMP_DO private(ipart) shared(var)
#define GATHER_OMP_START_IV START_OMP_DO private(ipart) shared(var_int)
#define GATHER_OMP_END !$omp end parallel do
#else
#define GATHER_OMP_START_RV
#define GATHER_OMP_START_IV
#define GATHER_OMP_END
#endif

