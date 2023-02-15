#ifdef USE_OMP
#define START_OMP_DO !$omp parallel do
#define END_OMP_DO !$omp end parallel do
#else
#define START_OMP_DO !just a comment in case this is followed by more parameters
#define END_OMP_DO 
#endif