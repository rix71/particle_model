!--------------
! zax_style
! Regular depth axis:
#define STATIC_DEPTH_VALUES 0
! Adaptive coordinates:
#define DEPTH_VALUES 1
#define LAYER_THICKNESS 2
!--------------
! density_method
#define RHO_DEFAULT 0
#define RHO_VARIABLE 1
#define RHO_CALC 2
!--------------
! viscosity_method
#define VISC_DEFAULT 0
#define VISC_VARIABLE 1
#define VISC_CALC 2
!--------------
! seamask
#define DOM_LAND 1
#define DOM_SEA 2
#define DOM_BEACH 3
#define DOM_BOUNDARY 4