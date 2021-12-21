#include "cppdefs.h"
module postprocessing
  use precdefs
  use particle_vars, only: particles, runparts
  use domain_vars, only: nx, ny, x0, y0, dx, dy, seamask, lons, lats
  use fields, only: get_indices2d
  use output, only: outDir
  use run_params, only: runid
  use nc_manager, only: nc_write, nc_initialise, &
                        nc_add_dimension, nc_add_variable, &
                        nc_add_attr, FILLVALUE_BIG
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: postprocess
  !---------------------------------------------
  character(len=512) :: nc_fileout_post
  integer            :: nc_x_dimid, nc_y_dimid
  !===================================================
contains
  !===========================================
  subroutine postprocess

    integer  :: ipart, i, j
    integer  :: counts(nx, ny)
    real(rk) :: mean_age(nx, ny), mean_distance(nx, ny)

    dbghead(postprocess)

    counts = 0
    mean_age = 0
    mean_distance = 0
    do ipart = 1, runparts
      call get_indices2d(particles(ipart)%xPos, particles(ipart)%yPos, x0, y0, dx, dy, i, j)
      counts(i, j) = counts(i, j) + 1
      mean_age(i, j) = mean_age(i, j) + particles(ipart)%age
      mean_distance(i, j) = mean_distance(i, j) + particles(ipart)%trajLen
    end do

    where (seamask == 1) counts = int(FILLVALUE_BIG)
    where (counts > 0)
      mean_age = mean_age / counts
      mean_distance = mean_distance / counts
    end where
    where (seamask == 1)
      mean_age = FILLVALUE_BIG
      mean_distance = FILLVALUE_BIG
    end where

    nc_fileout_post = trim(outDir)//"/"//trim(runid)//".post.nc"
    call nc_initialise(trim(nc_fileout_post))
    call nc_add_dimension(trim(nc_fileout_post), "lon", nc_x_dimid, nx)
    call nc_add_dimension(trim(nc_fileout_post), "lat", nc_y_dimid, ny)

    call nc_add_variable(trim(nc_fileout_post), "lon", "float", 1, [nc_x_dimid])
    call nc_write(trim(nc_fileout_post), lons, "lon", nx)

    call nc_add_variable(trim(nc_fileout_post), "lat", "float", 1, [nc_y_dimid])
    call nc_write(trim(nc_fileout_post), lats, "lat", ny)

    call nc_add_variable(trim(nc_fileout_post), "counts", "int", 2, [nc_x_dimid, nc_y_dimid], FILLVALUE_BIG)
    call nc_write(trim(nc_fileout_post), counts, "counts", nx, ny)

    call nc_add_variable(trim(nc_fileout_post), "mean_age", "float", 2, [nc_x_dimid, nc_y_dimid], FILLVALUE_BIG)
    call nc_write(trim(nc_fileout_post), mean_age, "mean_age", nx, ny)

    call nc_add_variable(trim(nc_fileout_post), "mean_distance", "float", 2, [nc_x_dimid, nc_y_dimid], FILLVALUE_BIG)
    call nc_write(trim(nc_fileout_post), mean_distance, "mean_distance", nx, ny)

    dbgtail(postprocess)
  end subroutine postprocess
end module postprocessing
