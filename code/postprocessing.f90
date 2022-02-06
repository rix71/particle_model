! #define DEBUG
#include "cppdefs.h"
module postprocessing
  use precdefs
  use errors
  use particle_vars, only: particles, runparts
  use domain_vars, only: nx, ny, x0, y0, dx, dy, seamask, lons, lats
  use fields, only: get_indices2d
  use output, only: outDir, write_active_particles
  use run_params, only: runid
  use nc_manager
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: postprocess
  !---------------------------------------------
  character(len=512) :: nc_fileout_post
  integer            :: nc_x_dimid, nc_y_dimid, nc_t_dimid
  !===================================================
contains
  !===========================================
  subroutine process_all_particles(counts, mean_age, mean_distance)
    !---------------------------------------------
    ! Process the final state of all particles
    !---------------------------------------------
    integer               :: ipart, i, j
    integer, intent(out)  :: counts(nx, ny)
    real(rk), intent(out) :: mean_age(nx, ny), mean_distance(nx, ny)

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

  end subroutine process_all_particles
  !===========================================
  subroutine process_active_file(counts, ntimes, timevals, timeunit, mean_age_time)

    character(len=512)                   :: nc_filein_active
    character(len=512), intent(out)      :: timeunit
    integer                              :: nparticles, itime, ipart, i, j
    real(rk), allocatable                :: x(:, :), y(:, :), age(:, :)
    real(rk), allocatable, intent(inout) :: timevals(:) !, mean_age_time(:, :, :)
    integer, intent(out)                 :: counts(nx, ny), ntimes
    real(rk), intent(out)                :: mean_age_time(nx, ny)
    ! integer, allocatable, intent(inout)  :: counts(:, :, :)
    integer                              :: ierr

    dbghead(process_active_file)

    nc_filein_active = trim(outDir)//'/'//trim(runid)//'.active.nc'

    call nc_get_dim(trim(nc_filein_active), "time", ntimes)
    call nc_get_dim(trim(nc_filein_active), "particle", nparticles)
    debug(ntimes); debug(nparticles)

    allocate (timevals(ntimes), x(ntimes, nparticles), y(ntimes, nparticles), age(ntimes, nparticles), stat=ierr)
    ! counts(nx, ny), &
    ! mean_age_time(nx, ny), &
    
    if (ierr .ne. 0) call throw_error("process_active_file", "Could not allocate", ierr)
    call nc_read_real_1d(trim(nc_filein_active), "time", ntimes, timevals)
    call nc_get_timeunit(trim(nc_filein_active), timeunit)
    call nc_read_real_2d(trim(nc_filein_active), "x", nparticles, ntimes, x)
    call nc_read_real_2d(trim(nc_filein_active), "y", nparticles, ntimes, y)
    call nc_read_real_2d(trim(nc_filein_active), "age", nparticles, ntimes, age)

    counts = 0
    mean_age_time = 0.
    do itime = 1, ntimes
      do ipart = 1, nparticles
        if ((x(itime, ipart) > lons(1)) .and. (x(itime, ipart) < lons(nx)) &
            .and. (y(itime, ipart) > lats(1)) .and. (y(itime, ipart) < lats(ny))) then
          debug(x(itime, ipart)); debug(y(itime, ipart))
          call get_indices2d(x(itime, ipart), y(itime, ipart), x0, y0, dx, dy, i, j)
          debug(i); debug(j)

          counts(i, j) = counts(i, j) + 1
          mean_age_time(i, j) = mean_age_time(i, j) + age(itime, ipart)

        end if
      end do
    end do
    where (counts > 0) mean_age_time = mean_age_time / counts
    ! counts = sum(counts_time, dim=3)

    where (seamask == LAND)
      counts = int(FILLVALUE_BIG)
      ! counts_time(i, j) = FILLVALUE_BIG
      mean_age_time = FILLVALUE_BIG
    end where
    ! do i = 1, nx
    !   do j = 1, ny
    !     if (seamask(i, j) ==) then
    !       counts_time(i, j, :) = FILLVALUE_BIG
    !       mean_age_time(i, j, :) = FILLVALUE_BIG
    !     end if
    !   end do
    ! end do

    dbgtail(process_active_file)
  end subroutine process_active_file
  !===========================================
  subroutine postprocess

    dbghead(postprocess)

    nc_fileout_post = trim(outDir)//"/"//trim(runid)//".post.nc"
    call nc_initialise(trim(nc_fileout_post))
    call nc_add_dimension(trim(nc_fileout_post), "lon", nc_x_dimid, nx)
    call nc_add_dimension(trim(nc_fileout_post), "lat", nc_y_dimid, ny)

    final_state: block
      integer  :: counts(nx, ny)
      real(rk) :: mean_age(nx, ny), mean_distance(nx, ny)

      call process_all_particles(counts, mean_age, mean_distance)

      call nc_add_variable(trim(nc_fileout_post), "lon", "float", 1, [nc_x_dimid])
      call nc_add_attr(trim(nc_fileout_post), "lon", "units", "degrees east")
      call nc_write(trim(nc_fileout_post), lons, "lon", nx)

      call nc_add_variable(trim(nc_fileout_post), "lat", "float", 1, [nc_y_dimid])
      call nc_add_attr(trim(nc_fileout_post), "lat", "units", "degrees north")
      call nc_write(trim(nc_fileout_post), lats, "lat", ny)

      call nc_add_variable(trim(nc_fileout_post), "counts", "int", 2, [nc_x_dimid, nc_y_dimid], FILLVALUE_BIG)
      call nc_add_attr(trim(nc_fileout_post), "counts", "units", "particles")
      call nc_write(trim(nc_fileout_post), counts, "counts", nx, ny)

      call nc_add_variable(trim(nc_fileout_post), "mean_age", "float", 2, [nc_x_dimid, nc_y_dimid], FILLVALUE_BIG)
      call nc_add_attr(trim(nc_fileout_post), "mean_age", "units", "s")
      call nc_write(trim(nc_fileout_post), mean_age, "mean_age", nx, ny)

      call nc_add_variable(trim(nc_fileout_post), "mean_distance", "float", 2, [nc_x_dimid, nc_y_dimid], FILLVALUE_BIG)
      call nc_add_attr(trim(nc_fileout_post), "mean_distance", "units", "m")
      call nc_write(trim(nc_fileout_post), mean_distance, "mean_distance", nx, ny)
    end block final_state

    if (write_active_particles) then
      active: block
        character(len=512)    :: timeunit
        integer               :: counts(nx, ny), ntimes
        real(rk)              :: mean_age_time(nx, ny)
        ! integer, allocatable  :: counts_time(:, :, :)
        real(rk), allocatable :: timevals(:) !, mean_age_time(:, :, :)

        call process_active_file(counts, ntimes, timevals, timeunit, mean_age_time)

        call nc_add_dimension(trim(nc_fileout_post), "time", nc_t_dimid, ntimes)

        call nc_add_variable(trim(nc_fileout_post), "time", "float", 1, [nc_t_dimid])
        call nc_add_attr(trim(nc_fileout_post), "time", "units", timeunit)
        call nc_write(trim(nc_fileout_post), timevals, "time", ntimes)

        call nc_add_variable(trim(nc_fileout_post), "counts_active", "int", 2, [nc_x_dimid, nc_y_dimid], FILLVALUE_BIG)
        call nc_add_attr(trim(nc_fileout_post), "counts_active", "units", "particles")
        call nc_write(trim(nc_fileout_post), counts, "counts_active", nx, ny)

        ! call nc_add_variable(trim(nc_fileout_post), "counts_time_active", "int", 3, &
        !                      [nc_x_dimid, nc_y_dimid, nc_t_dimid], FILLVALUE_BIG)
        ! call nc_add_attr(trim(nc_fileout_post), "counts_time_active", "units", "particles")
        ! call nc_write(trim(nc_fileout_post), counts_time, "counts_time_active", nx, ny, ntimes)

        call nc_add_variable(trim(nc_fileout_post), "mean_age_active", "float", 2, &
                             [nc_x_dimid, nc_y_dimid], FILLVALUE_BIG)
        call nc_add_attr(trim(nc_fileout_post), "mean_age_active", "units", "s")
        call nc_write(trim(nc_fileout_post), mean_age_time, "mean_age_active", nx, ny)

      end block active
    end if

    dbgtail(postprocess)
  end subroutine postprocess
end module postprocessing
