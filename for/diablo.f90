!******************************************************************************|
! diablo.f90                                                           VERSION 3
! v2.0 JRT, 2012
! v3.0 AFW, 2021
!
! This Fortran code computes incompressible flow in a channel.
!
! Primative variables (u,v,w,p) are used, and continuity is enforced with a
! fractional step algorithm.
!
! SPATIAL DERIVATIVES:
!   The 1 & 3 (X/Z) directions are periodic and handled spectrally
!   The 2 (Y) direction is taken to be bounded by walls and handled with
!   momentum- and energy-conserving second-order central finite differences.
!
! TIME ADVANCEMENT
!   Two main approaches are implemented:
!     1. RKW3 on nonlinear terms and CN on viscous terms over each RK substep.
!     2. RKW3 on y-derivative terms and CN on other terms over each RK substep.
!
! A few simple high-performance programming constructs are used:
!   -> The inner 2 loops are broken out in such a way as to enable out-of-order
!      execution of these loops as much as possible, thereby leveraging
!      vector and superscalar CPU architectures.
!   -> The outer loops are fairly long (including as many operations as
!      possible inside on a single J plane of data) in order to make effective
!      use of cache.
!
!******************************************************************************|
!
! This code is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation; either version 2 of the License, or (at your
! option) any later version. This code is distributed in the hope that it
! will be useful, but WITHOUT any WARRANTY; without even the implied
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details. You should have received a
! copy of the GNU General Public License along with this code; if not,
! write to the Free Software Foundation, Inc., 59 Temple Place - Suite
! 330, Boston, MA 02111-1307, USA.
!
!******************************************************************************|

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
program diablo
  use parameters
  use domain
  use flow
  use ics
  use advance
  use control
  use statistics
  use les
  implicit none 

  integer :: n
  logical :: flag

  call read_inputs
  call set_parameters
  call init_mpi 

  if (rank == 0) then
    write (*, *)
    write (*, *) '             ****** WELCOME TO DIABLO ******'
    write (*, *)
    call log_input_parameters
  end if

  ! Initialize channel geometry
  call create_grid_chan
  call init_chan_mpi
  if (save_movie_dt /= 0) then
    call init_chan_movie
  end if

  call init_flow
  call set_flow
  call pre_first_step(compute_pressure)

  ! Initialize start_wall_time for run timing
  call wall_time(start_wall_time)
  call save_stats_chan(save_movie_dt/=0,.false.)
  if (use_LES) call save_stats_LES_OOL(.true.)
  call wall_time(previous_wall_time)
  if (rank == 0) then
    write (*,'("Elapsed Wall Time to Save Stats: ", ES13.3)') (previous_wall_time - start_wall_time)
  end if
  

  if (rank == 0) then
    write (*, *)
    write (*, *) '             ****** Done Initialising ******'
    write (*, *)
  end if

  ! main loop
  do
    time_step = time_step + 1

    if (verbosity > 2 .and. rank == 0) &
      write (*, '("Now beginning time step ", I10)') time_step

    do rk_step = 1, 3
      if (time_ad_meth == 1) call rk_chan_1
      if (time_ad_meth == 2) call rk_chan_2
    end do
    time = time + delta_t

    ! Apply filters to the scalar fields if on
    call filter_scalars
    
    ! Check if we meet any stop conditions
    call end_run_mpi(flag)

    ! Save statistics to an output file
    if (time >= save_stats_time) then
      call save_stats
      if (flag .and. use_LES) then
        ! if stopping we must save LES stats now
        call save_stats_LES_OOL(.false.)
      end if 
    end if

    ! Save entire flow to a file 
    if (time >= save_flow_time) then
      call save_flow(.false.)
    end if

    ! Check if we're done
    if (flag) then 
      exit
    end if

  end do

  ! Calculate and display the runtime for the simulation
  call wall_time(end_wall_time)
  if (rank == 0) then
    write (*, '("Elapsed Wall Time (sec): ", ES11.3)') end_wall_time - start_wall_time
    write (*, '("Wall Seconds per Iteration: ", ES12.5)') (end_wall_time - start_wall_time) / time_step
  end if

  ! create the end file
  call save_flow(.true.)

  ! tidy up
  !call deallocate_all
  call mpi_finalize(ierror)

end
