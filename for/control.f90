! subroutines and timing parameters for saving stats and ending run
! called from the program main loop, defined here mostly for tidiness
module control
  use parameters
  use domain
  use flow
  use phdf5 
  use statistics, only: save_stats_chan
  implicit none 

  real(rkind) :: start_wall_time, previous_wall_time, end_wall_time

contains 

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine wall_time(wt)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    !
    ! Return wall-clock time as seconds after Jan. 1, 2016.
    ! Support for leap year is not included anymore.
    !
    ! By using a 'save' statement, the wall-time after the first
    ! call to the subroutine could be computed, but that is not
    ! intended with the present subroutine (e.g. the history file)
    !
    implicit none

    real(rkind), intent(out) :: wt
    integer :: val(8), i, shift, day

    integer mon(12, 2)
    data mon/ &
      31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
      31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
    !
    ! Get current date and time
    ! val(1) : year
    ! val(2) : month
    ! val(3) : day
    ! val(4) : difference to GMT
    ! val(5) : hour
    ! val(6) : minute
    ! val(7) : second
    ! val(8) : 1/1000 second
    !
    call date_and_time(values=val)
    !
    ! Determine leap year
    !
    if (mod(val(1), 4) == 0) then
      if (mod(val(1), 100) == 0) then
        if (mod(val(1), 400) == 0) then
          shift = 2
        else
          shift = 1
        end if
      else
        shift = 2
      end if
    else
      shift = 1
    end if
    !
    ! Construct day of the year
    !
    day = val(3) - 1
    do i = 1, val(2) - 1
      day = day + mon(i, shift)
    end do
    !
    ! And compute wall-clock time
    !
    wt = (val(1) - 2016) * 365 * 86400 + &
        day * 86400 + val(5) * 3600 + val(6) * 60 + val(7) + dble(val(8) / 1000.d0)

  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine save_stats
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    real(rkind) :: wall_begin

    call wall_time(wall_begin)

    
    if (time >= save_movie_time) then
      call save_stats_chan(.true.,.false.)
      save_movie_time = save_stats_time + save_movie_dt - save_stats_dt*1.d-5 ! Increment from stats_time
    else
      call save_stats_chan(.false.,.false.)
    end if

    save_stats_time = save_stats_time + save_stats_dt
    flag_save_LES = .true. ! save LES stats on next RK step

    call wall_time(end_wall_time)

    ! Timing Diagnostics

    if (rank == 0) then
      write (*,'("Elapsed Wall Time to Save Stats: ", ES13.3)') (end_wall_time - wall_begin)

      write (*,'("Wall Seconds per Stats Output: ", ES15.3)') &
                (end_wall_time - previous_wall_time)
      write (*,'("Wall Seconds per Iteration: ", ES18.3)') &
                (end_wall_time - previous_wall_time) / float(time_step - previous_time_step)
      write (*,'("Wall Seconds per Simulation Time: ", ES12.3)') &
                (end_wall_time - previous_wall_time) / save_stats_dt

      write (*, *)

    end if
    call wall_time(previous_wall_time)
    previous_time_step = time_step
    
  end 

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine save_flow(final)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    logical, intent(in) :: final 
    real(rkind) :: wall_begin

    call wall_time(wall_begin)

    if (rank == 0) &
      write (*, '("Saving Output File")')
    call write_flow(final)
    call wall_time(end_wall_time)
    if (rank == 0) &
      write (*,'("Elapsed Wall Time to Save Flow: ", ES12.5)') (end_wall_time - wall_begin)

    save_flow_time = save_flow_time + save_flow_dt

  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine end_run(flag)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    logical, intent(out) :: flag 
    logical :: file_exists

    flag = .false.
    ! Check for the time
    call wall_time(end_wall_time)
    if (end_wall_time - start_wall_time > wall_time_limit) then
      if (rank == 0) &
        write (*, '("STOP because of wall-time hit!")')
      flag = .true.
    end if

    if (time >= time_limit) then
      if (rank == 0) &
        write (*, '("STOP because of simulation end-time hit!")')
      flag = .true.
    end if

    inquire (file="stop.now", exist=file_exists)
    if (file_exists) then
      if (rank == 0) &
        write (*, '("STOP because of stop.now file!")')
      flag = .true.
    end if

    return
  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine end_run_mpi(flag)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    logical, intent(out) :: flag

    if (rank == 0) then
      call end_run(flag)
    end if
    call mpi_bcast(flag, 1, mpi_logical, 0, mpi_comm_world, ierror)

  end

end module control