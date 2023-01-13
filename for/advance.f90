! advance a timestep
module advance 
  use parameters
  use domain 
  use flow 
  use fft
  use boundary
  use les
  implicit none 

contains 

  include 'channel.f90'
  include 'forcing.f90'

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine pre_first_step(compute_pressure)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! This subroutine ensures the initial velocity is divergence free,
    ! sets the boundary conditions, ghost cells and if necesary computes 
    ! the pressure
    logical, intent(in) :: compute_pressure 

    ! Apply Boundary conditions to velocity field
    call apply_BC_vel_mpi_post
    call ghost_chan_mpi


    ! Remove the divergence of the velocity field
    call rem_div_chan
    call ghost_chan_mpi

    ! Get the pressure from the poisson equation
    if (compute_pressure) then
      call poisson_p_chan
      ! Fix for the pressure
      call ghost_chan_mpi
    end if

    if (variable_dt) then
      call courant
    end if

  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine courant
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! This subroutine sets the timestep based on the specified CFL number
    ! The subroutine should be called with the velocity in physical space

    real(rkind) dt_x, dt_y, dt_z
    real(rkind) Nmax
    integer i, j, k, n
    integer imin, jmin, kmin
    real(rkind) r_max ! Maximum fractional change in dt

    ! Set the initial dt to some arbitrary large number
    dt = max_dt

    ! Set the timestep based on viscosity and diffusivity
    ! dt = min(dt, 0.5d0 * min(dx(1), dz(1))**(2.d0 * beta) / nu)
    ! do n = 1, N_th
    !   dt = min(dt, dt * nu / (nu / Pr(n)))
    ! end do
    ! Make sure that we capture the inertial period (for rotating flows)
    if (Ro_inv /= 0.d0) then
      dt = min(dt, 2.d0 * pi / abs((Ro_inv / delta)) / 20.d0)
    end if
    ! Make sure that we capture the buoyancy period (for stratified flows)
    do n = 1, N_th
      Nmax = sqrt(abs(th_BC_Ymin_c1(n)))
      dt = min(dt, 0.1 * 2.d0 * pi / Nmax)
    end do

    ! Make sure we capture the Geostrophic Flow-through time
    !   because the split advection component may still be at most Vg !
    !dt=min(dt,min(dx(1),dy(1))/(LY*delta/Ro_inv*dTHdX))

    ! Use the model velocity to calculate the CFL number

    if (flavor == 'Front') then
      ! Add thermal wind to velocity when calculating the CFL number
      do n = 1, N_th
        do j = jstart, jend
          do k = 0, Nzp - 1
            do i = 0, Nxm1
              dt_x = CFL * dx(i) / abs(u1(i, k, j) - dTHdZ(n) &
                                      * gyf(j) / (Ro_inv / delta))
              dt_y = CFL * dy(j) / abs(u2(i, k, j))
              dt_z = CFL * dz(k) / abs(u3(i, k, j) + (1.d0 / (Ro_inv / delta)) &
                                      * dTHdX(n) * (gyf(j) - 0.5d0*Ly))
              dt = min(dt, dt_x, dt_y, dt_z)
            end do
          end do
        end do
      end do
    else
      do j = 1, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            dt_x = CFL * dx(i) / abs(u1(i, k, j))
            dt_y = CFL * dy(j) / abs(u2(i, k, j))
            dt_z = CFL * dz(k) / abs(u3(i, k, j))
            dt = min(dt, dt_x, dt_y, dt_z)
          end do
        end do
      end do
    end if

    call get_minimum_mpi(dt)

    if (dt <= 0) then
      if (rank == 0) &
        write (*, '("Error: dt <= 0 in courant")')
      ! Set DELTA_T to some small default value
      delta_t = 0.00001d0
    else if (dt >= 1000.) then
      write (*, '("Warning: dt > 1000, value capped at 1000")')
      delta_t = 1000.d0
    else
      delta_t = dt
    end if

    ! Find the next event we need to meet at
    delta_t_next_event = min(time_limit, &
                            save_stats_time - time, &
                            save_flow_time - time,  &
                            time_limit - time)

    ! Try to schedule the next 2 time-steps to not be too small (keep CFL ~ 0.5)
    r_max = (1. - 0.025)
    if (delta_t_next_event < dt) then ! Directly go there (collect $200)
      delta_t = delta_t_next_event + 1.d-14
    else if (delta_t_next_event < (1.+r_max)*dt) then ! Two steps this time
      delta_t = 1./(1.+r_max)*delta_t_next_event
    else
      delta_t = dt
    end if

    ! if (time + delta_t > save_stats_time) then
    !   delta_t = save_stats_time - time + 1.d-14
    ! else if (time + delta_t > save_flow_time) then
    !   delta_t = save_flow_time - time + 1.d-14
    ! else if (time + delta_t > time_limit) then
    !   delta_t = time_limit - time
    ! end if

    if (rank == 0 .and. delta_t < 1.d-10) &
      write (*, '("Warning: dt < 1e-10, so something may be going wrong!")')


    h_bar(1) = delta_t * (8.d0 / 15.d0)
    h_bar(2) = delta_t * (2.d0 / 15.d0)
    h_bar(3) = delta_t * (5.d0 / 15.d0)

    return
  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine filter_scalars
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    integer :: n

    do n = 1, N_th
      if (filter_th(n) &
          .and. (mod(time_step, filter_int(n)) == 0)) then
        call filter_chan(n)
      end if
    end do

  end

end module advance