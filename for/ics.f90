! set the initial conditions or read in from a file
module ics 
  use parameters 
  use domain 
  use fft
  use flow
  use phdf5
  implicit none

  integer :: num_read_th
  integer :: read_th_index(1:N_th)
  logical :: compute_pressure

contains

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine set_flow 
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    integer :: n
    ! Initialize values for reading of scalars
    num_read_th = 0
    do n = 1, N_th
      if (create_new_th(n)) then
        num_read_th = num_read_th
      else
        num_read_th = num_read_th + 1
        read_th_index(num_read_th) = n
      end if
    end do
    call create_th_chan

    ! Create flow.
    if (create_new_flow) then
      call create_flow_chan
      compute_pressure = .true.
      if (rank == 0) &
        write (*, '("A new flowfield has been created.")')
    else
      if (rank == 0) &
        write (*, '("Reading flow...")')
      call read_flow
      if (rank == 0) &
        write (*, '("Starting at time step ", I10)') time_step

    end if 
  
    return
  end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine create_flow_chan
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|


    integer i, j, k, n
    real(rkind) rnum1, rnum2, rnum3
    integer, dimension(:), allocatable :: seed

    ! Initialize the random number generator
    call random_seed(size=k)
    allocate (seed(1:k))
    do i = 1, k
      seed(i) = rank*k + i + 999
    end do
    call random_seed(put=seed)

    ! UBULK0 and kick should be set in input.dat

    ! IC_type is set in input_chan.dat and can be used to easily
    ! control which initial condition is used.  A few examples
    ! are given here. These can be modified, or new types can be
    ! added

    if (IC_type == 0) then
      ! Parabolic profile for laminar closed channel flow
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = (3./2.) * ubulk0 * (1.d0 - gyf(j)**2.)
            u2(i, k, j) = 0.
            u3(i, k, j) = 0.
          end do
        end do
      end do
    else if (IC_type == 1) then
      ! Laminar profile for open channel flow
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          do j = 1, Nyp
            u1(i, k, j) = -(3./2.) * ubulk0 * gyf(j)**2.+3.*ubulk0 * gyf(j)
            u2(i, k, j) = 0.
            u3(i, k, j) = 0.
          end do
          u1(i, k, 0) = 0.
          u3(i, k, 0) = 0.
          u1(i, k, Nyp + 1) = 0.
          u3(i, k, Nyp + 1) = 0.
        end do
      end do
    else if (IC_type == 2) then
      ! Linear profile for laminar Couette flow
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = gyf(j)
            u2(i, k, j) = 0.
            u3(i, k, j) = 0.
          end do
        end do
      end do
    else if (IC_type == 3) then
      ! Tanh shear layer
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = tanh(gyf(j))
            u2(i, k, j) = 0.d0
            u3(i, k, j) = 0.d0
          end do
        end do
      end do
    else if (IC_type == 4) then
      ! Infinite Front
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = 0.d0
            u2(i, k, j) = 0.d0
            u3(i, k, j) = 0.d0 ! TWS exactly balanced in solver
          end do
        end do
      end do
    else if (IC_type == 5) then
      ! Balanced Infinite Front + Barotropic Background vorticity (W)
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = 0.d0
            u2(i, k, j) = 0.d0
            u3(i, k, j) = (tanh((gx(i) - 0.5 * Lx) / delta) - &
                              (tanh(0.5*Lx/delta) - tanh(-0.5*Lx/delta)) / Lx * (gx(i) - 0.5 * Lx)) & ! Make periodic
                            / (delta*(1.d0/delta - (tanh(0.5*Lx/delta)-tanh(-0.5*Lx/delta))/Lx)) ! Scale so middle gradient is exactly 1/delta
          end do
        end do
      end do
    else if (IC_type == 6) then
      ! Finite Sloping Front
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = 0.d0
            u2(i, k, j) = 0.d0
            if (Ri(1) == 0) then
              u3(i, k, j) = (gyf(j) - 0.5 * Ly) / &
                            cosh((gx(i) - 0.5 * Lx) / delta)**2.d0 &
                            / tanh(0.5d0 * Lx / delta) &
                            - 2.d0 * delta / Lx * (gyf(j) - 0.5 * Ly)
            else
              u3(i, k, j) = delta * sqrt(Ri(1)**2.d0 + Ro_inv**2.d0) / Ri(1) * &
                            tanh(((gx(i) - 0.5 * Lx) * Ro_inv + &
                                  (gyf(j) - 0.5 * Ly) * Ri(1)) / &
                                 (delta * sqrt(Ri(1)**2.d0 + Ro_inv**2.d0))) &
                            - 2.d0 * delta / Lx * (gyf(j) - 0.5 * Ly)
            end if
          end do
        end do
      end do
    else if (IC_type == 7) then
      ! Finite Front, with u3 = 0 at bottom (for no-slip BC)
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = 0.d0
            u2(i, k, j) = 0.d0
            u3(i, k, j) = (gyf(j)) / &
                            cosh((gx(i) - 0.5 * Lx) / delta)**2.d0 &
                            / tanh(0.5d0 * Lx / delta) &
                            - 2.d0 * delta / Lx * (gyf(j))
          end do
        end do
      end do
    else if (IC_type == 8) then
      ! Finite Front Unbalanced
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = 0.d0
            u2(i, k, j) = 0.d0
            u3(i, k, j) = -2.d0 * delta / Lx * (gyf(j) - 0.5 * Ly) ! Remove the TWS from solver
          end do
        end do
      end do
    else
      write (*, '("Warning, unsupported IC_type in create_flow.")')
    end if


    !call add_SI_Mode


    if (physical_noise) then
      ! Add random noise in physical space
      call random_number(rnum1)
      call random_number(rnum1)
      call random_number(rnum1)
      do j = 0, Nyp + 1
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            call random_number(rnum1)
            u1(i, k, j) = u1(i, k, j) + kick * (rnum1 - 0.5d0)
            call random_number(rnum1)
            u2(i, k, j) = u2(i, k, j) + kick * (rnum1 - 0.5d0)
            call random_number(rnum1)
            u3(i, k, j) = u3(i, k, j) + kick * (rnum1 - 0.5d0)
          end do
        end do
      end do


      ! Convert to Fourier space
      call fft_xz_to_fourier(u1, cu1)
      call fft_xz_to_fourier(u2, cu2)
      call fft_xz_to_fourier(u3, cu3)
      call fft_xz_to_fourier(p, cp)

    else ! Optionally, add random noise in Fourier space instead

      ! Convert to Fourier space
      call fft_xz_to_fourier(u1, cu1)
      call fft_xz_to_fourier(u2, cu2)
      call fft_xz_to_fourier(u3, cu3)
      call fft_xz_to_fourier(p, cp)


      do i = 0, Nxp - 1
        do j = 1, Nyp
          do k = 0, twoNkz
            ! Now, give the velocity field a random perturbation
            call random_number(rnum1)
            call random_number(rnum2)
            cu1(i,k,j) = cu1(i,k,j) + cmplx((rnum1 - 0.5d0), (rnum2 - 0.5d0)) * kick
            call random_number(rnum1)
            call random_number(rnum2)
            cu2(i,k,j) = cu2(i,k,j) + cmplx((rnum1 - 0.5d0), (rnum2 - 0.5d0)) * kick
            call random_number(rnum1)
            call random_number(rnum2)
            cu3(i,k,j) = cu3(i,k,j) + cmplx((rnum1 - 0.5d0), (rnum2 - 0.5d0)) * kick
          end do
          if (twoNkz == 0) then
            ! Here, In the 2d case we want to add a kick to the mean in z
            k = 0
            call random_number(rnum1)
            call random_number(rnum2)
            call random_number(rnum3)

            if (IC_type == 3) then
              cu1(i,k,j) = cu1(i,k,j) + (rnum1 - 0.5) * kick * exp(-(gyf(j) * 20.d0)**2.d0)
              cu2(i,k,j) = cu2(i,k,j) + (rnum1 - 0.5) * kick * exp(-(gyf(j) * 20.d0)**2.d0)
              cu3(i,k,j) = cu3(i,k,j) + (rnum1 - 0.5) * kick * exp(-(gyf(j) * 20.d0)**2.d0)
            else
              cu1(i,k,j) = cu1(i,k,j) + (rnum1 - 0.5) * kick
              cu2(i,k,j) = cu2(i,k,j) + (rnum2 - 0.5) * kick
              cu3(i,k,j) = cu3(i,k,j) + (rnum3 - 0.5) * kick
            end if
          end if
        end do
      end do

    end if

    return
  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine Add_SI_Mode
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! Add exact SI Mode Perturbations


    integer i, j, k, n
    real(rkind) lambda1, lambda2, k_max, l_max, sigma
    complex(rkind) amplitude_scaled

    if (Ro_inv == 1/1.d0) then
      ! Ro = 1
      lambda1 = -19.943848405400693
      lambda2 = -13.660663098433014
      k_max = 26.922065852762007
      l_max = 17.093431565855951
      sigma = 0.765902833138420
    else if (Ro_inv == 1/10.d0) then
      ! Ro = 10
      lambda1 = -17.235267629576086
      lambda2 = -10.952082328361845
      k_max = 14.461814598675044
      l_max = 14.439573333787491
      sigma = 0.300145893642542
    else if (Ro_inv == 1/0.1d0) then
      ! Ro = 0.1
      lambda1 = -11.590816435595643
      lambda2 = -5.307631128443722
      k_max = 85.214180782654623
      l_max = 9.014376678159723
      sigma = 0.850936873842555
    else
      stop 'Exact Mode not defined for this Ro!'
    end if


    do j = 0, Nyp + 1
      do k = 0, Nzp - 1
        do i = 0, Nxm1

          amplitude_scaled = (0, 1) * (lambda1 * exp((0, 1) * lambda1 * gyf(j)) - &
                           lambda2 * exp((0, 1) * lambda2 * gyf(j))) * &
                 exp((0, -1) * k_max * gx(i))
          u1(i, k, j) = u1(i, k, j) + real(amplitude_scaled) * kick * 100

          amplitude_scaled = (0, -1) * k_max * &
                 (exp((0, 1) * lambda1 * gyf(j)) - exp((0, 1) * lambda2 * gyf(j))) * &
                 exp((0, -1) * k_max * gx(i))
          u2(i, k, j) = u2(i, k, j) + real(amplitude_scaled) * kick * 100

          amplitude_scaled = -(u1(i, k, j) * Ro_inv + u2(i, k, j)) / &
                 (sigma + nu * (k_max * k_max + l_max * l_max))
          u3(i, k, j) = u3(i, k, j) + real(amplitude_scaled) * kick * 100
        end do
      end do
    end do

    Return
  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine create_th_chan
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! Initialize the scalar fields
    ! In this subroutine, you should initialize each scalar field for the
    ! particular problem of interest


    integer i, j, k, n

    do n = 1, N_th
      if (create_new_th(n)) then

        if (IC_type == 0) then
          ! As an example, initialize TH1 with a sine in x
          do k = 0, Nzp - 1
            do i = 0, Nxm1
              do j = 1, Nyp
                th(i, k, j, n) = Ri(n) * sin(2.d0 * pi * gx(i) / Lx) / (4.d0 * pi**2.d0)
              end do
            end do
          end do
        else if ((IC_type == 1) .or. (IC_type == 2)) then
          ! Initialize with a linear profile using the BCs
          do k = 0, Nzp - 1
            do i = 0, Nxm1
              if ((th_BC_Ymin(n) == 0) .and. (th_BC_Ymax(n) == 0)) then
                do j = 1, Nyp
                  if (gyf(j) <= 2.0) then
                    th(i, k, j, n) = (th_BC_Ymax_c1(n) - th_BC_Ymin_c1(n)) &
                                     * (gyf(j) + 1.) / 2.0 + th_BC_Ymin_c1(n)
                  else
                    th(i, k, j, n) = th_BC_Ymax_c1(n)
                  end if
                end do
              else if ((th_BC_Ymin(n) == 1) &
                       .and. (th_BC_Ymax(n) == 1)) then
                do j = 1, Nyp
                  ! Linear profile with slope corresponding to lower value
                  th(i, k, j, n) = th_BC_Ymin_c1(n) * gyf(j)
                end do
              else
                if (rank == 0) then
                  write (*, '("Warning, THETA INITIALIZED TO ZERO ...")')
                  write (*, '("Create an initial value in create_flow_chan")')
                end if
              end if
            end do
          end do
        else if (IC_type == 3) then
          ! Tanh vertical profile (for shear layer)
          do k = 0, Nzp - 1
            do i = 0, Nxm1
              do j = 1, Nyp
                th(i, k, j, n) = Ri(n) * tanh(gyf(j))
              end do
            end do
          end do
        else if (IC_type == 4 .or. IC_type == 5) then
          ! Infinite Front (5 includes background sinusoidal lateral W shear)
          do k = 0, Nzp - 1
            do i = 0, Nxm1
              th(i, k, 0, n) = 0.d0
              do j = 1, Nyp
                th(i, k, j, n) = Ri(n) * (gyf(j) - 0.5 * Ly)
              end do
            end do
          end do
        else if (IC_type == 6 .or. IC_type == 7 .or. IC_type == 8) then
          ! Finite Sloped Front
          do k = 0, Nzp - 1
            do i = 0, Nxm1
              do j = 1, Nyp
                if (Ri(n) == 0) then
                  th(i, k, j, n) = Ro_inv &
                                   * tanh((gx(i) - 0.5 * Lx) / delta) &
                                   / tanh(0.5d0 * Lx / delta) &
                                   + Ro_inv * (1.d0 - 2.d0 / Lx * gx(i))
                else
                  th(i, k, j, n) = sqrt(Ri(n)**2.d0 + Ro_inv**2.d0) &
                                   * tanh(((gx(i) - 0.5 * Lx) * Ro_inv + &
                                           (gyf(j) - 0.5 * Ly) * Ri(n)) / &
                                          (delta * sqrt(Ri(n)**2.d0 + Ro_inv**2.d0))) &
                                   + Ro_inv * (1.d0 - 2.d0 / Lx * gx(i))
                end if
              end do
            end do
          end do
        else
          write (*, '("Warning, unsupported IC_type in create_flow.")')
        end if

        call fft_xz_to_fourier(th(:, :, :, n), cth(:, :, :, n))

      end if
    end do

    return
  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine read_flow
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    character(len=55) :: fname
    integer :: i, j, k, n
    logical :: read_pressure

    fname = 'start.h5'
    if (rank == 0) &
      write (*, '("Reading flow from " A10)') fname

    call mpi_barrier(mpi_comm_world, ierror)
    call ReadHDF5(fname,read_pressure)
    compute_pressure = (.not. read_pressure)
    
    return
  end

end module ics