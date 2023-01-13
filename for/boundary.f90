! routines for setting boundary conditions and ghost cells
! these routines are used by the advance, les, and statistics modules
module boundary 
    use parameters
    use domain
    use fft 
    use flow 
    implicit none 
    
contains 

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine ghost_chan_mpi
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! This subroutine is part of the MPI package for the channel flow
    ! Diablo package.
    ! Here, we define a set of ghost cells on each process
    ! the ghost cells contain information from the neighboring nodes
    ! and allow us to compute finite differences over the local gridpoints.
    ! We need to update the contents of the ghost cells at the start of
    ! each Runge-Kutta substep
  
    integer i, j, k, n
  
    ! Define the arrays that will be used for data packing.  This makes the
    ! communication between processes more efficient by only requiring one
    ! send and recieve.
    ! The communication will be done in Fourier space, so these arrays should
    ! be complex arrays to match the velocity
    ! The size of the buffer array is 0:Nkx,0:twoNkz,# of variables
    complex(rkind) ocpack(0:Nxp - 1, 0:twoNkz, 4 + N_th)
    complex(rkind) icpack(0:Nxp - 1, 0:twoNkz, 4 + N_th)
  
    ! If we are using more than one processor, then we need to pass data
  
    if (NprocY > 1) then
  
      ! First, Pass data up the chain to higher ranked processes
  
      if (rankY == 0) then
        ! If we are the lowest ranked process, then we don't need to recieve
        ! data at the lower ghost cells, these will be filled with boundary
        ! condition information
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            ocpack(i, k, 1) = cu1(i, k, Nyp)
            ocpack(i, k, 2) = cu2(i, k, Nyp)
            ocpack(i, k, 3) = cu3(i, k, Nyp)
            ocpack(i, k, 4) = cp(i, k, Nyp)
            do n = 1, N_th
              ocpack(i, k, 4 + n) = cth(i, k, Nyp, n)
            end do
          end do
        end do
        ! Now, we have packed the data into a compact array, pass the data up
        call mpi_send(ocpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                      , mpi_double_complex &
                      , rankY + 1, 1, mpi_comm_y, ierror)
  
        ! End if RANK=0
      else if (rankY < NprocY - 1) then
        ! Here, we are one of the middle processes and we need to pass data
        ! up and recieve data from below
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            ocpack(i, k, 1) = cu1(i, k, Nyp)
            ocpack(i, k, 2) = cu2(i, k, Nyp)
            ocpack(i, k, 3) = cu3(i, k, Nyp)
            ocpack(i, k, 4) = cp(i, k, Nyp)
            do n = 1, N_th
              ocpack(i, k, 4 + n) = cth(i, k, Nyp, n)
            end do
          end do
        end do
        ! Use MPI_SENDRECV since we need to recieve and send data
        call mpi_send(ocpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                      , mpi_double_complex &
                      , rankY + 1, 1, mpi_comm_y, ierror)
  
        call mpi_recv(icpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                      , mpi_double_complex &
                      , rankY - 1, 1, mpi_comm_y, status, ierror)
        ! Now, unpack the data that we have recieved
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cu1(i, k, 1) = icpack(i, k, 1)
            cu2(i, k, 1) = icpack(i, k, 2)
            cu3(i, k, 1) = icpack(i, k, 3)
            cp(i, k, 1) = icpack(i, k, 4)
            do n = 1, N_th
              cth(i, k, 1, n) = icpack(i, k, 4 + n)
            end do
          end do
        end do
  
      else
        ! Otherwise, we must be the uppermost process with RANK=Nprocs-1
        ! Here, we need to recieve data from below, but don't need to send data up
        call mpi_recv(icpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                      , mpi_double_complex &
                      , rankY - 1, 1, mpi_comm_y, status, ierror)
        ! Unpack the data that we have recieved
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cu1(i, k, 1) = icpack(i, k, 1)
            cu2(i, k, 1) = icpack(i, k, 2)
            cu3(i, k, 1) = icpack(i, k, 3)
            cp(i, k, 1) = icpack(i, k, 4)
            do n = 1, N_th
              cth(i, k, 1, n) = icpack(i, k, 4 + n)
            end do
          end do
        end do
      end if
  
      ! AT this point we have passed data up the chain
      if (rankY == NprocY - 1) then
        ! If we are the higest ranked process, then we don't need to recieve
        ! data at the upper ghost cells, these will be filled with boundary
        ! condition information
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            ocpack(i, k, 1) = cu1(i, k, 2)
            ocpack(i, k, 2) = cu2(i, k, 2)
            ocpack(i, k, 3) = cu3(i, k, 2)
            ocpack(i, k, 4) = cp(i, k, 2)
            do n = 1, N_th
              ocpack(i, k, 4 + n) = cth(i, k, 2, n)
            end do
          end do
        end do
        ! Now, we have packed the data into a compact array, pass the data up
        call mpi_send(ocpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                      , mpi_double_complex &
                      , rankY - 1, 3, mpi_comm_y, ierror)
      else if (rankY > 0) then
        ! Here, we are one of the middle processes and we need to pass data
        ! down and recieve data from above us
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            ocpack(i, k, 1) = cu1(i, k, 2)
            ocpack(i, k, 2) = cu2(i, k, 2)
            ocpack(i, k, 3) = cu3(i, k, 2)
            ocpack(i, k, 4) = cp(i, k, 2)
            do n = 1, N_th
              ocpack(i, k, 4 + n) = cth(i, k, 2, n)
            end do
          end do
        end do
  
        call mpi_send(ocpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                      , mpi_double_complex &
                      , rankY - 1, 3, mpi_comm_y, ierror)
  
        call mpi_recv(icpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                      , mpi_double_complex &
                      , rankY + 1, 3, mpi_comm_y, status, ierror)
        ! Now, unpack the data that we have recieved
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cu1(i, k, Nyp + 1) = icpack(i, k, 1)
            cu2(i, k, Nyp + 1) = icpack(i, k, 2)
            cu3(i, k, Nyp + 1) = icpack(i, k, 3)
            cp(i, k, Nyp + 1) = icpack(i, k, 4)
            do n = 1, N_th
              cth(i, k, Nyp + 1, n) = icpack(i, k, 4 + n)
            end do
          end do
        end do
      else
        ! Here, we must be the lowest process (RANK=0) and we need to recieve
        ! data from above
        call mpi_recv(icpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                      , mpi_double_complex &
                      , rankY + 1, 3, mpi_comm_y, status, ierror)
        ! Unpack the data that we have recieved
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cu1(i, k, Nyp + 1) = icpack(i, k, 1)
            cu2(i, k, Nyp + 1) = icpack(i, k, 2)
            cu3(i, k, Nyp + 1) = icpack(i, k, 3)
            cp(i, k, Nyp + 1) = icpack(i, k, 4)
            do n = 1, N_th
              cth(i, k, Nyp + 1, n) = icpack(i, k, 4 + n)
            end do
          end do
        end do
      end if
  
    end if
  
    return
  end
  
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine ghost_chan_mpi_j0
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! ghost_chan_mpi subroutine to update the j = 0 cell
    ! Needed for les output and epsilon calculation
  
    integer i, j, k, n
  
    ! Define the arrays that will be used for data packing.  This makes the
    ! communication between processes more efficient by only requiring one
    ! send and recieve.
    ! The communication will be done in Fourier space, so these arrays should
    ! be complex arrays to match the velocity
    ! The size of the buffer array is 0:Nkx,0:twoNkz,# of variables
    complex(rkind) ocpack(0:Nxp - 1, 0:twoNkz, 4 + N_th)
    complex(rkind) icpack(0:Nxp - 1, 0:twoNkz, 4 + N_th)
  
    ! If we are using more than one processor, then we need to pass data
  
    if (NprocY > 1) then
  
      if (rankY == 0) then
        ! If we are the lowest ranked process, then we don't need to recieve
        ! data at the lower ghost cells, these will be filled with boundary
        ! condition information
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            ocpack(i, k, 1) = cu1(i, k, Nyp - 1)
            ocpack(i, k, 2) = cu2(i, k, Nyp - 1)
            ocpack(i, k, 3) = cu3(i, k, Nyp - 1)
            ocpack(i, k, 4) = cp(i, k, Nyp - 1)
            do n = 1, N_th
              ocpack(i, k, 4 + n) = cth(i, k, Nyp - 1, n)
            end do
          end do
        end do
        ! Now, we have packed the data into a compact array, pass the data up
        call mpi_send(ocpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                      , mpi_double_complex &
                      , rankY + 1, 1, mpi_comm_y, ierror)
  
        ! End if RANK=0
      else if (rankY < NprocY - 1) then
        ! Here, we are one of the middle processes and we need to pass data
        ! up and recieve data from below
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            ocpack(i, k, 1) = cu1(i, k, Nyp - 1)
            ocpack(i, k, 2) = cu2(i, k, Nyp - 1)
            ocpack(i, k, 3) = cu3(i, k, Nyp - 1)
            ocpack(i, k, 4) = cp(i, k, Nyp - 1)
            do n = 1, N_th
              ocpack(i, k, 4 + n) = cth(i, k, Nyp - 1, n)
            end do
          end do
        end do
        ! Use MPI_SENDRECV since we need to recieve and send data
        call mpi_send(ocpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                      , mpi_double_complex &
                      , rankY + 1, 1, mpi_comm_y, ierror)
  
        call mpi_recv(icpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                      , mpi_double_complex &
                      , rankY - 1, 1, mpi_comm_y, status, ierror)
        ! Now, unpack the data that we have recieved
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cu1(i, k, 0) = icpack(i, k, 1)
            cu2(i, k, 0) = icpack(i, k, 2)
            cu3(i, k, 0) = icpack(i, k, 3)
            cp(i, k, 0) = icpack(i, k, 4)
            do n = 1, N_th
              cth(i, k, 0, n) = icpack(i, k, 4 + n)
            end do
          end do
        end do
  
      else
        ! Otherwise, we must be the uppermost process with RANK=Nprocs-1
        ! Here, we need to recieve data from below, but don't need to send data up
        call mpi_recv(icpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                      , mpi_double_complex &
                      , rankY - 1, 1, mpi_comm_y, status, ierror)
        ! Unpack the data that we have recieved
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cu1(i, k, 0) = icpack(i, k, 1)
            cu2(i, k, 0) = icpack(i, k, 2)
            cu3(i, k, 0) = icpack(i, k, 3)
            cp(i, k, 0) = icpack(i, k, 4)
            do n = 1, N_th
              cth(i, k, 0, n) = icpack(i, k, 4 + n)
            end do
          end do
        end do
  
      end if
  
    end if
  
    return
  end
  
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  subroutine apply_BC_vel_lower_post
    !----*|--.---------.---------.---------.---------.---------.---------.-|--
    ! This subroutine is called after initializing the flow
    ! It sets the appropriate boundary conditions including ghost cell values
    !  on the velocity field in Fourier space
  
    integer i, k
  
    ! Now, apply the boundary conditions depending on the type specified
    if (u_BC_Ymin == 0) then
      ! Dirichlet
      ! Start with zero
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu1(i, k, 0) = 0.d0
          cu1(i, k, 1) = 0.d0
        end do
      end do
      ! Now, set only the mean
      if (rankZ == 0) then
        cu1(0, 0, 1) = u_BC_Ymin_c1
        cu1(0, 0, 0) = u_BC_Ymin_c1
      end if
    else if (u_BC_Ymin == 1) then
      ! Neumann
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu1(i, k, 1) = cu1(i, k, 2)
          cu1(i, k, 0) = cu1(i, k, 1)
        end do
      end do
      ! Now, Apply BC to mean
      if (rankZ == 0) then
        cu1(0, 0, 1) = cu1(0, 0, 2) - dy(2) * u_BC_Ymin_c1
        cu1(0, 0, 0) = cu1(0, 0, 1) - dy(1) * u_BC_Ymin_c1
      end if
    else
      stop 'Error: u_BC_Zmin must be 0, or 1'
    end if
  
    if (w_BC_Ymin == 0) then
      ! Dirichlet
      ! Start with zero
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu3(i, k, 0) = 0.d0
          cu3(i, k, 1) = 0.d0
        end do
      end do
      ! Now, set only the mean
      if (rankZ == 0) then
        cu3(0, 0, 1) = w_BC_Ymin_c1
        cu3(0, 0, 0) = w_BC_Ymin_c1
      end if
    else if (w_BC_Ymin == 1) then
      ! Neumann
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu3(i, k, 1) = cu3(i, k, 2)
          cu3(i, k, 0) = cu3(i, k, 1)
        end do
      end do
      ! Now, Apply BC to mean
      if (rankZ == 0) then
        cu3(0, 0, 1) = cu3(0, 0, 2) - dy(2) * w_BC_Ymin_c1
        cu3(0, 0, 0) = cu3(0, 0, 1) - dy(1) * w_BC_Ymin_c1
      end if
    else
      stop 'Error: v_BC_Zmin must be 0, 1, or 2'
    end if
  
    if (v_BC_Ymin == 0) then
      ! Dirichlet
      ! Set the vertical velocity at GYF(1) (halfway between GY(2) and GY(1))
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu2(i, k, 1) = 2.d0 * v_BC_Ymin_c1 - cu2(i, k, 2)
          cu2(i, k, 0) = cu2(i, k, 1)
        end do
      end do
    else if (v_BC_Ymin == 1) then
      ! Neumann
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu2(i, k, 1) = cu2(i, k, 2)
          cu2(i, k, 0) = cu2(i, k, 1)
        end do
      end do
      if (rankZ == 0) then
        cu2(0, 0, 1) = cu2(0, 0, 2) - dyf(1) * v_BC_Ymin_c1
        cu2(0, 0, 0) = cu2(0, 0, 1) - dyf(1) * v_BC_Ymin_c1
      end if
    else
      stop 'Error: w_BC_Zmin must be 0 or 1'
    end if
  
    return
  end
  
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  subroutine apply_BC_vel_upper_post
    !----*|--.---------.---------.---------.---------.---------.---------.-|--
    ! This subroutine is called after initializing the flow
    ! It sets the appropriate boundary conditions including ghost cell values
    !  on the velocity field in Fourier space
  
    integer i, k
  
    ! Now, apply boundary conditions to the top of the domain
    if (u_BC_Ymax == 0) then
      ! Dirichlet
      ! Start with zero
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu1(i, k, Nyp) = 0.d0
          cu1(i, k, Nyp + 1) = 0.d0
        end do
      end do
      ! Now, set only the mean
      if (rankZ == 0) then
        cu1(0, 0, Nyp) = u_BC_Ymax_c1
        cu1(0, 0, Nyp + 1) = u_BC_Ymax_c1
      end if
    else if (u_BC_Ymax == 1) then
      ! Neumann
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu1(i, k, Nyp) = cu1(i, k, Nyp - 1)
          cu1(i, k, Nyp + 1) = cu1(i, k, Nyp)
        end do
      end do
      ! Now, Apply BC to mean
      if (rankZ == 0) then
        cu1(0, 0, Nyp) = cu1(0, 0, Nyp - 1) + dy(Nyp) * u_BC_Ymax_c1
        cu1(0, 0, Nyp + 1) = cu1(0, 0, Nyp) + dy(Nyp) * u_BC_Ymax_c1
      end if
    else
      stop 'Error: u_BC_Zmax must be 0 or 1'
    end if
  
    if (w_BC_Ymax == 0) then
      ! Dirichlet
      ! Start with zero
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu3(i, k, Nyp) = 0.d0
          cu3(i, k, Nyp + 1) = 0.d0
        end do
      end do
      ! Now, set only the mean
      if (rankZ == 0) then
        cu3(0, 0, Nyp) = w_BC_Ymax_c1
        cu3(0, 0, Nyp + 1) = w_BC_Ymax_c1
      end if
      ! Ghost cell not used
      cu3(0, 0, Nyp + 1) = 0.d0
    else if (w_BC_Ymax == 1) then
      ! Neumann
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu3(i, k, Nyp) = cu3(i, k, Nyp - 1)
          cu3(i, k, Nyp + 1) = cu3(i, k, Nyp)
        end do
      end do
      ! Now, Apply BC to mean
      if (rankZ == 0) then
        if (f_type == 4)  w_BC_Ymax_c1_transient = w_BC_Ymax_c1 + amp_omega0 * sin(omega0 * (Ro_inv/delta * time - force_start) ) ! force_start is the phase to start _down-front_ forcing
        cu3(0, 0, Nyp) = cu3(0, 0, Nyp - 1) + dy(Nyp) * w_BC_Ymax_c1_transient
        cu3(0, 0, Nyp + 1) = cu3(0, 0, Nyp) + dy(Nyp) * w_BC_Ymax_c1_transient
      end if
    else
      stop 'Error: v_BC_Zmax must be 0 or 1'
    end if
  
    if (v_BC_Ymax == 0) then
      ! Dirichlet
      ! Set the vertical velocity at GYF(Nyp) (halfway between GY(Nyp) and GY(Nyp+1))
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu2(0, 0, Nyp + 1) = 2.d0 * v_BC_Ymax_c1 - cu2(0, 0, Nyp)
        end do
      end do
    else if (v_BC_Ymax == 1) then
      ! Neumann
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu2(i, k, Nyp) = cu2(i, k, Nyp - 1)
          cu2(i, k, Nyp + 1) = cu2(i, k, Nyp)
        end do
      end do
      ! Now, Apply BC to mean
      if (rankZ == 0) then
        cu2(0, 0, Nyp + 1) = cu2(0, 0, Nyp) + dy(Nyp) * v_BC_Ymax_c1
      end if
    else
      stop 'Error: w_BC_Zmax must be 0 or 1'
    end if
  
    return
  end
  
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  subroutine apply_BC_th_phys_lower
    !----*|--.---------.---------.---------.---------.---------.---------.-|--
    ! This subroutine is called after initializing the flow
    ! It sets the appropriate boundary conditions including ghost cell values
    !  on the velocity field in Physical space
  
    integer i, k, n
  
    do n = 1, N_th
  
      ! Now, apply the boundary conditions depending on the type specified
      if (th_BC_Ymin(n) == 0) then
        ! Dirichlet
        ! Start with zero
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            th(i, k, 0, n) = th_BC_Ymin_c1(n)
            th(i, k, 1, n) = th_BC_Ymin_c1(n)
          end do
        end do
      else if (th_BC_Ymin(n) == 1) then
        ! Neumann
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            th(i, k, 1, n) = th(i, k, 2, n) - dy(2) * th_BC_Ymin_c1(n)
            th(i, k, 0, n) = th(i, k, 1, n) - dy(1) * th_BC_Ymin_c1(n)
          end do
        end do
      else
        stop 'Error: TH_BC_Zmin must be 0 or 1'
      end if
  
    end do
  
    return
  end
  
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  subroutine apply_BC_th_phys_upper
    !----*|--.---------.---------.---------.---------.---------.---------.-|--
    ! This subroutine is called after initializing the flow
    ! It sets the appropriate boundary conditions including ghost cell values
    !  on the velocity field in Fourier space
  
    integer i, k, n
  
    do n = 1, N_th
  
      ! Now, apply boundary conditions to the top of the domain
      if (th_BC_Ymax(n) == 0) then
        ! Dirichlet
        ! Start with zero
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            th(i, k, Nyp, n) = th_BC_Ymax_c1(n)
            th(i, k, Nyp + 1, n) = th_BC_Ymax_c1(n)
          end do
        end do
      else if (th_BC_Ymax(n) == 1) then
        ! Neumann
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            th(i, k, Nyp, n) = th(i, k, Nyp - 1, n) + dy(Nyp) * th_BC_Ymax_c1(n)
            th(i, k, Nyp + 1, n) = th(i, k, Nyp, n) + dy(Nyp) * th_BC_Ymax_c1(n)
          end do
        end do
      else
        stop 'Error: TH_BC_Zmax must be 0 or 1'
      end if
  
    end do
  
    return
  end
  
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  subroutine apply_BC_vel_phys_lower
    !----*|--.---------.---------.---------.---------.---------.---------.-|--
    ! This subroutine is called after initializing the flow
    ! It sets the appropriate boundary conditions including ghost cell values
    !  on the velocity field in Physical space
  
    integer i, k
  
    ! Now, apply the boundary conditions depending on the type specified
    if (u_BC_Ymin == 0) then
      ! Dirichlet
      ! Start with zero
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          u1(i, k, 0) = u_BC_Ymin_c1
          u1(i, k, 1) = u_BC_Ymin_c1
        end do
      end do
    else if (u_BC_Ymin == 1) then
      ! Neumann
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          u1(i, k, 1) = u1(i, k, 2) - dy(2) * u_BC_Ymin_c1
          u1(i, k, 0) = u1(i, k, 1) - dy(1) * u_BC_Ymin_c1
        end do
      end do
    else
      stop 'Error: u_BC_Zmin must be 0 or 1'
    end if
  
    if (w_BC_Ymin == 0) then
      ! Dirichlet
      ! Start with zero
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          u3(i, k, 0) = w_BC_Ymin_c1
          u3(i, k, 1) = w_BC_Ymin_c1
        end do
      end do
    else if (w_BC_Ymin == 1) then
      ! Neumann
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          u3(i, k, 1) = u3(i, k, 2) - dy(2) * w_BC_Ymin_c1
          u3(i, k, 0) = u3(i, k, 1) - dy(1) * w_BC_Ymin_c1
        end do
      end do
    else
      stop 'Error: v_BC_Zmin must be 0 or 1'
    end if
  
    if (v_BC_Ymin == 0) then
      ! Dirichlet
      ! Set the vertical velocity at GYF(1) (halfway between GY(2) and GY(1))
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          u2(i, k, 1) = 2.d0 * v_BC_Ymin_c1 - u2(i, k, 2)
          u2(i, k, 0) = u2(i, k, 1)
        end do
      end do
    else if (v_BC_Ymin == 1) then
      ! Neumann
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          u2(i, k, 1) = u2(i, k, 2) - dyf(1) * v_BC_Ymin_c1
          u2(i, k, 0) = u2(i, k, 1) - dyf(1) * v_BC_Ymin_c1
        end do
      end do
    else
      stop 'Error: w_BC_Zmin must be 0 1'
    end if
  
    return
  end
  
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  subroutine apply_BC_vel_phys_upper
    !----*|--.---------.---------.---------.---------.---------.---------.-|--
    ! This subroutine is called after initializing the flow
    ! It sets the appropriate boundary conditions including ghost cell values
    !  on the velocity field in Fourier space
  
    integer i, k
  
    ! Now, apply boundary conditions to the top of the domain
    if (u_BC_Ymax == 0) then
      ! Dirichlet
      ! Start with zero
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          u1(i, k, Nyp) = u_BC_Ymax_c1
          u1(i, k, Nyp + 1) = u_BC_Ymax_c1
        end do
      end do
    else if (u_BC_Ymax == 1) then
      ! Neumann
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          u1(i, k, Nyp) = u1(i, k, Nyp - 1) + dy(Nyp) * u_BC_Ymax_c1
          u1(i, k, Nyp + 1) = u1(i, k, Nyp) + dy(Nyp) * u_BC_Ymax_c1
        end do
      end do
    else
      stop 'Error: u_BC_Zmax must be 0 or 1'
    end if
  
    if (w_BC_Ymax == 0) then
      ! Dirichlet
      ! Start with zero
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          u3(i, k, Nyp) = w_BC_Ymax_c1
          u3(i, k, Nyp + 1) = w_BC_Ymax_c1
        end do
      end do
    else if (w_BC_Ymax == 1) then
      ! Neumann
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          if (f_type == 4)  w_BC_Ymax_c1_transient = w_BC_Ymax_c1 + amp_omega0 * sin(omega0 * (Ro_inv/delta * time - force_start) ) ! force_start is the phase to start _down-front_ forcing
          u3(i, k, Nyp) = u3(i, k, Nyp - 1) + dy(Nyp) * w_BC_Ymax_c1_transient
          u3(i, k, Nyp + 1) = u3(i, k, Nyp) + dy(Nyp) * w_BC_Ymax_c1_transient
        end do
      end do
    else
      stop 'Error: v_BC_Zmax must be 0 or 1'
    end if
  
    if (v_BC_Ymax == 0) then
      ! Dirichlet
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          u2(i, k, Nyp + 1) = v_BC_Ymax_c1
          u2(i, k, Nyp) = v_BC_Ymax_c1
        end do
      end do
    else if (v_BC_Ymax == 1) then
      ! Neumann
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          u2(i, k, Nyp + 1) = u2(i, k, Nyp) + dy(Nyp) * v_BC_Ymax_c1
        end do
      end do
    else
      stop 'Error: w_BC_Zmax must be 0 1'
    end if
  
    return
  end
  
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  subroutine apply_BC_th_lower_post
    !----*|--.---------.---------.---------.---------.---------.---------.-|--
    ! This subroutine is called after initializing the flow
    ! It sets the appropriate boundary conditions including ghost cell values
    !  on the scalar fields in Fourier space
  
    integer i, k, n
  
    do n = 1, N_th
  
      ! Now, apply the boundary conditions depending on the type specified
      if (th_BC_Ymin(n) == 0) then
        ! Dirichlet
        ! Start with zero
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cth(i, k, 0, n) = 0.d0
            cth(i, k, 1, n) = 0.d0
          end do
        end do
        ! Now, set only the mean
        if (rankZ == 0) then
          cth(0, 0, 1, n) = th_BC_Ymin_c1(n)
          cth(0, 0, 0, n) = th_BC_Ymin_c1(n)
        end if
      else if (th_BC_Ymin(n) == 1) then
        ! Neumann
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cth(i, k, 1, n) = cth(i, k, 2, n)
            cth(i, k, 0, n) = cth(i, k, 1, n)
          end do
        end do
        ! Now, Apply BC to mean
        if (rankZ == 0) then
          cth(0, 0, 1, n) = cth(0, 0, 2, n) - dy(2) * th_BC_Ymin_c1(n)
          cth(0, 0, 0, n) = cth(0, 0, 1, n) - dy(1) * th_BC_Ymin_c1(n)
        end if
      else
        stop 'Error: TH_BC_Zmin must be 0, or 1'
      end if
  
    end do
  
    return
  end
  
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  subroutine apply_BC_th_upper_post
    !----*|--.---------.---------.---------.---------.---------.---------.-|--
    ! This subroutine is called after initializing the flow
    ! It sets the appropriate boundary conditions including ghost cell values
    !  on the scalar fields in Fourier space
  
    integer i, k, n
  
    do n = 1, N_th
  
      ! Now, apply boundary conditions to the top of the domain
      if (th_BC_Ymax(n) == 0) then
        ! Dirichlet
        ! Start with zero
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cth(i, k, Nyp, n) = 0.d0
            cth(i, k, Nyp + 1, n) = 0.d0
          end do
        end do
        ! Now, set only the mean
        if (rankZ == 0) then
          cth(0, 0, Nyp, n) = th_BC_Ymax_c1(n)
          cth(0, 0, Nyp + 1, n) = th_BC_Ymax_c1(n)
        end if
      else if (th_BC_Ymax(n) == 1) then
        ! Neumann
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cth(i, k, Nyp, n) = cth(i, k, Nyp - 1, n)
            cth(i, k, Nyp + 1, n) = cth(i, k, Nyp, n)
          end do
        end do
        ! Now, Apply B! to mean
        if (rankZ == 0) then
          cth(0, 0, Nyp, n) = cth(0, 0, Nyp - 1, n) + dy(Nyp) * th_BC_Ymax_c1(n)
          cth(0, 0, Nyp + 1, n) = cth(0, 0, Nyp, n) + dy(Nyp) * th_BC_Ymax_c1(n)
        end if
      else
        stop 'Error: TH_BC_Zmax must be 0 or 1'
      end if
  
    end do
  
    return
  end
  
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine apply_BC_vel_mpi_post
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! This subroutine applies the boundary conditions for the Poisson Eq.
    ! Note, MATL, MATD, etc. are dimensioned in header
  
    ! Apply Boundary conditions to velocity field
    if (rankY == 0) then
      call apply_BC_vel_lower_post
    end if
    if (rankY == NprocY - 1) then
      call apply_BC_vel_upper_post
    end if
  
    return
  end
  
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine apply_BC_th_mpi_post
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! This subroutine applies the boundary conditions for the Poisson Eq.
    ! Note, MATL, MATD, etc. are dimensioned in header
  
    ! Apply Boundary conditions to scalar field
    if (rankY == 0) then
      call apply_BC_th_lower_post
    end if
    if (rankY == NprocY - 1) then
      call apply_BC_th_upper_post
    end if
  
    return
  end
  
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine apply_BC_vel_phys_mpi
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! This subroutine applies the boundary conditions for the Poisson Eq.
    ! Note, MATL, MATD, etc. are dimensioned in header
  
  
    ! Apply Boundary conditions to velocity field
    if (rankY == 0) then
      call apply_BC_vel_phys_lower
    end if
    if (rankY == NprocY - 1) then
      call apply_BC_vel_phys_upper
    end if
  
    return
  end

end module boundary