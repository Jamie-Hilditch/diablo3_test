! user defined forcing
! included in advance.f90

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine user_rhs_chan_physical
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

  ! Here, you can add terms to the right hand side
  ! of the momentum and scalar equations.
  ! The right hand side forcing arrays, CF1, CF2, CF3, CFTH
  ! are in Fourier space.  The velocity and scalars are available
  ! in physical space.
  ! S1 is available as a working variable

  ! integer i, j, k, n
  !
  ! real(rkind) alpha
  !
  ! call fft_xz_to_physical(cf3, f3)
  !
  ! do j = 1, Nyp
  !   do k = 0, Nzp - 1
  !     do i = 0, Nxm1
  !       f3(i, k, j) = f3(i, k, j) + 0.001 * exp(-((gx(i) - Lx/2)**2 + (gyf(j) - Ly/2)**2)/(2.d0*0.01**2))*sin(0.7*time);
  !     end do
  !   end do
  ! end do
  !
  ! call fft_xz_to_fourier(f3, cf3)

  ! call slip_vel

  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine user_rhs_chan_fourier
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

  ! Here, you can add terms to the right hand side
  ! of the momentum and scalar equations.
  ! The right hand side forcing arrays, CF1, CF2, CF3, CFTH
  ! are in Fourier space.  The velocity and scalars are available
  ! in Fourier space.
  ! S1 is available as a working variable

  integer i, j, k, n

  ! Advection owing to thermal wind
  if ((flavor == 'Front') .and. (Ro_inv /= 0.d0)) then
    do n = 1, N_th
      ! Loop over all scalars

      ! Add thermal wind advection to the momentum equations
      do j = jstart, jend
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cf1(i, k, j) = cf1(i, k, j) &
                           - (dTHdX(n) * (gyf(j) - 0.5d0*Ly) * delta / Ro_inv) &
                           * cikz(k) * cu1(i, k, j) &
                           - (-1.d0 * dTHdZ(n) * (gyf(j) - 0.5d0*Ly) * delta / Ro_inv) &
                           * cikx(i) * cu1(i, k, j) &
                           - (-1.d0 * dTHdZ(n) * delta / Ro_inv) &
                           * 0.5d0 * (cu2(i, k, j) + cu2(i, k, j + 1))
            cf3(i, k, j) = cf3(i, k, j) &
                           - (dTHdX(n) * (gyf(j) - 0.5d0*Ly) * delta / Ro_inv) &
                           * cikz(k) * cu3(i, k, j) &
                           - (-1.d0 * dTHdZ(n) * (gyf(j) - 0.5d0*Ly) * delta / Ro_inv) &
                           * cikx(i) * cu3(i, k, j) &
                           - (dTHdX(n) * delta / Ro_inv) &
                           * 0.5d0 * (cu2(i, k, j) + cu2(i, k, j + 1))
          end do
        end do
      end do

      do j = 2, Nyp
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cf2(i, k, j) = cf2(i, k, j) &
                           - (dTHdX(n) * (gy(j) - 0.5d0*Ly) * delta / Ro_inv) &
                           * cikz(k) * cu2(i, k, j) &
                           - (-1.d0 * dTHdZ(n) * (gy(j) - 0.5d0*Ly) * delta / Ro_inv) &
                           * cikx(i) * cu2(i, k, j)
          end do
        end do
      end do

      ! Add advection by thermal wind to the scalar equations
      do j = jstart_th(n), jend_th(n)
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cfth(i, k, j, n) = cfth(i, k, j, n) &
                               - (delta / Ro_inv) * dTHdX(n) * (gyf(j) - 0.5d0*Ly) &
                               * cikz(k) * cth(i, k, j, n) &
                               - (delta / Ro_inv) * (-1.d0) * dTHdZ(n) * (gyf(j) - 0.5d0*Ly) &
                               * cikx(i) * cth(i, k, j, n)
          end do
        end do
      end do

      ! End do N_th
    end do

  end if

  ! Add sponge layer forcing
  ! do n = 1, N_th
  !   call sponge_th(n)
  ! end do
  ! call sponge_vel

  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine sponge_th(n)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
  ! specified background state for the temperature field
  ! The intention is to allow an open boundary

  integer i, j, k, n
  real(rkind) L_sponge, L_bottom
  real(rkind) sponge_amp

  ! The following variables will store the background state
  real(rkind) th_0(-1:Nyp + 1)

  real(rkind) ri_b(0:Nyp + 1)

  ! This variable will hold the forcing rate
  real(rkind) sponge_sigma(0:Nyp + 1)

  ! Set the amplitude of the sponge
  sponge_amp = 0.005d0
  ! Set the top of the sponge layer in physical units
  L_sponge = -120.d0
  ! Set the bottom of the computational domain in physical units
  L_bottom = -140.d0
  do j = 0, Nyp + 1
    ! Quadratic damping at lower wall
    if (gyf(j) < L_sponge) then
      sponge_sigma(j) = sponge_amp * ((L_sponge - gyf(j)) &
                                      / (L_sponge - L_bottom))**2.d0
    else
      sponge_sigma(j) = 0.d0
    end if
  end do

  ! Set the profile for relaxing the mean TH
  do j = 0, Nyp + 1
    th_0(j) = th_BC_Ymin_c1(n) * gyf(j)
  end do

  ! For MLI latmix
  if (n == 1) then
    th_0(0) = 0.d0
    do j = 1, Nyp + 1
      ri_b(j) = 20.d0
      th_0(j) = th_0(j - 1) &
                + dy(j) * ri_b(j) * (dTHdX(n))**2.d0 &
                * (delta / Ro_inv)**2.d0
    end do
  else
    do j = 0, Nyp + 1
      th_0(j) = 0.d0
    end do
  end if

  ! Add damping to R-K terms
  ! Damp the perturbations towards 0
  do k = 0, twoNkz
    do i = 0, Nxp - 1
      if ((rankZ /= 0) .or. (i /= 0) .or. (k /= 0)) then
        do j = jstart_th(n), jend_th(n)
          cfth(i, k, j, n) = cfth(i, k, j, n) &
                             - sponge_sigma(j) * (cth(i, k, j, n) - 0.)
        end do
      end if
    end do
  end do
  ! Damp the mean gradient towards TH_0
  do j = jstart_th(n), jend_th(n)
    cfth(0, 0, j, n) = cfth(0, 0, j, n) - sponge_sigma(j) &
                       * (cth(0, 0, j, n) - th_0(j))
  end do

  return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine slip_vel
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine adds advection by a slip velocity to some scalars

  integer i, j, k, n, j1_th(1:N_th), j2_th(1:N_th)
  real(rkind) w_s(0:Nyp + 1, 1:N_th)
  real(rkind), dimension(5) :: slip_def = (/ 0.d0, 0.d0, 0.00005d0, 0.0005d0, 0.005d0/)

  ! Set indices corresponding to start and end of GYF grid
  do n = 1, N_th
    if (rankY == NprocY - 1) then
      ! We are at the upper wall
      j1_th(n) = jstart_th(n)
      j2_th(n) = Nyp - 1
    else if (rankY == 0) then
      ! We are at the lower wall
      j1_th(n) = 2
      j2_th(n) = jend_th(n)
    else
      ! We are on a middle process
      j1_th(n) = jstart_th(n)
      j2_th(n) = jend_th(n)
    end if
  end do

  ! First, set the slip velocity
  do j = 0, Nyp + 1
    do n = 1, N_th
      w_s(j, n) = slip_def(n)
    end do
  end do

  if (rankY == NprocY - 1) then
    ! We are on a process at the top boundary
    ! Set the slip velocity to zero at GY(Nyp) (and ghost cells)
    do n = 1, N_th
      w_s(Nyp, n) = 0.d0
      w_s(Nyp + 1, n) = 0.d0
    end do
  else if (rankY == 0) then
    ! We are on a process at the bottom boundary
    ! Set the slip velocity to zero at GY(2) (and ghost cells)
    do n = 1, N_th
      w_s(0, n) = 0.d0
      w_s(1, n) = 0.d0
      w_s(2, n) = 0.d0
    end do
  end if

  do n = 1, N_th
    do j = j1_th(n), j2_th(n)
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          ! Central differencing
          !              S1(I,K,J)=
          !     &     ((TH(I,K,J+1,N)*W_S(J+1,N) + TH(I,K,J,N)*W_S(J+1,N)
          !     &    -TH(I,K,J,N)*W_S(J,N)-TH(I,K,J-1,N)*W_S(J,n))/(2.d0*DYF(J)))
          ! Second order Upwinding
          !              S1(I,K,J)=(W_S(J+1,N)*TH(I,K,J,N)
          !     &               -W_S(J,N)*(TH(I,K,J,N)+TH(I,K,J-1,N))/2.d0)
          !     &                 /(GYF(j)-GY(j))
          ! First order upwinding
          s1(i, k, j) = (w_s(j + 1, n) * th(i, k, j, n) &
                        - w_s(j, n) * th(i, k, j - 1, n)) &
                        / (gyf(j) - gyf(j - 1))

          !              S1(I,K,J)=0.5d0*(W_S(J+1,N)+W_S(J,N))
          !     &              *(TH(I,K,J,N)-TH(I,K,J-1,N))/(GYF(J)-GYF(J-1))
        end do
      end do
    end do
    call fft_xz_to_fourier(s1, cs1)
    do j = j1_th(n), j2_th(n)
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cfth(i, k, j, n) = cfth(i, k, j, n) - cs1(i, k, j)
        end do
      end do
    end do
  end do

  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine sponge_vel
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
  ! specified background state
  ! The intention is to allow an open boundary

  integer i, j, k

  real(rkind) L_sponge, L_bottom
  real(rkind) sponge_amp

  ! The following variables will store the background state
  real(rkind) u1_0(-1:Nyp + 1), u2_0(0:Nyp + 1), u3_0(-1:Nyp + 1)

  ! This variable will hold the forcing rate
  real(rkind) sponge_sigma(0:Nyp + 1)

  ! Set the amplitude of the sponge
  sponge_amp = 0.0001d0
  ! Set the top of the sponge layer in physical units
  L_sponge = -120.d0
  ! Set the bottom of the computational domain in physical units
  L_bottom = -140.d0
  do j = 0, Nyp + 1
    ! Quadratic damping at lower wall
    if (gyf(j) < L_sponge) then
      sponge_sigma(j) = sponge_amp * ((L_sponge - gyf(j)) &
                                      / (L_sponge - L_bottom))**2.d0
    else
      sponge_sigma(j) = 0.d0
    end if
  end do

  ! Set the background state
  ! Here, set the background to be geostrophic, with a linear temperature profile
  do j = 0, Nyp + 1
    u1_0(j) = 0.d0
    u3_0(j) = 0.d0
  end do
  do j = 0, Nyp + 1
    u2_0(j) = 0.d0
  end do

  ! Add damping function to explicit R-K
  do k = 0, twoNkz
    do i = 0, Nxp - 1 ! Nkx
      if ((i /= 0) .or. (k /= 0)) then
        do j = jstart, jend
          cf1(i, k, j) = cf1(i, k, j) - sponge_sigma(j) * (cu1(i, k, j) - 0.d0)
          cf3(i, k, j) = cf3(i, k, j) - sponge_sigma(j) * (cu3(i, k, j) - 0.d0)
        end do
        do j = 1, Nyp
          cf2(i, k, j) = cf2(i, k, j) - &
                        0.5 * (sponge_sigma(j) + sponge_sigma(j + 1)) * (cu2(i, k, j) - 0.d0)
        end do
      end if
    end do
  end do
  ! Damp mean flow
  do j = jstart, jend
    cf1(0, 0, j) = cf1(0, 0, j) - sponge_sigma(j) * (cu1(0, 0, j) - u1_0(j))
    cf3(0, 0, j) = cf3(0, 0, j) - sponge_sigma(j) * (cu3(0, 0, j) - u3_0(j))
  end do
  do j = 1, Nyp
    cf2(0, 0, j) = cf2(0, 0, j) - sponge_sigma(j) * (cu2(0, 0, j) - u2_0(j))
  end do

  return
end










