! define and initialise flow variables
! define routines for fixing boundary conditions
module flow
  use parameters
  use domain
  use fft
  implicit none

  ! RK
  real(rkind) :: h_bar(3), beta_bar(3), zeta_bar(3)
 
  ! time
  real(rkind) :: time
  integer :: time_step, rk_step, previous_time_step
  real(rkind) :: dt, delta_t_next_event
  real(rkind) :: save_flow_time, save_stats_time, save_movie_time
  

  ! 3D arrays
  real(rkind), pointer, contiguous, dimension(:,:,:) :: u1,u2,u3,p,r1,r2,r3,f1,f2,f3,s1
  complex(rkind), pointer, contiguous, dimension(:,:,:) :: cu1,cu2,cu3,cp,cr1,cr2,cr3, &
                                                           cf1,cf2,cf3,cs1
  ! 4D arrays
  real(rkind), pointer, contiguous, dimension(:,:,:,:) :: th,fth,rth
  complex(rkind), pointer, contiguous, dimension(:,:,:,:) :: cth,cfth,crth

  ! LES !
  real(rkind), pointer, contiguous, dimension(:,:,:) :: Sij1,Sij2,Sij3,Sij4,Sij5,Sij6,temp
  complex(rkind), pointer, contiguous, dimension(:,:,:) :: CSij1,CSij2,CSij3,CSij4,CSij5,CSij6,ctemp
  ! AMD !
  real(rkind), pointer, contiguous, dimension(:,:,:) :: temp_th, s1_th
  complex(rkind), pointer, contiguous, dimension(:,:,:) :: Ctemp_th, cs1_th
  real(rkind), pointer, contiguous, dimension(:,:,:) :: Oij4,Oij5,Oij6
  complex(rkind), pointer, contiguous, dimension(:,:,:) :: COij4,COij5,COij6
  real(rkind), pointer, contiguous, dimension(:,:,:) :: SIij1,SIij2,SIij3,SIij4,SIij5,SIij6
  complex(rkind), pointer, contiguous, dimension(:,:,:) :: CSIij1,CSIij2,CSIij3,CSIij4,CSIij5,CSIij6
  real(rkind), pointer, contiguous, dimension(:,:,:) :: du1dx,du2dy,du3dz, &
                                                           du2dx,du3dx, &
                                                           du1dz,du2dz, &
                                                           du1dy,du3dy
  complex(rkind), pointer, contiguous, dimension(:,:,:) :: Cdu1dx,Cdu2dy,Cdu3dz, &
                                                           Cdu2dx,Cdu3dx, &
                                                           Cdu1dz,Cdu2dz, &
                                                           Cdu1dy,Cdu3dy
  real(rkind), pointer, contiguous, dimension(:,:,:,:) :: dthetadx,dthetady,dthetadz
  complex(rkind), pointer, contiguous, dimension(:,:,:,:) :: Cdthetadx,Cdthetady,Cdthetadz

  
  logical flag_save_LES
  real(rkind) :: nu_t(0:Nx+1, 0:Nzp+1, 0:Nyp+1) = 0.d0
  real(rkind) :: kappa_t(0:Nx+1, 0:Nzp+1, 0:Nyp+1, 1:N_th) = 0.d0
  integer j1, j2


contains

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine init_flow
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    call allocate_all
    if (use_LES) call allocate_all_les

    ! Initialize FFT package (includes defining the wavenumber vectors).
    call init_fft(u1,cu1)

    ! Initialize time
    if (reset_time .or. create_new_flow) then
      save_flow_time = save_flow_dt
      save_stats_time = save_stats_dt
      save_movie_time = save_movie_dt
      time_step = 0
      time = 0    
    end if
    previous_time_step = time_step

    ! Initialize RKW3 parameters.
    h_bar(1) = delta_t * (8.d0 / 15.d0)
    h_bar(2) = delta_t * (2.d0 / 15.d0)
    h_bar(3) = delta_t * (5.d0 / 15.d0)
    beta_bar(1) = 1.d0
    beta_bar(2) = 25.d0 / 8.d0
    beta_bar(3) = 9.d0 / 4.d0
    zeta_bar(1) = 0.d0
    zeta_bar(2) = -17.d0 / 8.d0
    zeta_bar(3) = -5.d0 / 4.d0

  end 

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine allocate_all
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    call alloc_array3D(u1,cu1)
    call alloc_array3D(u2,cu2)
    call alloc_array3D(u3,cu3)
    call alloc_array3D(p ,cp)
    call alloc_array3D(r1,cr1)
    call alloc_array3D(r2,cr2)
    call alloc_array3D(r3,cr3)
    call alloc_array3D(f1,cf1)
    call alloc_array3D(f2,cf2)
    call alloc_array3D(f3,cf3)
    call alloc_array3D(s1,cs1)

    call alloc_array4D(th, cth)   ! Not using the same memory!
    call alloc_array4D(fth,cfth)
    call alloc_array4D(rth,crth)

  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine allocate_all_les
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    call alloc_array3D(Sij1,CSij1)
    call alloc_array3D(Sij2,CSij2)
    call alloc_array3D(Sij3,CSij3)
    call alloc_array3D(Sij4,CSij4)
    call alloc_array3D(Sij5,CSij5)
    call alloc_array3D(Sij6,CSij6)
    call alloc_array3D(temp,ctemp)


    ! Need Extra AMD Arrays
    call alloc_array3D(temp_th,ctemp_th)
    call alloc_array3D(s1_th,cs1_th)

    call alloc_array3D(Oij4,COij4)
    call alloc_array3D(Oij5,COij5)
    call alloc_array3D(Oij6,COij6)

    call alloc_array3D(SIij1,CSIij1)
    call alloc_array3D(SIij2,CSIij2)
    call alloc_array3D(SIij3,CSIij3)
    call alloc_array3D(SIij4,CSIij4)
    call alloc_array3D(SIij5,CSIij5)
    call alloc_array3D(SIij6,CSIij6)

    ! Can we actually point some of these to the SIij ?
    call alloc_array3D(du1dx,Cdu1dx)
    call alloc_array3D(du2dy,Cdu2dy)
    call alloc_array3D(du3dz,Cdu3dz)

    call alloc_array3D(du2dx,Cdu2dx)
    call alloc_array3D(du3dx,Cdu3dx)

    call alloc_array3D(du1dz,Cdu1dz)
    call alloc_array3D(du2dz,Cdu2dz)

    call alloc_array3D(du1dy,Cdu1dy)
    call alloc_array3D(du3dy,Cdu3dy)

    call alloc_array4D(dthetadx,Cdthetadx)   ! Not using the same memory!
    call alloc_array4D(dthetady,Cdthetady)
    call alloc_array4D(dthetadz,Cdthetadz)

  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine deallocate_all
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    ! call fftw_free(p) ! at pointer p...

  end

end module flow
