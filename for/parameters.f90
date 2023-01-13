! Define parameters used throughout diablo
! Contains subroutines to read input files and set parameters
! Option to use toml-f to read in the inputs as a toml file

module parameters
  implicit none
  save

  ! current version - update if code is edited to invalidate old input files
  character(len=4) :: current_version = "3.5"

  ! Specify data-types
  integer, parameter :: single_kind = kind(0.0)
  integer, parameter :: double_kind = kind(0.d0)
  integer, parameter :: rkind = double_kind

  ! Grid
  ! (We hardwire these into the code so that the compiler may perform
  !  optimizations based on the grid size at compile time).
  include 'grid_def' ! integer, parameter :: Nx, Ny, Nz, N_th

  ! Parameters set in inputs
  character(len=4) :: version
  ! scheme
  character(len=35) :: flavor
  logical :: use_LES
  integer :: time_ad_meth ! must be 1
  integer :: les_model_type
  real(rkind) :: beta
  ! physical
  real(rkind) :: Lx, Ly, Lz
  real(rkind) :: Re, nu 
  real(rkind) :: nu_v_scale
  real(rkind) :: Ro, Ro_inv 
  real(rkind) :: delta
  real(rkind) :: grav_x, grav_y, grav_z
  ! timestepping
  real(rkind) :: wall_time_limit, time_limit
  real(rkind) :: delta_t
  real(rkind) :: max_dt
  logical :: variable_dt
  real(rkind) :: CFL
  integer :: update_dt
  ! output
  integer :: verbosity
  real(rkind) :: save_flow_dt, save_stats_dt, save_movie_dt
  real(rkind) :: XcMovie, YcMovie, ZcMovie
  ! initial conditions 
  logical :: create_new_flow, reset_time
  integer :: IC_Type
  real(rkind) :: kick 
  logical :: physical_noise
  ! forcing
  integer :: f_type
  real(rkind) :: ubulk0, px0, omega0, amp_omega0, force_start
  ! velocity bcs
  integer :: u_BC_Xmin, v_BC_Xmin, w_BC_Xmin
  integer :: u_BC_Xmax, v_BC_Xmax, w_BC_Xmax
  integer :: u_BC_Ymin, v_BC_Ymin, w_BC_Ymin
  integer :: u_BC_Ymax, v_BC_Ymax, w_BC_Ymax
  integer :: u_BC_Zmin, v_BC_Zmin, w_BC_Zmin
  integer :: u_BC_Zmax, v_BC_Zmax, w_BC_Zmax
  real(rkind) :: u_BC_Xmin_c1, v_BC_Xmin_c1, w_BC_Xmin_c1
  real(rkind) :: u_BC_Ymin_c1, v_BC_Ymin_c1, w_BC_Ymin_c1
  real(rkind) :: u_BC_Zmin_c1, v_BC_Zmin_c1, w_BC_Zmin_c1
  real(rkind) :: u_BC_Xmax_c1, v_BC_Xmax_c1, w_BC_Xmax_c1
  real(rkind) :: u_BC_Ymax_c1, v_BC_Ymax_c1, w_BC_Ymax_c1
  real(rkind) :: u_BC_Zmax_c1, v_BC_Zmax_c1, w_BC_Zmax_c1
  ! scalars 
  logical :: create_new_th(1:N_th)
  logical :: filter_th(1:N_th)
  integer :: filter_int(1:N_th)
  real(rkind) :: Ri(1:N_th), Pr(1:N_th)
  integer :: th_BC_Xmin(1:N_th), th_BC_Ymin(1:N_th), th_BC_Zmin(1:N_th)
  integer :: th_BC_Xmax(1:N_th), th_BC_Zmax(1:N_th), th_BC_Ymax(1:N_th)
  real(rkind) :: th_BC_Xmin_c1(1:N_th), th_BC_Ymin_c1(1:N_th), th_BC_Zmin_c1(1:N_th)
  real(rkind) :: th_BC_Xmax_c1(1:N_th), th_BC_Ymax_c1(1:N_th), th_BC_Zmax_c1(1:N_th)

  ! parameters defined in set_parameters
  logical :: homogeneousX
  real(rkind) :: w_BC_Ymax_c1_transient
  real(rkind) :: dTHdX(1:N_th), dTHdZ(1:N_th)

contains


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine read_inputs
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    ! read in input parameters
#ifdef TOML_INPUT
    ! read parameters from input.toml
    call read_input_toml
#else
    ! read in parameters from input.dat
    call read_input_dat
    call read_input_chan
#endif

    nu = 1.d0 / Re
    Ro_inv = 1.d0 / Ro
    max_dt = delta_t

  end


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine log_input_parameters
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! write out the input parameters
    ! should only be called by the rank 0 process

    integer n  
    
    write (*, '("Flavor: ", A35)') flavor
    write (*, '("Nx = ", I6)') Nx
    write (*, '("Ny = ", I6)') Ny
    write (*, '("Nz = ", I6)') Nz
    do n = 1, N_th
      write (*, '("Scalar Number: ", I2)') n
      write (*, '("  Richardson number = ", ES12.5)') Ri(n)
      write (*, '("  Prandtl number    = ", ES12.5)') Pr(n)
    end do
    write (*, '("Use LES: " L1)') use_LES
    if (use_LES) then 
      write(*, '("LES model type: ", I2.1)') les_model_type
    endif
    write (*, '("Nu   = ", ES12.5)') nu
    write (*, '("Beta = ", ES12.5)') beta
    write (*, '("Wall time limit = ", ES12.5)') wall_time_limit
    write (*, '("Sim time limit  = ", ES12.5)') time_limit

  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine read_input_dat
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    integer n
    
    open (11, file='input.dat', form='formatted', status='old')

    ! Read input file.
    !   (Note - if you change the following section of code, update the
    !    current_version number to make obsolete previous input files !)

    read (11, *)
    read (11, *)
    read (11, *)
    read (11, *)
    read (11, *) flavor, version
    if (version /= current_version) stop 'Wrong input data format.'
    read (11, *)
    read (11, *) use_LES
    read (11, *)
    read (11, *) Re, beta, Lx, Ly, Lz
    read (11, *)
    read (11, *) nu_v_scale
    read (11, *)
    read (11, *) create_new_flow
    read (11, *)
    read (11, *) wall_time_limit, time_limit, delta_t, reset_time, &
      variable_dt, CFL, update_dt
    read (11, *)
    read (11, *) verbosity, save_flow_dt, save_stats_dt, save_movie_dt, XcMovie, YcMovie, ZcMovie
    read (11, *)
    ! Read in the parameters for the N_th scalars
    do n = 1, N_th
      read (11, *)
      read (11, *) create_new_th(n)
      read (11, *)
      read (11, *) filter_th(n), filter_int(n)
      read (11, *)
      read (11, *) Ri(n), Pr(n)
    end do
  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine read_input_chan
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    
    integer n
    

    ! Read in input parameters specific for channel flow case
    open (11, file='input_chan.dat', form='formatted', status='old')
    ! Read input file.

    read (11, *)
    read (11, *)
    read (11, *)
    read (11, *)
    read (11, *) version
    if (version /= current_version) &
      stop 'Wrong input data format input_chan'
    read (11, *)
    read (11, *) time_ad_meth
    read (11, *)
    read (11, *) les_model_type
    read (11, *)
    read (11, *) IC_Type, kick, physical_noise
    read (11, *)
    read (11, *) Ro
    read (11, *)
    read (11, *) delta
    read (11, *)
    read (11, *) grav_x, grav_y, grav_z
    read (11, *)
    read (11, *) f_type, ubulk0, px0, omega0, amp_omega0, force_start
    read (11, *)
    read (11, *)
    read (11, *) u_BC_Ymin, u_BC_Ymin_c1
    read (11, *)
    read (11, *) w_BC_Ymin, w_BC_Ymin_c1
    read (11, *)
    read (11, *) v_BC_Ymin, v_BC_Ymin_c1
    read (11, *)
    read (11, *) u_BC_Ymax, u_BC_Ymax_c1
    read (11, *)
    read (11, *) w_BC_Ymax, w_BC_Ymax_c1
    read (11, *)
    read (11, *) v_BC_Ymax, v_BC_Ymax_c1
    read (11, *)
    ! Read in boundary conditions and background gradients for the N_th scalars
    do n = 1, N_th
      read (11, *)
      read (11, *) th_BC_Ymin(n), th_BC_Ymin_c1(n)
      read (11, *)
      read (11, *) th_BC_Ymax(n), th_BC_Ymax_c1(n)
    end do

    


    return
  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine set_parameters
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    ! Compensate no-slip BC in the GS flow direction due to dTHdx
    !   AND also define dTHdx & dTHdz
    if (IC_Type == 4 .or. IC_Type == 5) then ! Infinite Front
      if (w_BC_Ymin == 1) then
        w_BC_Ymin_c1 = w_BC_Ymin_c1 - 1.d0
      end if
      if (w_BC_Ymax == 1) then
        w_BC_Ymax_c1 = w_BC_Ymax_c1 - 1.d0
      end if
      dTHdX(1) = Ro_inv / delta
      dTHdZ(1) = 0.d0

    else if (IC_Type == 6 .or. IC_Type == 7 .or. IC_Type == 8) then ! Finite Front
      if (w_BC_Ymin == 1) then
        w_BC_Ymin_c1 = w_BC_Ymin_c1 - 2.d0 * delta / Lx
      end if
      if (w_BC_Ymax == 1) then
        w_BC_Ymax_c1 = w_BC_Ymax_c1 - 2.d0 * delta / Lx
      end if
      dTHdX(1) = 2.d0 / Lx * Ro_inv
      dTHdZ(1) = 0.d0

    else 
      dTHdX(1) = 0.d0
      dTHdZ(1) = 0.d0

    end if

    w_BC_Ymax_c1_transient = w_BC_Ymax_c1 ! Mean of surface forcing (compensating for GS flow)


    ! Set the valid averaging directions depending on the IC
    if (IC_Type == 5 .or. IC_Type == 6 .or. IC_Type == 7 .or. IC_Type == 8) then
      homogeneousX = .false.
    else ! Infinite, homogeneous front (or other IC...)
      homogeneousX = .true. ! Assume the x-direction is a valid averaging dimension
    endif
  end

  ! only define toml input option if using the tomlf library
#ifdef TOML_INPUT 

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine read_input_toml
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    use tomlf 
    type(toml_table), allocatable :: table ! root table
    type(toml_error), allocatable :: error
    type(toml_table), pointer :: child ! subtables
    type(toml_array), pointer :: array ! array of scalar tables
    real(rkind) :: vector(3) ! for reading GRAV and MOVIE
    integer :: number_of_scalars
    integer :: n

    ! read in the root table
    call toml_load(table,"input.toml",error=error)
    if (allocated(error)) then
      write (*,'("Error reading input.toml: ", A)') error%message
      stop
    end if

    ! check the version number
    call get_string_from_table(table,"VERSION", version)
    if (version /= current_version) &
      stop 'Wrong input version'

    ! set the scheme parameters
    call get_value(table,"SCHEME",child)
    call get_string_from_table(child,"FLAVOUR",flavor)
    call get_bool_from_table(child,"USE_LES",use_LES)
    call get_int_from_table(child,"LES_MODEL_TYPE",les_model_type)
    call get_float_from_table(child,"BETA",beta)
    time_ad_meth = 1

    ! set the physical parameters
    call get_value(table,"PHYSICAL",child)
    call get_float_from_table(child,"LX",LX)
    call get_float_from_table(child,"LY",LY)
    call get_float_from_table(child,"LZ",LZ)
    call get_float_from_table(child,"RE",Re)
    call get_float_from_table(child,"NU_V_SCALE",nu_v_scale)
    call get_float_from_table(child,"RO",Ro)
    call get_float_from_table(child,"DELTA",delta)
    call get_floats_from_table(child,"GRAV",vector)
    grav_x = vector(1)
    grav_y = vector(2)
    grav_z = vector(3)

    ! set the timestepping parameters
    call get_value(table,"TIMESTEPPING",child)
    call get_float_from_table(child,"WALL_TIME_LIMIT",wall_time_limit)
    call get_float_from_table(child,"TIME_LIMIT",time_limit)
    call get_float_from_table(child,"DELTA_T",delta_t)
    call get_bool_from_table(child,"VARIABLE_DT",variable_dt)
    call get_float_from_table(child,"CFL",CFL)
    call get_int_from_table(child,"UPDATE_DT",update_dt)

    ! set output parameters
    call get_value(table,"OUTPUT",child)
    call get_int_from_table(child,"VERBOSITY",verbosity)
    call get_float_from_table(child,"SAVE_FLOW_DT",save_flow_dt)
    call get_float_from_table(child,"SAVE_STATS_DT",save_stats_dt)
    call get_float_from_table(child,"SAVE_MOVIE_DT",save_movie_dt)
    call get_floats_from_table(child,"MOVIE",vector)
    XcMovie = vector(1)
    YcMovie = vector(2)
    ZcMovie = vector(3)

    ! set initial conditions parameters
    call get_value(table,"INITIAL_CONDITIONS",child)
    call get_bool_from_table(child,"CREATE_NEW_FLOW",create_new_flow)
    call get_bool_from_table(child,"RESET_TIME",reset_time)
    call get_int_from_table(child,"IC_TYPE",IC_Type)
    call get_float_from_table(child,"KICK",kick)
    call get_bool_from_table(child,"PHYSICAL_NOISE",physical_noise)

    ! set forcing parameters
    call get_value(table,"FORCING",child)
    call get_int_from_table(child,"F_TYPE",f_type)
    call get_float_from_table(child,"UBULK0",ubulk0)
    call get_float_from_table(child,"PX0",px0)
    call get_float_from_table(child,"OMEGA0",omega0)
    call get_float_from_table(child,"AMP_OMEGA0",amp_omega0)
    call get_float_from_table(child,"FORCE_START",force_start)

    ! set velocity bc parameters
    call get_value(table,"VELOCITY_BCS",child)
    call get_int_from_table(child,"U_BC_YMIN",u_BC_Ymin)
    call get_float_from_table(child,"U_BC_YMIN_C1",u_BC_Ymin_c1)
    call get_int_from_table(child,"V_BC_YMIN",w_BC_Ymin)
    call get_float_from_table(child,"V_BC_YMIN_C1",w_BC_Ymin_c1)
    call get_int_from_table(child,"W_BC_YMIN",v_BC_Ymin)
    call get_float_from_table(child,"W_BC_YMIN_C1",v_BC_Ymin_c1)
    call get_int_from_table(child,"U_BC_YMAX",u_BC_Ymax)
    call get_float_from_table(child,"U_BC_YMAX_C1",u_BC_Ymax_c1)
    call get_int_from_table(child,"V_BC_YMAX",w_BC_Ymax)
    call get_float_from_table(child,"V_BC_YMAX_C1",w_BC_Ymax_c1)
    call get_int_from_table(child,"W_BC_YMAX",v_BC_Ymax)
    call get_float_from_table(child,"W_BC_YMAX_C1",v_BC_Ymax_c1)

    ! set scalar parameters (array of tables)
    call get_value(table, "SCALARS", array)
    number_of_scalars = len(array)
    if (number_of_scalars /= N_th) then 
      write (*,'("Error: ", I2.1, " scalars defined but ", &
        I2.1, " scalars found in input.toml")') N_th, number_of_scalars
      stop
    end if
    do n = 1, N_th
      call get_value(array,n,child)
      call get_bool_from_table(child,"CREATE_FLOW_TH",create_new_th(n))
      call get_bool_from_table(child,"FILTER_TH",filter_th(n))
      call get_int_from_table(child,"FILTER_INT",filter_int(n))
      call get_float_from_table(child,"RI",Ri(n))
      call get_float_from_table(child,"PR",Pr(n))
      call get_int_from_table(child,"TH_BC_YMIN",th_BC_Ymin(n))
      call get_float_from_table(child,"TH_BC_YMIN_C1",th_BC_Ymin_c1(n))
      call get_int_from_table(child,"TH_BC_YMAX",th_BC_Ymax(n))
      call get_float_from_table(child,"TH_BC_YMAX_C1",th_BC_Ymax_c1(n))
    end do

  contains
    ! define a series of wrappers for get_value that get check the status 
    ! and print error messages
    
    subroutine get_float_from_table(table,varname,float)
      type(toml_table), intent(inout) :: table 
      character(len=*), intent(in) :: varname
      real(rkind), intent(out) :: float
      integer :: stat
      
      call get_value(table,varname,float,stat=stat)
      call check_stat(stat,varname,"float")
    end

    subroutine get_floats_from_table(table,varname,floats)
      ! read an array of floats
      type(toml_table), intent(inout) :: table
      character(len=*), intent(in) :: varname
      real(rkind), intent(out) :: floats(:)
      type(toml_array), pointer :: arr
      integer :: stat
      integer :: i, arr_len
      character(len=30) :: varname_index
      
      call get_value(table,varname,arr,stat=stat)
      call check_stat(stat,varname,"array")
      arr_len = len(arr)
      if (arr_len /= size(floats)) then 
        write(*,'("Error: Read ", I2.1, "floats from ", A, &
          " expected " ,I2.1)') arr_len, varname, size(floats)
        stop 
      endif
      do i = 1,arr_len
        call get_value(arr, i, floats(i),stat=stat)
        if (stat /= 0) then 
          write (varname_index,'(A,"-",I2.1)') varname, i
          call check_stat(stat,varname_index,"float")
        endif
      end do
    end

    subroutine get_int_from_table(table,varname,int)
      type(toml_table), intent(inout) :: table 
      character(len=*), intent(in) :: varname
      integer, intent(out) :: int
      integer :: stat
      
      call get_value(table,varname,int,stat=stat)
      call check_stat(stat,varname,"integer")
    end

    subroutine get_bool_from_table(table,varname,bool)
      type(toml_table), intent(inout) :: table 
      character(len=*), intent(in) :: varname
      logical, intent(out) :: bool
      integer :: stat
      
      call get_value(table,varname,bool,stat=stat)
      call check_stat(stat,varname,"boolean")
    end

    subroutine get_string_from_table(table,varname,string)
      type(toml_table), intent(inout) :: table 
      character(len=*), intent(in) :: varname
      character(len=*), intent(out) :: string
      character(len=:), allocatable :: temp
      integer :: stat
      
      call get_value(table,varname,temp,stat=stat)
      call check_stat(stat,varname,"string")
      string = temp
    end
    
    subroutine check_stat(stat,varname,datatype)
      ! catch all errors because otherwise toml-f will silently convert some types
      ! also display useful error message for fixing the input file
      integer, intent(in) :: stat
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: datatype

      if (stat == toml_stat%success) then
        return
      else if (stat == toml_stat%duplicate_key) then
        write(*,'("Error: Duplicate key when reading", A, ". Check your input file")') varname
        stop 
      else if (stat == toml_stat%type_mismatch) then
        write(*,'("Error: Wrong type when reading ", A, ". Expected ", A)') varname, datatype
        stop 
      else if (stat == toml_stat%conversion_error) then
        write(*,'("Error: Error converting", A, " to ", A)') varname, datatype
        stop 
      else if (stat == toml_stat%fatal) then
        write(*,'("Error: Fatal error in toml-f when reading", A, " :(")') varname
        stop 
      else 
        write(*,'("Error: Undefined error in toml-f when reading", A, " :(")') varname
        stop 
      end if
    end
    
  end
  
#endif 
! TOML_INPUT


end module parameters
