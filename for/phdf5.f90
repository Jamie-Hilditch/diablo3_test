!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Writing out Stats & Planes
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
module phdf5
  use parameters 
  use domain
  use fft 
  use flow
  use hdf5
  implicit none

contains
  
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine write_flow(final)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    logical, intent(in) :: final
    character(len=55) fname
    integer i, j, k, n
    logical save_pressure

    save_pressure = .false.
    if (final) then
      fname = 'end.h5'
      save_pressure = .true.
    else
      write (fname,'(A4, I0.6, A3)') 'out.', time_step, '.h5'
    end if

    call mpi_barrier(mpi_comm_world, ierror)
    call WriteHDF5(fname, save_pressure)

    return
  end subroutine write_flow

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine WriteHDF5_real(fname, gname, Diag)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! Writes a single real scalar to the specified file

    character(len=12), intent(in) :: fname
    character(len=20), intent(in) :: gname
    real(rkind), intent(in) :: Diag
    character(len=20) :: dname

    ! Identifiers
    integer(hid_t) :: file_id, dset_id
    integer(hid_t) :: filspace_id, memspace_id

    ! Identifiers
    integer(hid_t) :: gid, selspace_id
    integer(hid_t) :: plist_id_d
    integer(hsize_t), dimension(1) :: dimsf

    ! Dimensions in the memory and in the file
    integer(hsize_t), dimension(1) :: dims

    integer :: rHDF5 = 1
    integer(hsize_t), dimension(1) :: adims
    integer(hid_t) :: dspace_id, aid, tspace

    
    integer :: nsamp
    logical :: flage

    integer(hsize_t), dimension(1) :: count, offset
    integer(hsize_t), dimension(1) :: stride, block1, offset_m

    ! integer(hsize_t)  ::  my_dim
    integer :: Error, i, j

    ! We only need one process to write to file
    if (rank == 0) then

      ! Write one value at a time
      nsamp = 1
      dims = 1
      stride = 1
      block1 = 1
      count = 1
      dimsf = 1

      ! Initialize interface
      call h5open_f(Error)

      inquire (file=trim(fname), exist=flage)
      if (.not. flage) then
        ! Create the file
        call h5fcreate_f(trim(fname), h5f_acc_trunc_f, &
                        file_id, Error)
        call h5fclose_f(file_id, Error)
      end if

      ! Open the file to get information
      call h5fopen_f(trim(fname), h5f_acc_rdwr_f, &
                    file_id, Error)

      call h5screate_f(h5s_scalar_f, tspace, Error)

      ! Check to see if the group exists
      ! Open the right group or create if it does not exist
      adims = 1
      call h5lexists_f(file_id, "/"//trim(gname), flage, Error)
      if (.not. flage) then
        call h5gcreate_f(file_id, gname, gid, Error)

        call h5acreate_f(gid, 'SAMPLES', h5t_std_i32le, &
                        tspace, aid, Error)
        nsamp = 0;
        call h5awrite_f(aid, h5t_native_integer, nsamp, adims, Error)
        call h5aclose_f(aid, Error)
      else
        call h5gopen_f(file_id, "/"//trim(gname), gid, Error)
        call h5aopen_f(gid, 'SAMPLES', aid, Error)
        ! Read in the number of samples in the group
        call h5aread_f(aid, h5t_native_integer, nsamp, adims, Error)
        call h5aclose_f(aid, Error)
      end if

      ! Increase the sample index by 1
      nsamp = nsamp + 1

      ! The dataset name will be the sample number
      write (dname, '(1I0.4)') nsamp

      call h5screate_simple_f(rHDF5, dimsf, filspace_id, Error)
      call h5sselect_hyperslab_f(filspace_id, h5s_select_set_f, &
                                offset, count, Error, stride, block1)

      ! Write the new number of samples to the file
      call h5aopen_f(gid, 'SAMPLES', aid, Error)
      call h5awrite_f(aid, h5t_native_integer, nsamp, adims, Error)
      call h5aclose_f(aid, Error)

      call h5dcreate_f(gid, dname, h5t_ieee_f64le, filspace_id &
                      , dset_id, Error)
      call h5dwrite_f(dset_id, h5t_ieee_f64le, Diag, dims, Error)

      call h5dclose_f(dset_id, Error)
      call h5sclose_f(filspace_id, Error)
      call h5gclose_f(gid, Error)
      call h5fclose_f(file_id, Error)

      call h5close_f(Error)

    end if 

    ! Sync the cores as a precaution (probably not necessary)
    call mpi_barrier(mpi_comm_world, ierror)

    return
  end subroutine WriteHDF5_real

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine WriteStatH5_Y(fname, gname, Diag)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! Writes a vector of values in Y (Vertical)

    character(len=12), intent(in) :: fname
    real(rkind), intent(in) :: Diag(1:Nyp)
    character(len=20), intent(in) :: gname

    character(len=20) :: dname

    ! Identifiers
    integer(hid_t) :: file_id, dset_id
    integer(hid_t) :: filspace_id, memspace_id

    ! Identifiers
    integer(hid_t) :: gid, selspace_id
    integer(hid_t) :: plist_id_d

    ! Dimensions in the memory and in the file
    integer(hsize_t), dimension(1) :: dimsm, dimsf

    integer :: rHDF5 = 1
    integer(hsize_t), dimension(1)       :: adims
    integer(hid_t)                      :: aid, tspace

    
    integer nsamp
    logical flage

    integer(hsize_t), dimension(1) :: count, offset
    integer(hsize_t), dimension(1) :: stride, block1, offset_m

    ! integer(hsize_t)  ::  my_dim
    integer Error, i, j

    ! *********************
    ! START DEFINITION
    ! *********************

    dimsm = Nyp
    dimsf = (Nyp - 1) * NprocY + 1

    ! Stride and count for number of rows and columns in each dimension
    stride = 1
    count = 1

    ! Offset determined by the rank of a processor
    !  offset(1) = 0

    !  offset_m(1:2)=0
    if (rankY == 0) then
      block1 = Nyp
      offset = 0
      offset_m = 0
    else
      block1 = (Nyp - 1)
      offset = rankY * (Nyp - 1) + 1
      offset_m = 1
    end if

    ! *********************
    ! FINISH DEFINITION
    ! *********************

    ! Initialize interface
    call h5open_f(Error)

    ! Setup file access property list with parallel I/O access
    call h5pcreate_f(h5p_file_access_f, plist_id_d, Error)
    call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_y%mpi_val, &
                            mpi_info_null%mpi_val, Error)

    inquire (file=trim(fname), exist=flage)
    if (.not. flage) then
      ! Create the file collectively
      call h5fcreate_f(trim(fname), h5f_acc_trunc_f, &
                      file_id, Error, access_prp=plist_id_d)
      call h5fclose_f(file_id, Error)
    end if

    ! Create the file collectively
    call h5fopen_f(trim(fname), h5f_acc_rdwr_f, &
                  file_id, Error, access_prp=plist_id_d)
    call h5pclose_f(plist_id_d, Error)

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id_d, Error)
    call h5pset_dxpl_mpio_f(plist_id_d, h5fd_mpio_collective_f, &
                            Error)

    call h5screate_f(h5s_scalar_f, tspace, Error)

    ! Open the right group or create if it does not exist
    adims = 1
    call h5lexists_f(file_id, "/"//trim(gname), flage, Error)
    if (.not. flage) then
      call h5gcreate_f(file_id, gname, gid, Error)

      call h5acreate_f(gid, 'SAMPLES', h5t_std_i32le, &
                      tspace, aid, Error)
      nsamp = 0;
      call h5awrite_f(aid, h5t_native_integer, nsamp, adims, Error)
      call h5aclose_f(aid, Error)
    else
      call h5gopen_f(file_id, "/"//trim(gname), gid, Error)
      call h5aopen_f(gid, 'SAMPLES', aid, Error)
      call h5aread_f(aid, h5t_native_integer, nsamp, adims, Error)
      call h5aclose_f(aid, Error)
    end if

    nsamp = nsamp + 1

    write (dname, '(1I0.4)') nsamp

    call h5screate_simple_f(rHDF5, dimsf, filspace_id, Error)
    call h5screate_simple_f(rHDF5, dimsm, memspace_id, Error)

    call h5sselect_hyperslab_f(filspace_id, h5s_select_set_f, &
                              offset, count, Error, stride, block1)
    call h5sselect_hyperslab_f(memspace_id, h5s_select_set_f, &
                              offset_m, count, Error, stride, block1)

    call h5aopen_f(gid, 'SAMPLES', aid, Error)
    call h5awrite_f(aid, h5t_native_integer, nsamp, adims, Error)
    call h5aclose_f(aid, Error)

    call h5dcreate_f(gid, dname, h5t_ieee_f64le, &
                    filspace_id, dset_id, Error)
    ! Write the dataset collectively
    call h5dwrite_f(dset_id, h5t_native_double, &
                    Diag, &
                    dimsm, Error, file_space_id=filspace_id, &
                    mem_space_id=memspace_id, xfer_prp=plist_id_d)

    call h5acreate_f(dset_id, 'Time', h5t_ieee_f64le, tspace, &
                    aid, Error)
    call h5awrite_f(aid, h5t_ieee_f64le, time, adims, Error)
    call h5aclose_f(aid, Error)

    call h5dclose_f(dset_id, Error)

    call h5sclose_f(filspace_id, Error)
    call h5sclose_f(memspace_id, Error)
    call h5pclose_f(plist_id_d, Error)

    call h5gclose_f(gid, Error)
    call h5fclose_f(file_id, Error)
    call h5close_f(Error)

  end subroutine WriteStatH5_Y

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine WriteStatH5_X(fname, gname, Diag, NperProc)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! Writes a vector of values in X

    character(len=12) fname

    ! Dataset names
    character(len=20) :: gname, dname

    integer NperProc ! How many points is each procZ sending?

    ! Identifiers
    integer(hid_t) :: file_id, dset_id
    integer(hid_t) :: filspace_id, memspace_id

    ! Identifiers
    integer(hid_t) :: gid, selspace_id
    integer(hid_t) :: plist_id_d

    ! Dimensions in the memory and in the file
    integer(hsize_t), dimension(1) :: dimsm, dimsf

    integer :: rHDF5 = 1
    integer(hsize_t), dimension(1)       :: adims
    integer(hid_t)                      :: aid, tspace

    real(rkind) Diag(1:NperProc)
    integer nsamp
    logical flage

    integer(hsize_t), dimension(1) :: count, offset
    integer(hsize_t), dimension(1) :: stride, block1, offset_m

    ! integer(hsize_t)  ::  my_dim
    integer Error, i, j

    ! *********************
    ! START DEFINITION
    ! *********************

    dimsm = NperProc
    dimsf = NperProc * NprocZ

    ! Stride and count for number of rows and columns in each dimension
    stride = 1
    count = 1

    ! Offset determined by the rank of a processor
    !  offset(1) = 0
    !  offset_m(1:2)=0

    block1 = NperProc
    offset = rankZ * NperProc
    offset_m = 0


    ! *********************
    ! FINISH DEFINITION
    ! *********************

    ! Initialize interface
    call h5open_f(Error)

    ! Setup file access property list with parallel I/O access
    call h5pcreate_f(h5p_file_access_f, plist_id_d, Error)
    call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_z%mpi_val, &
                            mpi_info_null%mpi_val, Error)

    inquire (file=trim(fname), exist=flage)
    if (.not. flage) then
      ! Create the file collectively
      call h5fcreate_f(trim(fname), h5f_acc_trunc_f, &
                      file_id, Error, access_prp=plist_id_d)
      call h5fclose_f(file_id, Error)
    end if

    ! Create the file collectively
    call h5fopen_f(trim(fname), h5f_acc_rdwr_f, &
                  file_id, Error, access_prp=plist_id_d)
    call h5pclose_f(plist_id_d, Error)

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id_d, Error)
    call h5pset_dxpl_mpio_f(plist_id_d, h5fd_mpio_collective_f, &
                            Error)

    call h5screate_f(h5s_scalar_f, tspace, Error)

    ! Open the right group or create if it does not exist
    adims = 1
    call h5lexists_f(file_id, "/"//trim(gname), flage, Error)
    if (.not. flage) then
      call h5gcreate_f(file_id, gname, gid, Error)

      call h5acreate_f(gid, 'SAMPLES', h5t_std_i32le, &
                      tspace, aid, Error)
      nsamp = 0;
      call h5awrite_f(aid, h5t_native_integer, nsamp, adims, Error)
      call h5aclose_f(aid, Error)
    else
      call h5gopen_f(file_id, "/"//trim(gname), gid, Error)
      call h5aopen_f(gid, 'SAMPLES', aid, Error)
      call h5aread_f(aid, h5t_native_integer, nsamp, adims, Error)
      call h5aclose_f(aid, Error)
    end if

    nsamp = nsamp + 1

    write (dname, '(1I0.4)') nsamp

    call h5screate_simple_f(rHDF5, dimsf, filspace_id, Error)
    call h5screate_simple_f(rHDF5, dimsm, memspace_id, Error)

    call h5sselect_hyperslab_f(filspace_id, h5s_select_set_f, &
                              offset, count, Error, stride, block1)
    call h5sselect_hyperslab_f(memspace_id, h5s_select_set_f, &
                              offset_m, count, Error, stride, block1)

    call h5aopen_f(gid, 'SAMPLES', aid, Error)
    call h5awrite_f(aid, h5t_native_integer, nsamp, adims, Error)
    call h5aclose_f(aid, Error)

    call h5dcreate_f(gid, dname, h5t_ieee_f64le, &
                    filspace_id, dset_id, Error)
    ! Write the dataset collectively
    call h5dwrite_f(dset_id, h5t_native_double, &
                    Diag, &
                    dimsm, Error, file_space_id=filspace_id, &
                    mem_space_id=memspace_id, xfer_prp=plist_id_d)

    call h5acreate_f(dset_id, 'Time', h5t_ieee_f64le, tspace, &
                    aid, Error)
    call h5awrite_f(aid, h5t_ieee_f64le, time, adims, Error)
    call h5aclose_f(aid, Error)

    call h5dclose_f(dset_id, Error)

    call h5sclose_f(filspace_id, Error)
    call h5sclose_f(memspace_id, Error)
    call h5pclose_f(plist_id_d, Error)

    call h5gclose_f(gid, Error)
    call h5fclose_f(file_id, Error)
    call h5close_f(Error)

  end subroutine WriteStatH5_X

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine WriteHDF5_XYplane(fname, gname, var2d)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! Writes out entire X-Y plane

    character(len=35) fname

    ! Dataset names
    character(len=20) :: gname, dname

    ! Identifiers
    integer(hid_t) :: file_id, dset_id
    integer(hid_t) :: filspace_id, memspace_id

    ! Identifiers
    integer(hid_t) :: gid, selspace_id
    integer(hid_t) :: plist_id_d

    ! Dimensions in the memory and in the file
    integer(hsize_t), dimension(2) :: dimsm, dimsf

    integer :: rHDF5 = 2
    integer(hsize_t), dimension(1)       :: adims
    integer(hid_t)                      :: aid, tspace

    real(rkind) var2d(Nx, Nyp)
    integer nsamp
    logical flage

    integer(hsize_t), dimension(2) :: count, offset
    integer(hsize_t), dimension(2) :: stride, block1, offset_m

    ! integer(hsize_t)  ::  my_dim
    integer Error, i, j

    ! *********************
    ! START DEFINITION
    ! *********************

    dimsm(1:2) = (/Nx, Nyp/)
    dimsf(1:2) = (/Nx, (Nyp - 1) * NprocY + 1/)

    block1(1) = Nx

    ! Stride and count for number of rows and columns in each dimension
    stride = 1
    count = 1

    ! Offset determined by the rank of a processor
    offset(1) = 0

    offset_m(1:2) = 0
    if (rankY == 0) then
      block1(2) = Nyp
      offset(2) = 0
      offset_m(2) = 0
    else
      block1(2) = (Nyp - 1)
      offset(2) = rankY * (Nyp - 1) + 1
      offset_m(2) = 1
    end if

    ! *********************
    ! FINISH DEFINITION
    ! *********************

    ! Initialize interface
    call h5open_f(Error)

    ! Setup file access property list with parallel I/O access
    call h5pcreate_f(h5p_file_access_f, plist_id_d, Error)
    call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_y%mpi_val, &
                            mpi_info_null%mpi_val, Error)

    inquire (file=trim(fname), exist=flage)
    if (.not. flage) then
      ! Create the file collectively
      call h5fcreate_f(trim(fname), h5f_acc_trunc_f, &
                      file_id, Error, access_prp=plist_id_d)
      call h5fclose_f(file_id, Error)
    end if

    ! Create the file collectively
    call h5fopen_f(trim(fname), h5f_acc_rdwr_f, &
                  file_id, Error, access_prp=plist_id_d)
    call h5pclose_f(plist_id_d, Error)

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id_d, Error)
    call h5pset_dxpl_mpio_f(plist_id_d, h5fd_mpio_collective_f, &
                            Error)

    call h5screate_f(h5s_scalar_f, tspace, Error)

    ! Open the right group or create if it does not exist
    adims = 1
    call h5lexists_f(file_id, "/"//trim(gname), flage, Error)
    if (.not. flage) then
      call h5gcreate_f(file_id, gname, gid, Error)

      call h5acreate_f(gid, 'SAMPLES', h5t_std_i32le, &
                      tspace, aid, Error)
      nsamp = 0;
      call h5awrite_f(aid, h5t_native_integer, nsamp, adims, Error)
      call h5aclose_f(aid, Error)
    else
      call h5gopen_f(file_id, "/"//trim(gname), gid, Error)
      call h5aopen_f(gid, 'SAMPLES', aid, Error)
      call h5aread_f(aid, h5t_native_integer, nsamp, adims, Error)
      call h5aclose_f(aid, Error)
    end if

    nsamp = nsamp + 1

    write (dname, '(1I0.4)') nsamp

    call h5screate_simple_f(rHDF5, dimsf, filspace_id, Error)
    call h5screate_simple_f(rHDF5, dimsm, memspace_id, Error)

    call h5sselect_hyperslab_f(filspace_id, h5s_select_set_f, &
                              offset, count, Error, stride, block1)
    call h5sselect_hyperslab_f(memspace_id, h5s_select_set_f, &
                              offset_m, count, Error, stride, block1)

    call h5aopen_f(gid, 'SAMPLES', aid, Error)
    call h5awrite_f(aid, h5t_native_integer, nsamp, adims, Error)
    call h5aclose_f(aid, Error)

    call h5dcreate_f(gid, dname, h5t_ieee_f64le, &
                    filspace_id, dset_id, Error)
    ! Write the dataset collectively
    call h5dwrite_f(dset_id, h5t_native_double, &
                    var2d, &
                    dimsm, Error, file_space_id=filspace_id, &
                    mem_space_id=memspace_id, xfer_prp=plist_id_d)

    call h5acreate_f(dset_id, 'Time', h5t_ieee_f64le, tspace, &
                    aid, Error)
    call h5awrite_f(aid, h5t_ieee_f64le, time, adims, Error)
    call h5aclose_f(aid, Error)

    call h5dclose_f(dset_id, Error)

    call h5sclose_f(filspace_id, Error)
    call h5sclose_f(memspace_id, Error)
    call h5pclose_f(plist_id_d, Error)

    call h5gclose_f(gid, Error)
    call h5fclose_f(file_id, Error)
    call h5close_f(Error)

  end subroutine WriteHDF5_XYplane

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine WriteHDF5_XZplane(fname, gname, var2d)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! Writes out entire X-Z plane

    character(len=35) fname

    ! Dataset names
    character(len=20) :: gname, dname

    ! Identifiers
    integer(hid_t) :: file_id, dset_id
    integer(hid_t) :: filspace_id, memspace_id

    ! Identifiers
    integer(hid_t) :: gid, selspace_id
    integer(hid_t) :: plist_id_d

    ! Dimensions in the memory and in the file
    integer(hsize_t), dimension(2) :: dimsm, dimsf

    integer :: rHDF5 = 2
    integer(hsize_t), dimension(1)       :: adims
    integer(hid_t)                      :: aid, tspace

    real(rkind) var2d(Nx, Nzp)
    integer nsamp
    logical flage

    integer(hsize_t), dimension(2) :: count, offset
    integer(hsize_t), dimension(2) :: stride, block1, offset_m

    ! integer(hsize_t)  ::  my_dim
    integer Error, i, j

    ! *********************
    ! START DEFINITION
    ! *********************

    dimsm(1:2) = (/Nx, Nzp/)
    dimsf(1:2) = (/Nx, Nz/)

    block1(1) = Nx
    block1(2) = Nzp

    ! Stride and count for number of rows and columns in each dimension
    stride = 1
    count = 1

    ! Offset determined by the rank of a processor
    offset(1) = 0
    offset(2) = Nzp * rankZ

    offset_m(1:2) = 0

    ! *********************
    ! FINISH DEFINITION
    ! *********************

    ! Initialize interface
    call h5open_f(Error)

    ! Setup file access property list with parallel I/O access
    call h5pcreate_f(h5p_file_access_f, plist_id_d, Error)
    call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_z%mpi_val, &
                            mpi_info_null%mpi_val, Error)

    inquire (file=trim(fname), exist=flage)
    if (.not. flage) then
      ! Create the file collectively
      call h5fcreate_f(trim(fname), h5f_acc_trunc_f, &
                      file_id, Error, access_prp=plist_id_d)
      call h5fclose_f(file_id, Error)
    end if

    ! Create the file collectively
    call h5fopen_f(trim(fname), h5f_acc_rdwr_f, &
                  file_id, Error, access_prp=plist_id_d)
    call h5pclose_f(plist_id_d, Error)

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id_d, Error)
    call h5pset_dxpl_mpio_f(plist_id_d, h5fd_mpio_collective_f, &
                            Error)

    call h5screate_f(h5s_scalar_f, tspace, Error)

    ! Open the right group or create if it does not exist
    adims = 1
    call h5lexists_f(file_id, "/"//trim(gname), flage, Error)
    if (.not. flage) then
      call h5gcreate_f(file_id, gname, gid, Error)

      call h5acreate_f(gid, 'SAMPLES', h5t_std_i32le, &
                      tspace, aid, Error)
      nsamp = 0;
      call h5awrite_f(aid, h5t_native_integer, nsamp, adims, Error)
      call h5aclose_f(aid, Error)
    else
      call h5gopen_f(file_id, "/"//trim(gname), gid, Error)
      call h5aopen_f(gid, 'SAMPLES', aid, Error)
      call h5aread_f(aid, h5t_native_integer, nsamp, adims, Error)
      call h5aclose_f(aid, Error)
    end if

    nsamp = nsamp + 1

    write (dname, '(1I0.4)') nsamp

    call h5screate_simple_f(rHDF5, dimsf, filspace_id, Error)
    call h5screate_simple_f(rHDF5, dimsm, memspace_id, Error)

    call h5sselect_hyperslab_f(filspace_id, h5s_select_set_f, &
                              offset, count, Error, stride, block1)
    call h5sselect_hyperslab_f(memspace_id, h5s_select_set_f, &
                              offset_m, count, Error, stride, block1)

    call h5aopen_f(gid, 'SAMPLES', aid, Error)
    call h5awrite_f(aid, h5t_native_integer, nsamp, adims, Error)
    call h5aclose_f(aid, Error)

    call h5dcreate_f(gid, dname, h5t_ieee_f64le, &
                    filspace_id, dset_id, Error)
    ! Write the dataset collectively
    call h5dwrite_f(dset_id, h5t_native_double, &
                    var2d, &
                    dimsm, Error, file_space_id=filspace_id, &
                    mem_space_id=memspace_id, xfer_prp=plist_id_d)

    call h5acreate_f(dset_id, 'Time', h5t_ieee_f64le, tspace, &
                    aid, Error)
    call h5awrite_f(aid, h5t_ieee_f64le, time, adims, Error)
    call h5aclose_f(aid, Error)

    call h5dclose_f(dset_id, Error)

    call h5sclose_f(filspace_id, Error)
    call h5sclose_f(memspace_id, Error)
    call h5pclose_f(plist_id_d, Error)

    call h5gclose_f(gid, Error)
    call h5fclose_f(file_id, Error)
    call h5close_f(Error)

  end subroutine WriteHDF5_XZplane

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine WriteHDF5_ZYplane(fname, gname, var2d)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! Writes out entire Z-Y plane

    character(len=35) fname

    ! Dataset names
    character(len=20) :: gname, dname

    ! Identifiers
    integer(hid_t) :: file_id, dset_id
    integer(hid_t) :: filspace_id, memspace_id

    ! Identifiers
    integer(hid_t) :: gid, selspace_id
    integer(hid_t) :: plist_id_d

    ! Dimensions in the memory and in the file
    integer(hsize_t), dimension(2) :: dimsm, dimsf

    integer :: rHDF5 = 2
    integer(hsize_t), dimension(1)       :: adims
    integer(hid_t)                      :: aid, tspace

    real(rkind) var2d(Nzp, Nyp)
    integer nsamp
    logical flage

    integer(hsize_t), dimension(2) :: count, offset
    integer(hsize_t), dimension(2) :: stride, block1, offset_m

    ! integer(hsize_t)  ::  my_dim
    integer Error, i, j

    ! *********************
    ! START DEFINITION
    ! *********************

    dimsm(1:2) = (/Nzp, Nyp/)
    dimsf(1:2) = (/Nz, (Nyp - 1) * NprocY + 1/)

    block1(1) = Nzp

    ! Stride and count for number of rows and columns in each dimension
    stride = 1
    count = 1

    ! Offset determined by the rank of a processor
    offset(1) = Nzp * rankZ
    offset_m(1:2) = 0

    if (rankY == 0) then
      block1(2) = Nyp
      offset(2) = 0
      offset_m(2) = 0
    else
      block1(2) = (Nyp - 1)
      offset(2) = rankY * (Nyp - 1) + 1
      offset_m(2) = 1
    end if

    ! *********************
    ! FINISH DEFINITION
    ! *********************

    ! Initialize interface
    call h5open_f(Error)

    ! Setup file access property list with parallel I/O access
    call h5pcreate_f(h5p_file_access_f, plist_id_d, Error)
    call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_world%mpi_val, &
                            mpi_info_null%mpi_val, Error)

    inquire (file=trim(fname), exist=flage)
    if (.not. flage) then
      ! Create the file collectively
      call h5fcreate_f(trim(fname), h5f_acc_trunc_f, &
                      file_id, Error, access_prp=plist_id_d)
      call h5fclose_f(file_id, Error)
    end if

    ! Create the file collectively
    call h5fopen_f(trim(fname), h5f_acc_rdwr_f, &
                  file_id, Error, access_prp=plist_id_d)
    call h5pclose_f(plist_id_d, Error)

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id_d, Error)
    call h5pset_dxpl_mpio_f(plist_id_d, h5fd_mpio_collective_f, &
                            Error)

    call h5screate_f(h5s_scalar_f, tspace, Error)

    ! Open the right group or create if it does not exist
    adims = 1
    call h5lexists_f(file_id, "/"//trim(gname), flage, Error)
    if (.not. flage) then
      call h5gcreate_f(file_id, gname, gid, Error)

      call h5acreate_f(gid, 'SAMPLES', h5t_std_i32le, &
                      tspace, aid, Error)
      nsamp = 0;
      call h5awrite_f(aid, h5t_native_integer, nsamp, adims, Error)
      call h5aclose_f(aid, Error)
    else
      call h5gopen_f(file_id, "/"//trim(gname), gid, Error)
      call h5aopen_f(gid, 'SAMPLES', aid, Error)
      call h5aread_f(aid, h5t_native_integer, nsamp, adims, Error)
      call h5aclose_f(aid, Error)
    end if

    nsamp = nsamp + 1

    write (dname, '(1I0.4)') nsamp

    call h5screate_simple_f(rHDF5, dimsf, filspace_id, Error)
    call h5screate_simple_f(rHDF5, dimsm, memspace_id, Error)

    call h5sselect_hyperslab_f(filspace_id, h5s_select_set_f, &
                              offset, count, Error, stride, block1)
    call h5sselect_hyperslab_f(memspace_id, h5s_select_set_f, &
                              offset_m, count, Error, stride, block1)

    call h5aopen_f(gid, 'SAMPLES', aid, Error)
    call h5awrite_f(aid, h5t_native_integer, nsamp, adims, Error)
    call h5aclose_f(aid, Error)

    call h5dcreate_f(gid, dname, h5t_ieee_f64le, &
                    filspace_id, dset_id, Error)
    ! Write the dataset collectively
    call h5dwrite_f(dset_id, h5t_native_double, &
                    var2d, &
                    dimsm, Error, file_space_id=filspace_id, &
                    mem_space_id=memspace_id, xfer_prp=plist_id_d)

    call h5acreate_f(dset_id, 'Time', h5t_ieee_f64le, tspace, &
                    aid, Error)
    call h5awrite_f(aid, h5t_ieee_f64le, time, adims, Error)
    call h5aclose_f(aid, Error)

    call h5dclose_f(dset_id, Error)

    call h5sclose_f(filspace_id, Error)
    call h5sclose_f(memspace_id, Error)
    call h5pclose_f(plist_id_d, Error)

    call h5gclose_f(gid, Error)
    call h5fclose_f(file_id, Error)
    call h5close_f(Error)

  end subroutine WriteHDF5_ZYplane

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Reading / Writing out Checkpoints & Cubes
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine WriteHDF5(fname, save_pressure)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! Writes full 3D output/checkpoint file

    character(len=55) fname
    logical save_pressure

    real(rkind) tmp(Nx, Nyp, Nzp)

    ! Dataset names
    character(len=10) :: dname

    ! Identifiers
    integer(hid_t) :: file_id, dset_id
    integer(hid_t) :: filspace_id, memspace_id

    ! Identifiers
    integer(hid_t) :: gid, selspace_id
    integer(hid_t) :: plist_id_w, plist_id_d

    ! Dimensions in the memory and in the file
    integer(hsize_t), dimension(3) :: dimsm, dimsf

    integer(hsize_t), dimension(3) :: chunk_dims, count, offset
    integer(hsize_t), dimension(3) :: stride, block1, offset_m

    integer :: rHDF5 = 3, arank = 1

    integer(hsize_t), dimension(1)       :: adims
    integer(hid_t)                      :: aid, tspace
    real, dimension(1)               :: treal(3)
    integer, dimension(1)               :: tint(3)
    character(len=80)                        :: namnbuf
    character(len=20)                        :: sttimec

    integer Error, ith

    double precision En(4)

    dimsm(1:3) = (/Nx, Nyp, Nzp/)
    dimsf(1:3) = (/Nx, (Nyp - 1) * NprocY + 1, Nz/)

    ! Flow fields are saved in the fractional grid. We use a basic
    ! interpolation
    !
    ! u_j+1/2=u_j+u_j+1
    !
    ! on the way back we just invert the relation in order to have
    ! exact values.
    !
    ! NOTE that inverting the formula above requires a solution of a
    ! tridiagonal system in the vertical direction. This is similar
    ! to the Thomas algoritm in the implemantation, in the sensa that
    ! a pipeline strategy must be used. Thus, the parallel performance
    ! is rather poor

    if (rank == 0) &
      write (*, '("Writing flow to " A55)') fname

    chunk_dims(1) = Nx
    chunk_dims(2) = 1
    chunk_dims(3) = Nzp

    block1(1) = Nx
    block1(3) = Nzp

    ! Stride and count for number of rows and columns in each dimension
    stride = 1
    count = 1

    ! Offset determined by the rank of a processor
    offset(1) = 0
    offset(3) = Nzp * rankZ

    offset_m(1:3) = 0
    if (rankY == 0) then
      block1(2) = Nyp
      offset(2) = 0
      offset_m(2) = 0
    else
      block1(2) = (Nyp - 1)
      offset(2) = rankY * (Nyp - 1) + 1
      offset_m(2) = 1
    end if

    ! Initialize interface
    call h5open_f(Error)

    ! Setup file access property list with parallel I/O access
    call h5pcreate_f(h5p_file_access_f, plist_id_d, Error)
    call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_world%mpi_val, &
                            mpi_info_null%mpi_val, Error)

    ! Create the file collectively
    call h5fcreate_f(trim(fname), h5f_acc_trunc_f, &
                    file_id, Error, access_prp=plist_id_d)
    call h5pclose_f(plist_id_d, Error)


    ! -----------------------------
    ! Resolution
    arank = 1
    adims = 3
    call h5screate_simple_f(arank, adims, tspace, Error)

    call h5acreate_f(file_id, 'Resolution', h5t_std_i32le, tspace, &
                    aid, Error)
    tint(1) = Nx
    tint(2) = NprocY * (Nyp - 1) + 1
    tint(3) = Nz
    call h5awrite_f(aid, h5t_native_integer, tint, adims, Error)
    call h5aclose_f(aid, Error)

    call h5sclose_f(tspace, Error)
    ! -----------------------------


    ! -----------------------------
    ! Date
    adims = 20
    call h5screate_simple_f(arank, adims, tspace, Error)

    call h5acreate_f(file_id, 'Date', h5t_c_s1, tspace, aid, &
                    Error)
    call time_string(sttimec)
    call h5awrite_f(aid, h5t_c_s1, sttimec, adims, Error)
    call h5aclose_f(aid, Error)
    call h5sclose_f(tspace, Error)
    ! -----------------------------


    ! -----------------------------
    ! Timey Wimey Stuff

    ! Create the timestep space and group
    call h5screate_f(h5s_scalar_f, tspace, Error)
    call h5gcreate_f(file_id, 'Timestep', gid, Error)

    call h5acreate_f(gid, 'Time', h5t_ieee_f64le, tspace, aid, Error)
    call h5awrite_f(aid, h5t_native_double, time, adims, Error)
    call h5aclose_f(aid, Error)

    call h5acreate_f(gid, 'Save_Flow_Time', h5t_ieee_f64le, tspace, aid, Error)
    call h5awrite_f(aid, h5t_native_double, save_flow_time, adims, Error)
    call h5aclose_f(aid, Error)

    call h5acreate_f(gid, 'Save_Stats_Time', h5t_ieee_f64le, tspace, aid, Error)
    call h5awrite_f(aid, h5t_native_double, save_stats_time, adims, Error)
    call h5aclose_f(aid, Error)

    call h5acreate_f(gid, 'Save_Movie_Time', h5t_ieee_f64le, tspace, aid, Error)
    call h5awrite_f(aid, h5t_native_double, save_movie_time, adims, Error)
    call h5aclose_f(aid, Error)

    call h5acreate_f(gid, 'Time_Step', h5t_std_i32le, tspace, aid, Error)
    call h5awrite_f(aid, h5t_native_integer, time_step, adims, Error)
    call h5aclose_f(aid, Error)

    call h5sclose_f(tspace, Error)
    ! ----------------------------



    ! Convert to physical space
    call fft_xz_to_physical(cu1, u1)
    call fft_xz_to_physical(cu2, u2)
    call fft_xz_to_physical(cu3, u3)
    do ith = 1, N_th
      call fft_xz_to_physical(cth(:, :, :, ith), th(:, :, :, ith))
    end do

    ! Create property list for the chunked dataset creation
    call h5pcreate_f(h5p_dataset_create_f, plist_id_d, Error)
    call h5pset_chunk_f(plist_id_d, rHDF5, chunk_dims, Error)

    ! Create the dataspace for ur
    call h5screate_simple_f(rHDF5, dimsf, filspace_id, Error)
    call h5screate_simple_f(rHDF5, dimsm, memspace_id, &
                            Error)

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id_w, Error)
    call h5pset_dxpl_mpio_f(plist_id_w, h5fd_mpio_collective_f, &
                            Error)

    do ith = 1, 3 + N_th

      select case (ith)
      case (1)
        call swapzy(u1, tmp)
        dname = "U"
      case (2)
        ! Interpolation to the fractional grid
        call g2gf(u2)
        call swapzy(u2, tmp)
        call gf2g(u2)
        dname = "V"
      case (3)
        call swapzy(u3, tmp)
        dname = "W"
      case (4:)
        call swapzy(th(:, :, :, ith - 3), tmp)
        dname = "TH"//char(ith + 45)
      end select

      call h5dcreate_f(gid, trim(dname), h5t_ieee_f64le, &
                      filspace_id, dset_id, Error, dcpl_id=plist_id_d)

      ! Select hyperslab in the file.
      ! call h5dget_space_f(dsetur_id, selspace_id, Error)
      call h5sselect_hyperslab_f(filspace_id, h5s_select_set_f, &
                                offset, count, Error, stride, block1)

      call h5sselect_hyperslab_f(memspace_id, h5s_select_set_f, &
                                offset_m, count, Error, stride, block1)

      ! Write the dataset collectively
      call h5dwrite_f(dset_id, h5t_native_double, &
                      tmp, &
                      dimsm, Error, file_space_id=filspace_id, &
                      mem_space_id=memspace_id, xfer_prp=plist_id_w)

      ! Close dateset
      call h5dclose_f(dset_id, Error)
    end do

    ! If want to use this as a checkpoint file...
    if (save_pressure) then
      call fft_xz_to_physical(cp, p)

      call swapzy(p, tmp)
      dname = "P"

      call h5dcreate_f(gid, trim(dname), h5t_ieee_f64le, &
                      filspace_id, dset_id, Error, dcpl_id=plist_id_d)

      ! Select hyperslab in the file.
      ! call h5dget_space_f(dsetur_id, selspace_id, Error)
      call h5sselect_hyperslab_f(filspace_id, h5s_select_set_f, &
                                offset, count, Error, stride, block1)

      call h5sselect_hyperslab_f(memspace_id, h5s_select_set_f, &
                                offset_m, count, Error, stride, block1)

      ! Write the dataset collectively
      call h5dwrite_f(dset_id, h5t_native_double, &
                      tmp, &
                      dimsm, Error, file_space_id=filspace_id, &
                      mem_space_id=memspace_id, xfer_prp=plist_id_w)

      ! Close dateset
      call h5dclose_f(dset_id, Error)
      call fft_xz_to_fourier(p, cp)
    end if

    ! Close the dataspace for the memory and for the file
    call h5sclose_f(filspace_id, Error)
    call h5sclose_f(memspace_id, Error)

    ! Close the properties for the dataspace creation and the writing
    call h5pclose_f(plist_id_d, Error)
    call h5pclose_f(plist_id_w, Error)

    ! Close groups
    call h5gclose_f(gid, Error)
    call h5fclose_f(file_id, Error)
    call h5close_f(Error)

    call fft_xz_to_fourier(u1, cu1)
    call fft_xz_to_fourier(u2, cu2)
    call fft_xz_to_fourier(u3, cu3)
    do ith = 1, N_th
      call fft_xz_to_fourier(th(:, :, :, ith), cth(:, :, :, ith))
    end do

    ! call mpi_finalize(ierror)
    ! stop

  end subroutine WriteHDF5


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine ReadHDF5(fname, read_pressure)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! Reads full 3D output/checkpoint file

    character(len=55), intent(in) :: fname
    logical, intent(out) :: read_pressure

    real(rkind) tmp(Nx, Nyp, Nzp)
    !
    ! HDF5 ------------------------------------------------------
    !
    ! Dataset names
    character(len=10) :: dname

    ! Identifiers
    integer(hid_t) :: file_id, dset_id
    integer(hid_t) :: filspace_id, memspace_id

    ! Identifiers
    integer(hid_t) :: gid, selspace_id
    integer(hid_t) :: plist_id_w, plist_id_d

    ! Dimensions in the memory and in the file
    integer(hsize_t), dimension(3) :: dimsm, dimsf

    integer(hsize_t), dimension(3) :: count, offset
    integer(hsize_t), dimension(3) :: stride, block1, offset_m

    integer :: rHDF5 = 3

    integer(hsize_t), dimension(1)       :: adims
    integer(hid_t)                      :: aid, tspace
    real, dimension(1)               :: treal(3)
    integer, dimension(1)               :: tint(3)
    character(len=80)                        :: namnbuf
    character(len=20)                        :: sttimec

    integer Error, ith

    double precision En(4)

    dimsm(1:3) = (/Nx, Nyp, Nzp/)
    dimsf(1:3) = (/Nx, (Nyp - 1) * NprocY + 1, Nz/)

    ! Flow fields are saved in the fractional grid. We use a basic
    ! interpolation
    !
    ! u_j+1/2=u_j+u_j+1
    !
    ! on the way back we just invert the relation in order to have
    ! exact values.

    block1(1) = Nx
    block1(3) = Nzp

    ! Stride and count for number of rows and columns in each dimension
    stride = 1
    count = 1

    ! Offset determined by the rank of a processor
    offset(1) = 0
    offset(3) = rankZ * Nzp

    block1(2) = Nyp
    offset(2) = rankY * (Nyp - 1)

    ! Initialize interface
    call h5open_f(Error)

    ! Setup file access property list with parallel I/O access
    call h5pcreate_f(h5p_file_access_f, plist_id_d, Error)
    call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_world%mpi_val, &
                            mpi_info_null%mpi_val, Error)

    ! Create the file collectively
    call h5fopen_f(trim(fname), h5f_acc_rdonly_f, &
                  file_id, Error, access_prp=plist_id_d)
    call h5pclose_f(plist_id_d, Error)

    adims = 3

    ! -----------------------------
    ! Resolution
    ! -----------------------------
    call h5aopen_by_name_f(file_id, '.', 'Resolution', aid, Error)
    call h5aread_f(aid, h5t_native_integer, tint, adims, Error)
    call h5aclose_f(aid, Error)
    ! Check that the resolution is of the same kind
    if (verbosity > 4 .and. rank == 0) &
      write (*, '("Reading resolution")')
    if ((tint(1) /= Nx) .or. &
        (tint(2) /= (Nyp - 1) * NprocY + 1) .or. &
        (tint(3) /= Nz)) then
      if (rank == 0) then
        write (*, '("Error. File and program have different resolution.")')
        write (*, '("Program: ", 3I10 )') Nx, (Nyp - 1) * NprocY + 1, Nz
        write (*, '("File: ", 3I10 )') tint(1:3)
      end if
      call mpi_finalize(ierror)
      stop
    end if


    ! -----------------------------
    ! Timey Wimey Stuff
    ! -----------------------------
    if (verbosity > 4 .and. rank == 0) &
      write (*, '("Reading time and output parameters")')

    call h5gopen_f(file_id, "/Timestep", gid, Error)
    
    call h5aopen_by_name_f(gid, '.', 'Time', aid, Error)
    call h5aread_f(aid, h5t_native_double, time, adims, Error)
    call h5aclose_f(aid, Error)

    call h5aopen_by_name_f(gid, '.', 'Save_Flow_Time', aid, Error)
    call h5aread_f(aid, h5t_native_double, save_flow_time, adims, Error)
    call h5aclose_f(aid, Error)

    call h5aopen_by_name_f(gid, '.', 'Save_Stats_Time', aid, Error)
    call h5aread_f(aid, h5t_native_double, save_stats_time, adims, Error)
    call h5aclose_f(aid, Error)

    call h5aopen_by_name_f(gid, '.', 'Save_Movie_Time', aid, Error)
    call h5aread_f(aid, h5t_native_double, save_movie_time, adims, Error)
    call h5aclose_f(aid, Error)

    call h5aopen_by_name_f(gid, '.', 'Time_Step', aid, Error)
    call h5aread_f(aid, h5t_native_integer, time_step, adims, Error)
    call h5aclose_f(aid, Error)

    ! -----------------------------

    if (verbosity > 4 .and. rank == 0) &
      write (*, '("Reading variables ...")')

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id_w, Error)
    call h5pset_dxpl_mpio_f(plist_id_w, h5fd_mpio_collective_f, &
                            Error)

    ! Dataspace in memory
    call h5screate_simple_f(rHDF5, dimsm, memspace_id, &
                            Error)

    do ith = 1, 3 + N_th
      ! Here it starts the loop--->
      select case (ith)
      case (1)
        dname = "U"
      case (2)
        dname = "V"
      case (3)
        dname = "W"
      case (4:)
        dname = "TH"//char(ith + 45)
      end select

      if ((ith <= 3) .or. &
          ((ith > 3) .and. (.not. create_new_th(max(1, ith - 3))))) then

        ! Check to make sure that we should read in this scalar
          if (verbosity > 4 .and. rank == 0) &
            write (*, '("Reading ",A)') dname

        call h5dopen_f(gid, trim(dname), dset_id, Error)
        call h5dget_space_f(dset_id, filspace_id, Error)

        ! Select hyperslab in the file
        call h5sselect_hyperslab_f(filspace_id, h5s_select_set_f, &
                                  offset, count, Error, stride, block1)

        offset_m(1:3) = 0
        call h5sselect_hyperslab_f(memspace_id, h5s_select_set_f, &
                                  offset_m, count, Error, stride, block1)

        ! Write the dataset collectively
        call h5dread_f(dset_id, h5t_native_double, &
                      tmp, &
                      dimsm, Error, file_space_id=filspace_id, &
                      mem_space_id=memspace_id) !, xfer_prp = plist_id_w)

        select case (ith)
        case (1)
          call swapyz(tmp, u1)
        case (2)
          call swapyz(tmp, u2)
          ! Interpolation to the collocated grid
          call gf2g(u2)
        case (3)
          call swapyz(tmp, u3)
        case (4:)
          call swapyz(tmp, th(:, :, :, ith - 3))
        end select

        ! Close dateset
        call h5sclose_f(filspace_id, Error)
        call h5dclose_f(dset_id, Error)

      end if

    end do

    if (verbosity > 4 .and. rank == 0) &
      write (*, '("Checking if pressure exists")')

    ! Decide whether to compute the pressure or to read
    call h5lexists_f(gid, 'P', read_pressure, Error)
    if (read_pressure) then
      if (verbosity > 4 .and. rank == 0) &
      write (*, '("Reading pressure")')

      dname = "P"
      call h5dopen_f(gid, trim(dname), dset_id, Error)
      call h5dget_space_f(dset_id, filspace_id, Error)

      ! Select hyperslab in the file
      call h5sselect_hyperslab_f(filspace_id, h5s_select_set_f, &
                                offset, count, Error, stride, block1)

      offset_m(1:3) = 0
      call h5sselect_hyperslab_f(memspace_id, h5s_select_set_f, &
                                offset_m, count, Error, stride, block1)

      ! Write the dataset collectively
      call h5dread_f(dset_id, h5t_native_double, &
                    tmp, &
                    dimsm, Error, file_space_id=filspace_id, &
                    mem_space_id=memspace_id) !, xfer_prp = plist_id_w)

      call swapyz(tmp, p)
      ! Close dateset
      call h5sclose_f(filspace_id, Error)
      call h5dclose_f(dset_id, Error)
      call fft_xz_to_fourier(p, cp)
    end if

    ! Close the dataspace for the memory
    call h5sclose_f(memspace_id, Error)

    ! Close the properties for the reading
    call h5pclose_f(plist_id_w, Error)

    ! Close groups
    call h5gclose_f(gid, Error)
    call h5fclose_f(file_id, Error)
    call h5close_f(Error)

    if (verbosity > 4 .and. rank == 0) &
      write (*, '("Converting to Fourier space")')

    ! Convert to physical space
    call fft_xz_to_fourier(u1, cu1)
    call fft_xz_to_fourier(u2, cu2)
    call fft_xz_to_fourier(u3, cu3)
    do ith = 1, N_th
      if (.not. create_new_th(ith)) then
        call fft_xz_to_fourier(th(:, :, :, ith), cth(:, :, :, ith))
      end if
    end do

    if (verbosity > 4 .and. rank == 0) &
      write (*, '("Done reading from HDF5")')

  end subroutine ReadHDF5

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine swapzy(in, out)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    real(rkind) in(0:Nx + 1, 0:Nzp + 1, 0:Nyp + 1)
    real(rkind) out(1:Nx, 1:Nyp, 1:Nzp)
    integer x, z, y

    out = 0.d0
    do x = 0, Nx - 1
      do y = 1, Nyp
        do z = 0, Nzp - 1
          out(x + 1, y, z + 1) = in(x, z, y)
        end do
      end do
    end do

  end subroutine

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine swapyz(in, out)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    real(rkind) out(0:Nx + 1, 0:Nzp + 1, 0:Nyp + 1)
    real(rkind) in(1:Nx, 1:Nyp, 1:Nzp)
    integer x, z, y

    out = 0.d0
    do x = 0, Nx - 1
      do y = 1, Nyp
        do z = 0, Nzp - 1
          out(x, z, y) = in(x + 1, y, z + 1)
        end do
      end do
    end do

  end subroutine

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine gf2g(var)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    real(rkind) var(0:Nx + 1, 0:Nzp + 1, 0:Nyp + 1)
    integer x, z, y
    type(mpi_datatype) :: xzbox

    ! Define new data type
    call mpi_type_vector(Nzp, Nx, Nx + 2, mpi_double_precision, &
                         xzbox, ierror)
    call mpi_type_commit(xzbox, ierror)

    if (rankY /= NprocY - 1) then
      call mpi_recv(var(0, 0, Nyp + 1), 1, xzbox, rankY + 1, 101 + rankY, &
                    mpi_comm_y, status, ierror)
    else
      do x = 0, Nx - 1
        do z = 0, Nzp - 1
          var(x, z, Nyp + 1) = var(x, z, Nyp)
        end do
      end do
    end if

    do x = 0, Nx - 1
      do z = 0, Nzp - 1
        do y = Nyp, 1, -1
          var(x, z, y) = 2 * var(x, z, y) - var(x, z, y + 1)
        end do
      end do
    end do

    if (rankY /= 0) call mpi_send(var(0, 0, 2), 1, &
                                    xzbox, rankY - 1, 100 + rankY, &
                                    mpi_comm_y, ierror)

    ! Impose the values at the boundary as prescribed in the
    ! code in order to have zero mass flux
    if (rankY == 0) var(:, :, 1) = -var(:, :, 2)
    if (rankY == NprocY - 1) var(:, :, Nyp + 1) = -var(:, :, Nyp)

  end subroutine

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine g2gf(var)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    real(rkind) var(0:Nx + 1, 0:Nzp + 1, 0:Nyp + 1)
    integer x, z, y

    do x = 0, Nx - 1
      do z = 0, Nzp - 1
        if (rankY == 0) var(x, z, 1) = var(x, z, 2)
        if (rankY == NprocY - 1) var(x, z, Nyp + 1) = var(x, z, Nyp)
        do y = 1, Nyp
          var(x, z, y) = 0.5 * (var(x, z, y) + var(x, z, y + 1))
        end do
      end do
    end do

  end subroutine

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine time_string(cdt)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    !
    ! Construct string in the format '19-DEC-2005 22:47:06'
    !
    implicit none

    integer :: i

    integer :: val(8)
    character(len=20), intent(out) :: cdt
    character(len=3) :: monc

    call date_and_time(values=val)

    if (val(2) == 1) then
      monc = 'JAN'
    else if (val(2) == 2) then
      monc = 'FEB'
    else if (val(2) == 3) then
      monc = 'MAR'
    else if (val(2) == 4) then
      monc = 'APR'
    else if (val(2) == 5) then
      monc = 'MAY'
    else if (val(2) == 6) then
      monc = 'JUN'
    else if (val(2) == 7) then
      monc = 'JUL'
    else if (val(2) == 8) then
      monc = 'AUG'
    else if (val(2) == 9) then
      monc = 'SEP'
    else if (val(2) == 10) then
      monc = 'OCT'
    else if (val(2) == 11) then
      monc = 'NOV'
    else if (val(2) == 12) then
      monc = 'DEC'
    else
      monc = 'XXX'
    end if

    write (cdt, '(i2,a1,a3,a1,i4,a1,i2,a1,i2,a1,i2)') &
      val(3), '-', monc, '-', val(1), ' ', val(5), ':', val(6), ':', val(7)
    do i = 1, 2
      if (cdt(i:i) == ' ') then
        cdt(i:i) = '0'
      end if
    end do
    do i = 13, 20
      if (cdt(i:i) == ' ') then
        cdt(i:i) = '0'
      end if
    end do

  end subroutine time_string

end module phdf5