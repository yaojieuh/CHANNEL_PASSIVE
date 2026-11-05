!--------------------------------------------------------------------
! This module contains variables and subroutines 
! o be used when performing DFT and global data transpose 
! with 2D decomposition (a.k.a. pencil decomposition)
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------

module ppdfftw
  use, intrinsic :: iso_c_binding 
  use omp_lib
  use timer
  use runge_kutta
  use global_variable
  implicit none
  include 'mpif.h'
  include 'fftw3-mpi.f03'


  !--------------------------------------------------
  ! Global MPI variables
  !--------------------------------------------------
  integer :: ierr,nprocs,myrank
  !--------------------------------------------------


  !--------------------------------------------------
  ! Variables for 2 dimensional MPI topology
  !--------------------------------------------------
  integer :: nprocs2D, rank2D !Looks like not necessary
  integer :: ndims     ! dimension of MPI topology. Always 2
  integer :: dimsize(2)! size of MPI topology in 2 direction
  integer :: nprocs1DA ! =dimsize(1)
  integer :: nprocs1DB ! =dimsize(2)
  integer :: myrankA   ! rank in A direction
  integer :: myrankB   ! rank in B direction
  !--------------------------------------------------


  !--------------------------------------------------
  ! =0 : does not show detailed topology, show otherwise
  !--------------------------------------------------
  integer :: show_mpi_topology
  !--------------------------------------------------


  !--------------------------------------------------
  ! Communicators
  !--------------------------------------------------
  integer :: COMM2D, COMM1DA, COMM1DB
  !-------------------------------------------------- 


  !-----------------------------------------------------------------------
  ! Variables related to size of fields
  ! First capital character : pencil direction
  ! Second small character : direction of data
  ! Third and Fourth characters : sz(=size) st(=start) en(=end,redundant)
  ! Example : Yxsz --> x direction size of Y-pencil 
  !-----------------------------------------------------------------------
  integer(C_INTPTR_T) :: nx,ny,nz,nxh,nxhpad,nz32,nx32,nx32hpad
  integer(C_INTPTR_T) :: Xysz,Xyst,Xyen,Xzsz,Xzst,Xxen
  integer(C_INTPTR_T) :: Yxsz,Yxst,Yxen,Yzsz,Yzst,Yzen
  integer(C_INTPTR_T) :: Zxsz,Zxst,Zxen,Zysz,Zyst,Zyen
  integer(C_INTPTR_T) :: alloc_size1,alloc_size2,alloc_size,no_of_field
  integer(C_INTPTR_T) :: omp_x, omp_y, omp_z
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  ! Wavenumbers which will be used in wavespace
  !-----------------------------------------------------------------------
  integer, allocatable :: kx(:),kz(:)
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  ! Buf line to calculate statistics
  !-----------------------------------------------------------------------
  real(C_DOUBLE),allocatable :: stat_buf_line(:), stat_buf_line_dummy(:)
  real(C_DOUBLE),allocatable :: stat_PI1_uu(:), stat_PI1_vv(:), stat_PI1_ww(:), stat_PI1_uv(:)
  real(C_DOUBLE),allocatable :: stat_uu(:),stat_vv(:),stat_ww(:),stat_uv(:),stat_vw(:),stat_uw(:)
  real(C_DOUBLE),allocatable :: stat_tt(:),stat_ut(:),stat_vt(:),stat_wt(:)
  real(C_DOUBLE),allocatable :: stat_uuv(:),stat_vww(:),stat_vvv(:),stat_uvv(:)
  real(C_DOUBLE),allocatable :: stat_U(:), stat_W(:), stat_dUdy(:),stat_dWdy(:), stat_T(:), stat_dTdy(:)   !<---- need to allocate
  real(C_DOUBLE),allocatable :: stat_buf_plane_x(:,:), stat_buf_plane_z(:,:)
  real(C_DOUBLE),allocatable :: stat_buf_plane_x_dummy(:,:), stat_buf_plane_z_dummy(:,:)
  real(C_DOUBLE),allocatable :: stat_skewness(:,:), stat_flatness(:,:)
  real(C_DOUBLE),allocatable :: stat_hs6(:,:), stat_hs8(:,:), stat_hs10(:,:)
  real(C_DOUBLE),allocatable :: stat_real_omp_buf(:,:,:)
  complex(C_DOUBLE_COMPLEX),allocatable :: stat_PPE_boundary(:,:,:)
  complex(C_DOUBLE_COMPLEX),allocatable :: stat_buf_line_dummy1(:),stat_buf_line_dummy2(:)
  real(C_DOUBLE), pointer :: stat_buf_field(:,:)             ! buffer field for statistics calculation
  type(C_PTR) :: p_stat_buf_field
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  ! Memory space to hold actual data and buffer for transpose/dft
  !-----------------------------------------------------------------------
  real(C_DOUBLE), pointer :: Upad(:,:)          ! Field padded in z direction 
  real(C_DOUBLE), pointer :: U(:,:)             ! Field in original size
  real(C_DOUBLE), pointer :: RK_Store(:,:)      ! Field for Nonlinear Storage with RK method
  real(C_DOUBLE), pointer :: buf(:),buf_pad(:)  ! Buffer field for transpose
  real(C_DOUBLE), pointer :: buf1(:),buf_pad1(:)  ! Buffer field for transpose
  real(C_DOUBLE), pointer :: buf_line_x_r(:,:,:)  ! Buffer line for DFT in X (real)
  complex(C_DOUBLE_COMPLEX), pointer :: buf_line_x_c(:,:,:)  ! Buffer line for DFT in X (complex)
  complex(C_DOUBLE_COMPLEX), pointer :: buf_line_y_c1(:,:)   ! Buffer line for wavespace in Y (complex)
  complex(C_DOUBLE_COMPLEX), pointer :: buf_line_y_c2(:,:)   ! Buffer line for wavespace in Y (complex)
  complex(C_DOUBLE_COMPLEX), pointer :: buf_line_z(:,:)      ! Buffer line for DFT in Z (complex)
  type(C_PTR) :: p_U,p_Upad,p_RK_Store
  type(C_PTR) :: p_buf,p_buf_pad 
  type(C_PTR) :: p_buf1,p_buf_pad1 
  type(C_PTR) :: p_buf_line_x_c
  type(C_PTR) :: p_buf_line_y_c1
  type(C_PTR) :: p_buf_line_y_c2
  type(C_PTR) :: p_buf_line_x_r
  type(C_PTR) :: p_buf_line_z
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  ! parameters to make a fftw_transpose plan
  !-----------------------------------------------------------------------
  integer(C_INTPTR_T), parameter :: ONE=1
  integer(C_INTPTR_T), parameter :: TWO=2
  integer(C_INTPTR_T), parameter :: THREE=3
  integer(C_INTPTR_T), parameter :: FOUR=4
  integer(C_INTPTR_T), parameter :: FIVE=5
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  ! FFTW_FLAGs for transpose / dft
  ! Note : need to choose only one
  !-----------------------------------------------------------------------
  !  integer, parameter :: FFTW_FLAG_DFT       = FFTW_ESTIMATE
  !  integer, parameter :: FFTW_FLAG_TRANSPOSE = FFTW_ESTIMATE
    integer, parameter :: FFTW_FLAG_DFT       = FFTW_MEASURE
    integer, parameter :: FFTW_FLAG_TRANSPOSE = FFTW_MEASURE
  !  integer, parameter :: FFTW_FLAG_DFT       = FFTW_PATIENT
  !  integer, parameter :: FFTW_FLAG_TRANSPOSE = FFTW_PATIENT
  !  integer, parameter :: FFTW_FLAG_DFT       = FFTW_EXHAUSTIVE
  !  integer, parameter :: FFTW_FLAG_TRANSPOSE = FFTW_EXHAUSTIVE
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  ! FFTW_PLANs for transpose / dft
  !-----------------------------------------------------------------------
  type(C_PTR) :: plan_transpose_x_to_z
  type(C_PTR) :: plan_transpose_z_to_x
  type(C_PTR) :: plan_transpose_z_to_y
  type(C_PTR) :: plan_transpose_y_to_z
  type(C_PTR) :: plan_transform_forward_z
  type(C_PTR) :: plan_transform_inverse_z
  type(C_PTR) :: plan_transform_forward_x
  type(C_PTR) :: plan_transform_inverse_x
  !-----------------------------------------------------------------------  


  !-----------------------------------------------------------------------
  ! List of Public Subroutines
  !-----------------------------------------------------------------------
  public :: ppdfftw_x_to_z,ppdfftw_z_to_y,ppdfftw_y_to_z,ppdfftw_z_to_x
  public :: ppdfftw_dft_z_forward, ppdfftw_dft_z_inverse
  public :: ppdfftw_initialize, ppdfftw_finalize
  !-----------------------------------------------------------------------


contains


  !-----------------------------------------------------------------------------------------
  ! this subroutine will calculate the pencil size and start / end point 
  !-----------------------------------------------------------------------------------------
  subroutine calculate_memory_size
    implicit none
    integer :: dummy1, dummy2, dummy_max1, dummy_max2, i

    alloc_size1 = fftw_mpi_local_size_2d_transposed(nxh,nz32,COMM1DA,Zxsz,Zxst,Xzsz,Xzst)
    alloc_size2 = fftw_mpi_local_size_2d_transposed(nz  ,ny ,COMM1DB,Yzsz,Yzst,Zysz,Zyst)

    alloc_size1 = alloc_size1 * ny * 2 / dimsize(2)
    alloc_size2 = alloc_size2 * nx / dimsize(1)

    dummy1 = alloc_size1
    dummy2 = alloc_size2

    Xysz = Zysz
    Xyst = Zyst
    Xyen = Zyen
    Yxsz = Zxsz
    Yxst = Zxst
    Yxen = Zxen


    ! Using dummies looks very stupid... but necessary to use with XL compiler.
    call mpi_reduce(dummy1, dummy_max1, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(dummy2, dummy_max2, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    call mpi_bcast(dummy_max1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(dummy_max2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    alloc_size1 = dummy_max1
    alloc_size2 = dummy_max2

    ! print out the topology of mpi_ranks
    if (show_mpi_topology .ne. 0) then
       if (myrank .eq. 0) then
          write(*,'(11a10,2a14)') 'RankGlobal','RankA','RankB','Xzsz','Xzst','Zysz ','Zyst ','Zxsz ','Zxst ','Yzsz','Yzst',&
               'Alloc_size_1','Alloc_size_2'
          write(*,'(11a10)') '','','','','','(=Xysz)','(Xyst)','(=Yxsz)','(=Yxst)','',''
       endif
       call mpi_barrier(MPI_COMM_WORLD,ierr)

       do i=0,(nprocs-1)
          call mpi_barrier(MPI_COMM_WORLD,ierr)
          if (myrank.eq.i) then
             write(*,'(11i10,2i14)') myrank,myrankA,myrankB,Xzsz,Xzst,Zysz,Zyst,Zxsz,Zxst,Yzsz,Yzst,alloc_size1, alloc_size2
          endif
          call mpi_barrier(MPI_COMM_WORLD,ierr)
       enddo
       call mpi_barrier(MPI_COMM_WORLD,ierr)
       if (myrank.eq.0) print*, " "
    endif

    return

  end subroutine calculate_memory_size
  !-----------------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------------
  ! This Subroutine initiate MPI environment 
  ! Also create a 2D cartesian topology and
  ! Create the Communicator in each direction
  !-----------------------------------------------------------------------------------------
  subroutine set_mpi
    implicit none
    logical :: periods(2), belongs(2), reorder

    ! Decomposed by 2 direction (a.k.a. pencil decomposition)
    ndims = 2
    ! no periodicity
    periods(1) = .false.                                  
    periods(2) = .false.

    ! Create Catesian Topology 2D 
    reorder = .false. ! This is better to be false when using BG/Q
    call mpi_cart_create(MPI_COMM_WORLD,ndims,dimsize,periods,reorder,COMM2D,ierr)
    call mpi_comm_rank(COMM2D,rank2D,ierr)

    ! Create Communicator in 1st direction
    belongs(1) = .true.
    belongs(2) = .false.
    call mpi_cart_sub(COMM2D,belongs,COMM1DA,ierr)

    ! Create Communicator in 2nd direction
    belongs(1) = .false.
    belongs(2) = .true.
    call mpi_cart_sub(COMM2D,belongs,COMM1DB,ierr)

    ! Number of processes in each direction
    nprocs1DA = dimsize(1)
    nprocs1DB = dimsize(2)

    call mpi_comm_rank(COMM1DA,myrankA,ierr)
    call mpi_comm_rank(COMM1DB,myrankB,ierr)

    if (myrank.eq.0) then
       print*, "MPI Topology Information"
       print*, "nprocs    = ", nprocs
       print*, "nprocs1DA = ", nprocs1DA
       print*, "nprocs1DB = ", nprocs1DB
       print*, "-------------------------------------------------"
       print*, " " 
    endif

    return
  end subroutine set_mpi
  !-----------------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------------
  ! This module allocate aligned memory using FFTW
  !-----------------------------------------------------------------------------------------
  subroutine allocate_memory_field  
    implicit none

    omp_x = omp_compute
    omp_y = omp_compute
    omp_z = omp_compute

    p_Upad = fftw_alloc_real(int(alloc_size1*no_of_field,C_SIZE_T))
    call c_f_pointer(p_Upad,Upad,[alloc_size1,no_of_field])
 
    p_U = fftw_alloc_real(int(alloc_size2*no_of_field,C_SIZE_T))
    call c_f_pointer(p_U,U,[alloc_size2,no_of_field])

    p_RK_Store = fftw_alloc_real(int(alloc_size2*3,C_SIZE_T))
    call c_f_pointer(p_RK_Store,RK_Store,[alloc_size2,three])

    p_buf_pad = fftw_alloc_real(int(alloc_size1,C_SIZE_T))
    call c_f_pointer(p_buf_pad,buf_pad,[alloc_size1])

    p_buf = fftw_alloc_real(int(alloc_size2,C_SIZE_T))
    call c_f_pointer(p_buf,buf,[alloc_size2])

    p_buf_pad1 = fftw_alloc_real(int(alloc_size1,C_SIZE_T))
    call c_f_pointer(p_buf_pad1,buf_pad1,[alloc_size1])

    p_buf1 = fftw_alloc_real(int(alloc_size2,C_SIZE_T))
    call c_f_pointer(p_buf1,buf1,[alloc_size2])

    p_buf_line_x_r = fftw_alloc_real(int(nx32*no_of_field*omp_x,C_SIZE_T))
    call c_f_pointer(p_buf_line_x_r,buf_line_x_r,[nx32,no_of_field,omp_x])

    p_buf_line_x_c = fftw_alloc_complex(int((nx32/2+1)*no_of_field*omp_x,C_SIZE_T))
    call c_f_pointer(p_buf_line_x_c,buf_line_x_c,[(nx32/2+1),no_of_field,omp_x])

    p_buf_line_y_c1 = fftw_alloc_complex(int(ny*omp_y,C_SIZE_T))
    call c_f_pointer(p_buf_line_y_c1,buf_line_y_c1,[ny,omp_y])

    p_buf_line_y_c2 = fftw_alloc_complex(int(ny*omp_y,C_SIZE_T))
    call c_f_pointer(p_buf_line_y_c2,buf_line_y_c2,[ny,omp_y])

    p_buf_line_z = fftw_alloc_complex(int(nz32*omp_z,C_SIZE_T))
    call c_f_pointer(p_buf_line_z,buf_line_z,[nz32,omp_z])

    p_stat_buf_field = fftw_alloc_real(int(alloc_size2*three,C_SIZE_T))
    call c_f_pointer(p_stat_buf_field,stat_buf_field,[alloc_size2,three])



    return
  end subroutine allocate_memory_field
  !-----------------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------------
  !This subroutine assigns wave numbers which belongs to each MPI tasks
  !-----------------------------------------------------------------------------------------
  subroutine set_waveno
    implicit none
    integer, allocatable :: kx_total(:)
    integer, allocatable :: kz_total(:)
    integer :: i,nzh

    allocate(kx_total(nxhpad))
    allocate(kz_total(nz))
    allocate(kx(Yxsz))
    allocate(kz(Yzsz))

    nzh = nz/2

    do i=1,nxhpad
       kx_total(i) = i-1
    enddo

    do i=1,nzh+1
       kz_total(i) = i-1
    enddo

    do i=nzh+2,nz
       kz_total(i) = i-2*nzh-1
    enddo

    do i=1,Yxsz
       kx(i) = kx_total(i+Yxst)
    enddo

    do i=1,Yzsz
       kz(i) = kz_total(i+Yzst)
    enddo

    deallocate(kx_total)
    deallocate(kz_total)

    return
  end subroutine set_waveno
  !-----------------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------------
  ! This subroutine will create plan to use FFTW
  !
  ! This subroutine creat the plan for discrete Fourier tranform in X and Z direction. (forward / inverse)
  !
  ! This subroutine create the plan for transpose from z pencil to x pencil.
  ! The number 2 indicate the complex number.
  ! Data re-ordering is required in 2D MPI topology regardless of FFTW_MPI_TRANSPOSED_OUT/IN flag.
  ! Hence, FFTW_MPI_TRANSPOSED_OUT flag is added to turn off the underneath data reordering in FFTW_TRANSPOSE calls.
  !                                           
  ! IMPORTANT : Need to perform data re-order manually
  !-----------------------------------------------------------------------------------------
  subroutine create_plan
    use timer
    implicit none
    integer(C_INT) :: nx32_fortran, nz32_fortran
    COMPLEX(C_DOUBLE_COMPLEX), POINTER :: DUMMY(:)      ! DUMMY BUFFER LINE TO CREATE PLAN FOR DFT IN Z DIRECTION
    TYPE(C_PTR) :: P_DUMMY

    ! THIS IS NEEDED TO ENSURE COMPATIBILITY OF ISO_C_BINDING WITH SOME COMPILERS (E.G. XL COMPILER)
    NX32_FORTRAN = NX32*1
    NZ32_FORTRAN = NZ32*1

    P_DUMMY = FFTW_ALLOC_COMPLEX(INT(NZ32,C_SIZE_T))
    CALL C_F_POINTER(P_DUMMY,DUMMY,[NZ32])

    IF (MYRANK.EQ.0) time_plan_start = MPI_WTIME()


    !------------------------------------------------------------------------------------------------------------------------------
    ! Create Transpose Plan
    !------------------------------------------------------------------------------------------------------------------------------
    ! PLAN FOR TRANSPOSE FROM Y PENCIL TO Z PENCIL
    IF (MYRANK.EQ.0) WRITE(*,'(A)',ADVANCE='NO')  " CREATING TRANSPOSE PLAN Y TO Z........"
    PLAN_TRANSPOSE_Y_TO_Z = &
         FFTW_MPI_PLAN_MANY_TRANSPOSE(NZ/YZSZ,NY/ZYSZ,2*ZYSZ*YZSZ*YXSZ,ONE,ONE,U,BUF,COMM1DB,FFTW_FLAG_TRANSPOSE)
    IF (MYRANK.EQ.0) WRITE(*,'(A)') "DONE"
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

    ! PLAN FOR TRANSPOSE FROM Z PENCIL TO Y PENCIL
    IF (MYRANK.EQ.0) WRITE(*,'(A)',ADVANCE='NO')  " CREATING TRANSPOSE PLAN Z TO Y........"
    PLAN_TRANSPOSE_Z_TO_Y = &
         FFTW_MPI_PLAN_MANY_TRANSPOSE(NY/ZYSZ,NZ/YZSZ,2*YZSZ*ZYSZ*YXSZ,ONE,ONE,U,BUF,COMM1DB,FFTW_FLAG_TRANSPOSE)
    IF (MYRANK.EQ.0) WRITE(*,'(A)') "DONE"
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

    ! PLAN FOR TRANSPOSE FROM Z PENCIL TO X PENCIL 
    IF (MYRANK.EQ.0) WRITE(*,'(A)',ADVANCE='NO')  " CREATING TRANSPOSE PLAN Z TO X........"
    PLAN_TRANSPOSE_Z_TO_X = &
         FFTW_MPI_PLAN_MANY_TRANSPOSE(NXH/ZXSZ,NZ32/XZSZ,2*XZSZ*ZXSZ*XYSZ,ONE,ONE,UPAD,BUF_PAD,COMM1DA,FFTW_FLAG_TRANSPOSE)
    IF (MYRANK.EQ.0) WRITE(*,'(A)') "DONE"
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

    ! PLAN FOR TRANSPOSE FROM X PENCIL TO Z PENCIL
    IF (MYRANK.EQ.0) WRITE(*,'(A)',ADVANCE='NO')  " CREATING TRANSPOSE PLAN X TO Z........"
    PLAN_TRANSPOSE_X_TO_Z = &
         FFTW_MPI_PLAN_MANY_TRANSPOSE(NZ32/XZSZ,NXH/ZXSZ,2*ZXSZ*XYSZ*XZSZ,ONE,ONE,UPAD,BUF_PAD,COMM1DA,FFTW_FLAG_TRANSPOSE)
    IF (MYRANK.EQ.0) WRITE(*,'(A)') "DONE"
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    !------------------------------------------------------------------------------------------------------------------------------





    !------------------------------------------------------------------------------------------------------------------------------
    ! Create Transform Plan
    !------------------------------------------------------------------------------------------------------------------------------
    ! PLAN FOR FFT COMPLEX TO COMPLEX IN Z PENCIL (INVERSE)
    IF (MYRANK.EQ.0) WRITE(*,'(A)',ADVANCE='NO')  " CREATING TRANSFORM PLAN Z INVERSE....."
    PLAN_TRANSFORM_INVERSE_Z = FFTW_PLAN_DFT_1D(NZ32_FORTRAN,buf_line_z(:,1),dummy,FFTW_BACKWARD,FFTW_FLAG_DFT)
    IF (MYRANK.EQ.0) WRITE(*,'(A)') "DONE"
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

    ! PLAN FOR FFT COMPLEX TO COMPLEX IN Z PENCIL (FORWARD)
    IF (MYRANK.EQ.0) WRITE(*,'(A)',ADVANCE='NO')  " CREATING TRANSFORM PLAN Z FORWARD....."
    PLAN_TRANSFORM_FORWARD_Z = FFTW_PLAN_DFT_1D(NZ32_FORTRAN,dummy,buf_line_z(:,1),FFTW_FORWARD,FFTW_FLAG_DFT)
    if (myrank.eq.0) write(*,'(a)') "DONE"
    call mpi_barrier(MPI_COMM_WORLD,ierr)

    ! PLAN FOR FFT REAL TO COMPLEX IN X PENCIL (FORWARD)
    IF (MYRANK.EQ.0) WRITE(*,'(A)',ADVANCE='NO')  " CREATING TRANSFORM PLAN X FORWARD....."
    PLAN_TRANSFORM_FORWARD_X = FFTW_PLAN_DFT_R2C_1D(NX32_FORTRAN,BUF_LINE_X_R(:,1,1),BUF_LINE_X_C(:,1,1),FFTW_FLAG_DFT)
    IF (MYRANK.EQ.0) WRITE(*,'(A)') "DONE"
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)


    ! PLAN FOR FFT COMPLEX TO REAL IN X PENCIL (INVERSE)
    IF (MYRANK.EQ.0) WRITE(*,'(A)',ADVANCE='NO')  " CREATING TRANSFORM PLAN X INVERSE....."
    PLAN_TRANSFORM_INVERSE_X = FFTW_PLAN_DFT_C2R_1D(NX32_FORTRAN,BUF_LINE_X_C(:,1,1),BUF_LINE_X_R(:,1,1),FFTW_FLAG_DFT)
    IF (MYRANK.EQ.0) WRITE(*,'(A)') "DONE"
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

    if (myrank.eq.0) then 
    IF (MYRANK.EQ.0) time_plan_end = MPI_WTIME()
       write(*,'(a35,es10.2,a4)') " Elapsed Time for Creating Plans = ", &
            time_plan_end - time_plan_start, " sec"
       print*, " "
    endif
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call fftw_free(p_dummy)

    return

  end subroutine create_plan
  !-----------------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------------
  ! This subroutine set the Environment to use ppdfftw 
  !-----------------------------------------------------------------------------------------
  subroutine ppdfftw_initialize
    use global_variable, only: dyr
    use labs, only: LABS_Y
    implicit none
    integer :: i

    ! Setting MPI Topology / Sub communicator
    call set_mpi
    call mpi_barrier(MPI_COMM_WORLD,ierr)

    ! Calculate Memory Size for each MPI tasks
    call calculate_memory_size
    call mpi_barrier(MPI_COMM_WORLD,ierr)

    ! Allocate Memory
    call allocate_memory_field
    call mpi_barrier(MPI_COMM_WORLD,ierr)


    ! Set Wavenumber in Wavespace for each MPI Tasks
    call set_waveno
    call mpi_barrier(MPI_COMM_WORLD,ierr)


    ! allocate buf line to calculate statistics
    allocate(stat_uu(ny))
    allocate(stat_vv(ny))
    allocate(stat_ww(ny))
    allocate(stat_tt(ny))
    allocate(stat_uv(ny))
    allocate(stat_vw(ny))
    allocate(stat_uw(ny))
    allocate(stat_ut(ny))
    allocate(stat_vt(ny))
    allocate(stat_wt(ny))
    allocate(stat_uuv(ny))
    allocate(stat_vww(ny))
    allocate(stat_vvv(ny))
    allocate(stat_uvv(ny))
    allocate(stat_U(ny))
    allocate(stat_W(ny))
    allocate(stat_T(ny))
    allocate(stat_dUdy(ny))
    allocate(stat_dWdy(ny))
    allocate(stat_dTdy(ny))
    allocate(stat_buf_line(ny))
    allocate(stat_buf_line_dummy(ny))
    allocate(stat_buf_line_dummy1(ny))
    allocate(stat_buf_line_dummy2(ny))
    allocate(stat_skewness(ny,3))
    allocate(stat_flatness(ny,3))
    allocate(stat_hs6(ny,3))
    allocate(stat_hs8(ny,3))
    allocate(stat_hs10(ny,3))
    allocate(stat_buf_plane_x(ny,yxsz))
    allocate(stat_buf_plane_z(ny,yzsz))
    allocate(stat_real_omp_buf(21,ny,omp_compute))
    allocate(stat_buf_plane_x_dummy(ny,yxsz))
    allocate(stat_buf_plane_z_dummy(ny,yzsz))
    allocate(stat_PPE_boundary(2,yxsz,yzsz))
    allocate(stat_PI1_uu(ny))
    allocate(stat_PI1_vv(ny))
    allocate(stat_PI1_ww(ny))
    allocate(stat_PI1_uv(ny))


    ! Create FFTW Plans to use for global data transpose
    call create_plan
    call mpi_barrier(MPI_COMM_WORLD,ierr)




    ! dy to be used for calculating time step
    allocate(dyr(Xysz))
    ! This part set the smallest delta y from nearest grid 
    do i=1+Xyst,Xysz+Xyst
       if (i.eq.1) then 
          dyr(i-Xyst)=1d0/(dy_factor*dabs(LABS_Y(2)-LABS_Y(1)))
       else if (i.eq.ny) then
          dyr(i-Xyst) = 1d0/(dy_factor*dabs(LABS_Y(ny) - LABS_Y(ny-1)))
       else
          dyr(i-Xyst) = 1d0/(dy_factor*min(dabs(LABS_Y(i)-LABS_Y(i-1)),dabs(LABS_Y(i)-LABS_Y(i+1))))
       endif
    enddo



    return 
  end subroutine ppdfftw_initialize
  !-----------------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------------
  ! This subroutine deallocate all fields and destroy plans to use ppdfftw
  !-----------------------------------------------------------------------------------------
  subroutine ppdfftw_finalize
    use global_variable, only: dyr
    implicit none

    deallocate(dyr)
    deallocate(kx)
    deallocate(kz)
    deallocate(stat_buf_line)
    deallocate(stat_buf_line_dummy)
    deallocate(stat_buf_line_dummy1)
    deallocate(stat_buf_line_dummy2)
    deallocate(stat_buf_plane_x)
    deallocate(stat_buf_plane_z)
    deallocate(stat_buf_plane_x_dummy)
    deallocate(stat_buf_plane_z_dummy)
    deallocate(stat_uu)
    deallocate(stat_vv)
    deallocate(stat_ww)
    deallocate(stat_uv)
    deallocate(stat_vw)
    deallocate(stat_uw)
    deallocate(stat_uuv)
    deallocate(stat_vww)
    deallocate(stat_vvv)
    deallocate(stat_uvv)
    deallocate(stat_U)
    deallocate(stat_W)
    deallocate(stat_dUdy)
    deallocate(stat_dWdy)
    deallocate(stat_skewness)
    deallocate(stat_flatness)
    deallocate(stat_real_omp_buf)
    deallocate(stat_PPE_boundary)
    deallocate(stat_PI1_uu)
    deallocate(stat_PI1_vv)
    deallocate(stat_PI1_ww)
    deallocate(stat_PI1_uv)
    call fftw_free(p_stat_buf_field)









    call mpi_comm_free(COMM1DA,ierr)
    call mpi_comm_free(COMM1DB,ierr)
    call mpi_comm_free(COMM2D ,ierr)

    call fftw_destroy_plan(plan_transpose_y_to_z)
    call fftw_destroy_plan(plan_transpose_z_to_y)
    call fftw_destroy_plan(plan_transpose_x_to_z)
    call fftw_destroy_plan(plan_transpose_z_to_x)

    call fftw_destroy_plan(plan_transform_forward_z)
    call fftw_destroy_plan(plan_transform_inverse_z)
    call fftw_destroy_plan(plan_transform_inverse_x)
    call fftw_destroy_plan(plan_transform_forward_x)

    call fftw_free(p_U)
    call fftw_free(p_buf)
    call fftw_free(p_buf_pad)
    call fftw_free(p_buf1)
    call fftw_free(p_buf_pad1)
    call fftw_free(p_buf_line_y_c1)
    call fftw_free(p_buf_line_y_c2)
    call fftw_free(p_buf_line_x_c)
    call fftw_free(p_buf_line_x_r)
    call fftw_free(p_buf_line_z)
    call fftw_free(p_RK_Store)
    call fftw_free(p_Upad)


    return
  end subroutine ppdfftw_finalize
  !-----------------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------------
  ! This subroutine transpose the data from x pencil to z pencil
  !-----------------------------------------------------------------------------------------
  subroutine ppdfftw_x_to_z(U_BIG)
    implicit none
    real(C_DOUBLE) :: U_BIG(alloc_size1)

!!$    if(myrank.eq.0) time_transpose_start = MPI_WTIME()
    call on_node_preorder_x_to_z(nxh,nz32,ny,nprocs1DA,nprocs1DB,U_BIG,BUF_PAD1,omp_reorder)
!!$    if(myrank.eq.0) time_transpose_end = MPI_WTIME()
!!$    if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "x to z (pre-shuffle)", &
!!$         timestep,dt,time,time_transpose_end-time_transpose_start

!!$    if(myrank.eq.0) time_transpose_start = MPI_WTIME()
    call fftw_mpi_execute_r2r(plan_transpose_x_to_z,BUF_PAD1,BUF_PAD)
!!$    if(myrank.eq.0) time_transpose_end = MPI_WTIME()
!!$    if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "x to z (communication)", &
!!$         timestep,dt,time,time_transpose_end-time_transpose_start

!!$    if(myrank.eq.0) time_transpose_start = MPI_WTIME()
    call on_node_reorder_x_to_z(nxh,nz32,ny,nprocs1DA,nprocs1DB,BUF_PAD,U_BIG,omp_reorder)
!!$    if(myrank.eq.0) time_transpose_end = MPI_WTIME()
!!$    if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "x to z (re-shuffle)", &
!!$         timestep,dt,time,time_transpose_end-time_transpose_start

    return
  end subroutine ppdfftw_x_to_z
  !-----------------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------------
  ! This subroutine transpose the data from z pencil to y pencil
  !-----------------------------------------------------------------------------------------
  subroutine ppdfftw_z_to_y(U_SMALL)
    implicit none
    real(C_DOUBLE) :: U_SMALL(alloc_size2)
!!$    if(myrank.eq.0) time_transpose_start = MPI_WTIME()
    call on_node_preorder_z_to_y(nxh,nz,ny,nprocs1DA,nprocs1DB,U_SMALL,BUF1,omp_reorder)
!!$    if(myrank.eq.0) time_transpose_end = MPI_WTIME()
!!$    if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "z to y (pre-shuffle)", &
!!$         timestep,dt,time,time_transpose_end-time_transpose_start

!!$    if(myrank.eq.0) time_transpose_start = MPI_WTIME()
    call fftw_mpi_execute_r2r(plan_transpose_z_to_y,BUF1,BUF)
!!$    if(myrank.eq.0) time_transpose_end = MPI_WTIME()
!!$    if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "z to y (communication)", &
!!$         timestep,dt,time,time_transpose_end-time_transpose_start

!!$    if(myrank.eq.0) time_transpose_start = MPI_WTIME()
    call on_node_reorder_z_to_y(nxh,nz,ny,nprocs1DA,nprocs1DB,BUF,U_SMALL,omp_reorder)
!!$    if(myrank.eq.0) time_transpose_end = MPI_WTIME()
!!$    if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "z to y (re-shuffle)", &
!!$         timestep,dt,time,time_transpose_end-time_transpose_start

    return
  end subroutine ppdfftw_z_to_y
  !-----------------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------------
  ! This subroutine transpose the data from y pencil to z pencil
  !-----------------------------------------------------------------------------------------
  subroutine ppdfftw_y_to_z(U_SMALL)
    implicit none
    real(C_DOUBLE) :: U_SMALL(alloc_size2)
!!$    if(myrank.eq.0) time_transpose_start = MPI_WTIME()
    call on_node_preorder_y_to_z(nxh,nz,ny,nprocs1DA,nprocs1DB,U_SMALL,BUF1,omp_reorder)
!!$    if(myrank.eq.0) time_transpose_end = MPI_WTIME()
!!$    if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "y to z (pre-shuffle)", &
!!$         timestep,dt,time,time_transpose_end-time_transpose_start

!!$    if(myrank.eq.0) time_transpose_start = MPI_WTIME()
    call fftw_mpi_execute_r2r(plan_transpose_y_to_z,BUF1,BUF)
!!$    if(myrank.eq.0) time_transpose_end = MPI_WTIME()
!!$    if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "y to z (communication)", &
!!$         timestep,dt,time,time_transpose_end-time_transpose_start

!!$    if(myrank.eq.0) time_transpose_start = MPI_WTIME()
    call on_node_reorder_y_to_z(nxh,nz,ny,nprocs1DA,nprocs1DB,BUF,U_SMALL,omp_reorder)
!!$    if(myrank.eq.0) time_transpose_end = MPI_WTIME()
!!$    if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "y to z (re-shuffle)", &
!!$         timestep,dt,time,time_transpose_end-time_transpose_start

    return
  end subroutine ppdfftw_y_to_z
  !-----------------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------------
  ! This subroutine transpose the data from z pencil to x pencil
  !-----------------------------------------------------------------------------------------
  subroutine ppdfftw_z_to_x(U_BIG)
    implicit none
    real(C_DOUBLE) :: U_BIG(alloc_size1)
!!$    if(myrank.eq.0) time_transpose_start = MPI_WTIME()
    call on_node_preorder_z_to_x(nxh,nz32,ny,nprocs1DA,nprocs1DB,U_BIG,BUF_PAD1,omp_reorder)
!!$    if(myrank.eq.0) time_transpose_end = MPI_WTIME()
!!$    if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "z to x (pre-shuffle)", &
!!$         timestep,dt,time,time_transpose_end-time_transpose_start

!!$    if(myrank.eq.0) time_transpose_start = MPI_WTIME()
    call fftw_mpi_execute_r2r(plan_transpose_z_to_x,BUF_PAD1,BUF_PAD)
!!$    if(myrank.eq.0) time_transpose_end = MPI_WTIME()
!!$    if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "z to x (communication)", &
!!$         timestep,dt,time,time_transpose_end-time_transpose_start

!!$    if(myrank.eq.0) time_transpose_start = MPI_WTIME()
    call on_node_reorder_z_to_x(nxh,nz32,ny,nprocs1DA,nprocs1DB,BUF_PAD,U_BIG,omp_reorder)
!!$    if(myrank.eq.0) time_transpose_end = MPI_WTIME()
!!$    if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "z to x (re-shuffle)", &
!!$         timestep,dt,time,time_transpose_end-time_transpose_start
    return
  end subroutine ppdfftw_z_to_x
  !-----------------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------------
  ! This subroutine performs forward DFT in Z direction and packs the data from Upad to U
  !-----------------------------------------------------------------------------------------
  subroutine ppdfftw_dft_z_forward(U_BIG,U_SMALL,buf)
    implicit none
    complex(C_DOUBLE_COMPLEX) :: U_BIG(nz32,Zxsz*Zysz)
    complex(C_DOUBLE_COMPLEX) :: U_SMALL(nz,Zxsz*Zysz)
    complex(C_DOUBLE_COMPLEX) :: buf(nz32,omp_z)
    integer :: i,j,ntotal,nh,l

    nh = nz/2
    ntotal = Zxsz*Zysz

    call omp_set_num_threads(omp_compute)

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,l)
    l = omp_get_thread_num() + 1
    !$OMP DO SCHEDULE(STATIC)
    do i = 1,ntotal
       call fftw_execute_dft(plan_transform_forward_z,U_BIG(:,i),buf(:,l))
       do j=1,nh
          U_SMALL(j,i) = buf(j,l)
       enddo
       U_SMALL(nh+1,i) = (0d0,0d0)
       do j=nh*2+2,nz32
          U_SMALL(j-nh,i) = buf(j,l)
       enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    return
  end subroutine ppdfftw_dft_z_forward
  !-----------------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------------
  ! This subroutine performs forward DFT in Z direction and pads the data from U to Upad
  !-----------------------------------------------------------------------------------------
  subroutine ppdfftw_dft_z_inverse(U_SMALL,U_BIG,buf)
    implicit none
    complex(C_DOUBLE_COMPLEX) :: U_BIG(nz32,Zxsz*Zysz)
    complex(C_DOUBLE_COMPLEX) :: U_SMALL(nz,Zxsz*Zysz)
    complex(C_DOUBLE_COMPLEX) :: buf(nz32,omp_z)
    integer :: i,j,ntotal,nh,l

    nh = nz/2
    ntotal = Zxsz*Zysz

    call omp_set_num_threads(omp_compute)

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,l)
    l = omp_get_thread_num() + 1
    !$OMP DO SCHEDULE(STATIC)
    do i = 1,ntotal
       do j=1,nh
          buf(j,l) = U_SMALL(j,i)
       enddo
       do j=nh+1,nh*2+1
          buf(j,l) = (0d0,0d0)
       enddo
       do j=nh*2+2,nz32
          buf(j,l) = U_SMALL(j-nh,i)
       enddo
       call fftw_execute_dft(plan_transform_inverse_z,buf(:,l),U_BIG(:,i))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    return
  end subroutine ppdfftw_dft_z_inverse
  !-----------------------------------------------------------------------------------------


end module ppdfftw
!-----------------------------End of Module-----------------------------------------------



!--------------------------------------------------------------------------------------------------
! ADDITIONAL NECESSARY SUBROUINES
!--------------------------------------------------------------------------------------------------
! Following Subroutines are neccesary after global data transpose using FFTW-MPI
! Those are not included in module ppdfft because of the convenience of type/shape conversion.
! Conceptually, those are part of ppdfftw module and can't be seperated with ppdfftw module.
! 
! Situation : After FFTW_MPI_TRANSPOSE, data is safely moved to the each tasks                    
!             However, due to multidimensional MPI Topology, Data is shuffled in contiguous memory
! Goal : Re-arrange the data in x/z pencil directions
! Original Data Type : 1D contiguous Memory - Complex Type
! Data structure is declared in 4D array for convenience
!--------------------------------------------------------------------------------------------------
subroutine on_node_reorder_x_to_z(n1,n2,n3,nproc1,nproc2,A,B,omp_thread)
  use omp_lib
  use, intrinsic :: iso_c_binding 
  implicit none
  integer(C_INTPTR_T) :: n1,n2,n3
  integer :: i,j,k,l,na,nb,nc,nd,nproc1,nproc2,omp_thread
  complex(C_DOUBLE_COMPLEX) :: A(n1/nproc1, n2/nproc1, n3/nproc2,    nproc1)
  complex(C_DOUBLE_COMPLEX) :: B(n2/nproc1,    nproc1, n3/nproc2, n1/nproc1)

  ! This is needed to ensure compatibility of ISO_C_BINDING with some compilers (e.g. XL compiler)
  na = n2/nproc1
  nb = nproc1
  nc = n3/nproc2
  nd = n1/nproc1

  call omp_set_num_threads(omp_thread)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l)
  !$OMP DO SCHEDULE(STATIC) 
  do j = 1,nb
     do k = 1,nc
        do i = 1,na
           do l = 1,nd
              B(i,j,k,l) = A(l,i,k,j)
           enddo
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  return
end subroutine on_node_reorder_x_to_z

!---------------------------------------------------------------------------------------------------
subroutine on_node_preorder_x_to_z(n1,n2,n3,nproc1,nproc2,A,B,omp_thread)
  use omp_lib
  use, intrinsic :: iso_c_binding 
  implicit none
  integer(C_INTPTR_T) :: n1,n2,n3
  integer :: i,j,k,na,nb,nc,nproc1,nproc2, omp_thread
  complex(C_DOUBLE_COMPLEX) :: A(n1/nproc1, nproc1, n2/nproc1, n3/nproc2)
  complex(C_DOUBLE_COMPLEX) :: B(n1/nproc1, n2/nproc1, n3/nproc2, nproc1)

  ! This is needed to ensure compatibility of ISO_C_BINDING with some compilers (e.g. XL compiler)
  na = nproc1
  nb = n2/nproc1
  nc = n3/nproc2

  call omp_set_num_threads(omp_thread)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
  !$OMP DO SCHEDULE(STATIC) 
  do j = 1,nb
     do k = 1,nc
        do i = 1,na
           B(:,j,k,i) = A(:,i,j,k)
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  return
end subroutine on_node_preorder_x_to_z

!---------------------------------------------------------------------------------------------------
subroutine on_node_reorder_z_to_y(n1,n2,n3,nproc1,nproc2,A,B,omp_thread)
  use omp_lib
  use, intrinsic :: iso_c_binding 
  implicit none

  integer(C_INTPTR_T) :: n1,n2,n3
  integer :: i,j,k,l,na,nb,nc,nd,nproc1,nproc2,omp_thread
  complex(C_DOUBLE_COMPLEX) :: A(n2/nproc2, n3/nproc2, n1/nproc1,    nproc2)
  complex(C_DOUBLE_COMPLEX) :: B(n3/nproc2,    nproc2, n1/nproc1, n2/nproc2)

  ! This is needed to ensure compatibility of ISO_C_BINDING with some compilers (e.g. XL compiler)
  na = n3/nproc2
  nb = nproc2
  nc = n1/nproc1
  nd = n2/nproc2

  call omp_set_num_threads(omp_thread)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l)
  !$OMP DO SCHEDULE(STATIC) 
  do j = 1,nb
     do k = 1,nc
        do i = 1,na
           do l = 1,nd
              B(i,j,k,l) = A(l,i,k,j)
           enddo
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  return
end subroutine on_node_reorder_z_to_y

!---------------------------------------------------------------------------------------------------
subroutine on_node_preorder_z_to_y(n1,n2,n3,nproc1,nproc2,A,B,omp_thread)
  use omp_lib
  use, intrinsic :: iso_c_binding 
  implicit none

  integer(C_INTPTR_T) :: n1,n2,n3
  integer :: i,j,k,na,nb,nc,nproc1,nproc2,omp_thread
  complex(C_DOUBLE_COMPLEX) :: A(n2/nproc2,    nproc2, n3/nproc2, n1/nproc1)
  complex(C_DOUBLE_COMPLEX) :: B(n2/nproc2, n3/nproc2, n1/nproc1,    nproc2)

  ! This is needed to ensure compatibility of ISO_C_BINDING with some compilers (e.g. XL compiler)
  na = nproc2
  nb = n3/nproc2
  nc = n1/nproc1

  call omp_set_num_threads(omp_thread)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
  !$OMP DO SCHEDULE(STATIC) 
  do k = 1,nc
     do j = 1,nb
        do i = 1,na
           B(:,j,k,i) = A(:,i,j,k)
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  return
end subroutine on_node_preorder_z_to_y

!---------------------------------------------------------------------------------------------------
subroutine on_node_reorder_y_to_z(n1,n2,n3,nproc1,nproc2,A,B,omp_thread)
  use omp_lib
  use, intrinsic :: iso_c_binding 
  implicit none
  integer(C_INTPTR_T) :: n1, n2, n3
  integer :: na,nb,nc,i,j,k,l, nproc1, nproc2,omp_thread
  complex(C_DOUBLE_COMPLEX) :: A(n3/nproc2,n1/nproc1,n2)
  complex(C_DOUBLE_COMPLEX) :: B(n2,n3/nproc2,n1/nproc1)
  
  ! This is needed to ensure compatibility of ISO_C_BINDING with some compilers (e.g. XL compiler)
  na = n2
  nb = n3/nproc2
  nc = n1/nproc1

  call omp_set_num_threads(omp_thread)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
  !$OMP DO SCHEDULE(STATIC) 
  do i=1,na
     do k=1,nc
        do j=1,nb
           B(i,j,k) = A(j,k,i)
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  return
end subroutine on_node_reorder_y_to_z

!---------------------------------------------------------------------------------------------------
subroutine on_node_preorder_y_to_z(n1,n2,n3,nproc1,nproc2,A,B,omp_thread)
  use omp_lib
  use, intrinsic :: iso_c_binding 
  implicit none
  integer(C_INTPTR_T) :: n1, n2, n3
  integer :: na,nb,nc,i,j,k,l, nproc1, nproc2,omp_thread
  complex(C_DOUBLE_COMPLEX) :: A(n3/nproc2,nproc2,n1/nproc1,n2/nproc2)
  complex(C_DOUBLE_COMPLEX) :: B(n3/nproc2,n1/nproc1,n2/nproc2,nproc2)
  
  ! This is needed to ensure compatibility of ISO_C_BINDING with some compilers (e.g. XL compiler)
  na = nproc2
  nb = n1/nproc1
  nc = n2/nproc2

  call omp_set_num_threads(omp_thread)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
  !$OMP DO SCHEDULE(STATIC) 
     do j=1,nb
  do k=1,nc
        do i=1,na
           B(:,j,k,i) = A(:,i,j,k)
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  return
end subroutine on_node_preorder_y_to_z

!---------------------------------------------------------------------------------------------------
subroutine on_node_reorder_z_to_x(n1,n2,n3,nproc1,nproc2,A,B,omp_thread)
  use omp_lib
  use, intrinsic :: iso_c_binding 
  implicit none

  integer(C_INTPTR_T) :: n1, n2, n3
  integer :: na,nb,nc,i,j,k,l, nproc1, nproc2,omp_thread
  complex(C_DOUBLE_COMPLEX) :: A(n2/nproc1, n3/nproc2, n1)
  complex(C_DOUBLE_COMPLEX) :: B(n1       , n2/nproc1, n3/nproc2)

  ! This is needed to ensure compatibility of ISO_C_BINDING with some compilers (e.g. XL compiler)
  nc = n3/nproc2
  nb = n2/nproc1
  na = n1

  call omp_set_num_threads(omp_thread)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
  !$OMP DO SCHEDULE(STATIC) 
  do i=1,na
     do k=1,nc
        do j=1,nb
           B(i,j,k) = A(j,k,i)
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  return
end subroutine on_node_reorder_z_to_x

!---------------------------------------------------------------------------------------------------
subroutine on_node_preorder_z_to_x(n1,n2,n3,nproc1,nproc2,A,B,omp_thread)
  use omp_lib
  use, intrinsic :: iso_c_binding 
  implicit none

  integer(C_INTPTR_T) :: n1, n2, n3
  integer :: na,nb,nc,i,j,k,l, nproc1, nproc2,omp_thread
  complex(C_DOUBLE_COMPLEX) :: A(n2/nproc1, nproc1, n3/nproc2, n1/nproc1)
  complex(C_DOUBLE_COMPLEX) :: B(n2/nproc1, n3/nproc2, n1/nproc1, nproc1)

  ! This is needed to ensure compatibility of ISO_C_BINDING with some compilers (e.g. XL compiler)
  nc = n1/nproc1
  nb = n3/nproc2
  na = nproc1

  call omp_set_num_threads(omp_thread)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
  !$OMP DO SCHEDULE(STATIC) 
  do k=1,nc
     do j=1,nb
        do i=1,na
           B(:,j,k,i) = A(:,i,j,k)
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  return
end subroutine on_node_preorder_z_to_x

!---------------------------------------------------------------------------------------------------
