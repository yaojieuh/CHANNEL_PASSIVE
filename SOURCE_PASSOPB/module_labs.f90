!--------------------------------------------------------------------
! Linear Algebra for B-Spline (LABS)
!
! This module contains variables and subroutine to use B-Spline
!
! The matrices to calculate Zeroth/first/second order derivatives
! are calculated in CPP file called 'labs_fill_full.cpp
!
! The full size matrices are saved in modified banded matirx form
!
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------

module LABS
  use,intrinsic :: iso_c_binding
  implicit none

  integer :: LABS_K
  integer :: LABS_BW
  integer :: LABS_NY
  integer :: WRITE_LABS_INFO
  real(c_double) :: stretch
  real(c_double), allocatable, dimension(:)   :: LABS_Y
  real(c_double), allocatable, dimension(:)   :: LABS_Y_KNOT
  real(c_double), allocatable, dimension(:)   :: LABS_PRECONDITION
  real(c_double), allocatable, dimension(:)   :: LABS_CENTER_LINE
  real(c_double), allocatable, dimension(:)   :: LABS_INTEGRATE
  real(c_double), allocatable, dimension(:,:) :: LABS_D0_FULL
  real(c_double), allocatable, dimension(:,:) :: LABS_D1_FULL
  real(c_double), allocatable, dimension(:,:) :: LABS_D2_FULL
  real(c_double), allocatable, dimension(:,:) :: LABS_D3_FULL
  real(c_double), allocatable, dimension(:,:) :: LABS_D4_FULL
  real(c_double), allocatable, dimension(:,:,:) :: LABS_D0_B
  real(c_double), allocatable, dimension(:,:,:) :: LABS_D1_B
  real(c_double), allocatable, dimension(:,:,:) :: LABS_D2_B
  real(c_double), allocatable, dimension(:,:,:) :: LABS_BUF
  real(c_double), allocatable, dimension(:,:,:) :: LABS_D0_LU
  real(c_double), allocatable, dimension(:,:,:) :: LABS_D2_LU
  
!  public :: LABS_initialize, LABS_D1, LABS_D2, LABS_lu, LABS_bs, LABS_mat_mul, LABS_integration, LABS_finalize

  interface 
     subroutine LABS_fill_full(LABS_D0_FULL,LABS_D1_FULL,LABS_D2_FULL,LABS_D3_FULL,LABS_D4_FULL,&
          Y,Y_knot,LABS_INTEGRATE,LABS_CENTER_LINE,LABS_K,LABS_NY,stretch,Ly) bind(C, name="LABS_fill_full")
       use,intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: LABS_NY
       real(c_double) :: LABS_D0_FULL(LABS_NY,LABS_NY)
       real(c_double) :: LABS_D1_FULL(LABS_NY,LABS_NY)
       real(c_double) :: LABS_D2_FULL(LABS_NY,LABS_NY)
       real(c_double) :: LABS_D3_FULL(LABS_NY,LABS_NY)
       real(c_double) :: LABS_D4_FULL(LABS_NY,LABS_NY)
       real(c_double) :: Y(LABS_NY)
       real(c_double) :: Y_knot(LABS_NY-LABS_K+1)
       real(c_double) :: LABS_INTEGRATE(LABS_NY)
       real(c_double) :: LABS_CENTER_LINE(LABS_NY)
       integer(c_int) :: LABS_K
       real(c_double) :: stretch
       real(c_double) :: Ly
     end subroutine LABS_fill_full
  end interface

contains


  !----------------------------------------------------------
  ! This subroutine builds an environment to use LABS
  ! In detail This subroutine does
  ! - allocate necessary data space
  ! - calculate full size matrices for 0/1/2 order derivative
  ! - change full size matrices to modified banded matrices
  !----------------------------------------------------------
  subroutine LABS_initialize(ny,stretch,Ly,omp_compute)
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_intptr_t) :: ny
    real(c_double) :: stretch
    real(c_double) :: Ly
    complex(c_double_complex),allocatable :: dummy(:)
    integer :: i,j,omp_compute



    LABS_NY = ny
    allocate(LABS_Y(LABS_NY))
    allocate(LABS_Y_KNOT(LABS_NY-LABS_K+1))
    allocate(LABS_PRECONDITION(LABS_NY))
    allocate(LABS_INTEGRATE(LABS_NY))
    allocate(LABS_CENTER_LINE(LABS_NY))
    allocate(dummy(LABS_NY))
    allocate(LABS_D0_FULL(LABS_NY,LABS_NY))
    allocate(LABS_D1_FULL(LABS_NY,LABS_NY))
    allocate(LABS_D2_FULL(LABS_NY,LABS_NY))
    allocate(LABS_D3_FULL(LABS_NY,LABS_NY))
    allocate(LABS_D4_FULL(LABS_NY,LABS_NY))

    ! get the full size matrix for B-Spline from GSL
    call LABS_fill_full(LABS_D0_FULL,LABS_D1_FULL,LABS_D2_FULL,LABS_D3_FULL,LABS_D4_FULL, &
         LABS_Y,LABS_Y_KNOT,LABS_INTEGRATE,LABS_CENTER_LINE,LABS_K+1,LABS_NY,stretch,Ly)

    ! Calculate the bandwidth of full matrix
    call LABS_eval_bw

    ! change the coefficient matrix from full size to banded form
    allocate(LABS_D0_B(2*LABS_BW+1,LABS_NY,OMP_COMPUTE))
    allocate(LABS_D1_B(2*LABS_BW+1,LABS_NY,OMP_COMPUTE))
    allocate(LABS_D2_B(2*LABS_BW+1,LABS_NY,OMP_COMPUTE))
    call LABS_full_to_band(LABS_D0_FULL,LABS_D0_B(:,:,1),LABS_BW,LABS_NY)
    call LABS_full_to_band(LABS_D1_FULL,LABS_D1_B(:,:,1),LABS_BW,LABS_NY)
    call LABS_full_to_band(LABS_D2_FULL,LABS_D2_B(:,:,1),LABS_BW,LABS_NY)


    if (WRITE_LABS_INFO.ne.0) then
       open(28,file='LABS_D0_FULL.DAT')
       do j=1,ny
          do i=1,ny-1
             write(28,'(ES30.16)',advance='no') LABS_D0_FULL(i,j)
          enddo
          write(28,'(ES30.16)') LABS_D0_FULL(ny,j)
       enddo
       close(28)    

       open(28,file='LABS_D1_FULL.DAT')
       do j=1,ny
          do i=1,ny-1
             write(28,'(ES30.16)',advance='no') LABS_D1_FULL(i,j)
          enddo
          write(28,'(ES30.16)') LABS_D1_FULL(ny,j)
       enddo
       close(28)    

       open(28,file='LABS_D2_FULL.DAT')
       do j=1,ny
          do i=1,ny-1
             write(28,'(ES30.16)',advance='no') LABS_D2_FULL(i,j)
          enddo
          write(28,'(ES30.16)') LABS_D2_FULL(ny,j)
       enddo
       close(28)    

       open(28,file='LABS_D3_FULL.DAT')
       do j=1,ny
          do i=1,ny-1
             write(28,'(ES30.16)',advance='no') LABS_D3_FULL(i,j)
          enddo
          write(28,'(ES30.16)') LABS_D3_FULL(ny,j)
       enddo
       close(28)    

       open(28,file='LABS_D4_FULL.DAT')
       do j=1,ny
          do i=1,ny-1
             write(28,'(ES30.16)',advance='no') LABS_D4_FULL(i,j)
          enddo
          write(28,'(ES30.16)') LABS_D4_FULL(ny,j)
       enddo
       close(28)    


       open(28,file='LABS_D0_BAND.DAT')
       do j=1,ny
          do i=1,2*LABS_BW
             write(28,'(ES30.16)',advance='no') LABS_D0_B(i,j,1)
          enddo
          write(28,'(ES30.16)') LABS_D0_B(2*LABS_BW+1,ny,1)
       enddo
       close(28)    

       open(28,file='LABS_D1_BAND.DAT')
       do j=1,ny
          do i=1,2*LABS_BW
             write(28,'(ES30.16)',advance='no') LABS_D1_B(i,j,1)
          enddo
          write(28,'(ES30.16)') LABS_D1_B(2*LABS_BW+1,ny,1)
       enddo
       close(28)    

       open(28,file='LABS_D2_BAND.DAT')
       do j=1,ny
          do i=1,2*LABS_BW
             write(28,'(ES30.16)',advance='no') LABS_D2_B(i,j,1)
          enddo
          write(28,'(ES30.16)') LABS_D2_B(2*LABS_BW+1,ny,1)
       enddo
       close(28)    

       open(28,file='LABS_Y_LOCATION.DAT')
       do j=1,ny
          write(28,'(ES30.16)') LABS_Y(j)
       enddo
       close(28)    

       open(28,file='LABS_INTEGRATE_COEF.DAT')
       do j=1,ny
          write(28,'(ES30.16)') LABS_INTEGRATE(j)
       enddo

       open(28,file='LABS_CENTER_LINE.DAT')
       do j=1,ny
          write(28,'(ES30.16)') LABS_CENTER_LINE(j)
       enddo

       close(28)    
    endif


    dummy(1) = (0.d0,0.d0)
    do i=2,LABS_NY-1
       dummy(i) = (1.d0,0.d0)
    enddo
    dummy(LABS_NY) = (0.d0,0.d0)

    allocate(LABS_D0_LU(2*LABS_BW+1,LABS_NY,OMP_COMPUTE))
    allocate(LABS_D2_LU(2*LABS_BW+1,LABS_NY,OMP_COMPUTE))
    allocate(LABS_BUF(2*LABS_BW+1,LABS_NY,OMP_COMPUTE))
    LABS_D0_LU(:,:,1) = LABS_D0_B(:,:,1)
    LABS_D2_LU(:,:,1) = LABS_D2_B(:,:,1)
    call LABS_lu(LABS_D0_LU(:,:,1))
    call LABS_lu(LABS_D2_LU(:,:,1))

    call LABS_bs(LABS_D0_LU(:,:,1),dummy)

    do i=1,OMP_COMPUTE
       LABS_D0_B(:,:,i) = LABS_D0_B(:,:,1)
       LABS_D1_B(:,:,i) = LABS_D1_B(:,:,1)
       LABS_D2_B(:,:,i) = LABS_D2_B(:,:,1)
       LABS_D0_LU(:,:,i) = LABS_D0_LU(:,:,1)
       LABS_D2_LU(:,:,i) = LABS_D2_LU(:,:,1)
    enddo

    deallocate(dummy)

    return
  end subroutine LABS_initialize


  !-------------------------------------------------------------------
  ! This subroutine calculate the bandwidth of full matrix for B-Spline
  ! The bandwidth will be used through entire LABS module 
  !-------------------------------------------------------------------
  subroutine LABS_eval_bw
    implicit none
    integer :: i,j,k
    real(c_double) :: error

    error =0.d0
    k=1

    do
       do j=1,k+1
          do i=2*k+2,LABS_NY
             error = max(error,abs(LABS_D0_FULL(i,j)))
             error = max(error,abs(LABS_D1_FULL(i,j)))
             error = max(error,abs(LABS_D2_FULL(i,j)))
          enddo
       enddo

       do j=k+2,LABS_NY-k-1
          do i=1,j-k-1
             error = max(error,abs(LABS_D0_FULL(i,j)))
             error = max(error,abs(LABS_D1_FULL(i,j)))
             error = max(error,abs(LABS_D2_FULL(i,j)))
          enddo

          do i=j+k+1,LABS_NY
             error = max(error,abs(LABS_D0_FULL(i,j)))
             error = max(error,abs(LABS_D1_FULL(i,j)))
             error = max(error,abs(LABS_D2_FULL(i,j)))
          enddo
       enddo

       do j=LABS_NY-k ,LABS_NY
          do i=1,LABS_NY-(2*k)-1
             error = max(error,abs(LABS_D0_FULL(i,j)))
             error = max(error,abs(LABS_D1_FULL(i,j)))
             error = max(error,abs(LABS_D2_FULL(i,j)))
          enddo
       enddo

       if (error.le.0d0) then
          LABS_BW = k
          exit
       else
          error = 0.d0
          k=k+1
       end if
    enddo
    return
  end subroutine LABS_eval_bw



  !-------------------------------------------------------------------
  ! This subroutine change the form of matrix for B-Spline
  ! from original densed matrix to near boundary packed banded matrix
  !-------------------------------------------------------------------
  subroutine LABS_full_to_band(FULL,BAND,BW,N)
    integer :: BW,N
    real(c_double) :: FULL(N,N)
    real(c_double) :: BAND(2*BW+1,N)
    integer :: i,j
    do j=1,BW
       do i=1,BW+j
          BAND(BW+1+i-j,j) = FULL(i,j)
       enddo
       do i=j+BW+1,2*BW+1
          BAND(2*BW+2-i,j) = FULL(i,j)
       enddo
    enddo
    do j=BW+1,N-BW
       do i=1,2*BW+1
          BAND(i,j) = FULL(j+i-BW-1,j)
       enddo
    enddo
    do j=N-BW+1,N
       do i=j-BW,N
          BAND(BW+1+i-j,j) = FULL(i,j)
       enddo
       do i=N-2*BW,j-BW-1
          BAND(N-i+1,j) = FULL(i,j)
       enddo
    enddo
    return
  end subroutine LABS_full_to_band


  !-------------------------------------------------------
  ! This subroutine calculate the first order derivative
  ! 
  ! input : dZ0 - complex vector with LABS_NY length
  ! output : dZ1 - complex vector with LABS_NY length
  !-------------------------------------------------------
  subroutine LABS_d1(dZ0,dZ1)
    complex(c_double_complex) :: dZ0(LABS_NY)
    complex(c_double_complex) :: dZ1(LABS_NY)
    complex(c_double_complex) :: c(LABS_NY)

    c = dZ0
    call LABS_bs(LABS_D0_lu(:,:,1),c)
    call LABS_mat_mul(LABS_D1_B(:,:,1),c,dZ1)

    return
  end subroutine LABS_d1

  subroutine LABS_d1_r(dZ0,dZ1)
    real(c_double) :: dZ0(LABS_NY)
    real(c_double) :: dZ1(LABS_NY)
    real(c_double) :: c(LABS_NY)

    c = dZ0
    call LABS_bs_r(LABS_D0_lu(:,:,1),c)
    call LABS_mat_mul_r(LABS_D1_B(:,:,1),c,dZ1)

    return
  end subroutine LABS_d1_r


  !-------------------------------------------------------
  ! This subroutine calculate the second order derivative
  ! 
  ! input : dZ0 - complex vector with LABS_NY length
  ! output : dZ2 - complex vector with LABS_NY length
  !-------------------------------------------------------
  subroutine LABS_d2(dZ0,dZ2)
    complex(c_double_complex) :: dZ0(LABS_NY)
    complex(c_double_complex) :: dZ2(LABS_NY)
    complex(c_double_complex) :: c(LABS_NY)

    c(:) = dZ0(:)
    call LABS_bs(LABS_D0_lu(:,:,1),c)
    call LABS_mat_mul(LABS_D2_B(:,:,1),c,dZ2)

    return
  end subroutine LABS_d2


  subroutine LABS_d2_r(dZ0,dZ2)
    real(c_double) :: dZ0(LABS_NY)
    real(c_double) :: dZ2(LABS_NY)
    real(c_double) :: c(LABS_NY)

    c(:) = dZ0(:)
    call LABS_bs_r(LABS_D0_lu(:,:,1),c)
    call LABS_mat_mul_r(LABS_D2_B(:,:,1),c,dZ2)

    return
  end subroutine LABS_d2_r



  !-------------------------------------------------------
  ! This subroutine calculate the integration of vector in given domain
  ! 
  ! input : src - complex vector with LABS_NY length
  ! output : tgt - complex number (result)
  !-------------------------------------------------------
  subroutine LABS_integration(src,tgt)
    use, intrinsic :: iso_c_binding
    implicit none
    complex(c_double_complex) :: src(LABS_NY)
    complex(c_double_complex) :: tgt
    integer :: i
    
    tgt = dcmplx(0.d0,0.d0)
    do i=1,LABS_NY
       tgt = tgt + LABS_INTEGRATE(i) * src(i) 
    enddo
    return
  end subroutine LABS_integration

  subroutine LABS_integration_r(src,tgt)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double) :: src(LABS_NY)
    real(c_double) :: tgt
    integer :: i
    
    tgt = 0d0
    do i=1,LABS_NY
       tgt = tgt + LABS_INTEGRATE(i) * src(i) 
    enddo
    return
  end subroutine LABS_integration_r


  !----------------------------------------------------------
  ! This subroutine destroy the environment to use LABS
  ! In detail This subroutine does
  ! - deallocate allocated data space
  !----------------------------------------------------------
  subroutine LABS_finalize
    deallocate(LABS_Y)
    deallocate(LABS_PRECONDITION)
    deallocate(LABS_INTEGRATE)
    deallocate(LABS_CENTER_LINE)
    deallocate(LABS_D0_B)
    deallocate(LABS_D1_B)
    deallocate(LABS_D2_B) 
    deallocate(LABS_BUF)
    deallocate(LABS_D0_LU)
    deallocate(LABS_D2_LU)
    deallocate(LABS_D0_FULL)
    deallocate(LABS_D1_FULL)
    deallocate(LABS_D2_FULL)
    deallocate(LABS_D3_FULL)
    deallocate(LABS_D4_FULL)
    deallocate(LABS_Y_KNOT)
    return
  end subroutine LABS_finalize




  !-------------------------------------------------------------------
  ! This subroutine perform the LU decomposition
  ! which is specialized for B-Spline Packed Banded Matrix
  !-------------------------------------------------------------------
  subroutine LABS_lu(A)
    implicit none
    real(c_double) :: A(2*LABS_BW+1,LABS_NY)

    select case (LABS_BW)
    case (1)
       call fixed_bandec1(A,LABS_NY)
    case (2)
       call fixed_bandec2(A,LABS_NY)
    case (3)
       call fixed_bandec3(A,LABS_NY)
    case (4)
       call fixed_bandec4(A,LABS_NY)
    case (5)
       call fixed_bandec5(A,LABS_NY)
    case (6)
       call fixed_bandec6(A,LABS_NY)
    case (7)
       call fixed_bandec7(A,LABS_NY)
    case (8)
       call fixed_bandec8(A,LABS_NY)
    case (9)
       call fixed_bandec9(A,LABS_NY)
    case (10)
       call fixed_bandec10(A,LABS_NY)
    case default
       call fixed_bandec(A,LABS_NY,LABS_BW)
    end select

    return
  end subroutine LABS_lu


  !-------------------------------------------------------------------
  ! This subroutine perform the Back substitution 
  ! from the Packed Banded Matrix which is specialized for B-Spline 
  !-------------------------------------------------------------------
  subroutine LABS_bs(A,b)
    implicit none
    real(c_double) :: A(2*LABS_BW+1,LABS_NY)
    complex(c_double_complex) :: b(LABS_NY)

    select case (LABS_BW)
    case (1)
       call fixed_banbks1(A,LABS_NY,b)
    case (2)
       call fixed_banbks2(A,LABS_NY,b)
    case (3)
       call fixed_banbks3(A,LABS_NY,b)
    case (4)
       call fixed_banbks4(A,LABS_NY,b)
    case (5)
       call fixed_banbks5(A,LABS_NY,b)
    case (6)
       call fixed_banbks6(A,LABS_NY,b)
    case (7)
       call fixed_banbks7(A,LABS_NY,b)
    case (8)
       call fixed_banbks8(A,LABS_NY,b)
    case (9)
       call fixed_banbks9(A,LABS_NY,b)
    case (10)
       call fixed_banbks10(A,LABS_NY,b)
    case default
       call fixed_banbks(A,LABS_NY,b,LABS_BW)
    end select
    return
  end subroutine LABS_bs

  subroutine LABS_bs_r(A,b)
    implicit none
    real(c_double) :: A(2*LABS_BW+1,LABS_NY)
    real(c_double) :: b(LABS_NY)

    select case (LABS_BW)
    case (1)
       call fixed_banbks1_r(A,LABS_NY,b)
    case (2)
       call fixed_banbks2_r(A,LABS_NY,b)
    case (3)
       call fixed_banbks3_r(A,LABS_NY,b)
    case (4)
       call fixed_banbks4_r(A,LABS_NY,b)
    case (5)
       call fixed_banbks5_r(A,LABS_NY,b)
    case (6)
       call fixed_banbks6_r(A,LABS_NY,b)
    case (7)
       call fixed_banbks7_r(A,LABS_NY,b)
    case (8)
       call fixed_banbks8_r(A,LABS_NY,b)
    case (9)
       call fixed_banbks9_r(A,LABS_NY,b)
    case (10)
       call fixed_banbks10_r(A,LABS_NY,b)
    case default
       call fixed_banbks_r(A,LABS_NY,b,LABS_BW)
    end select
    return
  end subroutine LABS_bs_r



  !-------------------------------------------------------------------
  ! This subroutine perform the Matrix-vector Multiplication
  ! which is specialized for B-Spline Packed Banded Matrix
  ! Caution : b = Ax   (x and b cannot be in same memory space) 
  !-------------------------------------------------------------------
  subroutine LABS_mat_mul(A,x,b)
    implicit none
    real(c_double) :: A(2*LABS_BW+1,LABS_NY)
    complex(c_double_complex) :: x(LABS_NY), b(LABS_NY)
    integer :: i,j,k

    select case (LABS_BW)
    case (1) 
       call LABS_mat_mul_1(A,x,b,LABS_NY)
    case (2) 
       call LABS_mat_mul_2(A,x,b,LABS_NY)
    case (3) 
       call LABS_mat_mul_3(A,x,b,LABS_NY)
    case (4) 
       call LABS_mat_mul_4(A,x,b,LABS_NY)
    case (5) 
       call LABS_mat_mul_5(A,x,b,LABS_NY)
    case (6) 
       call LABS_mat_mul_6(A,x,b,LABS_NY)
    case (7) 
       call LABS_mat_mul_7(A,x,b,LABS_NY)
    case (8) 
       call LABS_mat_mul_8(A,x,b,LABS_NY)
    case (9) 
       call LABS_mat_mul_9(A,x,b,LABS_NY)
    case (10) 
       call LABS_mat_mul_10(A,x,b,LABS_NY)
    case default
       call LABS_mat_mul_bw(A,x,b,LABS_NY,LABS_BW)
    end select
    return
  end subroutine LABS_mat_mul

  subroutine LABS_mat_mul_r(A,x,b)
    implicit none
    real(c_double) :: A(2*LABS_BW+1,LABS_NY)
    real(c_double) :: x(LABS_NY), b(LABS_NY)
    integer :: i,j,k

    select case (LABS_BW)
    case (1) 
       call LABS_mat_mul_1_r(A,x,b,LABS_NY)
    case (2) 
       call LABS_mat_mul_2_r(A,x,b,LABS_NY)
    case (3) 
       call LABS_mat_mul_3_r(A,x,b,LABS_NY)
    case (4) 
       call LABS_mat_mul_4_r(A,x,b,LABS_NY)
    case (5) 
       call LABS_mat_mul_5_r(A,x,b,LABS_NY)
    case (6) 
       call LABS_mat_mul_6_r(A,x,b,LABS_NY)
    case (7) 
       call LABS_mat_mul_7_r(A,x,b,LABS_NY)
    case (8) 
       call LABS_mat_mul_8_r(A,x,b,LABS_NY)
    case (9) 
       call LABS_mat_mul_9_r(A,x,b,LABS_NY)
    case (10) 
       call LABS_mat_mul_10_r(A,x,b,LABS_NY)
    case default
       call LABS_mat_mul_bw_r(A,x,b,LABS_NY,LABS_BW)
    end select
    return
  end subroutine LABS_mat_mul_r




!----------------------------------------------------------------------
! 
!  W A R N I N G !!!!!!
!  W A R N I N G !!!!!!
!  W A R N I N G !!!!!!
!
! From this point the routine is hand unrolled upto BW=7
! For BW=8 and after the routine is not unrolled manually.
! 
! Just, do not change from this point to the end of this module. 
!
!----------------------------------------------------------------------
  subroutine fixed_bandec1(A,n)
    implicit none
    integer :: n
    real(c_double) :: A(3,n)
    integer :: i
    integer :: m,bw
    real(c_double) :: dummy

    m = 2
    bw = 1

    dummy = 1.d0 / A(m,1)
    A(m+1,1) = dummy * A(m+1,1)
    A(m-1,1) = dummy * A(m-1,1)

    A(m  ,2) = A(m  ,2) -  A(m-1,2) * A(m+1,1)

    dummy = 1.d0 / A(m,2)
    A(m+1,2) = dummy * (A(m+1,2) -  A(m-1,2) * A(m-1,1)) 

    do i = 2*bw+1,n-2*bw-1
       A(m  ,i  ) = A(m  ,i  ) - A(m-1,i  ) * A(m+1,i-1)

       dummy = 1.d0 / A(m,i)
       A(m+1,i) = dummy * A(m+1,i)
    enddo

    A(m,n-2*bw)    = A(m,n-2*bw)   - A(m-1,n-2*bw)   * A(m+1,n-2*bw-1)

    dummy = 1.d0 / A(m,n-2*bw)
    A(m+1,n-2*bw) = dummy * A(m+1,n-2*bw)

    A(m,n-2*bw+1)   = A(m,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+1,n-2*bw)
    A(m-1,n-2*bw+2) = A(m-1,n-2*bw+2) - A(m+1,n-2*bw+2)*A(m+1,n-2*bw)  

    dummy = 1.d0 / A(m,n-2*bw+1)
    A(m+1,n-2*bw+1) = dummy * A(m+1,n-2*bw+1)

    A(m  ,n-2*bw+2)  = A(m  ,n-2*bw+2)  - A(m-1,n-2*bw+2)  * A(m+1,n-2*bw+1)

    return

  end subroutine fixed_bandec1


  subroutine fixed_bandec2(A,n)
    implicit none
    integer :: n
    real(c_double) :: A(5,n)
    integer :: i
    integer :: m,bw
    real(c_double) :: dummy

    m = 3
    bw = 2

    dummy = 1.d0 / A(m,1)
    A(m+1,1) = dummy * A(m+1,1)
    A(m+2,1) = dummy * A(m+2,1)
    A(m-1,1) = dummy * A(m-1,1)
    A(m-2,1) = dummy * A(m-2,1)

    A(m  ,2) = A(m  ,2) -  A(m-1,2) * A(m+1,1)
    A(m-1,3) = A(m-1,3) -  A(m-2,3) * A(m+1,1)

    dummy = 1.d0 / A(m,2)
    A(m+1,2) = dummy * (A(m+1,2) -  A(m-1,2) * A(m+2,1)) 
    A(m+2,2) = dummy * (A(m+2,2) -  A(m-1,2) * A(m-1,1)) 
    A(m-2,2) = dummy * (A(m-2,2) -  A(m-1,2) * A(m-2,1)) 

    A(m  , 3) = A(m  ,3) - A(m-1,3) * A(m+1,2) - A(m-2,3) * A(m+2,1)
    A(m-1, 4) = A(m-1,4) - A(m-2,4) * A(m+1,2) 

    dummy = 1.d0 / A(m,3)
    A(m+1,3) = dummy * (A(m+1,3) - A(m-1,3) * A(m+2,2) - A(m-2,3) * A(m-1,1))
    A(m+2,3) = dummy * (A(m+2,3) - A(m-1,3) * A(m-2,2) - A(m-2,3) * A(m-2,1))

    A(m  , 4) = A(m  , 4) - A(m-1, 4) * A(m+1,3) - A(m-2,4) * A(m+2,2)
    A(m-1, 5) = A(m-1, 5) - A(m-2, 5) * A(m+1,3)
    A(m-2, 6) = A(m-2, 6)

    dummy = 1.d0 / A(m,4)
    A(m+1,4) = dummy * (A(m+1,4) - A(m-1,4) * A(m+2,3) - A(m-2,4) * A(m-2,2))
    A(m+2,4) = dummy *  A(m+2,4)

    do i = 2*bw+1,n-2*bw-1
       A(m  ,i  ) = A(m  ,i  ) - A(m-1,i  ) * A(m+1,i-1) - A(m-2,i)   * A(m+2,i-2)
       A(m-1,i+1) = A(m-1,i+1) - A(m-2,i+1) * A(m+1,i-1)
       A(m-2,i+2) = A(m-2,i+2)

       dummy = 1.d0 / A(m,i)
       A(m+1,i) = dummy * (A(m+1,i) - A(m-1,i) * A(m+2,i-1))
       A(m+2,i) = dummy *  A(m+2,i)
    enddo

    A(m  ,n-2*bw)    = A(m  ,n-2*bw)   - A(m-1,n-2*bw)   * A(m+1,n-2*bw-1) - A(m-2,n-2*bw)   * A(m+2,n-2*bw-2)
    A(m-1,n-2*bw+1)  = A(m-1,n-2*bw+1) - A(m-2,n-2*bw+1) * A(m+1,n-2*bw-1)
    A(m-2,n-2*bw+2)  = A(m-2,n-2*bw+2)

    dummy = 1.d0 / A(m,n-2*bw)
    A(m+1,n-2*bw) = dummy * (A(m+1,n-2*bw) - A(m-1,n-2*bw) * A(m+2,n-2*bw-1))
    A(m+2,n-2*bw) = dummy * A(m+2,n-2*bw)

    A(m  ,n-2*bw+1)   = A(m  ,n-2*bw+1)  - A(m-1,n-2*bw+1)  * A(m+1,n-2*bw) - A(m-2,n-2*bw+1) * A(m+2,n-2*bw-1)
    A(m-1,n-2*bw+2)   = A(m-1,n-2*bw+2)  - A(m-2,n-2*bw+2)  * A(m+1,n-2*bw)
    A(m-2,n-2*bw+3)   = A(m-2,n-2*bw+3)  - A(m+2,n-2*bw+3)  * A(m+1,n-2*bw) 
    A(m+1,n-2*bw+4)   = A(m+1,n-2*bw+4)  - A(m+2,n-2*bw+4)  * A(m+1,n-2*bw) 

    dummy = 1.d0 / A(m,n-2*bw+1)
    A(m+1,n-2*bw+1) = dummy * (A(m+1,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+2,n-2*bw))
    A(m+2,n-2*bw+1) = dummy *  A(m+2,n-2*bw+1)

    A(m  ,n-2*bw+2)  = A(m  ,n-2*bw+2)  - A(m-1,n-2*bw+2)  * A(m+1,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+2,n-2*bw)
    A(m-1,n-2*bw+3)  = A(m-1,n-2*bw+3)  - A(m-2,n-2*bw+3)  * A(m+1,n-2*bw+1) - A(m+2,n-2*bw+3) * A(m+2,n-2*bw)
    A(m-2,n-2*bw+4)  = A(m-2,n-2*bw+4)  - A(m+1,n-2*bw+4)  * A(m+1,n-2*bw+1) - A(m+2,n-2*bw+4) * A(m+2,n-2*bw)

    dummy = 1.d0 / A(m,n-2*bw+2)
    A(m+1,n-2*bw+2) = dummy * (A(m+1,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+2,n-2*bw+1))
    A(m+2,n-2*bw+2) = dummy *  A(m+2,n-2*bw+2)

    A(m  ,n-2*bw+3) = A(m  ,n-2*bw+3)  - A(m-1,n-2*bw+3)  * A(m+1,n-2*bw+2) - A(m-2,n-2*bw+3)  * A(m+2,n-2*bw+1)
    A(m-1,n-2*bw+4) = A(m-1,n-2*bw+4)  - A(m-2,n-2*bw+4)  * A(m+1,n-2*bw+2) - A(m+1,n-2*bw+4)  * A(m+2,n-2*bw+1)


    dummy = 1.d0 / A(m,n-2*bw+3)
    A(m+1,n-2*bw+3) = dummy * (A(m+1,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+2,n-2*bw+2))

    A(m  ,n-2*bw+4) = A(m  ,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+1,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+2,n-2*bw+2)

    return

  end subroutine fixed_bandec2


  subroutine fixed_bandec3(A,n)
    implicit none
    integer :: n
    real(c_double) :: A(7,n)
    integer :: i
    integer :: m,bw
    real(c_double) :: dummy

    m = 4
    bw = 3

    dummy = 1.d0 / A(m,1)
    A(m+1,1) = dummy * A(m+1,1)
    A(m+2,1) = dummy * A(m+2,1)
    A(m+3,1) = dummy * A(m+3,1)
    A(m-1,1) = dummy * A(m-1,1)
    A(m-2,1) = dummy * A(m-2,1)
    A(m-3,1) = dummy * A(m-3,1)

    A(m  ,2) = A(m  ,2) -  A(m-1,2) * A(m+1,1)
    A(m-1,3) = A(m-1,3) -  A(m-2,3) * A(m+1,1)
    A(m-2,4) = A(m-2,4) -  A(m-3,4) * A(m+1,1)

    dummy = 1.d0 / A(m,2)
    A(m+1,2) = dummy * (A(m+1,2) -  A(m-1,2) * A(m+2,1)) 
    A(m+2,2) = dummy * (A(m+2,2) -  A(m-1,2) * A(m+3,1)) 
    A(m+3,2) = dummy * (A(m+3,2) -  A(m-1,2) * A(m-1,1)) 
    A(m-2,2) = dummy * (A(m-2,2) -  A(m-1,2) * A(m-2,1)) 
    A(m-3,2) = dummy * (A(m-3,2) -  A(m-1,2) * A(m-3,1)) 

    A(m  , 3) = A(m  ,3) - A(m-1,3) * A(m+1,2) - A(m-2,3) * A(m+2,1) 
    A(m-1, 4) = A(m-1,4) - A(m-2,4) * A(m+1,2) - A(m-3,4) * A(m+2,1) 
    A(m-2, 5) = A(m-2,5) - A(m-3,5) * A(m+1,2)

    dummy = 1.d0 / A(m,3)
    A(m+1,3) = dummy * (A(m+1,3) - A(m-1,3) * A(m+2,2) - A(m-2,3) * A(m+3,1))
    A(m+2,3) = dummy * (A(m+2,3) - A(m-1,3) * A(m+3,2) - A(m-2,3) * A(m-1,1))
    A(m+3,3) = dummy * (A(m+3,3) - A(m-1,3) * A(m-2,2) - A(m-2,3) * A(m-2,1))
    A(m-3,3) = dummy * (A(m-3,3) - A(m-1,3) * A(m-3,2) - A(m-2,3) * A(m-3,1))

    A(m  , 4) = A(m  , 4) - A(m-1, 4) * A(m+1,3) - A(m-2,4) * A(m+2,2) - A(m-3,4) * A(m+3,1)  
    A(m-1, 5) = A(m-1, 5) - A(m-2, 5) * A(m+1,3) - A(m-3,5) * A(m+2,2)
    A(m-2, 6) = A(m-2, 6) - A(m-3, 6) * A(m+1,3)

    dummy = 1.d0 / A(m,4)
    A(m+1,4) = dummy * (A(m+1,4) - A(m-1,4) * A(m+2,3) - A(m-2,4) * A(m+3,2) - A(m-3,4) * A(m-1,1))
    A(m+2,4) = dummy * (A(m+2,4) - A(m-1,4) * A(m+3,3) - A(m-2,4) * A(m-2,2) - A(m-3,4) * A(m-2,1))
    A(m+3,4) = dummy * (A(m+3,4) - A(m-1,4) * A(m-3,3) - A(m-2,4) * A(m-3,2) - A(m-3,4) * A(m-3,1))

    A(m  , 5) = A(m  , 5) - A(m-1, 5) * A(m+1,4) - A(m-2, 5) * A(m+2,3) - A(m-3,5) * A(m+3,2) 
    A(m-1, 6) = A(m-1, 6) - A(m-2, 6) * A(m+1,4) - A(m-3, 6) * A(m+2,3)
    A(m-2, 7) = A(m-2, 7) - A(m-3, 7) * A(m+1,4)

    dummy = 1.d0 / A(m,5)
    A(m+1,5) = dummy * (A(m+1,5) - A(m-1,5) * A(m+2,4) - A(m-2,5) * A(m+3,3) - A(m-3,5) * A(m-2,2))
    A(m+2,5) = dummy * (A(m+2,5) - A(m-1,5) * A(m+3,4) - A(m-2,5) * A(m-3,3) - A(m-3,5) * A(m-3,2))
    A(m+3,5) = dummy *  A(m+3,5)

    A(m  , 6) = A(m  , 6) - A(m-1, 6) * A(m+1,5) - A(m-2, 6) * A(m+2,4) - A(m-3, 6) * A(m+3,3)
    A(m-1, 7) = A(m-1, 7) - A(m-2, 7) * A(m+1,5) - A(m-3, 7) * A(m+2,4)
    A(m-2, 8) = A(m-2, 8) - A(m-3, 8) * A(m+1,5)

    dummy = 1.d0 / A(m,6)
    A(m+1,6) = dummy * (A(m+1,6) - A(m-1,6) * A(m+2,5) - A(m-2,6) * A(m+3,4) - A(m-3,6) * A(m-3,3))
    A(m+2,6) = dummy * (A(m+2,6) - A(m-1,6) * A(m+3,5))
    A(m+3,6) = dummy *  A(m+3,6)

    do i = 2*bw+1,n-2*bw-1
       A(m  ,i  ) = A(m  ,i  ) - A(m-1,i  ) * A(m+1,i-1) - A(m-2,i)   * A(m+2,i-2) - A(m-3,i  ) * A(m+3,i-3)
       A(m-1,i+1) = A(m-1,i+1) - A(m-2,i+1) * A(m+1,i-1) - A(m-3,i+1) * A(m+2,i-2)
       A(m-2,i+2) = A(m-2,i+2) - A(m-3,i+2) * A(m+1,i-1)
       A(m-3,i+3) = A(m-3,i+3)

       dummy = 1.d0 / A(m,i)
       A(m+1,i) = dummy * (A(m+1,i) - A(m-1,i) * A(m+2,i-1) - A(m-2,i) * A(m+3,i-2))
       A(m+2,i) = dummy * (A(m+2,i) - A(m-1,i) * A(m+3,i-1))
       A(m+3,i) = dummy *  A(m+3,i)
    enddo

    A(m  ,n-2*bw)    = A(m  ,n-2*bw)   - A(m-1,n-2*bw)   * A(m+1,n-2*bw-1) - A(m-2,n-2*bw)   * A(m+2,n-2*bw-2) &
         - A(m-3,n-2*bw)   * A(m+3,n-2*bw-3)
    A(m-1,n-2*bw+1)  = A(m-1,n-2*bw+1) - A(m-2,n-2*bw+1) * A(m+1,n-2*bw-1) - A(m-3,n-2*bw+1) * A(m+2,n-2*bw-2)
    A(m-2,n-2*bw+2)  = A(m-2,n-2*bw+2) - A(m-3,n-2*bw+2) * A(m+1,n-2*bw-1)
    A(m-3,n-2*bw+3)  = A(m-3,n-2*bw+3)

    dummy = 1.d0 / A(m,n-2*bw)
    A(m+1,n-2*bw) = dummy * (A(m+1,n-2*bw) - A(m-1,n-2*bw) * A(m+2,n-2*bw-1) - A(m-2,n-2*bw) * A(m+3,n-2*bw-2))
    A(m+2,n-2*bw) = dummy * (A(m+2,n-2*bw) - A(m-1,n-2*bw) * A(m+3,n-2*bw-1))
    A(m+3,n-2*bw) = dummy *  A(m+3,n-2*bw)

    A(m  ,n-2*bw+1)   = A(m  ,n-2*bw+1)  - A(m-1,n-2*bw+1)  * A(m+1,n-2*bw) - A(m-2,n-2*bw+1) * A(m+2,n-2*bw-1)&
         - A(m-3,n-2*bw+1) * A(m+3,n-2*bw-2)
    A(m-1,n-2*bw+2)   = A(m-1,n-2*bw+2)  - A(m-2,n-2*bw+2)  * A(m+1,n-2*bw) - A(m-3,n-2*bw+2) * A(m+2,n-2*bw-1)
    A(m-2,n-2*bw+3)   = A(m-2,n-2*bw+3)  - A(m-3,n-2*bw+3)  * A(m+1,n-2*bw) 
    A(m-3,n-2*bw+4)   = A(m-3,n-2*bw+4)  - A(m+3,n-2*bw+4)  * A(m+1,n-2*bw)
    A(m+2,n-2*bw+5)   = A(m+2,n-2*bw+5)  - A(m+3,n-2*bw+5)  * A(m+1,n-2*bw) 
    A(m+2,n-2*bw+6)   = A(m+2,n-2*bw+6)  - A(m+3,n-2*bw+6)  * A(m+1,n-2*bw) 

    dummy = 1.d0 / A(m,n-2*bw+1)
    A(m+1,n-2*bw+1) = dummy * (A(m+1,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+2,n-2*bw) - A(m-2,n-2*bw+1) * A(m+3,n-2*bw-1))
    A(m+2,n-2*bw+1) = dummy * (A(m+2,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+3,n-2*bw))
    A(m+3,n-2*bw+1) = dummy *  A(m+3,n-2*bw+1)

    A(m  ,n-2*bw+2)  = A(m  ,n-2*bw+2)  - A(m-1,n-2*bw+2)  * A(m+1,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+2,n-2*bw) &
         - A(m-3,n-2*bw+2) * A(m+3,n-2*bw-1)
    A(m-1,n-2*bw+3)  = A(m-1,n-2*bw+3)  - A(m-2,n-2*bw+3)  * A(m+1,n-2*bw+1) - A(m-3,n-2*bw+3) * A(m+2,n-2*bw)
    A(m-2,n-2*bw+4)  = A(m-2,n-2*bw+4)  - A(m-3,n-2*bw+4)  * A(m+1,n-2*bw+1) - A(m+3,n-2*bw+4) * A(m+2,n-2*bw)
    A(m-3,n-2*bw+5)  = A(m-3,n-2*bw+5)  - A(m+2,n-2*bw+5)  * A(m+1,n-2*bw+1) - A(m+3,n-2*bw+5) * A(m+2,n-2*bw)
    A(m+1,n-2*bw+6)  = A(m+1,n-2*bw+6)  - A(m+2,n-2*bw+6)  * A(m+1,n-2*bw+1) - A(m+3,n-2*bw+6) * A(m+2,n-2*bw)

    dummy = 1.d0 / A(m,n-2*bw+2)
    A(m+1,n-2*bw+2) = dummy * (A(m+1,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+2,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+3,n-2*bw))
    A(m+2,n-2*bw+2) = dummy * (A(m+2,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+3,n-2*bw+1))
    A(m+3,n-2*bw+2) = dummy *  A(m+3,n-2*bw+2)

    A(m  ,n-2*bw+3) = A(m  ,n-2*bw+3)  - A(m-1,n-2*bw+3)  * A(m+1,n-2*bw+2) - A(m-2,n-2*bw+3)  * A(m+2,n-2*bw+1) &
         - A(m-3,n-2*bw+3)  * A(m+3,n-2*bw)
    A(m-1,n-2*bw+4) = A(m-1,n-2*bw+4)  - A(m-2,n-2*bw+4)  * A(m+1,n-2*bw+2) - A(m-3,n-2*bw+4)  * A(m+2,n-2*bw+1) & 
         - A(m+3,n-2*bw+4)  * A(m+3,n-2*bw)
    A(m-2,n-2*bw+5) = A(m-2,n-2*bw+5)  - A(m-3,n-2*bw+5)  * A(m+1,n-2*bw+2) - A(m+2,n-2*bw+5)  * A(m+2,n-2*bw+1) &
         - A(m+3,n-2*bw+5)  * A(m+3,n-2*bw)
    A(m-3,n-2*bw+6) = A(m-3,n-2*bw+6)  - A(m+1,n-2*bw+6)  * A(m+1,n-2*bw+2) - A(m+2,n-2*bw+6)  * A(m+2,n-2*bw+1) &
         - A(m+3,n-2*bw+6)  * A(m+3,n-2*bw)


    dummy = 1.d0 / A(m,n-2*bw+3)
    A(m+1,n-2*bw+3) = dummy * (A(m+1,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+2,n-2*bw+2) - A(m-2,n-2*bw+3) * A(m+3,n-2*bw+1))
    A(m+2,n-2*bw+3) = dummy * (A(m+2,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+3,n-2*bw+2))
    A(m+3,n-2*bw+3) = dummy *  A(m+3,n-2*bw+3)

    A(m  ,n-2*bw+4) = A(m  ,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+1,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+2,n-2*bw+2) &
         - A(m-3,n-2*bw+4) * A(m+3,n-2*bw+1)
    A(m-1,n-2*bw+5) = A(m-1,n-2*bw+5) - A(m-2,n-2*bw+5) * A(m+1,n-2*bw+3) - A(m-3,n-2*bw+5) * A(m+2,n-2*bw+2) &
         - A(m+2,n-2*bw+5) * A(m+3,n-2*bw+1)
    A(m-2,n-2*bw+6) = A(m-2,n-2*bw+6) - A(m-3,n-2*bw+6) * A(m+1,n-2*bw+3) - A(m+1,n-2*bw+6) * A(m+2,n-2*bw+2) &
         - A(m+2,n-2*bw+6) * A(m+3,n-2*bw+1)


    dummy = 1.d0 / A(m,n-2*bw+4)
    A(m+1,n-2*bw+4) = dummy * (A(m+1,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+2,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+3,n-2*bw+2))
    A(m+2,n-2*bw+4) = dummy * (A(m+2,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+3,n-2*bw+3))

    A(m  ,n-2*bw+5)   = A(m  ,n-2*bw+5)  - A(m-1,n-2*bw+5)  * A(m+1,n-2*bw+4) - A(m-2,n-2*bw+5)  * A(m+2,n-2*bw+3) &
         - A(m-3,n-2*bw+5)  * A(m+3,n-2*bw+2)
    A(m-1,n-2*bw+6)   = A(m-1,n-2*bw+6)  - A(m-2,n-2*bw+6)  * A(m+1,n-2*bw+4) - A(m-3,n-2*bw+6)  * A(m+2,n-2*bw+3) &
         - A(m+1,n-2*bw+6)  * A(m+3,n-2*bw+2)

    dummy = 1.d0 / A(m,n-2*bw+5)
    A(m+1,n-2*bw+5) = dummy * (A(m+1,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+2,n-2*bw+4) - A(m-2,n-2*bw+5) * A(m+3,n-2*bw+3))

    A(m  ,n-2*bw+6)   = A(m  ,n-2*bw+6)  - A(m-1,n-2*bw+6)  * A(m+1,n-2*bw+5) - A(m-2,n-2*bw+6)  * A(m+2,n-2*bw+4) &
         - A(m-3,n-2*bw+6)  * A(m+3,n-2*bw+3)

    return

  end subroutine fixed_bandec3


  subroutine fixed_bandec4(A,n)
    implicit none
    integer :: n
    real(c_double) :: A(9,n)
    integer :: i
    integer :: m,bw
    real(c_double) :: dummy

    m = 5
    bw = 4

    dummy = 1.d0 / A(5,1)
    A(6,1) = dummy * A(6,1)
    A(7,1) = dummy * A(7,1)
    A(8,1) = dummy * A(8,1)
    A(9,1) = dummy * A(9,1)
    A(4,1) = dummy * A(4,1)
    A(3,1) = dummy * A(3,1)
    A(2,1) = dummy * A(2,1)
    A(1,1) = dummy * A(1,1)

    A(5,2) = A(5,2) -  A(4,2) * A(6,1)
    A(4,3) = A(4,3) -  A(3,3) * A(6,1)
    A(3,4) = A(3,4) -  A(2,4) * A(6,1)
    A(2,5) = A(2,5) -  A(1,5) * A(6,1)

    dummy = 1.d0 / A(5,2)
    A(6,2) = dummy * (A(6,2) -  A(4,2) * A(7,1)) 
    A(7,2) = dummy * (A(7,2) -  A(4,2) * A(8,1)) 
    A(8,2) = dummy * (A(8,2) -  A(4,2) * A(9,1)) 
    A(9,2) = dummy * (A(9,2) -  A(4,2) * A(4,1)) 
    A(3,2) = dummy * (A(3,2) -  A(4,2) * A(3,1)) 
    A(2,2) = dummy * (A(2,2) -  A(4,2) * A(2,1)) 
    A(1,2) = dummy * (A(1,2) -  A(4,2) * A(1,1)) 

    A(5, 3) = A(5, 3) - A(4,3) * A(6,2) - A(3,3) * A(7,1) 
    A(4, 4) = A(4, 4) - A(3,4) * A(6,2) - A(2,4) * A(7,1) 
    A(3, 5) = A(3, 5) - A(2,5) * A(6,2) - A(1,5) * A(7,1)
    A(2, 6) = A(2, 6) - A(1,6) * A(6,2) 

    dummy = 1.d0 / A(5,3)
    A(6,3) = dummy * (A(6,3) - A(4,3) * A(7,2) - A(3,3) * A(8,1))
    A(7,3) = dummy * (A(7,3) - A(4,3) * A(8,2) - A(3,3) * A(9,1))
    A(8,3) = dummy * (A(8,3) - A(4,3) * A(9,2) - A(3,3) * A(4,1))
    A(9,3) = dummy * (A(9,3) - A(4,3) * A(3,2) - A(3,3) * A(3,1))
    A(2,3) = dummy * (A(2,3) - A(4,3) * A(2,2) - A(3,3) * A(2,1))
    A(1,3) = dummy * (A(1,3) - A(4,3) * A(1,2) - A(3,3) * A(1,1))

    A(5, 4) = A(5, 4) - A(4, 4) * A(6,3) - A(3,4) * A(7,2) - A(2,4) * A(8,1)  
    A(4, 5) = A(4, 5) - A(3, 5) * A(6,3) - A(2,5) * A(7,2) - A(1,5) * A(8,1)
    A(3, 6) = A(3, 6) - A(2, 6) * A(6,3) - A(1,6) * A(7,2)
    A(2, 7) = A(2, 7) - A(1, 7) * A(6,3)

    dummy = 1.d0 / A(5,4)
    A(6,4) = dummy * (A(6,4) - A(4,4) * A(7,3) - A(3,4) * A(8,2) - A(2,4) * A(9,1))
    A(7,4) = dummy * (A(7,4) - A(4,4) * A(8,3) - A(3,4) * A(9,2) - A(2,4) * A(4,1))
    A(8,4) = dummy * (A(8,4) - A(4,4) * A(9,3) - A(3,4) * A(3,2) - A(2,4) * A(3,1))
    A(9,4) = dummy * (A(9,4) - A(4,4) * A(2,3) - A(3,4) * A(2,2) - A(2,4) * A(2,1))
    A(1,4) = dummy * (A(1,4) - A(4,4) * A(1,3) - A(3,4) * A(1,2) - A(2,4) * A(1,1))

    A(5, 5) = A(5, 5) - A(4, 5) * A(6,4) - A(3, 5) * A(7,3) - A(2,5) * A(8,2) - A(1,5) * A(9,1)    
    A(4, 6) = A(4, 6) - A(3, 6) * A(6,4) - A(2, 6) * A(7,3) - A(1,6) * A(8,2)
    A(3, 7) = A(3, 7) - A(2, 7) * A(6,4) - A(1, 7) * A(7,3)
    A(2, 8) = A(2, 8) - A(1, 8) * A(6,4)

    dummy = 1.d0 / A(5,5)
    A(6,5) = dummy * (A(6,5) - A(4,5) * A(7,4) - A(3,5) * A(8,3) - A(2,5) * A(9,2) - A(1,5) * A(4,1))
    A(7,5) = dummy * (A(7,5) - A(4,5) * A(8,4) - A(3,5) * A(9,3) - A(2,5) * A(3,2) - A(1,5) * A(3,1))
    A(8,5) = dummy * (A(8,5) - A(4,5) * A(9,4) - A(3,5) * A(2,3) - A(2,5) * A(2,2) - A(1,5) * A(2,1))
    A(9,5) = dummy * (A(9,5) - A(4,5) * A(1,4) - A(3,5) * A(1,3) - A(2,5) * A(1,2) - A(1,5) * A(1,1))

    A(5, 6) = A(5, 6) - A(4, 6) * A(6,5) - A(3, 6) * A(7,4) - A(2, 6) * A(8,3) - A(1,6) * A(9,2)
    A(4, 7) = A(4, 7) - A(3, 7) * A(6,5) - A(2, 7) * A(7,4) - A(1, 7) * A(8,3)
    A(3, 8) = A(3, 8) - A(2, 8) * A(6,5) - A(1, 8) * A(7,4)
    A(2, 9) = A(2, 9) - A(1, 9) * A(6,5) 

    dummy = 1.d0 / A(5,6)
    A(6,6) = dummy * (A(6,6) - A(4,6) * A(7,5) - A(3,6) * A(8,4) - A(2,6) * A(9,3) - A(1,6) * A(3,2))
    A(7,6) = dummy * (A(7,6) - A(4,6) * A(8,5) - A(3,6) * A(9,4) - A(2,6) * A(2,3) - A(1,6) * A(2,2))
    A(8,6) = dummy * (A(8,6) - A(4,6) * A(9,5) - A(3,6) * A(1,4) - A(2,6) * A(1,3) - A(1,6) * A(1,2))
    A(9,6) = dummy *  A(9,6)


    A(5, 7) = A(5, 7) - A(4, 7) * A(6,6) - A(3, 7) * A(7,5) - A(2, 7) * A(8,4) - A(1, 7) * A(9,3)
    A(4, 8) = A(4, 8) - A(3, 8) * A(6,6) - A(2, 8) * A(7,5) - A(1, 8) * A(8,4)
    A(3, 9) = A(3, 9) - A(2, 9) * A(6,6) - A(1, 9) * A(7,5)
    A(2,10) = A(2,10) - A(1,10) * A(6,6)

    dummy = 1.d0 / A(5,7)
    A(6,7) = dummy * (A(6,7) - A(4,7) * A(7,6) - A(3,7) * A(8,5) - A(2,7) * A(9,4) - A(1,7) * A(2,3))
    A(7,7) = dummy * (A(7,7) - A(4,7) * A(8,6) - A(3,7) * A(9,5) - A(2,7) * A(1,4) - A(1,7) * A(1,3))
    A(8,7) = dummy * (A(8,7) - A(4,7) * A(9,6))
    A(9,7) = dummy *  A(9,7)

    A(5, 8) = A(5, 8) - A(4, 8) * A(6,7) - A(3, 8) * A(7,6) - A(2, 8) * A(8,5) - A(1, 8) * A(9,4)
    A(4, 9) = A(4, 9) - A(3, 9) * A(6,7) - A(2, 9) * A(7,6) - A(1, 9) * A(8,5)
    A(3,10) = A(3,10) - A(2,10) * A(6,7) - A(1,10) * A(7,6)
    A(2,11) = A(2,11) - A(1,11) * A(6,7) 


    dummy = 1.d0 / A(5,8)
    A(6,8) = dummy * (A(6,8) - A(4,8) * A(7,7) - A(3,8) * A(8,6) - A(2,8) * A(9,5) - A(1,8) * A(1,4))
    A(7,8) = dummy * (A(7,8) - A(4,8) * A(8,7) - A(3,8) * A(9,6))
    A(8,8) = dummy * (A(8,8) - A(4,8) * A(9,7))
    A(9,8) = dummy *  A(9,8)


    do i = 2*bw+1,n-2*bw-1
       A(5,i  ) = A(5,i  ) - A(4,i  ) * A(6,i-1) - A(3,i)   * A(7,i-2) - A(2,i  ) * A(8,i-3) - A(1,i)  &
            * A(9,i-4)
       A(4,i+1) = A(4,i+1) - A(3,i+1) * A(6,i-1) - A(2,i+1) * A(7,i-2) - A(1,i+1) * A(8,i-3)
       A(3,i+2) = A(3,i+2) - A(2,i+2) * A(6,i-1) - A(1,i+2) * A(7,i-2)
       A(2,i+3) = A(2,i+3) - A(1,i+3) * A(6,i-1)

       dummy = 1.d0 / A(5,i)
       A(6,i) = dummy * (A(6,i) - A(4,i) * A(7,i-1) - A(3,i) * A(8,i-2) - A(2,i) * A(9,i-3))
       A(7,i) = dummy * (A(7,i) - A(4,i) * A(8,i-1) - A(3,i) * A(9,i-2))
       A(8,i) = dummy * (A(8,i) - A(4,i) * A(9,i-1))
       A(9,i) = dummy *  A(9,i)
    enddo

    A(5,n-2*bw)    = A(5,n-2*bw)    - A(4,n-2*bw)   * A(6,n-2*bw-1) - A(3,n-2*bw)   * A(7,n-2*bw-2) - A(2,n-2*bw) &
         * A(8,n-2*bw-3) - A(1,n-2*bw)   * A(9,n-2*bw-4)
    A(4,n-2*bw+1)  = A(4,n-2*bw+1)  - A(3,n-2*bw+1) * A(6,n-2*bw-1) - A(2,n-2*bw+1) * A(7,n-2*bw-2) - A(1,n-2*bw+1)&
         * A(8,n-2*bw-3)
    A(3,n-2*bw+2)  = A(3,n-2*bw+2)  - A(2,n-2*bw+2) * A(6,n-2*bw-1) - A(1,n-2*bw+2) * A(7,n-2*bw-2)
    A(2,n-2*bw+3)  = A(2,n-2*bw+3)  - A(1,n-2*bw+3) * A(6,n-2*bw-1)
    A(1,n-2*bw+4)  = A(1,n-2*bw+4)

    dummy = 1.d0 / A(5,n-2*bw)
    A(6,n-2*bw) = dummy * (A(6,n-2*bw) - A(4,n-2*bw) * A(7,n-2*bw-1) - A(3,n-2*bw) * A(8,n-2*bw-2) - A(2,n-2*bw)&
         * A(9,n-2*bw-3))
    A(7,n-2*bw) = dummy * (A(7,n-2*bw) - A(4,n-2*bw) * A(8,n-2*bw-1) - A(3,n-2*bw) * A(9,n-2*bw-2))
    A(8,n-2*bw) = dummy * (A(8,n-2*bw) - A(4,n-2*bw) * A(9,n-2*bw-1))
    A(9,n-2*bw) = dummy *  A(9,n-2*bw)

    A(5,n-2*bw+1)   = A(5,n-2*bw+1)  - A(4,n-2*bw+1)  * A(6,n-2*bw) - A(3,n-2*bw+1) * A(7,n-2*bw-1) - A(2,n-2*bw+1)&
         * A(8,n-2*bw-2) - A(1,n-2*bw+1) * A(9,n-2*bw-3)
    A(4,n-2*bw+2)   = A(4,n-2*bw+2)  - A(3,n-2*bw+2)  * A(6,n-2*bw) - A(2,n-2*bw+2) * A(7,n-2*bw-1) - A(1,n-2*bw+2)&
         * A(8,n-2*bw-2)
    A(3,n-2*bw+3)   = A(3,n-2*bw+3)  - A(2,n-2*bw+3)  * A(6,n-2*bw) - A(1,n-2*bw+3) * A(7,n-2*bw-1) 
    A(2,n-2*bw+4)   = A(2,n-2*bw+4)  - A(1,n-2*bw+4)  * A(6,n-2*bw)
    A(1,n-2*bw+5)   = A(1,n-2*bw+5)  - A(9,n-2*bw+5)  * A(6,n-2*bw)
    A(8,n-2*bw+6)   = A(8,n-2*bw+6)  - A(9,n-2*bw+6)  * A(6,n-2*bw) 
    A(8,n-2*bw+7)   = A(8,n-2*bw+7)  - A(9,n-2*bw+7)  * A(6,n-2*bw) 
    A(8,n-2*bw+8)   = A(8,n-2*bw+8)  - A(9,n-2*bw+8)  * A(6,n-2*bw) 

    dummy = 1.d0 / A(5,n-2*bw+1)
    A(6,n-2*bw+1) = dummy * (A(6,n-2*bw+1) - A(4,n-2*bw+1) * A(7,n-2*bw) - A(3,n-2*bw+1) * A(8,n-2*bw-1)&
         - A(2,n-2*bw+1) * A(9,n-2*bw-2))
    A(7,n-2*bw+1) = dummy * (A(7,n-2*bw+1) - A(4,n-2*bw+1) * A(8,n-2*bw) - A(3,n-2*bw+1) * A(9,n-2*bw-1))
    A(8,n-2*bw+1) = dummy * (A(8,n-2*bw+1) - A(4,n-2*bw+1) * A(9,n-2*bw))
    A(9,n-2*bw+1) = dummy *  A(9,n-2*bw+1)

    A(5,n-2*bw+2)  = A(5,n-2*bw+2)  - A(4,n-2*bw+2)  * A(6,n-2*bw+1) - A(3,n-2*bw+2) * A(7,n-2*bw) - A(2,n-2*bw+2)&
         * A(8,n-2*bw-1) - A(1,n-2*bw+2) * A(9,n-2*bw-2)
    A(4,n-2*bw+3)  = A(4,n-2*bw+3)  - A(3,n-2*bw+3)  * A(6,n-2*bw+1) - A(2,n-2*bw+3) * A(7,n-2*bw) - A(1,n-2*bw+3)&
         * A(8,n-2*bw-1)
    A(3,n-2*bw+4)  = A(3,n-2*bw+4)  - A(2,n-2*bw+4)  * A(6,n-2*bw+1) - A(1,n-2*bw+4) * A(7,n-2*bw)
    A(2,n-2*bw+5)  = A(2,n-2*bw+5)  - A(1,n-2*bw+5)  * A(6,n-2*bw+1) - A(9,n-2*bw+5) * A(7,n-2*bw)
    A(1,n-2*bw+6)  = A(1,n-2*bw+6)  - A(8,n-2*bw+6)  * A(6,n-2*bw+1) - A(9,n-2*bw+6) * A(7,n-2*bw)
    A(7,n-2*bw+7)  = A(7,n-2*bw+7)  - A(8,n-2*bw+7)  * A(6,n-2*bw+1) - A(9,n-2*bw+7) * A(7,n-2*bw)
    A(7,n-2*bw+8)  = A(7,n-2*bw+8)  - A(8,n-2*bw+8)  * A(6,n-2*bw+1) - A(9,n-2*bw+8) * A(7,n-2*bw)

    dummy = 1.d0 / A(5,n-2*bw+2)
    A(6,n-2*bw+2) = dummy * (A(6,n-2*bw+2) - A(4,n-2*bw+2) * A(7,n-2*bw+1) - A(3,n-2*bw+2) * A(8,n-2*bw)&
         - A(2,n-2*bw+2) * A(9,n-2*bw-1))
    A(7,n-2*bw+2) = dummy * (A(7,n-2*bw+2) - A(4,n-2*bw+2) * A(8,n-2*bw+1) - A(3,n-2*bw+2) * A(9,n-2*bw))
    A(8,n-2*bw+2) = dummy * (A(8,n-2*bw+2) - A(4,n-2*bw+2) * A(9,n-2*bw+1))
    A(9,n-2*bw+2) = dummy *  A(9,n-2*bw+2)

    A(5,n-2*bw+3) = A(5,n-2*bw+3)  - A(4,n-2*bw+3)  * A(6,n-2*bw+2) - A(3,n-2*bw+3)  * A(7,n-2*bw+1) - A(2,n-2*bw+3) &
         * A(8,n-2*bw) - A(1,n-2*bw+3) * A(9,n-2*bw-1)
    A(4,n-2*bw+4) = A(4,n-2*bw+4)  - A(3,n-2*bw+4)  * A(6,n-2*bw+2) - A(2,n-2*bw+4)  * A(7,n-2*bw+1) - A(1,n-2*bw+4)&
         * A(8,n-2*bw)
    A(3,n-2*bw+5) = A(3,n-2*bw+5)  - A(2,n-2*bw+5)  * A(6,n-2*bw+2) - A(1,n-2*bw+5)  * A(7,n-2*bw+1) - A(9,n-2*bw+5)&
         * A(8,n-2*bw)
    A(2,n-2*bw+6) = A(2,n-2*bw+6)  - A(1,n-2*bw+6)  * A(6,n-2*bw+2) - A(8,n-2*bw+6)  * A(7,n-2*bw+1) - A(9,n-2*bw+6)&
         * A(8,n-2*bw)
    A(1,n-2*bw+7) = A(1,n-2*bw+7)  - A(7,n-2*bw+7)  * A(6,n-2*bw+2) - A(8,n-2*bw+7)  * A(7,n-2*bw+1) - A(9,n-2*bw+7)&
         * A(8,n-2*bw)
    A(6,n-2*bw+8) = A(6,n-2*bw+8)  - A(7,n-2*bw+8)  * A(6,n-2*bw+2) - A(8,n-2*bw+8)  * A(7,n-2*bw+1) - A(9,n-2*bw+8)&
         * A(8,n-2*bw)


    dummy = 1.d0 / A(5,n-2*bw+3)
    A(6,n-2*bw+3) = dummy * (A(6,n-2*bw+3) - A(4,n-2*bw+3) * A(7,n-2*bw+2) - A(3,n-2*bw+3) * A(8,n-2*bw+1)&
         - A(2,n-2*bw+3) * A(9,n-2*bw))
    A(7,n-2*bw+3) = dummy * (A(7,n-2*bw+3) - A(4,n-2*bw+3) * A(8,n-2*bw+2) - A(3,n-2*bw+3) * A(9,n-2*bw+1))
    A(8,n-2*bw+3) = dummy * (A(8,n-2*bw+3) - A(4,n-2*bw+3) * A(9,n-2*bw+2))
    A(9,n-2*bw+3) = dummy *  A(9,n-2*bw+3)


    A(5,n-2*bw+4)  = A(5,n-2*bw+4)  - A(4,n-2*bw+4)  * A(6,n-2*bw+3) - A(3,n-2*bw+4)  * A(7,n-2*bw+2) - A(2,n-2*bw+4)&
         * A(8,n-2*bw+1) - A(1,n-2*bw+4) * A(9,n-2*bw)
    A(4,n-2*bw+5)  = A(4,n-2*bw+5)  - A(3,n-2*bw+5)  * A(6,n-2*bw+3) - A(2,n-2*bw+5)  * A(7,n-2*bw+2) - A(1,n-2*bw+5)&
         * A(8,n-2*bw+1) - A(9,n-2*bw+5) * A(9,n-2*bw)
    A(3,n-2*bw+6)  = A(3,n-2*bw+6)  - A(2,n-2*bw+6)  * A(6,n-2*bw+3) - A(1,n-2*bw+6)  * A(7,n-2*bw+2) - A(8,n-2*bw+6)&
         * A(8,n-2*bw+1) - A(9,n-2*bw+6) * A(9,n-2*bw)
    A(2,n-2*bw+7)  = A(2,n-2*bw+7)  - A(1,n-2*bw+7)  * A(6,n-2*bw+3) - A(7,n-2*bw+7)  * A(7,n-2*bw+2) - A(8,n-2*bw+7)&
         * A(8,n-2*bw+1) - A(9,n-2*bw+7) * A(9,n-2*bw) 
    A(1,n-2*bw+8)  = A(1,n-2*bw+8)  - A(6,n-2*bw+8)  * A(6,n-2*bw+3) - A(7,n-2*bw+8)  * A(7,n-2*bw+2) - A(8,n-2*bw+8)&
         * A(8,n-2*bw+1) - A(9,n-2*bw+8) * A(9,n-2*bw)


    dummy = 1.d0 / A(5,n-2*bw+4)
    A(6,n-2*bw+4) = dummy * (A(6,n-2*bw+4) - A(4,n-2*bw+4) * A(7,n-2*bw+3) - A(3,n-2*bw+4) * A(8,n-2*bw+2)&
         - A(2,n-2*bw+4) * A(9,n-2*bw+1))
    A(7,n-2*bw+4) = dummy * (A(7,n-2*bw+4) - A(4,n-2*bw+4) * A(8,n-2*bw+3) - A(3,n-2*bw+4) * A(9,n-2*bw+2))
    A(8,n-2*bw+4) = dummy * (A(8,n-2*bw+4) - A(4,n-2*bw+4) * A(9,n-2*bw+3))
    A(9,n-2*bw+4) = dummy *  A(9,n-2*bw+4)

    A(5,n-2*bw+5)   = A(5,n-2*bw+5)  - A(4,n-2*bw+5)  * A(6,n-2*bw+4) - A(3,n-2*bw+5)  * A(7,n-2*bw+3) - A(2,n-2*bw+5)&
         * A(8,n-2*bw+2) - A(1,n-2*bw+5)  * A(9,n-2*bw+1)
    A(4,n-2*bw+6)   = A(4,n-2*bw+6)  - A(3,n-2*bw+6)  * A(6,n-2*bw+4) - A(2,n-2*bw+6)  * A(7,n-2*bw+3) - A(1,n-2*bw+6)&
         * A(8,n-2*bw+2) - A(8,n-2*bw+6)  * A(9,n-2*bw+1)
    A(3,n-2*bw+7)   = A(3,n-2*bw+7)  - A(2,n-2*bw+7)  * A(6,n-2*bw+4) - A(1,n-2*bw+7)  * A(7,n-2*bw+3) - A(7,n-2*bw+7)&
         * A(8,n-2*bw+2) - A(8,n-2*bw+7)  * A(9,n-2*bw+1)
    A(2,n-2*bw+8)   = A(2,n-2*bw+8)  - A(1,n-2*bw+8)  * A(6,n-2*bw+4) - A(6,n-2*bw+8)  * A(7,n-2*bw+3) - A(7,n-2*bw+8)&
         * A(8,n-2*bw+2) - A(8,n-2*bw+8)  * A(9,n-2*bw+1)

    dummy = 1.d0 / A(5,n-2*bw+5)
    A(6,n-2*bw+5) = dummy * (A(6,n-2*bw+5) - A(4,n-2*bw+5) * A(7,n-2*bw+4) - A(3,n-2*bw+5) * A(8,n-2*bw+3)&
         - A(2,n-2*bw+5) * A(9,n-2*bw+2))
    A(7,n-2*bw+5) = dummy * (A(7,n-2*bw+5) - A(4,n-2*bw+5) * A(8,n-2*bw+4) - A(3,n-2*bw+5) * A(9,n-2*bw+3))
    A(8,n-2*bw+5) = dummy * (A(8,n-2*bw+5) - A(4,n-2*bw+5) * A(9,n-2*bw+4))

    A(5,n-2*bw+6)   = A(5,n-2*bw+6)  - A(4,n-2*bw+6)  * A(6,n-2*bw+5) - A(3,n-2*bw+6)  * A(7,n-2*bw+4) - A(2,n-2*bw+6) &
         * A(8,n-2*bw+3) - A(1,n-2*bw+6)  * A(9,n-2*bw+2)
    A(4,n-2*bw+7)   = A(4,n-2*bw+7)  - A(3,n-2*bw+7)  * A(6,n-2*bw+5) - A(2,n-2*bw+7)  * A(7,n-2*bw+4) - A(1,n-2*bw+7)&
         * A(8,n-2*bw+3) - A(7,n-2*bw+7)  * A(9,n-2*bw+2)
    A(3,n-2*bw+8)   = A(3,n-2*bw+8)  - A(2,n-2*bw+8)  * A(6,n-2*bw+5) - A(1,n-2*bw+8)  * A(7,n-2*bw+4) - A(6,n-2*bw+8)&
         * A(8,n-2*bw+3) - A(7,n-2*bw+8)  * A(9,n-2*bw+2)

    dummy = 1.d0 / A(5,n-2*bw+6)
    A(6,n-2*bw+6) = dummy * (A(6,n-2*bw+6) - A(4,n-2*bw+6) * A(7,n-2*bw+5) - A(3,n-2*bw+6) * A(8,n-2*bw+4)&
         - A(2,n-2*bw+6) * A(9,n-2*bw+3))
    A(7,n-2*bw+6) = dummy * (A(7,n-2*bw+6) - A(4,n-2*bw+6) * A(8,n-2*bw+5) - A(3,n-2*bw+6) * A(9,n-2*bw+4))


    A(5,n-2*bw+7)   = A(5,n-2*bw+7)  - A(4,n-2*bw+7)  * A(6,n-2*bw+6) - A(3,n-2*bw+7)  * A(7,n-2*bw+5) - A(2,n-2*bw+7) &
         * A(8,n-2*bw+4) - A(1,n-2*bw+7)  * A(9,n-2*bw+3)
    A(4,n-2*bw+8)   = A(4,n-2*bw+8)  - A(3,n-2*bw+8)  * A(6,n-2*bw+6) - A(2,n-2*bw+8)  * A(7,n-2*bw+5) - A(1,n-2*bw+8)&
         * A(8,n-2*bw+4) - A(6,n-2*bw+8)  * A(9,n-2*bw+3)


    dummy = 1.d0 / A(5,n-2*bw+7)
    A(6,n-2*bw+7) = dummy * (A(6,n-2*bw+7) - A(4,n-2*bw+7) * A(7,n-2*bw+6) - A(3,n-2*bw+7) * A(8,n-2*bw+5) - A(2,n-2*bw+7)&
         * A(9,n-2*bw+4))

    A(5,n-2*bw+8)   = A(5,n-2*bw+8)  - A(4,n-2*bw+8)  * A(6,n-2*bw+7) - A(3,n-2*bw+8)  * A(7,n-2*bw+6) - A(2,n-2*bw+8)&
         * A(8,n-2*bw+5) - A(1,n-2*bw+8)  * A(9,n-2*bw+4)

    return

  end subroutine fixed_bandec4


  subroutine fixed_bandec5(A,n)
    implicit none
    integer :: n
    real(c_double) :: A(11,n)
    integer :: i
    integer :: m,bw
    real(c_double) :: dummy

    m = 6
    bw = 5

    dummy = 1.d0 / A(m,1)
    A(m+1,1) = dummy * A(m+1,1)
    A(m+2,1) = dummy * A(m+2,1)
    A(m+3,1) = dummy * A(m+3,1)
    A(m+4,1) = dummy * A(m+4,1)
    A(m+5,1) = dummy * A(m+5,1)
    A(m-1,1) = dummy * A(m-1,1)
    A(m-2,1) = dummy * A(m-2,1)
    A(m-3,1) = dummy * A(m-3,1)
    A(m-4,1) = dummy * A(m-4,1)
    A(m-5,1) = dummy * A(m-5,1)

    A(m  ,2) = A(m  ,2) -  A(m-1,2) * A(m+1,1)
    A(m-1,3) = A(m-1,3) -  A(m-2,3) * A(m+1,1)
    A(m-2,4) = A(m-2,4) -  A(m-3,4) * A(m+1,1)
    A(m-3,5) = A(m-3,5) -  A(m-4,5) * A(m+1,1)
    A(m-4,6) = A(m-4,6) -  A(m-5,6) * A(m+1,1)

    dummy = 1.d0 / A(m,2)
    A(m+1,2) = dummy * (A(m+1,2) -  A(m-1,2) * A(m+2,1)) 
    A(m+2,2) = dummy * (A(m+2,2) -  A(m-1,2) * A(m+3,1)) 
    A(m+3,2) = dummy * (A(m+3,2) -  A(m-1,2) * A(m+4,1)) 
    A(m+4,2) = dummy * (A(m+4,2) -  A(m-1,2) * A(m+5,1)) 
    A(m+5,2) = dummy * (A(m+5,2) -  A(m-1,2) * A(m-1,1)) 
    A(m-2,2) = dummy * (A(m-2,2) -  A(m-1,2) * A(m-2,1)) 
    A(m-3,2) = dummy * (A(m-3,2) -  A(m-1,2) * A(m-3,1)) 
    A(m-4,2) = dummy * (A(m-4,2) -  A(m-1,2) * A(m-4,1)) 
    A(m-5,2) = dummy * (A(m-5,2) -  A(m-1,2) * A(m-5,1)) 

    A(m  , 3) = A(m  , 3) - A(m-1,3) * A(m+1,2) - A(m-2,3) * A(m+2,1) 
    A(m-1, 4) = A(m-1, 4) - A(m-2,4) * A(m+1,2) - A(m-3,4) * A(m+2,1) 
    A(m-2, 5) = A(m-2, 5) - A(m-3,5) * A(m+1,2) - A(m-4,5) * A(m+2,1)
    A(m-3, 6) = A(m-3, 6) - A(m-4,6) * A(m+1,2) - A(m-5,6) * A(m+2,1)
    A(m-4, 7) = A(m-4, 7) - A(m-5,7) * A(m+1,2)


    dummy = 1.d0 / A(m,3)
    A(m+1,3) = dummy * (A(m+1,3) - A(m-1,3) * A(m+2,2) - A(m-2,3) * A(m+3,1))
    A(m+2,3) = dummy * (A(m+2,3) - A(m-1,3) * A(m+3,2) - A(m-2,3) * A(m+4,1))
    A(m+3,3) = dummy * (A(m+3,3) - A(m-1,3) * A(m+4,2) - A(m-2,3) * A(m+5,1))
    A(m+4,3) = dummy * (A(m+4,3) - A(m-1,3) * A(m+5,2) - A(m-2,3) * A(m-1,1))
    A(m+5,3) = dummy * (A(m+5,3) - A(m-1,3) * A(m-2,2) - A(m-2,3) * A(m-2,1))
    A(m-3,3) = dummy * (A(m-3,3) - A(m-1,3) * A(m-3,2) - A(m-2,3) * A(m-3,1))
    A(m-4,3) = dummy * (A(m-4,3) - A(m-1,3) * A(m-4,2) - A(m-2,3) * A(m-4,1))
    A(m-5,3) = dummy * (A(m-5,3) - A(m-1,3) * A(m-5,2) - A(m-2,3) * A(m-5,1))

    A(m  , 4) = A(m  , 4) - A(m-1, 4) * A(m+1,3) - A(m-2,4) * A(m+2,2) - A(m-3,4) * A(m+3,1)  
    A(m-1, 5) = A(m-1, 5) - A(m-2, 5) * A(m+1,3) - A(m-3,5) * A(m+2,2) - A(m-4,5) * A(m+3,1)
    A(m-2, 6) = A(m-2, 6) - A(m-3, 6) * A(m+1,3) - A(m-4,6) * A(m+2,2) - A(m-5,6) * A(m+3,1)
    A(m-3, 7) = A(m-3, 7) - A(m-4, 7) * A(m+1,3) - A(m-5,7) * A(m+2,2)
    A(m-4, 8) = A(m-4, 8) - A(m-5, 8) * A(m+1,3)

    dummy = 1.d0 / A(m,4)
    A(m+1,4) = dummy * (A(m+1,4) - A(m-1,4) * A(m+2,3) - A(m-2,4) * A(m+3,2) - A(m-3,4) * A(m+4,1))
    A(m+2,4) = dummy * (A(m+2,4) - A(m-1,4) * A(m+3,3) - A(m-2,4) * A(m+4,2) - A(m-3,4) * A(m+5,1))
    A(m+3,4) = dummy * (A(m+3,4) - A(m-1,4) * A(m+4,3) - A(m-2,4) * A(m+5,2) - A(m-3,4) * A(m-1,1))
    A(m+4,4) = dummy * (A(m+4,4) - A(m-1,4) * A(m+5,3) - A(m-2,4) * A(m-2,2) - A(m-3,4) * A(m-2,1))
    A(m+5,4) = dummy * (A(m+5,4) - A(m-1,4) * A(m-3,3) - A(m-2,4) * A(m-3,2) - A(m-3,4) * A(m-3,1))
    A(m-4,4) = dummy * (A(m-4,4) - A(m-1,4) * A(m-4,3) - A(m-2,4) * A(m-4,2) - A(m-3,4) * A(m-4,1))
    A(m-5,4) = dummy * (A(m-5,4) - A(m-1,4) * A(m-5,3) - A(m-2,4) * A(m-5,2) - A(m-3,4) * A(m-5,1))

    A(m  , 5) = A(m  , 5) - A(m-1, 5) * A(m+1,4) - A(m-2, 5) * A(m+2,3) - A(m-3,5) * A(m+3,2) - A(m-4,5) * A(m+4,1)    
    A(m-1, 6) = A(m-1, 6) - A(m-2, 6) * A(m+1,4) - A(m-3, 6) * A(m+2,3) - A(m-4,6) * A(m+3,2) - A(m-5,6) * A(m+4,1)    
    A(m-2, 7) = A(m-2, 7) - A(m-3, 7) * A(m+1,4) - A(m-4, 7) * A(m+2,3) - A(m-5,7) * A(m+3,2)
    A(m-3, 8) = A(m-3, 8) - A(m-4, 8) * A(m+1,4) - A(m-5, 8) * A(m+2,3) 
    A(m-4, 9) = A(m-4, 9) - A(m-5, 9) * A(m+1,4)

    dummy = 1.d0 / A(m,5)
    A(m+1,5) = dummy * (A(m+1,5) - A(m-1,5) * A(m+2,4) - A(m-2,5) * A(m+3,3) - A(m-3,5) * A(m+4,2) - A(m-4,5) * A(m+5,1))
    A(m+2,5) = dummy * (A(m+2,5) - A(m-1,5) * A(m+3,4) - A(m-2,5) * A(m+4,3) - A(m-3,5) * A(m+5,2) - A(m-4,5) * A(m-1,1))
    A(m+3,5) = dummy * (A(m+3,5) - A(m-1,5) * A(m+4,4) - A(m-2,5) * A(m+5,3) - A(m-3,5) * A(m-2,2) - A(m-4,5) * A(m-2,1))
    A(m+4,5) = dummy * (A(m+4,5) - A(m-1,5) * A(m+5,4) - A(m-2,5) * A(m-3,3) - A(m-3,5) * A(m-3,2) - A(m-4,5) * A(m-3,1))
    A(m+5,5) = dummy * (A(m+5,5) - A(m-1,5) * A(m-4,4) - A(m-2,5) * A(m-4,3) - A(m-3,5) * A(m-4,2) - A(m-4,5) * A(m-4,1))
    A(m-5,5) = dummy * (A(m-5,5) - A(m-1,5) * A(m-5,4) - A(m-2,5) * A(m-5,3) - A(m-3,5) * A(m-5,2) - A(m-4,5) * A(m-5,1))

    A(m  , 6) = A(m  , 6) - A(m-1, 6) * A(m+1,5) - A(m-2, 6) * A(m+2,4) - A(m-3, 6) * A(m+3,3) - A(m-4,6) * A(m+4,2) &
         - A(m-5,6) * A(m+5,1) 
    A(m-1, 7) = A(m-1, 7) - A(m-2, 7) * A(m+1,5) - A(m-3, 7) * A(m+2,4) - A(m-4, 7) * A(m+3,3) - A(m-5,7) * A(m+4,2)
    A(m-2, 8) = A(m-2, 8) - A(m-3, 8) * A(m+1,5) - A(m-4, 8) * A(m+2,4) - A(m-5, 8) * A(m+3,3)
    A(m-3, 9) = A(m-3, 9) - A(m-4, 9) * A(m+1,5) - A(m-5, 9) * A(m+2,4) 
    A(m-4,10) = A(m-4,10) - A(m-5,10) * A(m+1,5)

    dummy = 1.d0 / A(m,6)
    A(m+1,6) = dummy * (A(m+1,6) - A(m-1,6) * A(m+2,5) - A(m-2,6) * A(m+3,4) - A(m-3,6) * A(m+4,3) - A(m-4,6) * A(m+5,2)&
         - A(m-5,6) * A(m-1,1))
    A(m+2,6) = dummy * (A(m+2,6) - A(m-1,6) * A(m+3,5) - A(m-2,6) * A(m+4,4) - A(m-3,6) * A(m+5,3) - A(m-4,6) * A(m-2,2) &
         - A(m-5,6) * A(m-2,1))
    A(m+3,6) = dummy * (A(m+3,6) - A(m-1,6) * A(m+4,5) - A(m-2,6) * A(m+5,4) - A(m-3,6) * A(m-3,3) - A(m-4,6) * A(m-3,2)&
         - A(m-5,6) * A(m-3,1))
    A(m+4,6) = dummy * (A(m+4,6) - A(m-1,6) * A(m+5,5) - A(m-2,6) * A(m-4,4) - A(m-3,6) * A(m-4,3) - A(m-4,6) * A(m-4,2)&
         - A(m-5,6) * A(m-4,1))
    A(m+5,6) = dummy * (A(m+5,6) - A(m-1,6) * A(m-5,5) - A(m-2,6) * A(m-5,4) - A(m-3,6) * A(m-5,3) - A(m-4,6) * A(m-5,2)&
         - A(m-5,6) * A(m-5,1))


    A(m  , 7) = A(m  , 7) - A(m-1, 7) * A(m+1,6) - A(m-2, 7) * A(m+2,5) - A(m-3, 7) * A(m+3,4) - A(m-4, 7) * A(m+4,3)&
         - A(m-5,7) * A(m+5,2)
    A(m-1, 8) = A(m-1, 8) - A(m-2, 8) * A(m+1,6) - A(m-3, 8) * A(m+2,5) - A(m-4, 8) * A(m+3,4) - A(m-5, 8) * A(m+4,3)
    A(m-2, 9) = A(m-2, 9) - A(m-3, 9) * A(m+1,6) - A(m-4, 9) * A(m+2,5) - A(m-5, 9) * A(m+3,4)
    A(m-3,10) = A(m-3,10) - A(m-4,10) * A(m+1,6) - A(m-5,10) * A(m+2,5)                        
    A(m-4,11) = A(m-4,11) - A(m-5,11) * A(m+1,6)

    dummy = 1.d0 / A(m,7)
    A(m+1,7) = dummy * (A(m+1,7) - A(m-1,7) * A(m+2,6) - A(m-2,7) * A(m+3,5) - A(m-3,7) * A(m+4,4) - A(m-4,7) * A(m+5,3) &
         - A(m-5,7) * A(m-2,2))
    A(m+2,7) = dummy * (A(m+2,7) - A(m-1,7) * A(m+3,6) - A(m-2,7) * A(m+4,5) - A(m-3,7) * A(m+5,4) - A(m-4,7) * A(m-3,3)&
         - A(m-5,7) * A(m-3,2))
    A(m+3,7) = dummy * (A(m+3,7) - A(m-1,7) * A(m+4,6) - A(m-2,7) * A(m+5,5) - A(m-3,7) * A(m-4,4) - A(m-4,7) * A(m-4,3)&
         - A(m-5,7) * A(m-4,2))
    A(m+4,7) = dummy * (A(m+4,7) - A(m-1,7) * A(m+5,6) - A(m-2,7) * A(m-5,5) - A(m-3,7) * A(m-5,4) - A(m-4,7) * A(m-5,3)&
         - A(m-5,7) * A(m-5,2))
    A(m+5,7) = dummy *  A(m+5,7) 


    A(m  , 8) = A(m  , 8) - A(m-1, 8) * A(m+1,7) - A(m-2, 8) * A(m+2,6) - A(m-3, 8) * A(m+3,5) - A(m-4, 8) * A(m+4,4)&
         - A(m-5, 8) * A(m+5,3)
    A(m-1, 9) = A(m-1, 9) - A(m-2, 9) * A(m+1,7) - A(m-3, 9) * A(m+2,6) - A(m-4, 9) * A(m+3,5) - A(m-5, 9) * A(m+4,4)
    A(m-2,10) = A(m-2,10) - A(m-3,10) * A(m+1,7) - A(m-4,10) * A(m+2,6) - A(m-5,10) * A(m+3,5)
    A(m-3,11) = A(m-3,11) - A(m-4,11) * A(m+1,7) - A(m-5,11) * A(m+2,6)
    A(m-4,12) = A(m-4,12) - A(m-5,12) * A(m+1,7) 


    dummy = 1.d0 / A(m,8)
    A(m+1,8) = dummy * (A(m+1,8) - A(m-1,8) * A(m+2,7) - A(m-2,8) * A(m+3,6) - A(m-3,8) * A(m+4,5) - A(m-4,8) * A(m+5,4)&
         - A(m-5,8) * A(m-3,3))
    A(m+2,8) = dummy * (A(m+2,8) - A(m-1,8) * A(m+3,7) - A(m-2,8) * A(m+4,6) - A(m-3,8) * A(m+5,5) - A(m-4,8) * A(m-4,4)&
         - A(m-5,8) * A(m-4,3))
    A(m+3,8) = dummy * (A(m+3,8) - A(m-1,8) * A(m+4,7) - A(m-2,8) * A(m+5,6) - A(m-3,8) * A(m-5,5) - A(m-4,8) * A(m-5,4)&
         - A(m-5,8) * A(m-5,3))
    A(m+4,8) = dummy * (A(m+4,8) - A(m-1,8) * A(m+5,7))
    A(m+5,8) = dummy *  A(m+5,8)  

    A(m  , 9) = A(m  , 9) - A(m-1, 9) * A(m+1,8) - A(m-2, 9) * A(m+2,7) - A(m-3, 9) * A(m+3,6) - A(m-4, 9) * A(m+4,5) &
         - A(m-5, 9) * A(m+5,4)
    A(m-1,10) = A(m-1,10) - A(m-2,10) * A(m+1,8) - A(m-3,10) * A(m+2,7) - A(m-4,10) * A(m+3,6) - A(m-5,10) * A(m+4,5)
    A(m-2,11) = A(m-2,11) - A(m-3,11) * A(m+1,8) - A(m-4,11) * A(m+2,7) - A(m-5,11) * A(m+3,6)
    A(m-3,12) = A(m-3,12) - A(m-4,12) * A(m+1,8) - A(m-5,12) * A(m+2,7)
    A(m-4,13) = A(m-4,13) - A(m-5,13) * A(m+1,8)


    dummy = 1.d0 / A(m,9)
    A(m+1,9) = dummy * (A(m+1,9) - A(m-1,9) * A(m+2,8) - A(m-2,9) * A(m+3,7) - A(m-3,9) * A(m+4,6) - A(m-4,9) * A(m+5,5)&
         - A(m-5,9) * A(m-4,4))
    A(m+2,9) = dummy * (A(m+2,9) - A(m-1,9) * A(m+3,8) - A(m-2,9) * A(m+4,7) - A(m-3,9) * A(m+5,6) - A(m-4,9) * A(m-5,5)&
         - A(m-5,9) * A(m-5,4))
    A(m+3,9) = dummy * (A(m+3,9) - A(m-1,9) * A(m+4,8) - A(m-2,9) * A(m+5,7))
    A(m+4,9) = dummy * (A(m+4,9) - A(m-1,9) * A(m+5,8))
    A(m+5,9) = dummy *  A(m+5,9)

    A(m  ,10) = A(m  ,10) - A(m-1,10) * A(m+1,9) - A(m-2,10) * A(m+2,8) - A(m-3,10) * A(m+3,7) - A(m-4,10) * A(m+4,6) &
         - A(m-5,10) * A(m+5,5)
    A(m-1,11) = A(m-1,11) - A(m-2,11) * A(m+1,9) - A(m-3,11) * A(m+2,8) - A(m-4,11) * A(m+3,7) - A(m-5,11) * A(m+4,6)
    A(m-2,12) = A(m-2,12) - A(m-3,12) * A(m+1,9) - A(m-4,12) * A(m+2,8) - A(m-5,12) * A(m+3,7)
    A(m-3,13) = A(m-3,13) - A(m-4,13) * A(m+1,9) - A(m-5,13) * A(m+2,8)
    A(m-4,14) = A(m-4,14) - A(m-5,14) * A(m+1,9)

    dummy = 1.d0 / A(m,10)
    A(m+1,10) = dummy * (A(m+1,10) - A(m-1,10) * A(m+2,9) - A(m-2,10) * A(m+3,8) - A(m-3,10) * A(m+4,7) - A(m-4,10) * A(m+5,6)&
         - A(m-5,10) * A(m-5,5))
    A(m+2,10) = dummy * (A(m+2,10) - A(m-1,10) * A(m+3,9) - A(m-2,10) * A(m+4,8) - A(m-3,10) * A(m+5,7))
    A(m+3,10) = dummy * (A(m+3,10) - A(m-1,10) * A(m+4,9) - A(m-2,10) * A(m+5,8))
    A(m+4,10) = dummy * (A(m+4,10) - A(m-1,10) * A(m+5,9))
    A(m+5,10) = dummy *  A(m+5,10)

    do i = 2*bw+1,n-2*bw-1
       A(m  ,i  ) = A(m  ,i  ) - A(m-1,i  ) * A(m+1,i-1) - A(m-2,i)   * A(m+2,i-2) - A(m-3,i  ) * A(m+3,i-3) - A(m-4,i)   &
            * A(m+4,i-4) - A(m-5,i)   * A(m+5,i-5)
       A(m-1,i+1) = A(m-1,i+1) - A(m-2,i+1) * A(m+1,i-1) - A(m-3,i+1) * A(m+2,i-2) - A(m-4,i+1) * A(m+3,i-3) - A(m-5,i+1)&
            * A(m+4,i-4)
       A(m-2,i+2) = A(m-2,i+2) - A(m-3,i+2) * A(m+1,i-1) - A(m-4,i+2) * A(m+2,i-2) - A(m-5,i+2) * A(m+3,i-3)
       A(m-3,i+3) = A(m-3,i+3) - A(m-4,i+3) * A(m+1,i-1) - A(m-5,i+3) * A(m+2,i-2)
       A(m-4,i+4) = A(m-4,i+4) - A(m-5,i+4) * A(m+1,i-1)

       dummy = 1.d0 / A(m,i)
       A(m+1,i) = dummy * (A(m+1,i) - A(m-1,i) * A(m+2,i-1) - A(m-2,i) * A(m+3,i-2) - A(m-3,i) * A(m+4,i-3) - A(m-4,i)&
            * A(m+5,i-4))
       A(m+2,i) = dummy * (A(m+2,i) - A(m-1,i) * A(m+3,i-1) - A(m-2,i) * A(m+4,i-2) - A(m-3,i) * A(m+5,i-3))
       A(m+3,i) = dummy * (A(m+3,i) - A(m-1,i) * A(m+4,i-1) - A(m-2,i) * A(m+5,i-2))
       A(m+4,i) = dummy * (A(m+4,i) - A(m-1,i) * A(m+5,i-1))
       A(m+5,i) = dummy *  A(m+5,i)
    enddo

    A(m  ,n-2*bw)    = A(m  ,n-2*bw)    - A(m-1,n-2*bw)   * A(m+1,n-2*bw-1) - A(m-2,n-2*bw)   * A(m+2,n-2*bw-2) - A(m-3,n-2*bw)  &
         * A(m+3,n-2*bw-3) - A(m-4,n-2*bw)   * A(m+4,n-2*bw-4) - A(m-5,n-2*bw)   * A(m+5,n-2*bw-5)
    A(m-1,n-2*bw+1)  = A(m-1,n-2*bw+1)  - A(m-2,n-2*bw+1) * A(m+1,n-2*bw-1) - A(m-3,n-2*bw+1) * A(m+2,n-2*bw-2) - A(m-4,n-2*bw+1)&
         * A(m+3,n-2*bw-3) - A(m-5,n-2*bw+1) * A(m+4,n-2*bw-4)
    A(m-2,n-2*bw+2)  = A(m-2,n-2*bw+2)  - A(m-3,n-2*bw+2) * A(m+1,n-2*bw-1) - A(m-4,n-2*bw+2) * A(m+2,n-2*bw-2) - A(m-5,n-2*bw+2)&
         * A(m+3,n-2*bw-3)
    A(m-3,n-2*bw+3)  = A(m-3,n-2*bw+3)  - A(m-4,n-2*bw+3) * A(m+1,n-2*bw-1) - A(m-5,n-2*bw+3) * A(m+2,n-2*bw-2)
    A(m-4,n-2*bw+4)  = A(m-4,n-2*bw+4)  - A(m-5,n-2*bw+4) * A(m+1,n-2*bw-1)

    dummy = 1.d0 / A(m,n-2*bw)
    A(m+1,n-2*bw) = dummy * (A(m+1,n-2*bw) - A(m-1,n-2*bw) * A(m+2,n-2*bw-1) - A(m-2,n-2*bw) * A(m+3,n-2*bw-2) - A(m-3,n-2*bw)&
         * A(m+4,n-2*bw-3) - A(m-4,n-2*bw) * A(m+5,n-2*bw-4))
    A(m+2,n-2*bw) = dummy * (A(m+2,n-2*bw) - A(m-1,n-2*bw) * A(m+3,n-2*bw-1) - A(m-2,n-2*bw) * A(m+4,n-2*bw-2) - A(m-3,n-2*bw)&
         * A(m+5,n-2*bw-3))
    A(m+3,n-2*bw) = dummy * (A(m+3,n-2*bw) - A(m-1,n-2*bw) * A(m+4,n-2*bw-1) - A(m-2,n-2*bw) * A(m+5,n-2*bw-2))
    A(m+4,n-2*bw) = dummy * (A(m+4,n-2*bw) - A(m-1,n-2*bw) * A(m+5,n-2*bw-1))
    A(m+5,n-2*bw) = dummy *  A(m+5,n-2*bw)

    A(m  ,n-2*bw+1)   = A(m  ,n-2*bw+1)  - A(m-1,n-2*bw+1)  * A(m+1,n-2*bw) - A(m-2,n-2*bw+1) * A(m+2,n-2*bw-1) - A(m-3,n-2*bw+1)&
         * A(m+3,n-2*bw-2) - A(m-4,n-2*bw+1) * A(m+4,n-2*bw-3) - A(m-5,n-2*bw+1) * A(m+5,n-2*bw-4)
    A(m-1,n-2*bw+2)   = A(m-1,n-2*bw+2)  - A(m-2,n-2*bw+2)  * A(m+1,n-2*bw) - A(m-3,n-2*bw+2) * A(m+2,n-2*bw-1) - A(m-4,n-2*bw+2) &
         * A(m+3,n-2*bw-2) - A(m-5,n-2*bw+2) * A(m+4,n-2*bw-3)
    A(m-2,n-2*bw+3)   = A(m-2,n-2*bw+3)  - A(m-3,n-2*bw+3)  * A(m+1,n-2*bw) - A(m-4,n-2*bw+3) * A(m+2,n-2*bw-1) - A(m-5,n-2*bw+3)&
         * A(m+3,n-2*bw-2) 
    A(m-3,n-2*bw+4)   = A(m-3,n-2*bw+4)  - A(m-4,n-2*bw+4)  * A(m+1,n-2*bw) - A(m-5,n-2*bw+4) * A(m+2,n-2*bw-1)
    A(m-4,n-2*bw+5)   = A(m-4,n-2*bw+5)  - A(m-5,n-2*bw+5)  * A(m+1,n-2*bw)
    A(m-5,n-2*bw+6)   = A(m-5,n-2*bw+6)  - A(m+5,n-2*bw+6)  * A(m+1,n-2*bw)
    A(m+4,n-2*bw+7)   = A(m+4,n-2*bw+7)  - A(m+5,n-2*bw+7)  * A(m+1,n-2*bw) 
    A(m+4,n-2*bw+8)   = A(m+4,n-2*bw+8)  - A(m+5,n-2*bw+8)  * A(m+1,n-2*bw) 
    A(m+4,n-2*bw+9)   = A(m+4,n-2*bw+9)  - A(m+5,n-2*bw+9)  * A(m+1,n-2*bw) 
    A(m+4,n-2*bw+10)  = A(m+4,n-2*bw+10) - A(m+5,n-2*bw+10) * A(m+1,n-2*bw) 

    dummy = 1.d0 / A(m,n-2*bw+1)
    A(m+1,n-2*bw+1) = dummy * (A(m+1,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+2,n-2*bw) - A(m-2,n-2*bw+1) * A(m+3,n-2*bw-1) &
    - A(m-3,n-2*bw+1) * A(m+4,n-2*bw-2) - A(m-4,n-2*bw+1) * A(m+5,n-2*bw-3))
    A(m+2,n-2*bw+1) = dummy * (A(m+2,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+3,n-2*bw) - A(m-2,n-2*bw+1) * A(m+4,n-2*bw-1) &
         - A(m-3,n-2*bw+1)  * A(m+5,n-2*bw-2))
    A(m+3,n-2*bw+1) = dummy * (A(m+3,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+4,n-2*bw) - A(m-2,n-2*bw+1) * A(m+5,n-2*bw-1))
    A(m+4,n-2*bw+1) = dummy * (A(m+4,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+5,n-2*bw))
    A(m+5,n-2*bw+1) = dummy *  A(m+5,n-2*bw+1)

    A(m  ,n-2*bw+2)  = A(m  ,n-2*bw+2)  - A(m-1,n-2*bw+2)  * A(m+1,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+2,n-2*bw) - A(m-3,n-2*bw+2)&
         * A(m+3,n-2*bw-1) - A(m-4,n-2*bw+2) * A(m+4,n-2*bw-2) - A(m-5,n-2*bw+2) * A(m+5,n-2*bw-3) 
    A(m-1,n-2*bw+3)  = A(m-1,n-2*bw+3)  - A(m-2,n-2*bw+3)  * A(m+1,n-2*bw+1) - A(m-3,n-2*bw+3) * A(m+2,n-2*bw) - A(m-4,n-2*bw+3)&
         * A(m+3,n-2*bw-1) - A(m-5,n-2*bw+3) * A(m+4,n-2*bw-2)
    A(m-2,n-2*bw+4)  = A(m-2,n-2*bw+4)  - A(m-3,n-2*bw+4)  * A(m+1,n-2*bw+1) - A(m-4,n-2*bw+4) * A(m+2,n-2*bw) - A(m-5,n-2*bw+4)&
         * A(m+3,n-2*bw-1)
    A(m-3,n-2*bw+5)  = A(m-3,n-2*bw+5)  - A(m-4,n-2*bw+5)  * A(m+1,n-2*bw+1) - A(m-5,n-2*bw+5) * A(m+2,n-2*bw)
    A(m-4,n-2*bw+6)  = A(m-4,n-2*bw+6)  - A(m-5,n-2*bw+6)  * A(m+1,n-2*bw+1) - A(m+5,n-2*bw+6) * A(m+2,n-2*bw)
    A(m-5,n-2*bw+7)  = A(m-5,n-2*bw+7)  - A(m+4,n-2*bw+7)  * A(m+1,n-2*bw+1) - A(m+5,n-2*bw+7) * A(m+2,n-2*bw)
    A(m+3,n-2*bw+8)  = A(m+3,n-2*bw+8)  - A(m+4,n-2*bw+8)  * A(m+1,n-2*bw+1) - A(m+5,n-2*bw+8) * A(m+2,n-2*bw)
    A(m+3,n-2*bw+9)  = A(m+3,n-2*bw+9)  - A(m+4,n-2*bw+9)  * A(m+1,n-2*bw+1) - A(m+5,n-2*bw+9) * A(m+2,n-2*bw)
    A(m+3,n-2*bw+10) = A(m+3,n-2*bw+10) - A(m+4,n-2*bw+10) * A(m+1,n-2*bw+1) - A(m+5,n-2*bw+10) * A(m+2,n-2*bw)

    dummy = 1.d0 / A(m,n-2*bw+2)
    A(m+1,n-2*bw+2) = dummy * (A(m+1,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+2,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+3,n-2*bw) &
         - A(m-3,n-2*bw+2) * A(m+4,n-2*bw-1) - A(m-4,n-2*bw+2) * A(m+5,n-2*bw-2))
    A(m+2,n-2*bw+2) = dummy * (A(m+2,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+3,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+4,n-2*bw) &
         - A(m-3,n-2*bw+2)  * A(m+5,n-2*bw-1))
    A(m+3,n-2*bw+2) = dummy * (A(m+3,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+4,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+5,n-2*bw))
    A(m+4,n-2*bw+2) = dummy * (A(m+4,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+5,n-2*bw+1))
    A(m+5,n-2*bw+2) = dummy *  A(m+5,n-2*bw+2)

    A(m  ,n-2*bw+3) = A(m  ,n-2*bw+3)  - A(m-1,n-2*bw+3)  * A(m+1,n-2*bw+2) - A(m-2,n-2*bw+3)  * A(m+2,n-2*bw+1) - A(m-3,n-2*bw+3) &
         * A(m+3,n-2*bw) - A(m-4,n-2*bw+3) * A(m+4,n-2*bw-1) - A(m-5,n-2*bw+3) * A(m+5,n-2*bw-2)
    A(m-1,n-2*bw+4) = A(m-1,n-2*bw+4)  - A(m-2,n-2*bw+4)  * A(m+1,n-2*bw+2) - A(m-3,n-2*bw+4)  * A(m+2,n-2*bw+1) - A(m-4,n-2*bw+4) &
         * A(m+3,n-2*bw) - A(m-5,n-2*bw+4) * A(m+4,n-2*bw-1)
    A(m-2,n-2*bw+5) = A(m-2,n-2*bw+5)  - A(m-3,n-2*bw+5)  * A(m+1,n-2*bw+2) - A(m-4,n-2*bw+5)  * A(m+2,n-2*bw+1) - A(m-5,n-2*bw+5) &
         * A(m+3,n-2*bw)
    A(m-3,n-2*bw+6) = A(m-3,n-2*bw+6)  - A(m-4,n-2*bw+6)  * A(m+1,n-2*bw+2) - A(m-5,n-2*bw+6)  * A(m+2,n-2*bw+1) - A(m+5,n-2*bw+6) &
         * A(m+3,n-2*bw)
    A(m-4,n-2*bw+7) = A(m-4,n-2*bw+7)  - A(m-5,n-2*bw+7)  * A(m+1,n-2*bw+2) - A(m+4,n-2*bw+7)  * A(m+2,n-2*bw+1) - A(m+5,n-2*bw+7) &
         * A(m+3,n-2*bw)
    A(m-5,n-2*bw+8) = A(m-5,n-2*bw+8)  - A(m+3,n-2*bw+8)  * A(m+1,n-2*bw+2) - A(m+4,n-2*bw+8)  * A(m+2,n-2*bw+1) - A(m+5,n-2*bw+8) &
         * A(m+3,n-2*bw)
    A(m+2,n-2*bw+9) = A(m+2,n-2*bw+9)  - A(m+3,n-2*bw+9)  * A(m+1,n-2*bw+2) - A(m+4,n-2*bw+9)  * A(m+2,n-2*bw+1) - A(m+5,n-2*bw+9) &
         * A(m+3,n-2*bw)
    A(m+2,n-2*bw+10) = A(m+2,n-2*bw+10) - A(m+3,n-2*bw+10) * A(m+1,n-2*bw+2) - A(m+4,n-2*bw+10) * A(m+2,n-2*bw+1) &
         - A(m+5,n-2*bw+10) * A(m+3,n-2*bw)


    dummy = 1.d0 / A(m,n-2*bw+3)
    A(m+1,n-2*bw+3) = dummy * (A(m+1,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+2,n-2*bw+2) - A(m-2,n-2*bw+3) * A(m+3,n-2*bw+1) &
         - A(m-3,n-2*bw+3) * A(m+4,n-2*bw) - A(m-4,n-2*bw+3) * A(m+5,n-2*bw-1))
    A(m+2,n-2*bw+3) = dummy * (A(m+2,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+3,n-2*bw+2) - A(m-2,n-2*bw+3) * A(m+4,n-2*bw+1) &
         - A(m-3,n-2*bw+3) * A(m+5,n-2*bw))
    A(m+3,n-2*bw+3) = dummy * (A(m+3,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+4,n-2*bw+2) - A(m-2,n-2*bw+3) * A(m+5,n-2*bw+1))
    A(m+4,n-2*bw+3) = dummy * (A(m+4,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+5,n-2*bw+2))
    A(m+5,n-2*bw+3) = dummy *  A(m+5,n-2*bw+3)


    A(m  ,n-2*bw+4)  = A(m  ,n-2*bw+4)  - A(m-1,n-2*bw+4)  * A(m+1,n-2*bw+3) - A(m-2,n-2*bw+4)  * A(m+2,n-2*bw+2)&
         - A(m-3,n-2*bw+4)  * A(m+3,n-2*bw+1) - A(m-4,n-2*bw+4) * A(m+4,n-2*bw) - A(m-5,n-2*bw+4) * A(m+5,n-2*bw-1)
    A(m-1,n-2*bw+5)  = A(m-1,n-2*bw+5)  - A(m-2,n-2*bw+5)  * A(m+1,n-2*bw+3) - A(m-3,n-2*bw+5)  * A(m+2,n-2*bw+2) &
         - A(m-4,n-2*bw+5)  * A(m+3,n-2*bw+1) - A(m-5,n-2*bw+5) * A(m+4,n-2*bw)
    A(m-2,n-2*bw+6)  = A(m-2,n-2*bw+6)  - A(m-3,n-2*bw+6)  * A(m+1,n-2*bw+3) - A(m-4,n-2*bw+6)  * A(m+2,n-2*bw+2) &
         - A(m-5,n-2*bw+6)  * A(m+3,n-2*bw+1) - A(m+5,n-2*bw+6) * A(m+4,n-2*bw)
    A(m-3,n-2*bw+7)  = A(m-3,n-2*bw+7)  - A(m-4,n-2*bw+7)  * A(m+1,n-2*bw+3) - A(m-5,n-2*bw+7)  * A(m+2,n-2*bw+2) &
         - A(m+4,n-2*bw+7)  * A(m+3,n-2*bw+1) - A(m+5,n-2*bw+7) * A(m+4,n-2*bw) 
    A(m-4,n-2*bw+8)  = A(m-4,n-2*bw+8)  - A(m-5,n-2*bw+8)  * A(m+1,n-2*bw+3) - A(m+3,n-2*bw+8)  * A(m+2,n-2*bw+2) &
         - A(m+4,n-2*bw+8)  * A(m+3,n-2*bw+1) - A(m+5,n-2*bw+8) * A(m+4,n-2*bw)
    A(m-5,n-2*bw+9)  = A(m-5,n-2*bw+9)  - A(m+2,n-2*bw+9)  * A(m+1,n-2*bw+3) - A(m+3,n-2*bw+9)  * A(m+2,n-2*bw+2) &
         - A(m+4,n-2*bw+9)  * A(m+3,n-2*bw+1) - A(m+5,n-2*bw+9) * A(m+4,n-2*bw)
    A(m+1,n-2*bw+10) = A(m+1,n-2*bw+10) - A(m+2,n-2*bw+10) * A(m+1,n-2*bw+3) - A(m+3,n-2*bw+10) * A(m+2,n-2*bw+2) &
         - A(m+4,n-2*bw+10) * A(m+3,n-2*bw+1) - A(m+5,n-2*bw+10) * A(m+4,n-2*bw)


    dummy = 1.d0 / A(m,n-2*bw+4)
    A(m+1,n-2*bw+4) = dummy * (A(m+1,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+2,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+3,n-2*bw+2)&
         - A(m-3,n-2*bw+4) * A(m+4,n-2*bw+1) - A(m-4,n-2*bw+4) * A(m+5,n-2*bw))
    A(m+2,n-2*bw+4) = dummy * (A(m+2,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+3,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+4,n-2*bw+2)&
         - A(m-3,n-2*bw+4) * A(m+5,n-2*bw+1))
    A(m+3,n-2*bw+4) = dummy * (A(m+3,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+4,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+5,n-2*bw+2))
    A(m+4,n-2*bw+4) = dummy * (A(m+4,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+5,n-2*bw+3))
    A(m+5,n-2*bw+4) = dummy *  A(m+5,n-2*bw+4)

    A(m  ,n-2*bw+5)   = A(m  ,n-2*bw+5)  - A(m-1,n-2*bw+5)  * A(m+1,n-2*bw+4) - A(m-2,n-2*bw+5)  * A(m+2,n-2*bw+3) &
         - A(m-3,n-2*bw+5)  * A(m+3,n-2*bw+2) - A(m-4,n-2*bw+5)  * A(m+4,n-2*bw+1) - A(m-5,n-2*bw+5)  * A(m+5,n-2*bw)
    A(m-1,n-2*bw+6)   = A(m-1,n-2*bw+6)  - A(m-2,n-2*bw+6)  * A(m+1,n-2*bw+4) - A(m-3,n-2*bw+6)  * A(m+2,n-2*bw+3)&
         - A(m-4,n-2*bw+6)  * A(m+3,n-2*bw+2) - A(m-5,n-2*bw+6)  * A(m+4,n-2*bw+1) - A(m+5,n-2*bw+6)  * A(m+5,n-2*bw)   
    A(m-2,n-2*bw+7)   = A(m-2,n-2*bw+7)  - A(m-3,n-2*bw+7)  * A(m+1,n-2*bw+4) - A(m-4,n-2*bw+7)  * A(m+2,n-2*bw+3) &
         - A(m-5,n-2*bw+7)  * A(m+3,n-2*bw+2) - A(m+4,n-2*bw+7)  * A(m+4,n-2*bw+1) - A(m+5,n-2*bw+7)  * A(m+5,n-2*bw)    
    A(m-3,n-2*bw+8)   = A(m-3,n-2*bw+8)  - A(m-4,n-2*bw+8)  * A(m+1,n-2*bw+4) - A(m-5,n-2*bw+8)  * A(m+2,n-2*bw+3) &
         - A(m+3,n-2*bw+8)  * A(m+3,n-2*bw+2) - A(m+4,n-2*bw+8)  * A(m+4,n-2*bw+1) - A(m+5,n-2*bw+8)  * A(m+5,n-2*bw)    
    A(m-4,n-2*bw+9)   = A(m-4,n-2*bw+9)  - A(m-5,n-2*bw+9)  * A(m+1,n-2*bw+4) - A(m+2,n-2*bw+9)  * A(m+2,n-2*bw+3)&
         - A(m+3,n-2*bw+9)  * A(m+3,n-2*bw+2) - A(m+4,n-2*bw+9)  * A(m+4,n-2*bw+1) - A(m+5,n-2*bw+9)  * A(m+5,n-2*bw)    
    A(m-5,n-2*bw+10)  = A(m-5,n-2*bw+10) - A(m+1,n-2*bw+10) * A(m+1,n-2*bw+4) - A(m+2,n-2*bw+10) * A(m+2,n-2*bw+3)&
         - A(m+3,n-2*bw+10) * A(m+3,n-2*bw+2) - A(m+4,n-2*bw+10) * A(m+4,n-2*bw+1) - A(m+5,n-2*bw+10) * A(m+5,n-2*bw)    

    dummy = 1.d0 / A(m,n-2*bw+5)
    A(m+1,n-2*bw+5) = dummy * (A(m+1,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+2,n-2*bw+4) - A(m-2,n-2*bw+5) * A(m+3,n-2*bw+3)&
         - A(m-3,n-2*bw+5) * A(m+4,n-2*bw+2) - A(m-4,n-2*bw+5) * A(m+5,n-2*bw+1))
    A(m+2,n-2*bw+5) = dummy * (A(m+2,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+3,n-2*bw+4) - A(m-2,n-2*bw+5) * A(m+4,n-2*bw+3)&
         - A(m-3,n-2*bw+5) * A(m+5,n-2*bw+2))
    A(m+3,n-2*bw+5) = dummy * (A(m+3,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+4,n-2*bw+4) - A(m-2,n-2*bw+5) * A(m+5,n-2*bw+3))
    A(m+4,n-2*bw+5) = dummy * (A(m+4,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+5,n-2*bw+4))
    A(m+5,n-2*bw+5) = dummy *  A(m+5,n-2*bw+5)

    A(m  ,n-2*bw+6)   = A(m  ,n-2*bw+6)  - A(m-1,n-2*bw+6)  * A(m+1,n-2*bw+5) - A(m-2,n-2*bw+6)  * A(m+2,n-2*bw+4) &
         - A(m-3,n-2*bw+6)  * A(m+3,n-2*bw+3) - A(m-4,n-2*bw+6)  * A(m+4,n-2*bw+2) - A(m-5,n-2*bw+6)  * A(m+5,n-2*bw+1)
    A(m-1,n-2*bw+7)   = A(m-1,n-2*bw+7)  - A(m-2,n-2*bw+7)  * A(m+1,n-2*bw+5) - A(m-3,n-2*bw+7)  * A(m+2,n-2*bw+4)&
         - A(m-4,n-2*bw+7)  * A(m+3,n-2*bw+3) - A(m-5,n-2*bw+7)  * A(m+4,n-2*bw+2) - A(m+4,n-2*bw+7)  * A(m+5,n-2*bw+1)  
    A(m-2,n-2*bw+8)   = A(m-2,n-2*bw+8)  - A(m-3,n-2*bw+8)  * A(m+1,n-2*bw+5) - A(m-4,n-2*bw+8)  * A(m+2,n-2*bw+4)&
         - A(m-5,n-2*bw+8)  * A(m+3,n-2*bw+3) - A(m+3,n-2*bw+8)  * A(m+4,n-2*bw+2) - A(m+4,n-2*bw+8)  * A(m+5,n-2*bw+1)   
    A(m-3,n-2*bw+9)   = A(m-3,n-2*bw+9)  - A(m-4,n-2*bw+9)  * A(m+1,n-2*bw+5) - A(m-5,n-2*bw+9)  * A(m+2,n-2*bw+4)&
         - A(m+2,n-2*bw+9)  * A(m+3,n-2*bw+3) - A(m+3,n-2*bw+9)  * A(m+4,n-2*bw+2) - A(m+4,n-2*bw+9)  * A(m+5,n-2*bw+1)    
    A(m-4,n-2*bw+10)  = A(m-4,n-2*bw+10) - A(m-5,n-2*bw+10) * A(m+1,n-2*bw+5) - A(m+1,n-2*bw+10) * A(m+2,n-2*bw+4)&
         - A(m+2,n-2*bw+10) * A(m+3,n-2*bw+3) - A(m+3,n-2*bw+10) * A(m+4,n-2*bw+2) - A(m+4,n-2*bw+10) * A(m+5,n-2*bw+1)

    dummy = 1.d0 / A(m,n-2*bw+6)
    A(m+1,n-2*bw+6) = dummy * (A(m+1,n-2*bw+6) - A(m-1,n-2*bw+6) * A(m+2,n-2*bw+5) - A(m-2,n-2*bw+6) * A(m+3,n-2*bw+4) &
         - A(m-3,n-2*bw+6) * A(m+4,n-2*bw+3) - A(m-4,n-2*bw+6) * A(m+5,n-2*bw+2))
    A(m+2,n-2*bw+6) = dummy * (A(m+2,n-2*bw+6) - A(m-1,n-2*bw+6) * A(m+3,n-2*bw+5) - A(m-2,n-2*bw+6) * A(m+4,n-2*bw+4)&
         - A(m-3,n-2*bw+6) * A(m+5,n-2*bw+3))
    A(m+3,n-2*bw+6) = dummy * (A(m+3,n-2*bw+6) - A(m-1,n-2*bw+6) * A(m+4,n-2*bw+5) - A(m-2,n-2*bw+6) * A(m+5,n-2*bw+4))
    A(m+4,n-2*bw+6) = dummy * (A(m+4,n-2*bw+6) - A(m-1,n-2*bw+6) * A(m+5,n-2*bw+5))


    A(m  ,n-2*bw+7)   = A(m  ,n-2*bw+7)  - A(m-1,n-2*bw+7)  * A(m+1,n-2*bw+6) - A(m-2,n-2*bw+7)  * A(m+2,n-2*bw+5)&
         - A(m-3,n-2*bw+7)  * A(m+3,n-2*bw+4) - A(m-4,n-2*bw+7)  * A(m+4,n-2*bw+3) - A(m-5,n-2*bw+7)  * A(m+5,n-2*bw+2)     
    A(m-1,n-2*bw+8)   = A(m-1,n-2*bw+8)  - A(m-2,n-2*bw+8)  * A(m+1,n-2*bw+6) - A(m-3,n-2*bw+8)  * A(m+2,n-2*bw+5) &
         - A(m-4,n-2*bw+8)  * A(m+3,n-2*bw+4) - A(m-5,n-2*bw+8)  * A(m+4,n-2*bw+3) - A(m+3,n-2*bw+8)  * A(m+5,n-2*bw+2)
    A(m-2,n-2*bw+9)   = A(m-2,n-2*bw+9)  - A(m-3,n-2*bw+9)  * A(m+1,n-2*bw+6) - A(m-4,n-2*bw+9)  * A(m+2,n-2*bw+5)&
         - A(m-5,n-2*bw+9)  * A(m+3,n-2*bw+4) - A(m+2,n-2*bw+9)  * A(m+4,n-2*bw+3) - A(m+3,n-2*bw+9)  * A(m+5,n-2*bw+2)
    A(m-3,n-2*bw+10)  = A(m-3,n-2*bw+10) - A(m-4,n-2*bw+10) * A(m+1,n-2*bw+6) - A(m-5,n-2*bw+10) * A(m+2,n-2*bw+5) &
         - A(m+1,n-2*bw+10) * A(m+3,n-2*bw+4) - A(m+2,n-2*bw+10) * A(m+4,n-2*bw+3) - A(m+3,n-2*bw+10) * A(m+5,n-2*bw+2)


    dummy = 1.d0 / A(m,n-2*bw+7)
    A(m+1,n-2*bw+7) = dummy * (A(m+1,n-2*bw+7) - A(m-1,n-2*bw+7) * A(m+2,n-2*bw+6) - A(m-2,n-2*bw+7) * A(m+3,n-2*bw+5) &
         - A(m-3,n-2*bw+7) * A(m+4,n-2*bw+4) - A(m-4,n-2*bw+7) * A(m+5,n-2*bw+3))
    A(m+2,n-2*bw+7) = dummy * (A(m+2,n-2*bw+7) - A(m-1,n-2*bw+7) * A(m+3,n-2*bw+6) - A(m-2,n-2*bw+7) * A(m+4,n-2*bw+5)&
         - A(m-3,n-2*bw+7) * A(m+5,n-2*bw+4))
    A(m+3,n-2*bw+7) = dummy * (A(m+3,n-2*bw+7) - A(m-1,n-2*bw+7) * A(m+4,n-2*bw+6) - A(m-2,n-2*bw+7) * A(m+5,n-2*bw+5))

    A(m  ,n-2*bw+8)   = A(m  ,n-2*bw+8)  - A(m-1,n-2*bw+8)  * A(m+1,n-2*bw+7) - A(m-2,n-2*bw+8)  * A(m+2,n-2*bw+6) &
         - A(m-3,n-2*bw+8)  * A(m+3,n-2*bw+5) - A(m-4,n-2*bw+8)  * A(m+4,n-2*bw+4) - A(m-5,n-2*bw+8)  * A(m+5,n-2*bw+3) 
    A(m-1,n-2*bw+9)   = A(m-1,n-2*bw+9)  - A(m-2,n-2*bw+9)  * A(m+1,n-2*bw+7) - A(m-3,n-2*bw+9)  * A(m+2,n-2*bw+6) &
         - A(m-4,n-2*bw+9)  * A(m+3,n-2*bw+5) - A(m-5,n-2*bw+9)  * A(m+4,n-2*bw+4) - A(m+2,n-2*bw+9)  * A(m+5,n-2*bw+3) 
    A(m-2,n-2*bw+10)  = A(m-2,n-2*bw+10) - A(m-3,n-2*bw+10) * A(m+1,n-2*bw+7) - A(m-4,n-2*bw+10) * A(m+2,n-2*bw+6) &
         - A(m-5,n-2*bw+10) * A(m+3,n-2*bw+5) - A(m+1,n-2*bw+10) * A(m+4,n-2*bw+4) - A(m+2,n-2*bw+10) * A(m+5,n-2*bw+3) 

    dummy = 1.d0 / A(m,n-2*bw+8)
    A(m+1,n-2*bw+8) = dummy * (A(m+1,n-2*bw+8) - A(m-1,n-2*bw+8) * A(m+2,n-2*bw+7) - A(m-2,n-2*bw+8) * A(m+3,n-2*bw+6)&
         - A(m-3,n-2*bw+8) * A(m+4,n-2*bw+5) - A(m-4,n-2*bw+8) * A(m+5,n-2*bw+4))
    A(m+2,n-2*bw+8) = dummy * (A(m+2,n-2*bw+8) - A(m-1,n-2*bw+8) * A(m+3,n-2*bw+7) - A(m-2,n-2*bw+8) * A(m+4,n-2*bw+6) &
         - A(m-3,n-2*bw+8) * A(m+5,n-2*bw+5))

    A(m  ,n-2*bw+9)   = A(m  ,n-2*bw+9)  - A(m-1,n-2*bw+9)  * A(m+1,n-2*bw+8) - A(m-2,n-2*bw+9)  * A(m+2,n-2*bw+7) &
         - A(m-3,n-2*bw+9)  * A(m+3,n-2*bw+6) - A(m-4,n-2*bw+9)  * A(m+4,n-2*bw+5) - A(m-5,n-2*bw+9)  * A(m+5,n-2*bw+4)
    A(m-1,n-2*bw+10)  = A(m-1,n-2*bw+10) - A(m-2,n-2*bw+10) * A(m+1,n-2*bw+8) - A(m-3,n-2*bw+10) * A(m+2,n-2*bw+7)&
         - A(m-4,n-2*bw+10) * A(m+3,n-2*bw+6) - A(m-5,n-2*bw+10) * A(m+4,n-2*bw+5) - A(m+1,n-2*bw+10) * A(m+5,n-2*bw+4)

    dummy = 1.d0 / A(m,n-2*bw+9)
    A(m+1,n-2*bw+9) = dummy * (A(m+1,n-2*bw+9) - A(m-1,n-2*bw+9) * A(m+2,n-2*bw+8) - A(m-2,n-2*bw+9) * A(m+3,n-2*bw+7)&
         - A(m-3,n-2*bw+9) * A(m+4,n-2*bw+6) - A(m-4,n-2*bw+9) * A(m+5,n-2*bw+5))

    A(m  ,n-2*bw+10)  = A(m  ,n-2*bw+10) - A(m-1,n-2*bw+10) * A(m+1,n-2*bw+9) - A(m-2,n-2*bw+10) * A(m+2,n-2*bw+8) &
         - A(m-3,n-2*bw+10) * A(m+3,n-2*bw+7) - A(m-4,n-2*bw+10) * A(m+4,n-2*bw+6) - A(m-5,n-2*bw+10) * A(m+5,n-2*bw+5)

    return

  end subroutine fixed_bandec5


  subroutine fixed_bandec6(A,n)
    implicit none
    integer :: n
    real(c_double) :: A(13,n)
    integer :: i
    integer :: m,bw
    real(c_double) :: dummy

    m = 7
    bw = 6

    dummy = 1.d0 / A(m,1)
    A(m+1,1) = dummy * A(m+1,1)
    A(m+2,1) = dummy * A(m+2,1)
    A(m+3,1) = dummy * A(m+3,1)
    A(m+4,1) = dummy * A(m+4,1)
    A(m+5,1) = dummy * A(m+5,1)
    A(m+6,1) = dummy * A(m+6,1)
    A(m-1,1) = dummy * A(m-1,1)
    A(m-2,1) = dummy * A(m-2,1)
    A(m-3,1) = dummy * A(m-3,1)
    A(m-4,1) = dummy * A(m-4,1)
    A(m-5,1) = dummy * A(m-5,1)
    A(m-6,1) = dummy * A(m-6,1)

    A(m  ,2) = A(m  ,2) -  A(m-1,2) * A(m+1,1)
    A(m-1,3) = A(m-1,3) -  A(m-2,3) * A(m+1,1)
    A(m-2,4) = A(m-2,4) -  A(m-3,4) * A(m+1,1)
    A(m-3,5) = A(m-3,5) -  A(m-4,5) * A(m+1,1)
    A(m-4,6) = A(m-4,6) -  A(m-5,6) * A(m+1,1)
    A(m-5,7) = A(m-5,7) -  A(m-6,7) * A(m+1,1)

    dummy = 1.d0 / A(m,2)
    A(m+1,2) = dummy * (A(m+1,2) -  A(m-1,2) * A(m+2,1)) 
    A(m+2,2) = dummy * (A(m+2,2) -  A(m-1,2) * A(m+3,1)) 
    A(m+3,2) = dummy * (A(m+3,2) -  A(m-1,2) * A(m+4,1)) 
    A(m+4,2) = dummy * (A(m+4,2) -  A(m-1,2) * A(m+5,1)) 
    A(m+5,2) = dummy * (A(m+5,2) -  A(m-1,2) * A(m+6,1)) 
    A(m+6,2) = dummy * (A(m+6,2) -  A(m-1,2) * A(m-1,1)) 
    A(m-2,2) = dummy * (A(m-2,2) -  A(m-1,2) * A(m-2,1)) 
    A(m-3,2) = dummy * (A(m-3,2) -  A(m-1,2) * A(m-3,1)) 
    A(m-4,2) = dummy * (A(m-4,2) -  A(m-1,2) * A(m-4,1)) 
    A(m-5,2) = dummy * (A(m-5,2) -  A(m-1,2) * A(m-5,1)) 
    A(m-6,2) = dummy * (A(m-6,2) -  A(m-1,2) * A(m-6,1)) 

    A(m  , 3) = A(m  , 3) - A(m-1,3) * A(m+1,2) - A(m-2,3) * A(m+2,1) 
    A(m-1, 4) = A(m-1, 4) - A(m-2,4) * A(m+1,2) - A(m-3,4) * A(m+2,1) 
    A(m-2, 5) = A(m-2, 5) - A(m-3,5) * A(m+1,2) - A(m-4,5) * A(m+2,1)
    A(m-3, 6) = A(m-3, 6) - A(m-4,6) * A(m+1,2) - A(m-5,6) * A(m+2,1)
    A(m-4, 7) = A(m-4, 7) - A(m-5,7) * A(m+1,2) - A(m-6,7) * A(m+2,1)
    A(m-5, 8) = A(m-5, 8) - A(m-6,8) * A(m+1,2)

    dummy = 1.d0 / A(m,3)
    A(m+1,3) = dummy * (A(m+1,3) - A(m-1,3) * A(m+2,2) - A(m-2,3) * A(m+3,1))
    A(m+2,3) = dummy * (A(m+2,3) - A(m-1,3) * A(m+3,2) - A(m-2,3) * A(m+4,1))
    A(m+3,3) = dummy * (A(m+3,3) - A(m-1,3) * A(m+4,2) - A(m-2,3) * A(m+5,1))
    A(m+4,3) = dummy * (A(m+4,3) - A(m-1,3) * A(m+5,2) - A(m-2,3) * A(m+6,1))
    A(m+5,3) = dummy * (A(m+5,3) - A(m-1,3) * A(m+6,2) - A(m-2,3) * A(m-1,1))
    A(m+6,3) = dummy * (A(m+6,3) - A(m-1,3) * A(m-2,2) - A(m-2,3) * A(m-2,1))
    A(m-3,3) = dummy * (A(m-3,3) - A(m-1,3) * A(m-3,2) - A(m-2,3) * A(m-3,1))
    A(m-4,3) = dummy * (A(m-4,3) - A(m-1,3) * A(m-4,2) - A(m-2,3) * A(m-4,1))
    A(m-5,3) = dummy * (A(m-5,3) - A(m-1,3) * A(m-5,2) - A(m-2,3) * A(m-5,1))
    A(m-6,3) = dummy * (A(m-6,3) - A(m-1,3) * A(m-6,2) - A(m-2,3) * A(m-6,1))

    A(m  , 4) = A(m  , 4) - A(m-1, 4) * A(m+1,3) - A(m-2,4) * A(m+2,2) - A(m-3,4) * A(m+3,1)  
    A(m-1, 5) = A(m-1, 5) - A(m-2, 5) * A(m+1,3) - A(m-3,5) * A(m+2,2) - A(m-4,5) * A(m+3,1)
    A(m-2, 6) = A(m-2, 6) - A(m-3, 6) * A(m+1,3) - A(m-4,6) * A(m+2,2) - A(m-5,6) * A(m+3,1)
    A(m-3, 7) = A(m-3, 7) - A(m-4, 7) * A(m+1,3) - A(m-5,7) * A(m+2,2) - A(m-6,7) * A(m+3,1)
    A(m-4, 8) = A(m-4, 8) - A(m-5, 8) * A(m+1,3) - A(m-6,8) * A(m+2,2)
    A(m-5, 9) = A(m-5, 9) - A(m-6, 9) * A(m+1,3) 

    dummy = 1.d0 / A(m,4)
    A(m+1,4) = dummy * (A(m+1,4) - A(m-1,4) * A(m+2,3) - A(m-2,4) * A(m+3,2) - A(m-3,4) * A(m+4,1))
    A(m+2,4) = dummy * (A(m+2,4) - A(m-1,4) * A(m+3,3) - A(m-2,4) * A(m+4,2) - A(m-3,4) * A(m+5,1))
    A(m+3,4) = dummy * (A(m+3,4) - A(m-1,4) * A(m+4,3) - A(m-2,4) * A(m+5,2) - A(m-3,4) * A(m+6,1))
    A(m+4,4) = dummy * (A(m+4,4) - A(m-1,4) * A(m+5,3) - A(m-2,4) * A(m+6,2) - A(m-3,4) * A(m-1,1))
    A(m+5,4) = dummy * (A(m+5,4) - A(m-1,4) * A(m+6,3) - A(m-2,4) * A(m-2,2) - A(m-3,4) * A(m-2,1))
    A(m+6,4) = dummy * (A(m+6,4) - A(m-1,4) * A(m-3,3) - A(m-2,4) * A(m-3,2) - A(m-3,4) * A(m-3,1))
    A(m-4,4) = dummy * (A(m-4,4) - A(m-1,4) * A(m-4,3) - A(m-2,4) * A(m-4,2) - A(m-3,4) * A(m-4,1))
    A(m-5,4) = dummy * (A(m-5,4) - A(m-1,4) * A(m-5,3) - A(m-2,4) * A(m-5,2) - A(m-3,4) * A(m-5,1))
    A(m-6,4) = dummy * (A(m-6,4) - A(m-1,4) * A(m-6,3) - A(m-2,4) * A(m-6,2) - A(m-3,4) * A(m-6,1))

    A(m  , 5) = A(m  , 5) - A(m-1, 5) * A(m+1,4) - A(m-2, 5) * A(m+2,3) - A(m-3,5) * A(m+3,2) - A(m-4,5) * A(m+4,1)    
    A(m-1, 6) = A(m-1, 6) - A(m-2, 6) * A(m+1,4) - A(m-3, 6) * A(m+2,3) - A(m-4,6) * A(m+3,2) - A(m-5,6) * A(m+4,1)    
    A(m-2, 7) = A(m-2, 7) - A(m-3, 7) * A(m+1,4) - A(m-4, 7) * A(m+2,3) - A(m-5,7) * A(m+3,2) - A(m-6,7) * A(m+4,1) 
    A(m-3, 8) = A(m-3, 8) - A(m-4, 8) * A(m+1,4) - A(m-5, 8) * A(m+2,3) - A(m-6,8) * A(m+3,2) 
    A(m-4, 9) = A(m-4, 9) - A(m-5, 9) * A(m+1,4) - A(m-6, 9) * A(m+2,3)
    A(m-5,10) = A(m-5,10) - A(m-6,10) * A(m+1,4)

    dummy = 1.d0 / A(m,5)
    A(m+1,5) = dummy * (A(m+1,5) - A(m-1,5) * A(m+2,4) - A(m-2,5) * A(m+3,3) - A(m-3,5) * A(m+4,2) - A(m-4,5) * A(m+5,1))
    A(m+2,5) = dummy * (A(m+2,5) - A(m-1,5) * A(m+3,4) - A(m-2,5) * A(m+4,3) - A(m-3,5) * A(m+5,2) - A(m-4,5) * A(m+6,1))
    A(m+3,5) = dummy * (A(m+3,5) - A(m-1,5) * A(m+4,4) - A(m-2,5) * A(m+5,3) - A(m-3,5) * A(m+6,2) - A(m-4,5) * A(m-1,1))
    A(m+4,5) = dummy * (A(m+4,5) - A(m-1,5) * A(m+5,4) - A(m-2,5) * A(m+6,3) - A(m-3,5) * A(m-2,2) - A(m-4,5) * A(m-2,1))
    A(m+5,5) = dummy * (A(m+5,5) - A(m-1,5) * A(m+6,4) - A(m-2,5) * A(m-3,3) - A(m-3,5) * A(m-3,2) - A(m-4,5) * A(m-3,1))
    A(m+6,5) = dummy * (A(m+6,5) - A(m-1,5) * A(m-4,4) - A(m-2,5) * A(m-4,3) - A(m-3,5) * A(m-4,2) - A(m-4,5) * A(m-4,1))
    A(m-5,5) = dummy * (A(m-5,5) - A(m-1,5) * A(m-5,4) - A(m-2,5) * A(m-5,3) - A(m-3,5) * A(m-5,2) - A(m-4,5) * A(m-5,1))
    A(m-6,5) = dummy * (A(m-6,5) - A(m-1,5) * A(m-6,4) - A(m-2,5) * A(m-6,3) - A(m-3,5) * A(m-6,2) - A(m-4,5) * A(m-6,1))

    A(m  , 6) = A(m  , 6) - A(m-1, 6) * A(m+1,5) - A(m-2, 6) * A(m+2,4) - A(m-3, 6) * A(m+3,3) - A(m-4,6) * A(m+4,2) &
         - A(m-5,6) * A(m+5,1) 
    A(m-1, 7) = A(m-1, 7) - A(m-2, 7) * A(m+1,5) - A(m-3, 7) * A(m+2,4) - A(m-4, 7) * A(m+3,3) - A(m-5,7) * A(m+4,2)&
         - A(m-6,7) * A(m+5,1) 
    A(m-2, 8) = A(m-2, 8) - A(m-3, 8) * A(m+1,5) - A(m-4, 8) * A(m+2,4) - A(m-5, 8) * A(m+3,3) - A(m-6,8) * A(m+4,2)
    A(m-3, 9) = A(m-3, 9) - A(m-4, 9) * A(m+1,5) - A(m-5, 9) * A(m+2,4) - A(m-6, 9) * A(m+3,3) 
    A(m-4,10) = A(m-4,10) - A(m-5,10) * A(m+1,5) - A(m-6,10) * A(m+2,4)
    A(m-5,11) = A(m-5,11) - A(m-6,11) * A(m+1,5)

    dummy = 1.d0 / A(m,6)
    A(m+1,6) = dummy * (A(m+1,6) - A(m-1,6) * A(m+2,5) - A(m-2,6) * A(m+3,4) - A(m-3,6) * A(m+4,3) - A(m-4,6) * A(m+5,2) &
         - A(m-5,6) * A(m+6,1))
    A(m+2,6) = dummy * (A(m+2,6) - A(m-1,6) * A(m+3,5) - A(m-2,6) * A(m+4,4) - A(m-3,6) * A(m+5,3) - A(m-4,6) * A(m+6,2) &
         - A(m-5,6) * A(m-1,1))
    A(m+3,6) = dummy * (A(m+3,6) - A(m-1,6) * A(m+4,5) - A(m-2,6) * A(m+5,4) - A(m-3,6) * A(m+6,3) - A(m-4,6) * A(m-2,2)&
         - A(m-5,6) * A(m-2,1))
    A(m+4,6) = dummy * (A(m+4,6) - A(m-1,6) * A(m+5,5) - A(m-2,6) * A(m+6,4) - A(m-3,6) * A(m-3,3) - A(m-4,6) * A(m-3,2)&
         - A(m-5,6) * A(m-3,1))
    A(m+5,6) = dummy * (A(m+5,6) - A(m-1,6) * A(m+6,5) - A(m-2,6) * A(m-4,4) - A(m-3,6) * A(m-4,3) - A(m-4,6) * A(m-4,2)&
         - A(m-5,6) * A(m-4,1))
    A(m+6,6) = dummy * (A(m+6,6) - A(m-1,6) * A(m-5,5) - A(m-2,6) * A(m-5,4) - A(m-3,6) * A(m-5,3) - A(m-4,6) * A(m-5,2)&
         - A(m-5,6) * A(m-5,1))
    A(m-6,6) = dummy * (A(m-6,6) - A(m-1,6) * A(m-6,5) - A(m-2,6) * A(m-6,4) - A(m-3,6) * A(m-6,3) - A(m-4,6) * A(m-6,2)&
         - A(m-5,6) * A(m-6,1))


    ! Subupper Region

    A(m  , 7) = A(m  , 7) - A(m-1, 7) * A(m+1,6) - A(m-2, 7) * A(m+2,5) - A(m-3, 7) * A(m+3,4) - A(m-4, 7) * A(m+4,3)&
         - A(m-5,7) * A(m+5,2) - A(m-6,7) * A(m+6,1) 
    A(m-1, 8) = A(m-1, 8) - A(m-2, 8) * A(m+1,6) - A(m-3, 8) * A(m+2,5) - A(m-4, 8) * A(m+3,4) - A(m-5, 8) * A(m+4,3)&
         - A(m-6,8) * A(m+5,2)
    A(m-2, 9) = A(m-2, 9) - A(m-3, 9) * A(m+1,6) - A(m-4, 9) * A(m+2,5) - A(m-5, 9) * A(m+3,4) - A(m-6, 9) * A(m+4,3)
    A(m-3,10) = A(m-3,10) - A(m-4,10) * A(m+1,6) - A(m-5,10) * A(m+2,5) - A(m-6,10) * A(m+3,4) 
    A(m-4,11) = A(m-4,11) - A(m-5,11) * A(m+1,6) - A(m-6,11) * A(m+2,5)
    A(m-5,12) = A(m-5,12) - A(m-6,12) * A(m+1,6)

    dummy = 1.d0 / A(m,7)
    A(m+1,7) = dummy * (A(m+1,7) - A(m-1,7) * A(m+2,6) - A(m-2,7) * A(m+3,5) - A(m-3,7) * A(m+4,4) - A(m-4,7) * A(m+5,3)&
         - A(m-5,7) * A(m+6,2) - A(m-6,7) * A(m-1,1))
    A(m+2,7) = dummy * (A(m+2,7) - A(m-1,7) * A(m+3,6) - A(m-2,7) * A(m+4,5) - A(m-3,7) * A(m+5,4) - A(m-4,7) * A(m+6,3)&
         - A(m-5,7) * A(m-2,2) - A(m-6,7) * A(m-2,1))
    A(m+3,7) = dummy * (A(m+3,7) - A(m-1,7) * A(m+4,6) - A(m-2,7) * A(m+5,5) - A(m-3,7) * A(m+6,4) - A(m-4,7) * A(m-3,3)&
         - A(m-5,7) * A(m-3,2) - A(m-6,7) * A(m-3,1))
    A(m+4,7) = dummy * (A(m+4,7) - A(m-1,7) * A(m+5,6) - A(m-2,7) * A(m+6,5) - A(m-3,7) * A(m-4,4) - A(m-4,7) * A(m-4,3) &
         - A(m-5,7) * A(m-4,2) - A(m-6,7) * A(m-4,1))
    A(m+5,7) = dummy * (A(m+5,7) - A(m-1,7) * A(m+6,6) - A(m-2,7) * A(m-5,5) - A(m-3,7) * A(m-5,4) - A(m-4,7) * A(m-5,3)&
         - A(m-5,7) * A(m-5,2) - A(m-6,7) * A(m-5,1))
    A(m+6,7) = dummy * (A(m+6,7) - A(m-1,7) * A(m-6,6) - A(m-2,7) * A(m-6,5) - A(m-3,7) * A(m-6,4) - A(m-4,7) * A(m-6,3)&
         - A(m-5,7) * A(m-6,2) - A(m-6,7) * A(m-6,1))

    A(m  , 8) = A(m  , 8) - A(m-1, 8) * A(m+1,7) - A(m-2, 8) * A(m+2,6) - A(m-3, 8) * A(m+3,5) - A(m-4, 8) * A(m+4,4) &
         - A(m-5, 8) * A(m+5,3) - A(m-6,8) * A(m+6,2)
    A(m-1, 9) = A(m-1, 9) - A(m-2, 9) * A(m+1,7) - A(m-3, 9) * A(m+2,6) - A(m-4, 9) * A(m+3,5) - A(m-5, 9) * A(m+4,4)&
         - A(m-6, 9) * A(m+5,3)
    A(m-2,10) = A(m-2,10) - A(m-3,10) * A(m+1,7) - A(m-4,10) * A(m+2,6) - A(m-5,10) * A(m+3,5) - A(m-6,10) * A(m+4,4)
    A(m-3,11) = A(m-3,11) - A(m-4,11) * A(m+1,7) - A(m-5,11) * A(m+2,6) - A(m-6,11) * A(m+3,5)
    A(m-4,12) = A(m-4,12) - A(m-5,12) * A(m+1,7) - A(m-6,12) * A(m+2,6)
    A(m-5,13) = A(m-5,13) - A(m-6,13) * A(m+1,7) 


    dummy = 1.d0 / A(m,8)
    A(m+1,8) = dummy * (A(m+1,8) - A(m-1,8) * A(m+2,7) - A(m-2,8) * A(m+3,6) - A(m-3,8) * A(m+4,5) - A(m-4,8) * A(m+5,4)&
         - A(m-5,8) * A(m+6,3) - A(m-6,8) * A(m-2,2))
    A(m+2,8) = dummy * (A(m+2,8) - A(m-1,8) * A(m+3,7) - A(m-2,8) * A(m+4,6) - A(m-3,8) * A(m+5,5) - A(m-4,8) * A(m+6,4)&
         - A(m-5,8) * A(m-3,3) - A(m-6,8) * A(m-3,2))
    A(m+3,8) = dummy * (A(m+3,8) - A(m-1,8) * A(m+4,7) - A(m-2,8) * A(m+5,6) - A(m-3,8) * A(m+6,5) - A(m-4,8) * A(m-4,4)&
         - A(m-5,8) * A(m-4,3) - A(m-6,8) * A(m-4,2))
    A(m+4,8) = dummy * (A(m+4,8) - A(m-1,8) * A(m+5,7) - A(m-2,8) * A(m+6,6) - A(m-3,8) * A(m-5,5) - A(m-4,8) * A(m-5,4)&
         - A(m-5,8) * A(m-5,3) - A(m-6,8) * A(m-5,2))
    A(m+5,8) = dummy * (A(m+5,8) - A(m-1,8) * A(m+6,7) - A(m-2,8) * A(m-6,6) - A(m-3,8) * A(m-6,5) - A(m-4,8) * A(m-6,4)&
         - A(m-5,8) * A(m-6,3) - A(m-6,8) * A(m-6,2)) 
    A(m+6,8) = dummy *  A(m+6,8) 

    A(m  , 9) = A(m  , 9) - A(m-1, 9) * A(m+1,8) - A(m-2, 9) * A(m+2,7) - A(m-3, 9) * A(m+3,6) - A(m-4, 9) * A(m+4,5) &
         - A(m-5, 9) * A(m+5,4) - A(m-6, 9) * A(m+6,3)
    A(m-1,10) = A(m-1,10) - A(m-2,10) * A(m+1,8) - A(m-3,10) * A(m+2,7) - A(m-4,10) * A(m+3,6) - A(m-5,10) * A(m+4,5) &
         - A(m-6,10) * A(m+5,4)
    A(m-2,11) = A(m-2,11) - A(m-3,11) * A(m+1,8) - A(m-4,11) * A(m+2,7) - A(m-5,11) * A(m+3,6) - A(m-6,11) * A(m+4,5)
    A(m-3,12) = A(m-3,12) - A(m-4,12) * A(m+1,8) - A(m-5,12) * A(m+2,7) - A(m-6,12) * A(m+3,6)
    A(m-4,13) = A(m-4,13) - A(m-5,13) * A(m+1,8) - A(m-6,13) * A(m+2,7)
    A(m-5,14) = A(m-5,14) - A(m-6,14) * A(m+1,8) 


    dummy = 1.d0 / A(m,9)
    A(m+1,9) = dummy * (A(m+1,9) - A(m-1,9) * A(m+2,8) - A(m-2,9) * A(m+3,7) - A(m-3,9) * A(m+4,6) - A(m-4,9) * A(m+5,5) &
         - A(m-5,9) * A(m+6,4) - A(m-6,9) * A(m-3,3))
    A(m+2,9) = dummy * (A(m+2,9) - A(m-1,9) * A(m+3,8) - A(m-2,9) * A(m+4,7) - A(m-3,9) * A(m+5,6) - A(m-4,9) * A(m+6,5) &
         - A(m-5,9) * A(m-4,4) - A(m-6,9) * A(m-4,3))
    A(m+3,9) = dummy * (A(m+3,9) - A(m-1,9) * A(m+4,8) - A(m-2,9) * A(m+5,7) - A(m-3,9) * A(m+6,6) - A(m-4,9) * A(m-5,5) &
         - A(m-5,9) * A(m-5,4) - A(m-6,9) * A(m-5,3))
    A(m+4,9) = dummy * (A(m+4,9) - A(m-1,9) * A(m+5,8) - A(m-2,9) * A(m+6,7) - A(m-3,9) * A(m-6,6) - A(m-4,9) * A(m-6,5) &
         - A(m-5,9) * A(m-6,4) - A(m-6,9) * A(m-6,3))
    A(m+5,9) = dummy * (A(m+5,9) - A(m-1,9) * A(m+6,8))
    A(m+6,9) = dummy *  A(m+6,9) 

    A(m  ,10) = A(m  ,10) - A(m-1,10) * A(m+1,9) - A(m-2,10) * A(m+2,8) - A(m-3,10) * A(m+3,7) - A(m-4,10) * A(m+4,6) &
         - A(m-5,10) * A(m+5,5) - A(m-6,10) * A(m+6,4)
    A(m-1,11) = A(m-1,11) - A(m-2,11) * A(m+1,9) - A(m-3,11) * A(m+2,8) - A(m-4,11) * A(m+3,7) - A(m-5,11) * A(m+4,6)&
         - A(m-6,11) * A(m+5,5)
    A(m-2,12) = A(m-2,12) - A(m-3,12) * A(m+1,9) - A(m-4,12) * A(m+2,8) - A(m-5,12) * A(m+3,7) - A(m-6,12) * A(m+4,6)
    A(m-3,13) = A(m-3,13) - A(m-4,13) * A(m+1,9) - A(m-5,13) * A(m+2,8) - A(m-6,13) * A(m+3,7)
    A(m-4,14) = A(m-4,14) - A(m-5,14) * A(m+1,9) - A(m-6,14) * A(m+2,8)
    A(m-5,15) = A(m-5,15) - A(m-6,15) * A(m+1,9) 
    A(m-6,16) = A(m-6,16) 

    dummy = 1.d0 / A(m,10)
    A(m+1,10) = dummy * (A(m+1,10) - A(m-1,10) * A(m+2,9) - A(m-2,10) * A(m+3,8) - A(m-3,10) * A(m+4,7) - A(m-4,10) * A(m+5,6)&
         - A(m-5,10) * A(m+6,5) - A(m-6,10) * A(m-4,4))
    A(m+2,10) = dummy * (A(m+2,10) - A(m-1,10) * A(m+3,9) - A(m-2,10) * A(m+4,8) - A(m-3,10) * A(m+5,7) - A(m-4,10) * A(m+6,6)&
         - A(m-5,10) * A(m-5,5) - A(m-6,10) * A(m-5,4))
    A(m+3,10) = dummy * (A(m+3,10) - A(m-1,10) * A(m+4,9) - A(m-2,10) * A(m+5,8) - A(m-3,10) * A(m+6,7) - A(m-4,10) * A(m-6,6)&
         - A(m-5,10) * A(m-6,5) - A(m-6,10) * A(m-6,4))
    A(m+4,10) = dummy * (A(m+4,10) - A(m-1,10) * A(m+5,9) - A(m-2,10) * A(m+6,8))
    A(m+5,10) = dummy * (A(m+5,10) - A(m-1,10) * A(m+6,9))
    A(m+6,10) = dummy *  A(m+6,10)

    A(m  ,11) = A(m  ,11) - A(m-1,11) * A(m+1,10) - A(m-2,11) * A(m+2,9) - A(m-3,11) * A(m+3,8) - A(m-4,11) * A(m+4,7)&
         - A(m-5,11) * A(m+5,6) - A(m-6,11) * A(m+6,5)
    A(m-1,12) = A(m-1,12) - A(m-2,12) * A(m+1,10) - A(m-3,12) * A(m+2,9) - A(m-4,12) * A(m+3,8) - A(m-5,12) * A(m+4,7)&
         - A(m-6,12) * A(m+5,6)
    A(m-2,13) = A(m-2,13) - A(m-3,13) * A(m+1,10) - A(m-4,13) * A(m+2,9) - A(m-5,13) * A(m+3,8) - A(m-6,13) * A(m+4,7)
    A(m-3,14) = A(m-3,14) - A(m-4,14) * A(m+1,10) - A(m-5,14) * A(m+2,9) - A(m-6,14) * A(m+3,8)
    A(m-4,15) = A(m-4,15) - A(m-5,15) * A(m+1,10) - A(m-6,15) * A(m+2,9)
    A(m-5,16) = A(m-5,16) - A(m-6,16) * A(m+1,10) 
    A(m-6,17) = A(m-6,17) 

    dummy = 1.d0 / A(m,11)
    A(m+1,11) = dummy * (A(m+1,11) - A(m-1,11) * A(m+2,10) - A(m-2,11) * A(m+3,9) - A(m-3,11) * A(m+4,8) - A(m-4,11) * A(m+5,7)&
         - A(m-5,11) * A(m+6,6) - A(m-6,11) * A(m-5,5))
    A(m+2,11) = dummy * (A(m+2,11) - A(m-1,11) * A(m+3,10) - A(m-2,11) * A(m+4,9) - A(m-3,11) * A(m+5,8) - A(m-4,11) * A(m+6,7)&
         - A(m-5,11) * A(m-6,6) - A(m-6,11) * A(m-6,5))
    A(m+3,11) = dummy * (A(m+3,11) - A(m-1,11) * A(m+4,10) - A(m-2,11) * A(m+5,9) - A(m-3,11) * A(m+6,8))
    A(m+4,11) = dummy * (A(m+4,11) - A(m-1,11) * A(m+5,10) - A(m-2,11) * A(m+6,9))
    A(m+5,11) = dummy * (A(m+5,11) - A(m-1,11) * A(m+6,10))
    A(m+6,11) = dummy *  A(m+6,11) 

    A(m  ,12) = A(m  ,12) - A(m-1,12) * A(m+1,11) - A(m-2,12) * A(m+2,10) - A(m-3,12) * A(m+3,9) - A(m-4,12) * A(m+4,8)&
         - A(m-5,12) * A(m+5,7) - A(m-6,12) * A(m+6,6)
    A(m-1,13) = A(m-1,13) - A(m-2,13) * A(m+1,11) - A(m-3,13) * A(m+2,10) - A(m-4,13) * A(m+3,9) - A(m-5,13) * A(m+4,8)&
         - A(m-6,13) * A(m+5,7)
    A(m-2,14) = A(m-2,14) - A(m-3,14) * A(m+1,11) - A(m-4,14) * A(m+2,10) - A(m-5,14) * A(m+3,9) - A(m-6,14) * A(m+4,8)
    A(m-3,15) = A(m-3,15) - A(m-4,15) * A(m+1,11) - A(m-5,15) * A(m+2,10) - A(m-6,15) * A(m+3,9)
    A(m-4,16) = A(m-4,16) - A(m-5,16) * A(m+1,11) - A(m-6,16) * A(m+2,10)
    A(m-5,17) = A(m-5,17) - A(m-6,17) * A(m+1,11) 


    dummy = 1.d0 / A(m,12)
    A(m+1,12) = dummy * (A(m+1,12) - A(m-1,12) * A(m+2,11) - A(m-2,12) * A(m+3,10) - A(m-3,12) * A(m+4,9) - A(m-4,12) * A(m+5,8)&
         - A(m-5,12) * A(m+6,7) - A(m-6,12) * A(m-6,6))
    A(m+2,12) = dummy * (A(m+2,12) - A(m-1,12) * A(m+3,11) - A(m-2,12) * A(m+4,10) - A(m-3,12) * A(m+5,9) - A(m-4,12) * A(m+6,8))
    A(m+3,12) = dummy * (A(m+3,12) - A(m-1,12) * A(m+4,11) - A(m-2,12) * A(m+5,10) - A(m-3,12) * A(m+6,9))
    A(m+4,12) = dummy * (A(m+4,12) - A(m-1,12) * A(m+5,11) - A(m-2,12) * A(m+6,10))
    A(m+5,12) = dummy * (A(m+5,12) - A(m-1,12) * A(m+6,11))
    A(m+6,12) = dummy *  A(m+6,12)


    !center region

    do i = 2*bw+1,n-2*bw-1
       A(m  ,i  ) = A(m  ,i  ) - A(m-1,i  ) * A(m+1,i-1) - A(m-2,i)   * A(m+2,i-2) - A(m-3,i  ) * A(m+3,i-3) &
            - A(m-4,i)   * A(m+4,i-4) - A(m-5,i)   * A(m+5,i-5) - A(m-6,i)   * A(m+6,i-6) 
       A(m-1,i+1) = A(m-1,i+1) - A(m-2,i+1) * A(m+1,i-1) - A(m-3,i+1) * A(m+2,i-2) - A(m-4,i+1) * A(m+3,i-3)&
            - A(m-5,i+1) * A(m+4,i-4) - A(m-6,i+1) * A(m+5,i-5)
       A(m-2,i+2) = A(m-2,i+2) - A(m-3,i+2) * A(m+1,i-1) - A(m-4,i+2) * A(m+2,i-2) - A(m-5,i+2) * A(m+3,i-3)&
            - A(m-6,i+2) * A(m+4,i-4)
       A(m-3,i+3) = A(m-3,i+3) - A(m-4,i+3) * A(m+1,i-1) - A(m-5,i+3) * A(m+2,i-2) - A(m-6,i+3) * A(m+3,i-3)
       A(m-4,i+4) = A(m-4,i+4) - A(m-5,i+4) * A(m+1,i-1) - A(m-6,i+4) * A(m+2,i-2)
       A(m-5,i+5) = A(m-5,i+5) - A(m-6,i+5) * A(m+1,i-1)

       dummy = 1.d0 / A(m,i)
       A(m+1,i) = dummy * (A(m+1,i) - A(m-1,i) * A(m+2,i-1) - A(m-2,i) * A(m+3,i-2) - A(m-3,i) * A(m+4,i-3) &
            - A(m-4,i) * A(m+5,i-4) - A(m-5,i) * A(m+6,i-5))
       A(m+2,i) = dummy * (A(m+2,i) - A(m-1,i) * A(m+3,i-1) - A(m-2,i) * A(m+4,i-2) - A(m-3,i) * A(m+5,i-3) &
            - A(m-4,i) * A(m+6,i-4))
       A(m+3,i) = dummy * (A(m+3,i) - A(m-1,i) * A(m+4,i-1) - A(m-2,i) * A(m+5,i-2) - A(m-3,i) * A(m+6,i-3))
       A(m+4,i) = dummy * (A(m+4,i) - A(m-1,i) * A(m+5,i-1) - A(m-2,i) * A(m+6,i-2))
       A(m+5,i) = dummy * (A(m+5,i) - A(m-1,i) * A(m+6,i-1))
       A(m+6,i) = dummy *  A(m+6,i) 
    enddo

    ! sublower region
    A(m  ,n-2*bw)    = A(m  ,n-2*bw)    - A(m-1,n-2*bw)   * A(m+1,n-2*bw-1) - A(m-2,n-2*bw)   * A(m+2,n-2*bw-2)&
         - A(m-3,n-2*bw)   * A(m+3,n-2*bw-3) - A(m-4,n-2*bw)   * A(m+4,n-2*bw-4) - A(m-5,n-2*bw)   * A(m+5,n-2*bw-5)&
         - A(m-6,n-2*bw)   * A(m+6,n-2*bw-6)
    A(m-1,n-2*bw+1)  = A(m-1,n-2*bw+1)  - A(m-2,n-2*bw+1) * A(m+1,n-2*bw-1) - A(m-3,n-2*bw+1) * A(m+2,n-2*bw-2) &
         - A(m-4,n-2*bw+1) * A(m+3,n-2*bw-3) - A(m-5,n-2*bw+1) * A(m+4,n-2*bw-4) - A(m-6,n-2*bw+1) * A(m+5,n-2*bw-5)  
    A(m-2,n-2*bw+2)  = A(m-2,n-2*bw+2)  - A(m-3,n-2*bw+2) * A(m+1,n-2*bw-1) - A(m-4,n-2*bw+2) * A(m+2,n-2*bw-2)&
         - A(m-5,n-2*bw+2) * A(m+3,n-2*bw-3) - A(m-6,n-2*bw+2) * A(m+4,n-2*bw-4)    
    A(m-3,n-2*bw+3)  = A(m-3,n-2*bw+3)  - A(m-4,n-2*bw+3) * A(m+1,n-2*bw-1) - A(m-5,n-2*bw+3) * A(m+2,n-2*bw-2)&
         - A(m-6,n-2*bw+3) * A(m+3,n-2*bw-3) 
    A(m-4,n-2*bw+4)  = A(m-4,n-2*bw+4)  - A(m-5,n-2*bw+4) * A(m+1,n-2*bw-1) - A(m-6,n-2*bw+4) * A(m+2,n-2*bw-2)
    A(m-5,n-2*bw+5)  = A(m-5,n-2*bw+5)  - A(m-6,n-2*bw+5) * A(m+1,n-2*bw-1)

    dummy = 1.d0 / A(m,n-2*bw)
    A(m+1,n-2*bw) = dummy * (A(m+1,n-2*bw) - A(m-1,n-2*bw) * A(m+2,n-2*bw-1) - A(m-2,n-2*bw) * A(m+3,n-2*bw-2)&
         - A(m-3,n-2*bw) * A(m+4,n-2*bw-3) - A(m-4,n-2*bw) * A(m+5,n-2*bw-4) - A(m-5,n-2*bw) * A(m+6,n-2*bw-5))
    A(m+2,n-2*bw) = dummy * (A(m+2,n-2*bw) - A(m-1,n-2*bw) * A(m+3,n-2*bw-1) - A(m-2,n-2*bw) * A(m+4,n-2*bw-2)&
         - A(m-3,n-2*bw) * A(m+5,n-2*bw-3) - A(m-4,n-2*bw) * A(m+6,n-2*bw-4))
    A(m+3,n-2*bw) = dummy * (A(m+3,n-2*bw) - A(m-1,n-2*bw) * A(m+4,n-2*bw-1) - A(m-2,n-2*bw) * A(m+5,n-2*bw-2)&
         - A(m-3,n-2*bw) * A(m+6,n-2*bw-3))
    A(m+4,n-2*bw) = dummy * (A(m+4,n-2*bw) - A(m-1,n-2*bw) * A(m+5,n-2*bw-1) - A(m-2,n-2*bw) * A(m+6,n-2*bw-2))
    A(m+5,n-2*bw) = dummy * (A(m+5,n-2*bw) - A(m-1,n-2*bw) * A(m+6,n-2*bw-1))
    A(m+6,n-2*bw) = dummy *  A(m+6,n-2*bw)

    A(m  ,n-2*bw+1)   = A(m  ,n-2*bw+1)  - A(m-1,n-2*bw+1)  * A(m+1,n-2*bw) - A(m-2,n-2*bw+1) * A(m+2,n-2*bw-1) &
         - A(m-3,n-2*bw+1) * A(m+3,n-2*bw-2) - A(m-4,n-2*bw+1) * A(m+4,n-2*bw-3) - A(m-5,n-2*bw+1) * A(m+5,n-2*bw-4)&
         - A(m-6,n-2*bw+1) * A(m+6,n-2*bw-5)
    A(m-1,n-2*bw+2)   = A(m-1,n-2*bw+2)  - A(m-2,n-2*bw+2)  * A(m+1,n-2*bw) - A(m-3,n-2*bw+2) * A(m+2,n-2*bw-1)&
         - A(m-4,n-2*bw+2) * A(m+3,n-2*bw-2) - A(m-5,n-2*bw+2) * A(m+4,n-2*bw-3) - A(m-6,n-2*bw+2) * A(m+5,n-2*bw-4)
    A(m-2,n-2*bw+3)   = A(m-2,n-2*bw+3)  - A(m-3,n-2*bw+3)  * A(m+1,n-2*bw) - A(m-4,n-2*bw+3) * A(m+2,n-2*bw-1)&
         - A(m-5,n-2*bw+3) * A(m+3,n-2*bw-2) - A(m-6,n-2*bw+3) * A(m+4,n-2*bw-3) 
    A(m-3,n-2*bw+4)   = A(m-3,n-2*bw+4)  - A(m-4,n-2*bw+4)  * A(m+1,n-2*bw) - A(m-5,n-2*bw+4) * A(m+2,n-2*bw-1)&
         - A(m-6,n-2*bw+4) * A(m+3,n-2*bw-2)
    A(m-4,n-2*bw+5)   = A(m-4,n-2*bw+5)  - A(m-5,n-2*bw+5)  * A(m+1,n-2*bw) - A(m-6,n-2*bw+5) * A(m+2,n-2*bw-1)
    A(m-5,n-2*bw+6)   = A(m-5,n-2*bw+6)  - A(m-6,n-2*bw+6)  * A(m+1,n-2*bw)
    A(m-6,n-2*bw+7)   = A(m-6,n-2*bw+7)  - A(m+6,n-2*bw+7)  * A(m+1,n-2*bw) 
    A(m+5,n-2*bw+8)   = A(m+5,n-2*bw+8)  - A(m+6,n-2*bw+8)  * A(m+1,n-2*bw) 
    A(m+5,n-2*bw+9)   = A(m+5,n-2*bw+9)  - A(m+6,n-2*bw+9)  * A(m+1,n-2*bw) 
    A(m+5,n-2*bw+10)  = A(m+5,n-2*bw+10) - A(m+6,n-2*bw+10) * A(m+1,n-2*bw) 
    A(m+5,n-2*bw+11)  = A(m+5,n-2*bw+11) - A(m+6,n-2*bw+11) * A(m+1,n-2*bw) 
    A(m+5,n-2*bw+12)  = A(m+5,n-2*bw+12) - A(m+6,n-2*bw+12) * A(m+1,n-2*bw) 

    dummy = 1.d0 / A(m,n-2*bw+1)
    A(m+1,n-2*bw+1) = dummy * (A(m+1,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+2,n-2*bw) - A(m-2,n-2*bw+1) * A(m+3,n-2*bw-1)&
         - A(m-3,n-2*bw+1) * A(m+4,n-2*bw-2) - A(m-4,n-2*bw+1) * A(m+5,n-2*bw-3) - A(m-5,n-2*bw+1) * A(m+6,n-2*bw-4))
    A(m+2,n-2*bw+1) = dummy * (A(m+2,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+3,n-2*bw) - A(m-2,n-2*bw+1) * A(m+4,n-2*bw-1)&
         - A(m-3,n-2*bw+1) * A(m+5,n-2*bw-2) - A(m-4,n-2*bw+1) * A(m+6,n-2*bw-3))
    A(m+3,n-2*bw+1) = dummy * (A(m+3,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+4,n-2*bw) - A(m-2,n-2*bw+1) * A(m+5,n-2*bw-1)&
         - A(m-3,n-2*bw+1) * A(m+6,n-2*bw-2))
    A(m+4,n-2*bw+1) = dummy * (A(m+4,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+5,n-2*bw) - A(m-2,n-2*bw+1) * A(m+6,n-2*bw-1))
    A(m+5,n-2*bw+1) = dummy * (A(m+5,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+6,n-2*bw))
    A(m+6,n-2*bw+1) = dummy *  A(m+6,n-2*bw+1)

    A(m  ,n-2*bw+2)   = A(m  ,n-2*bw+2)  - A(m-1,n-2*bw+2) * A(m+1,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+2,n-2*bw) &
         - A(m-3,n-2*bw+2) * A(m+3,n-2*bw-1) - A(m-4,n-2*bw+2) * A(m+4,n-2*bw-2) - A(m-5,n-2*bw+2) * A(m+5,n-2*bw-3) &
         - A(m-6,n-2*bw+2) * A(m+6,n-2*bw-4)
    A(m-1,n-2*bw+3)   = A(m-1,n-2*bw+3)  - A(m-2,n-2*bw+3) * A(m+1,n-2*bw+1) - A(m-3,n-2*bw+3) * A(m+2,n-2*bw)&
         - A(m-4,n-2*bw+3) * A(m+3,n-2*bw-1) - A(m-5,n-2*bw+3) * A(m+4,n-2*bw-2) - A(m-6,n-2*bw+3) * A(m+5,n-2*bw-3)
    A(m-2,n-2*bw+4)   = A(m-2,n-2*bw+4)  - A(m-3,n-2*bw+4) * A(m+1,n-2*bw+1) - A(m-4,n-2*bw+4) * A(m+2,n-2*bw) &
         - A(m-5,n-2*bw+4) * A(m+3,n-2*bw-1) - A(m-6,n-2*bw+4) * A(m+4,n-2*bw-2)
    A(m-3,n-2*bw+5)   = A(m-3,n-2*bw+5)  - A(m-4,n-2*bw+5) * A(m+1,n-2*bw+1) - A(m-5,n-2*bw+5) * A(m+2,n-2*bw)&
         - A(m-6,n-2*bw+5) * A(m+3,n-2*bw-1)
    A(m-4,n-2*bw+6)   = A(m-4,n-2*bw+6)  - A(m-5,n-2*bw+6) * A(m+1,n-2*bw+1) - A(m-6,n-2*bw+6) * A(m+2,n-2*bw)
    A(m-5,n-2*bw+7)   = A(m-5,n-2*bw+7)  - A(m-6,n-2*bw+7) * A(m+1,n-2*bw+1) - A(m+6,n-2*bw+7) * A(m+2,n-2*bw)
    A(m-6,n-2*bw+8)   = A(m-6,n-2*bw+8)  - A(m+5,n-2*bw+8) * A(m+1,n-2*bw+1) - A(m+6,n-2*bw+8) * A(m+2,n-2*bw)
    A(m+4,n-2*bw+9)   = A(m+4,n-2*bw+9)  - A(m+5,n-2*bw+9) * A(m+1,n-2*bw+1) - A(m+6,n-2*bw+9) * A(m+2,n-2*bw)
    A(m+4,n-2*bw+10)  = A(m+4,n-2*bw+10) - A(m+5,n-2*bw+10) * A(m+1,n-2*bw+1) - A(m+6,n-2*bw+10) * A(m+2,n-2*bw)
    A(m+4,n-2*bw+11)  = A(m+4,n-2*bw+11) - A(m+5,n-2*bw+11) * A(m+1,n-2*bw+1) - A(m+6,n-2*bw+11) * A(m+2,n-2*bw)
    A(m+4,n-2*bw+12)  = A(m+4,n-2*bw+12) - A(m+5,n-2*bw+12) * A(m+1,n-2*bw+1) - A(m+6,n-2*bw+12) * A(m+2,n-2*bw)

    dummy = 1.d0 / A(m,n-2*bw+2)
    A(m+1,n-2*bw+2) = dummy * (A(m+1,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+2,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+3,n-2*bw)&
         - A(m-3,n-2*bw+2) * A(m+4,n-2*bw-1) - A(m-4,n-2*bw+2) * A(m+5,n-2*bw-2) - A(m-5,n-2*bw+2) * A(m+6,n-2*bw-3))
    A(m+2,n-2*bw+2) = dummy * (A(m+2,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+3,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+4,n-2*bw)&
         - A(m-3,n-2*bw+2) * A(m+5,n-2*bw-1) - A(m-4,n-2*bw+2) * A(m+6,n-2*bw-2))
    A(m+3,n-2*bw+2) = dummy * (A(m+3,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+4,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+5,n-2*bw) &
         - A(m-3,n-2*bw+2) * A(m+6,n-2*bw-1))
    A(m+4,n-2*bw+2) = dummy * (A(m+4,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+5,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+6,n-2*bw))
    A(m+5,n-2*bw+2) = dummy * (A(m+5,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+6,n-2*bw+1))
    A(m+6,n-2*bw+2) = dummy * (A(m+6,n-2*bw+2))

    A(m  ,n-2*bw+3)   = A(m  ,n-2*bw+3)  - A(m-1,n-2*bw+3) * A(m+1,n-2*bw+2) - A(m-2,n-2*bw+3) * A(m+2,n-2*bw+1)&
         - A(m-3,n-2*bw+3) * A(m+3,n-2*bw) - A(m-4,n-2*bw+3) * A(m+4,n-2*bw-1) - A(m-5,n-2*bw+3) * A(m+5,n-2*bw-2)&
         - A(m-6,n-2*bw+3) * A(m+6,n-2*bw-3)
    A(m-1,n-2*bw+4)   = A(m-1,n-2*bw+4)  - A(m-2,n-2*bw+4) * A(m+1,n-2*bw+2) - A(m-3,n-2*bw+4) * A(m+2,n-2*bw+1)&
         - A(m-4,n-2*bw+4) * A(m+3,n-2*bw) - A(m-5,n-2*bw+4) * A(m+4,n-2*bw-1) - A(m-6,n-2*bw+4) * A(m+5,n-2*bw-2)
    A(m-2,n-2*bw+5)   = A(m-2,n-2*bw+5)  - A(m-3,n-2*bw+5) * A(m+1,n-2*bw+2) - A(m-4,n-2*bw+5) * A(m+2,n-2*bw+1) &
         - A(m-5,n-2*bw+5) * A(m+3,n-2*bw) - A(m-6,n-2*bw+5) * A(m+4,n-2*bw-1)
    A(m-3,n-2*bw+6)   = A(m-3,n-2*bw+6)  - A(m-4,n-2*bw+6) * A(m+1,n-2*bw+2) - A(m-5,n-2*bw+6) * A(m+2,n-2*bw+1) &
         - A(m-6,n-2*bw+6) * A(m+3,n-2*bw)
    A(m-4,n-2*bw+7)   = A(m-4,n-2*bw+7)  - A(m-5,n-2*bw+7) * A(m+1,n-2*bw+2) - A(m-6,n-2*bw+7) * A(m+2,n-2*bw+1)&
         - A(m+6,n-2*bw+7) * A(m+3,n-2*bw)
    A(m-5,n-2*bw+8)   = A(m-5,n-2*bw+8)  - A(m-6,n-2*bw+8)  * A(m+1,n-2*bw+2) - A(m+5,n-2*bw+8)  * A(m+2,n-2*bw+1)&
         - A(m+6,n-2*bw+8) * A(m+3,n-2*bw)
    A(m-6,n-2*bw+9)   = A(m-6,n-2*bw+9)  - A(m+4,n-2*bw+9)  * A(m+1,n-2*bw+2) - A(m+5,n-2*bw+9)  * A(m+2,n-2*bw+1)&
         - A(m+6,n-2*bw+9) * A(m+3,n-2*bw)
    A(m+3,n-2*bw+10)  = A(m+3,n-2*bw+10) - A(m+4,n-2*bw+10) * A(m+1,n-2*bw+2) - A(m+5,n-2*bw+10) * A(m+2,n-2*bw+1)&
         - A(m+6,n-2*bw+10) * A(m+3,n-2*bw)
    A(m+3,n-2*bw+11)  = A(m+3,n-2*bw+11) - A(m+4,n-2*bw+11) * A(m+1,n-2*bw+2) - A(m+5,n-2*bw+11) * A(m+2,n-2*bw+1)&
         - A(m+6,n-2*bw+11) * A(m+3,n-2*bw)
    A(m+3,n-2*bw+12)  = A(m+3,n-2*bw+12) - A(m+4,n-2*bw+12) * A(m+1,n-2*bw+2) - A(m+5,n-2*bw+12) * A(m+2,n-2*bw+1)&
         - A(m+6,n-2*bw+12) * A(m+3,n-2*bw)


    dummy = 1.d0 / A(m,n-2*bw+3)
    A(m+1,n-2*bw+3) = dummy * (A(m+1,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+2,n-2*bw+2) - A(m-2,n-2*bw+3) * A(m+3,n-2*bw+1)&
         - A(m-3,n-2*bw+3) * A(m+4,n-2*bw) - A(m-4,n-2*bw+3) * A(m+5,n-2*bw-1) - A(m-5,n-2*bw+3) * A(m+6,n-2*bw-2))
    A(m+2,n-2*bw+3) = dummy * (A(m+2,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+3,n-2*bw+2) - A(m-2,n-2*bw+3) * A(m+4,n-2*bw+1)&
         - A(m-3,n-2*bw+3) * A(m+5,n-2*bw) - A(m-4,n-2*bw+3) * A(m+6,n-2*bw-1))
    A(m+3,n-2*bw+3) = dummy * (A(m+3,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+4,n-2*bw+2) - A(m-2,n-2*bw+3) * A(m+5,n-2*bw+1)&
         - A(m-3,n-2*bw+3) * A(m+6,n-2*bw))
    A(m+4,n-2*bw+3) = dummy * (A(m+4,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+5,n-2*bw+2) - A(m-2,n-2*bw+3) * A(m+6,n-2*bw+1))
    A(m+5,n-2*bw+3) = dummy * (A(m+5,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+6,n-2*bw+2))
    A(m+6,n-2*bw+3) = dummy * (A(m+6,n-2*bw+3))

    A(m  ,n-2*bw+4)   = A(m  ,n-2*bw+4)  - A(m-1,n-2*bw+4) * A(m+1,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+2,n-2*bw+2)&
         - A(m-3,n-2*bw+4) * A(m+3,n-2*bw+1) - A(m-4,n-2*bw+4) * A(m+4,n-2*bw) - A(m-5,n-2*bw+4) * A(m+5,n-2*bw-1)&
         - A(m-6,n-2*bw+4) * A(m+6,n-2*bw-2)
    A(m-1,n-2*bw+5)   = A(m-1,n-2*bw+5)  - A(m-2,n-2*bw+5) * A(m+1,n-2*bw+3) - A(m-3,n-2*bw+5) * A(m+2,n-2*bw+2) &
         - A(m-4,n-2*bw+5) * A(m+3,n-2*bw+1) - A(m-5,n-2*bw+5) * A(m+4,n-2*bw) - A(m-6,n-2*bw+5) * A(m+5,n-2*bw-1)
    A(m-2,n-2*bw+6)   = A(m-2,n-2*bw+6)  - A(m-3,n-2*bw+6) * A(m+1,n-2*bw+3) - A(m-4,n-2*bw+6) * A(m+2,n-2*bw+2)&
         - A(m-5,n-2*bw+6) * A(m+3,n-2*bw+1) - A(m-6,n-2*bw+6) * A(m+4,n-2*bw)
    A(m-3,n-2*bw+7)   = A(m-3,n-2*bw+7)  - A(m-4,n-2*bw+7) * A(m+1,n-2*bw+3) - A(m-5,n-2*bw+7) * A(m+2,n-2*bw+2)&
         - A(m-6,n-2*bw+7) * A(m+3,n-2*bw+1) - A(m+6,n-2*bw+7) * A(m+4,n-2*bw) 
    A(m-4,n-2*bw+8)   = A(m-4,n-2*bw+8)  - A(m-5,n-2*bw+8) * A(m+1,n-2*bw+3) - A(m-6,n-2*bw+8) * A(m+2,n-2*bw+2)&
         - A(m+5,n-2*bw+8) * A(m+3,n-2*bw+1) - A(m+6,n-2*bw+8) * A(m+4,n-2*bw)
    A(m-5,n-2*bw+9)    = A(m-5,n-2*bw+9)   - A(m-6,n-2*bw+9)  * A(m+1,n-2*bw+3) - A(m+4,n-2*bw+9)  * A(m+2,n-2*bw+2)&
         - A(m+5,n-2*bw+9)  * A(m+3,n-2*bw+1) - A(m+6,n-2*bw+9) * A(m+4,n-2*bw)
    A(m-6,n-2*bw+10)   = A(m-6,n-2*bw+10)  - A(m+3,n-2*bw+10) * A(m+1,n-2*bw+3) - A(m+4,n-2*bw+10) * A(m+2,n-2*bw+2)&
         - A(m+5,n-2*bw+10) * A(m+3,n-2*bw+1)- A(m+6,n-2*bw+10) * A(m+4,n-2*bw)
    A(m+2,n-2*bw+11)   = A(m+2,n-2*bw+11)  - A(m+3,n-2*bw+11) * A(m+1,n-2*bw+3) - A(m+4,n-2*bw+11) * A(m+2,n-2*bw+2)&
         - A(m+5,n-2*bw+11) * A(m+3,n-2*bw+1)- A(m+6,n-2*bw+11) * A(m+4,n-2*bw)
    A(m+2,n-2*bw+12)   = A(m+2,n-2*bw+12)  - A(m+3,n-2*bw+12) * A(m+1,n-2*bw+3) - A(m+4,n-2*bw+12) * A(m+2,n-2*bw+2)&
         - A(m+5,n-2*bw+12) * A(m+3,n-2*bw+1)- A(m+6,n-2*bw+12) * A(m+4,n-2*bw)


    dummy = 1.d0 / A(m,n-2*bw+4)
    A(m+1,n-2*bw+4) = dummy * (A(m+1,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+2,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+3,n-2*bw+2) &
         - A(m-3,n-2*bw+4) * A(m+4,n-2*bw+1) - A(m-4,n-2*bw+4) * A(m+5,n-2*bw) - A(m-5,n-2*bw+4) * A(m+6,n-2*bw-1))
    A(m+2,n-2*bw+4) = dummy * (A(m+2,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+3,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+4,n-2*bw+2) &
         - A(m-3,n-2*bw+4) * A(m+5,n-2*bw+1) - A(m-4,n-2*bw+4) * A(m+6,n-2*bw))
    A(m+3,n-2*bw+4) = dummy * (A(m+3,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+4,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+5,n-2*bw+2) &
         - A(m-3,n-2*bw+4) * A(m+6,n-2*bw+1))
    A(m+4,n-2*bw+4) = dummy * (A(m+4,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+5,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+6,n-2*bw+2))
    A(m+5,n-2*bw+4) = dummy * (A(m+5,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+6,n-2*bw+3))
    A(m+6,n-2*bw+4) = dummy *  A(m+6,n-2*bw+4)

    A(m  ,n-2*bw+5)   = A(m  ,n-2*bw+5)  - A(m-1,n-2*bw+5)  * A(m+1,n-2*bw+4) - A(m-2,n-2*bw+5)  * A(m+2,n-2*bw+3) &
         - A(m-3,n-2*bw+5)  * A(m+3,n-2*bw+2) - A(m-4,n-2*bw+5)  * A(m+4,n-2*bw+1) - A(m-5,n-2*bw+5)  * A(m+5,n-2*bw) &
         - A(m-6,n-2*bw+5) * A(m+6,n-2*bw-1)     
    A(m-1,n-2*bw+6)   = A(m-1,n-2*bw+6)  - A(m-2,n-2*bw+6)  * A(m+1,n-2*bw+4) - A(m-3,n-2*bw+6)  * A(m+2,n-2*bw+3) &
         - A(m-4,n-2*bw+6)  * A(m+3,n-2*bw+2) - A(m-5,n-2*bw+6)  * A(m+4,n-2*bw+1) - A(m-6,n-2*bw+6)  * A(m+5,n-2*bw)   
    A(m-2,n-2*bw+7)   = A(m-2,n-2*bw+7)  - A(m-3,n-2*bw+7)  * A(m+1,n-2*bw+4) - A(m-4,n-2*bw+7)  * A(m+2,n-2*bw+3)&
         - A(m-5,n-2*bw+7)  * A(m+3,n-2*bw+2) - A(m-6,n-2*bw+7)  * A(m+4,n-2*bw+1) - A(m+6,n-2*bw+7)  * A(m+5,n-2*bw)    
    A(m-3,n-2*bw+8)   = A(m-3,n-2*bw+8)  - A(m-4,n-2*bw+8)  * A(m+1,n-2*bw+4) - A(m-5,n-2*bw+8)  * A(m+2,n-2*bw+3) &
         - A(m-6,n-2*bw+8)  * A(m+3,n-2*bw+2) - A(m+5,n-2*bw+8)  * A(m+4,n-2*bw+1) - A(m+6,n-2*bw+8)  * A(m+5,n-2*bw)    
    A(m-4,n-2*bw+9)   = A(m-4,n-2*bw+9)  - A(m-5,n-2*bw+9)  * A(m+1,n-2*bw+4) - A(m-6,n-2*bw+9)  * A(m+2,n-2*bw+3) &
         - A(m+4,n-2*bw+9)  * A(m+3,n-2*bw+2) - A(m+5,n-2*bw+9)  * A(m+4,n-2*bw+1) - A(m+6,n-2*bw+9)  * A(m+5,n-2*bw)    
    A(m-5,n-2*bw+10)  = A(m-5,n-2*bw+10) - A(m-6,n-2*bw+10) * A(m+1,n-2*bw+4) - A(m+3,n-2*bw+10) * A(m+2,n-2*bw+3) &
         - A(m+4,n-2*bw+10) * A(m+3,n-2*bw+2) - A(m+5,n-2*bw+10) * A(m+4,n-2*bw+1) - A(m+6,n-2*bw+10) * A(m+5,n-2*bw)    
    A(m-6,n-2*bw+11)  = A(m-6,n-2*bw+11) - A(m+2,n-2*bw+11) * A(m+1,n-2*bw+4) - A(m+3,n-2*bw+11) * A(m+2,n-2*bw+3) &
         - A(m+4,n-2*bw+11) * A(m+3,n-2*bw+2) - A(m+5,n-2*bw+11) * A(m+4,n-2*bw+1) - A(m+6,n-2*bw+11) * A(m+5,n-2*bw)    
    A(m+1,n-2*bw+12)  = A(m+1,n-2*bw+12) - A(m+2,n-2*bw+12) * A(m+1,n-2*bw+4) - A(m+3,n-2*bw+12) * A(m+2,n-2*bw+3)&
         - A(m+4,n-2*bw+12) * A(m+3,n-2*bw+2) - A(m+5,n-2*bw+12) * A(m+4,n-2*bw+1) - A(m+6,n-2*bw+12) * A(m+5,n-2*bw)    

    dummy = 1.d0 / A(m,n-2*bw+5)
    A(m+1,n-2*bw+5) = dummy * (A(m+1,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+2,n-2*bw+4) - A(m-2,n-2*bw+5) * A(m+3,n-2*bw+3)&
         - A(m-3,n-2*bw+5) * A(m+4,n-2*bw+2) - A(m-4,n-2*bw+5) * A(m+5,n-2*bw+1) - A(m-5,n-2*bw+5) * A(m+6,n-2*bw))
    A(m+2,n-2*bw+5) = dummy * (A(m+2,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+3,n-2*bw+4) - A(m-2,n-2*bw+5) * A(m+4,n-2*bw+3)&
         - A(m-3,n-2*bw+5) * A(m+5,n-2*bw+2) - A(m-4,n-2*bw+5) * A(m+6,n-2*bw+1))
    A(m+3,n-2*bw+5) = dummy * (A(m+3,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+4,n-2*bw+4) - A(m-2,n-2*bw+5) * A(m+5,n-2*bw+3)&
         - A(m-3,n-2*bw+5) * A(m+6,n-2*bw+2))
    A(m+4,n-2*bw+5) = dummy * (A(m+4,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+5,n-2*bw+4) - A(m-2,n-2*bw+5) * A(m+6,n-2*bw+3))
    A(m+5,n-2*bw+5) = dummy * (A(m+5,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+6,n-2*bw+4))
    A(m+6,n-2*bw+5) = dummy *  A(m+6,n-2*bw+5)

    A(m  ,n-2*bw+6)   = A(m  ,n-2*bw+6)  - A(m-1,n-2*bw+6)  * A(m+1,n-2*bw+5) - A(m-2,n-2*bw+6)  * A(m+2,n-2*bw+4)&
         - A(m-3,n-2*bw+6)  * A(m+3,n-2*bw+3) - A(m-4,n-2*bw+6)  * A(m+4,n-2*bw+2) - A(m-5,n-2*bw+6)  * A(m+5,n-2*bw+1)&
         - A(m-6,n-2*bw+6)  * A(m+6,n-2*bw)     
    A(m-1,n-2*bw+7)   = A(m-1,n-2*bw+7)  - A(m-2,n-2*bw+7)  * A(m+1,n-2*bw+5) - A(m-3,n-2*bw+7)  * A(m+2,n-2*bw+4)&
         - A(m-4,n-2*bw+7)  * A(m+3,n-2*bw+3) - A(m-5,n-2*bw+7)  * A(m+4,n-2*bw+2) - A(m-6,n-2*bw+7)  * A(m+5,n-2*bw+1)&
         - A(m+6,n-2*bw+7)  * A(m+6,n-2*bw)  
    A(m-2,n-2*bw+8)   = A(m-2,n-2*bw+8)  - A(m-3,n-2*bw+8)  * A(m+1,n-2*bw+5) - A(m-4,n-2*bw+8)  * A(m+2,n-2*bw+4)&
         - A(m-5,n-2*bw+8)  * A(m+3,n-2*bw+3) - A(m-6,n-2*bw+8)  * A(m+4,n-2*bw+2) - A(m+5,n-2*bw+8)  * A(m+5,n-2*bw+1)&
         - A(m+6,n-2*bw+8)  * A(m+6,n-2*bw)    
    A(m-3,n-2*bw+9)   = A(m-3,n-2*bw+9)  - A(m-4,n-2*bw+9)  * A(m+1,n-2*bw+5) - A(m-5,n-2*bw+9)  * A(m+2,n-2*bw+4)&
         - A(m-6,n-2*bw+9)  * A(m+3,n-2*bw+3) - A(m+4,n-2*bw+9)  * A(m+4,n-2*bw+2) - A(m+5,n-2*bw+9)  * A(m+5,n-2*bw+1)&
         - A(m+6,n-2*bw+9)  * A(m+6,n-2*bw)    
    A(m-4,n-2*bw+10)  = A(m-4,n-2*bw+10) - A(m-5,n-2*bw+10) * A(m+1,n-2*bw+5) - A(m-6,n-2*bw+10) * A(m+2,n-2*bw+4)&
         - A(m+3,n-2*bw+10) * A(m+3,n-2*bw+3) - A(m+4,n-2*bw+10) * A(m+4,n-2*bw+2) - A(m+5,n-2*bw+10) * A(m+5,n-2*bw+1)&
         - A(m+6,n-2*bw+10) * A(m+6,n-2*bw)    
    A(m-5,n-2*bw+11)  = A(m-5,n-2*bw+11) - A(m-6,n-2*bw+11) * A(m+1,n-2*bw+5) - A(m+2,n-2*bw+11) * A(m+2,n-2*bw+4)&
         - A(m+3,n-2*bw+11) * A(m+3,n-2*bw+3) - A(m+4,n-2*bw+11) * A(m+4,n-2*bw+2) - A(m+5,n-2*bw+11) * A(m+5,n-2*bw+1) &
         - A(m+6,n-2*bw+11) * A(m+6,n-2*bw)    
    A(m-6,n-2*bw+12)  = A(m-6,n-2*bw+12) - A(m+1,n-2*bw+12) * A(m+1,n-2*bw+5) - A(m+2,n-2*bw+12) * A(m+2,n-2*bw+4) &
         - A(m+3,n-2*bw+12) * A(m+3,n-2*bw+3) - A(m+4,n-2*bw+12) * A(m+4,n-2*bw+2) - A(m+5,n-2*bw+12) * A(m+5,n-2*bw+1)&
         - A(m+6,n-2*bw+12) * A(m+6,n-2*bw)    

    dummy = 1.d0 / A(m,n-2*bw+6)
    A(m+1,n-2*bw+6) = dummy * (A(m+1,n-2*bw+6) - A(m-1,n-2*bw+6) * A(m+2,n-2*bw+5) - A(m-2,n-2*bw+6) * A(m+3,n-2*bw+4) &
         - A(m-3,n-2*bw+6) * A(m+4,n-2*bw+3) - A(m-4,n-2*bw+6) * A(m+5,n-2*bw+2) - A(m-5,n-2*bw+6) * A(m+6,n-2*bw+1))
    A(m+2,n-2*bw+6) = dummy * (A(m+2,n-2*bw+6) - A(m-1,n-2*bw+6) * A(m+3,n-2*bw+5) - A(m-2,n-2*bw+6) * A(m+4,n-2*bw+4) &
         - A(m-3,n-2*bw+6) * A(m+5,n-2*bw+3) - A(m-4,n-2*bw+6) * A(m+6,n-2*bw+2))
    A(m+3,n-2*bw+6) = dummy * (A(m+3,n-2*bw+6) - A(m-1,n-2*bw+6) * A(m+4,n-2*bw+5) - A(m-2,n-2*bw+6) * A(m+5,n-2*bw+4) &
         - A(m-3,n-2*bw+6) * A(m+6,n-2*bw+3))
    A(m+4,n-2*bw+6) = dummy * (A(m+4,n-2*bw+6) - A(m-1,n-2*bw+6) * A(m+5,n-2*bw+5) - A(m-2,n-2*bw+6) * A(m+6,n-2*bw+4))
    A(m+5,n-2*bw+6) = dummy * (A(m+5,n-2*bw+6) - A(m-1,n-2*bw+6) * A(m+6,n-2*bw+5))
    A(m+6,n-2*bw+6) = dummy *  A(m+6,n-2*bw+6)


    A(m  ,n-2*bw+7)   = A(m  ,n-2*bw+7)  - A(m-1,n-2*bw+7)  * A(m+1,n-2*bw+6) - A(m-2,n-2*bw+7)  * A(m+2,n-2*bw+5)&
         - A(m-3,n-2*bw+7)  * A(m+3,n-2*bw+4) - A(m-4,n-2*bw+7)  * A(m+4,n-2*bw+3) - A(m-5,n-2*bw+7)  * A(m+5,n-2*bw+2)&
         - A(m-6,n-2*bw+7)  * A(m+6,n-2*bw+1)     
    A(m-1,n-2*bw+8)   = A(m-1,n-2*bw+8)  - A(m-2,n-2*bw+8)  * A(m+1,n-2*bw+6) - A(m-3,n-2*bw+8)  * A(m+2,n-2*bw+5)&
         - A(m-4,n-2*bw+8)  * A(m+3,n-2*bw+4) - A(m-5,n-2*bw+8)  * A(m+4,n-2*bw+3) - A(m-6,n-2*bw+8)  * A(m+5,n-2*bw+2)&
         - A(m+5,n-2*bw+8)  * A(m+6,n-2*bw+1) 
    A(m-2,n-2*bw+9)   = A(m-2,n-2*bw+9)  - A(m-3,n-2*bw+9)  * A(m+1,n-2*bw+6) - A(m-4,n-2*bw+9)  * A(m+2,n-2*bw+5)&
         - A(m-5,n-2*bw+9)  * A(m+3,n-2*bw+4) - A(m-6,n-2*bw+9)  * A(m+4,n-2*bw+3) - A(m+4,n-2*bw+9)  * A(m+5,n-2*bw+2)&
         - A(m+5,n-2*bw+9)  * A(m+6,n-2*bw+1) 
    A(m-3,n-2*bw+10)  = A(m-3,n-2*bw+10) - A(m-4,n-2*bw+10) * A(m+1,n-2*bw+6) - A(m-5,n-2*bw+10) * A(m+2,n-2*bw+5)&
         - A(m-6,n-2*bw+10) * A(m+3,n-2*bw+4) - A(m+3,n-2*bw+10) * A(m+4,n-2*bw+3) - A(m+4,n-2*bw+10) * A(m+5,n-2*bw+2)&
         - A(m+5,n-2*bw+10) * A(m+6,n-2*bw+1) 
    A(m-4,n-2*bw+11)  = A(m-4,n-2*bw+11) - A(m-5,n-2*bw+11) * A(m+1,n-2*bw+6) - A(m-6,n-2*bw+11) * A(m+2,n-2*bw+5)&
         - A(m+2,n-2*bw+11) * A(m+3,n-2*bw+4) - A(m+3,n-2*bw+11) * A(m+4,n-2*bw+3) - A(m+4,n-2*bw+11) * A(m+5,n-2*bw+2) &
         - A(m+5,n-2*bw+11) * A(m+6,n-2*bw+1) 
    A(m-5,n-2*bw+12)  = A(m-5,n-2*bw+12) - A(m-6,n-2*bw+12) * A(m+1,n-2*bw+6) - A(m+1,n-2*bw+12) * A(m+2,n-2*bw+5)&
         - A(m+2,n-2*bw+12) * A(m+3,n-2*bw+4) - A(m+3,n-2*bw+12) * A(m+4,n-2*bw+3) - A(m+4,n-2*bw+12) * A(m+5,n-2*bw+2)&
         - A(m+5,n-2*bw+12) * A(m+6,n-2*bw+1) 


    dummy = 1.d0 / A(m,n-2*bw+7)
    A(m+1,n-2*bw+7) = dummy * (A(m+1,n-2*bw+7) - A(m-1,n-2*bw+7) * A(m+2,n-2*bw+6) - A(m-2,n-2*bw+7) * A(m+3,n-2*bw+5) &
         - A(m-3,n-2*bw+7) * A(m+4,n-2*bw+4) - A(m-4,n-2*bw+7) * A(m+5,n-2*bw+3) - A(m-5,n-2*bw+7) * A(m+6,n-2*bw+2))
    A(m+2,n-2*bw+7) = dummy * (A(m+2,n-2*bw+7) - A(m-1,n-2*bw+7) * A(m+3,n-2*bw+6) - A(m-2,n-2*bw+7) * A(m+4,n-2*bw+5)&
         - A(m-3,n-2*bw+7) * A(m+5,n-2*bw+4) - A(m-4,n-2*bw+7) * A(m+6,n-2*bw+3))
    A(m+3,n-2*bw+7) = dummy * (A(m+3,n-2*bw+7) - A(m-1,n-2*bw+7) * A(m+4,n-2*bw+6) - A(m-2,n-2*bw+7) * A(m+5,n-2*bw+5)&
         - A(m-3,n-2*bw+7) * A(m+6,n-2*bw+4))
    A(m+4,n-2*bw+7) = dummy * (A(m+4,n-2*bw+7) - A(m-1,n-2*bw+7) * A(m+5,n-2*bw+6) - A(m-2,n-2*bw+7) * A(m+6,n-2*bw+5))
    A(m+5,n-2*bw+7) = dummy * (A(m+5,n-2*bw+7) - A(m-1,n-2*bw+7) * A(m+6,n-2*bw+6))

    A(m  ,n-2*bw+8)   = A(m  ,n-2*bw+8)  - A(m-1,n-2*bw+8)  * A(m+1,n-2*bw+7) - A(m-2,n-2*bw+8)  * A(m+2,n-2*bw+6)&
         - A(m-3,n-2*bw+8)  * A(m+3,n-2*bw+5) - A(m-4,n-2*bw+8)  * A(m+4,n-2*bw+4) - A(m-5,n-2*bw+8)  * A(m+5,n-2*bw+3)&
         - A(m-6,n-2*bw+8)  * A(m+6,n-2*bw+2) 
    A(m-1,n-2*bw+9)   = A(m-1,n-2*bw+9)  - A(m-2,n-2*bw+9)  * A(m+1,n-2*bw+7) - A(m-3,n-2*bw+9)  * A(m+2,n-2*bw+6)&
         - A(m-4,n-2*bw+9)  * A(m+3,n-2*bw+5) - A(m-5,n-2*bw+9)  * A(m+4,n-2*bw+4) - A(m-6,n-2*bw+9)  * A(m+5,n-2*bw+3)&
         - A(m+4,n-2*bw+9)  * A(m+6,n-2*bw+2) 
    A(m-2,n-2*bw+10)  = A(m-2,n-2*bw+10) - A(m-3,n-2*bw+10) * A(m+1,n-2*bw+7) - A(m-4,n-2*bw+10) * A(m+2,n-2*bw+6)&
         - A(m-5,n-2*bw+10) * A(m+3,n-2*bw+5) - A(m-6,n-2*bw+10) * A(m+4,n-2*bw+4) - A(m+3,n-2*bw+10) * A(m+5,n-2*bw+3)&
         - A(m+4,n-2*bw+10) * A(m+6,n-2*bw+2)
    A(m-3,n-2*bw+11)  = A(m-3,n-2*bw+11) - A(m-4,n-2*bw+11) * A(m+1,n-2*bw+7) - A(m-5,n-2*bw+11) * A(m+2,n-2*bw+6)&
         - A(m-6,n-2*bw+11) * A(m+3,n-2*bw+5) - A(m+2,n-2*bw+11) * A(m+4,n-2*bw+4) - A(m+3,n-2*bw+11) * A(m+5,n-2*bw+3)&
         - A(m+4,n-2*bw+11) * A(m+6,n-2*bw+2)
    A(m-4,n-2*bw+12)  = A(m-4,n-2*bw+12) - A(m-5,n-2*bw+12) * A(m+1,n-2*bw+7) - A(m-6,n-2*bw+12) * A(m+2,n-2*bw+6)&
         - A(m+1,n-2*bw+12) * A(m+3,n-2*bw+5) - A(m+2,n-2*bw+12) * A(m+4,n-2*bw+4) - A(m+3,n-2*bw+12) * A(m+5,n-2*bw+3)&
         - A(m+4,n-2*bw+12) * A(m+6,n-2*bw+2)

    dummy = 1.d0 / A(m,n-2*bw+8)
    A(m+1,n-2*bw+8) = dummy * (A(m+1,n-2*bw+8) - A(m-1,n-2*bw+8) * A(m+2,n-2*bw+7) - A(m-2,n-2*bw+8) * A(m+3,n-2*bw+6)&
         - A(m-3,n-2*bw+8) * A(m+4,n-2*bw+5) - A(m-4,n-2*bw+8) * A(m+5,n-2*bw+4) - A(m-5,n-2*bw+8) * A(m+6,n-2*bw+3))
    A(m+2,n-2*bw+8) = dummy * (A(m+2,n-2*bw+8) - A(m-1,n-2*bw+8) * A(m+3,n-2*bw+7) - A(m-2,n-2*bw+8) * A(m+4,n-2*bw+6)&
         - A(m-3,n-2*bw+8) * A(m+5,n-2*bw+5) - A(m-4,n-2*bw+8) * A(m+6,n-2*bw+4))
    A(m+3,n-2*bw+8) = dummy * (A(m+3,n-2*bw+8) - A(m-1,n-2*bw+8) * A(m+4,n-2*bw+7) - A(m-2,n-2*bw+8) * A(m+5,n-2*bw+6)&
         - A(m-3,n-2*bw+8) * A(m+6,n-2*bw+5))
    A(m+4,n-2*bw+8) = dummy * (A(m+4,n-2*bw+8) - A(m-1,n-2*bw+8) * A(m+5,n-2*bw+7) - A(m-2,n-2*bw+8) * A(m+6,n-2*bw+6))

    A(m  ,n-2*bw+9)   = A(m  ,n-2*bw+9)  - A(m-1,n-2*bw+9)  * A(m+1,n-2*bw+8) - A(m-2,n-2*bw+9)  * A(m+2,n-2*bw+7) &
         - A(m-3,n-2*bw+9)  * A(m+3,n-2*bw+6) - A(m-4,n-2*bw+9)  * A(m+4,n-2*bw+5) - A(m-5,n-2*bw+9)  * A(m+5,n-2*bw+4)&
         - A(m-6,n-2*bw+9)  * A(m+6,n-2*bw+3)    
    A(m-1,n-2*bw+10)  = A(m-1,n-2*bw+10) - A(m-2,n-2*bw+10) * A(m+1,n-2*bw+8) - A(m-3,n-2*bw+10) * A(m+2,n-2*bw+7) &
         - A(m-4,n-2*bw+10) * A(m+3,n-2*bw+6) - A(m-5,n-2*bw+10) * A(m+4,n-2*bw+5) - A(m-6,n-2*bw+10) * A(m+5,n-2*bw+4)&
         - A(m+3,n-2*bw+10) * A(m+6,n-2*bw+3)
    A(m-2,n-2*bw+11)  = A(m-2,n-2*bw+11) - A(m-3,n-2*bw+11) * A(m+1,n-2*bw+8) - A(m-4,n-2*bw+11) * A(m+2,n-2*bw+7) &
         - A(m-5,n-2*bw+11) * A(m+3,n-2*bw+6) - A(m-6,n-2*bw+11) * A(m+4,n-2*bw+5) - A(m+2,n-2*bw+11) * A(m+5,n-2*bw+4)&
         - A(m+3,n-2*bw+11) * A(m+6,n-2*bw+3)
    A(m-3,n-2*bw+12)  = A(m-3,n-2*bw+12) - A(m-4,n-2*bw+12) * A(m+1,n-2*bw+8) - A(m-5,n-2*bw+12) * A(m+2,n-2*bw+7)&
         - A(m-6,n-2*bw+12) * A(m+3,n-2*bw+6) - A(m+1,n-2*bw+12) * A(m+4,n-2*bw+5) - A(m+2,n-2*bw+12) * A(m+5,n-2*bw+4)&
         - A(m+3,n-2*bw+12) * A(m+6,n-2*bw+3)

    dummy = 1.d0 / A(m,n-2*bw+9)
    A(m+1,n-2*bw+9) = dummy * (A(m+1,n-2*bw+9) - A(m-1,n-2*bw+9) * A(m+2,n-2*bw+8) - A(m-2,n-2*bw+9) * A(m+3,n-2*bw+7)&
         - A(m-3,n-2*bw+9) * A(m+4,n-2*bw+6) - A(m-4,n-2*bw+9) * A(m+5,n-2*bw+5) - A(m-5,n-2*bw+9) * A(m+6,n-2*bw+4))
    A(m+2,n-2*bw+9) = dummy * (A(m+2,n-2*bw+9) - A(m-1,n-2*bw+9) * A(m+3,n-2*bw+8) - A(m-2,n-2*bw+9) * A(m+4,n-2*bw+7)&
         - A(m-3,n-2*bw+9) * A(m+5,n-2*bw+6) - A(m-4,n-2*bw+9) * A(m+6,n-2*bw+5))
    A(m+3,n-2*bw+9) = dummy * (A(m+3,n-2*bw+9) - A(m-1,n-2*bw+9) * A(m+4,n-2*bw+8) - A(m-2,n-2*bw+9) * A(m+5,n-2*bw+7)&
         - A(m-3,n-2*bw+9) * A(m+6,n-2*bw+6))

    A(m  ,n-2*bw+10)  = A(m  ,n-2*bw+10) - A(m-1,n-2*bw+10) * A(m+1,n-2*bw+9) - A(m-2,n-2*bw+10) * A(m+2,n-2*bw+8)&
         - A(m-3,n-2*bw+10) * A(m+3,n-2*bw+7) - A(m-4,n-2*bw+10) * A(m+4,n-2*bw+6) - A(m-5,n-2*bw+10) * A(m+5,n-2*bw+5) &
         - A(m-6,n-2*bw+10) * A(m+6,n-2*bw+4)
    A(m-1,n-2*bw+11)  = A(m-1,n-2*bw+11) - A(m-2,n-2*bw+11) * A(m+1,n-2*bw+9) - A(m-3,n-2*bw+11) * A(m+2,n-2*bw+8)&
         - A(m-4,n-2*bw+11) * A(m+3,n-2*bw+7) - A(m-5,n-2*bw+11) * A(m+4,n-2*bw+6) - A(m-6,n-2*bw+11) * A(m+5,n-2*bw+5)&
         - A(m+2,n-2*bw+11) * A(m+6,n-2*bw+4)
    A(m-2,n-2*bw+12)  = A(m-2,n-2*bw+12) - A(m-3,n-2*bw+12) * A(m+1,n-2*bw+9) - A(m-4,n-2*bw+12) * A(m+2,n-2*bw+8)&
         - A(m-5,n-2*bw+12) * A(m+3,n-2*bw+7) - A(m-6,n-2*bw+12) * A(m+4,n-2*bw+6) - A(m+1,n-2*bw+12) * A(m+5,n-2*bw+5) &
         - A(m+2,n-2*bw+12) * A(m+6,n-2*bw+4)

    dummy = 1.d0 / A(m,n-2*bw+10)
    A(m+1,n-2*bw+10) = dummy * (A(m+1,n-2*bw+10) - A(m-1,n-2*bw+10) * A(m+2,n-2*bw+9) - A(m-2,n-2*bw+10) * A(m+3,n-2*bw+8)&
         - A(m-3,n-2*bw+10) * A(m+4,n-2*bw+7) - A(m-4,n-2*bw+10) * A(m+5,n-2*bw+6) - A(m-5,n-2*bw+10) * A(m+6,n-2*bw+5))
    A(m+2,n-2*bw+10) = dummy * (A(m+2,n-2*bw+10) - A(m-1,n-2*bw+10) * A(m+3,n-2*bw+9) - A(m-2,n-2*bw+10) * A(m+4,n-2*bw+8)&
         - A(m-3,n-2*bw+10) * A(m+5,n-2*bw+7) - A(m-4,n-2*bw+10) * A(m+6,n-2*bw+6))

    A(m  ,n-2*bw+11)  = A(m  ,n-2*bw+11) - A(m-1,n-2*bw+11) * A(m+1,n-2*bw+10) - A(m-2,n-2*bw+11) * A(m+2,n-2*bw+9) &
         - A(m-3,n-2*bw+11) * A(m+3,n-2*bw+8) - A(m-4,n-2*bw+11) * A(m+4,n-2*bw+7) - A(m-5,n-2*bw+11) * A(m+5,n-2*bw+6)&
         - A(m-6,n-2*bw+11) * A(m+6,n-2*bw+5)    
    A(m-1,n-2*bw+12)  = A(m-1,n-2*bw+12) - A(m-2,n-2*bw+12) * A(m+1,n-2*bw+10) - A(m-3,n-2*bw+12) * A(m+2,n-2*bw+9)&
         - A(m-4,n-2*bw+12) * A(m+3,n-2*bw+8) - A(m-5,n-2*bw+12) * A(m+4,n-2*bw+7) - A(m-6,n-2*bw+12) * A(m+5,n-2*bw+6)&
         - A(m+1,n-2*bw+12) * A(m+6,n-2*bw+5) 

    dummy = 1.d0 / A(m,n-2*bw+11)
    A(m+1,n-2*bw+11) = dummy * (A(m+1,n-2*bw+11) - A(m-1,n-2*bw+11) * A(m+2,n-2*bw+10) - A(m-2,n-2*bw+11) * A(m+3,n-2*bw+9)&
         - A(m-3,n-2*bw+11) * A(m+4,n-2*bw+8) - A(m-4,n-2*bw+11) * A(m+5,n-2*bw+7) - A(m-5,n-2*bw+11) * A(m+6,n-2*bw+6))

    A(m  ,n-2*bw+12)  = A(m  ,n-2*bw+12) - A(m-1,n-2*bw+12) * A(m+1,n-2*bw+11) - A(m-2,n-2*bw+12) * A(m+2,n-2*bw+10) &
         - A(m-3,n-2*bw+12) * A(m+3,n-2*bw+9) - A(m-4,n-2*bw+12) * A(m+4,n-2*bw+8) - A(m-5,n-2*bw+12) * A(m+5,n-2*bw+7) &
         - A(m-6,n-2*bw+12) * A(m+6,n-2*bw+6)

    return

  end subroutine fixed_bandec6


  subroutine fixed_bandec7(A,n)
    implicit none
    integer :: n
    real(c_double) :: A(15,n)
    integer :: i
    integer :: m,bw
    real(c_double) :: dummy

    m = 8
    bw = 7

    dummy = 1.d0 / A(m,1)
    A(m+1,1) = dummy * A(m+1,1)
    A(m+2,1) = dummy * A(m+2,1)
    A(m+3,1) = dummy * A(m+3,1)
    A(m+4,1) = dummy * A(m+4,1)
    A(m+5,1) = dummy * A(m+5,1)
    A(m+6,1) = dummy * A(m+6,1)
    A(m+7,1) = dummy * A(m+7,1)
    A(m-1,1) = dummy * A(m-1,1)
    A(m-2,1) = dummy * A(m-2,1)
    A(m-3,1) = dummy * A(m-3,1)
    A(m-4,1) = dummy * A(m-4,1)
    A(m-5,1) = dummy * A(m-5,1)
    A(m-6,1) = dummy * A(m-6,1)
    A(m-7,1) = dummy * A(m-7,1)

    A(m  ,2) = A(m  ,2) -  A(m-1,2) * A(m+1,1)
    A(m-1,3) = A(m-1,3) -  A(m-2,3) * A(m+1,1)
    A(m-2,4) = A(m-2,4) -  A(m-3,4) * A(m+1,1)
    A(m-3,5) = A(m-3,5) -  A(m-4,5) * A(m+1,1)
    A(m-4,6) = A(m-4,6) -  A(m-5,6) * A(m+1,1)
    A(m-5,7) = A(m-5,7) -  A(m-6,7) * A(m+1,1)
    A(m-6,8) = A(m-6,8) -  A(m-7,8) * A(m+1,1)

    dummy = 1.d0 / A(m,2)
    A(m+1,2) = dummy * (A(m+1,2) -  A(m-1,2) * A(m+2,1)) 
    A(m+2,2) = dummy * (A(m+2,2) -  A(m-1,2) * A(m+3,1)) 
    A(m+3,2) = dummy * (A(m+3,2) -  A(m-1,2) * A(m+4,1)) 
    A(m+4,2) = dummy * (A(m+4,2) -  A(m-1,2) * A(m+5,1)) 
    A(m+5,2) = dummy * (A(m+5,2) -  A(m-1,2) * A(m+6,1)) 
    A(m+6,2) = dummy * (A(m+6,2) -  A(m-1,2) * A(m+7,1)) 
    A(m+7,2) = dummy * (A(m+7,2) -  A(m-1,2) * A(m-1,1)) 
    A(m-2,2) = dummy * (A(m-2,2) -  A(m-1,2) * A(m-2,1)) 
    A(m-3,2) = dummy * (A(m-3,2) -  A(m-1,2) * A(m-3,1)) 
    A(m-4,2) = dummy * (A(m-4,2) -  A(m-1,2) * A(m-4,1)) 
    A(m-5,2) = dummy * (A(m-5,2) -  A(m-1,2) * A(m-5,1)) 
    A(m-6,2) = dummy * (A(m-6,2) -  A(m-1,2) * A(m-6,1)) 
    A(m-7,2) = dummy * (A(m-7,2) -  A(m-1,2) * A(m-7,1)) 

    A(m  , 3) = A(m  , 3) - A(m-1,3) * A(m+1,2) - A(m-2,3) * A(m+2,1) 
    A(m-1, 4) = A(m-1, 4) - A(m-2,4) * A(m+1,2) - A(m-3,4) * A(m+2,1) 
    A(m-2, 5) = A(m-2, 5) - A(m-3,5) * A(m+1,2) - A(m-4,5) * A(m+2,1)
    A(m-3, 6) = A(m-3, 6) - A(m-4,6) * A(m+1,2) - A(m-5,6) * A(m+2,1)
    A(m-4, 7) = A(m-4, 7) - A(m-5,7) * A(m+1,2) - A(m-6,7) * A(m+2,1)
    A(m-5, 8) = A(m-5, 8) - A(m-6,8) * A(m+1,2) - A(m-7,8) * A(m+2,1)
    A(m-6, 9) = A(m-6, 9) - A(m-7,9) * A(m+1,2) 

    dummy = 1.d0 / A(m,3)
    A(m+1,3) = dummy * (A(m+1,3) - A(m-1,3) * A(m+2,2) - A(m-2,3) * A(m+3,1))
    A(m+2,3) = dummy * (A(m+2,3) - A(m-1,3) * A(m+3,2) - A(m-2,3) * A(m+4,1))
    A(m+3,3) = dummy * (A(m+3,3) - A(m-1,3) * A(m+4,2) - A(m-2,3) * A(m+5,1))
    A(m+4,3) = dummy * (A(m+4,3) - A(m-1,3) * A(m+5,2) - A(m-2,3) * A(m+6,1))
    A(m+5,3) = dummy * (A(m+5,3) - A(m-1,3) * A(m+6,2) - A(m-2,3) * A(m+7,1))
    A(m+6,3) = dummy * (A(m+6,3) - A(m-1,3) * A(m+7,2) - A(m-2,3) * A(m-1,1))
    A(m+7,3) = dummy * (A(m+7,3) - A(m-1,3) * A(m-2,2) - A(m-2,3) * A(m-2,1))
    A(m-3,3) = dummy * (A(m-3,3) - A(m-1,3) * A(m-3,2) - A(m-2,3) * A(m-3,1))
    A(m-4,3) = dummy * (A(m-4,3) - A(m-1,3) * A(m-4,2) - A(m-2,3) * A(m-4,1))
    A(m-5,3) = dummy * (A(m-5,3) - A(m-1,3) * A(m-5,2) - A(m-2,3) * A(m-5,1))
    A(m-6,3) = dummy * (A(m-6,3) - A(m-1,3) * A(m-6,2) - A(m-2,3) * A(m-6,1))
    A(m-7,3) = dummy * (A(m-7,3) - A(m-1,3) * A(m-7,2) - A(m-2,3) * A(m-7,1))

    A(m  , 4) = A(m  , 4) - A(m-1, 4) * A(m+1,3) - A(m-2,4) * A(m+2,2) - A(m-3,4) * A(m+3,1)  
    A(m-1, 5) = A(m-1, 5) - A(m-2, 5) * A(m+1,3) - A(m-3,5) * A(m+2,2) - A(m-4,5) * A(m+3,1)
    A(m-2, 6) = A(m-2, 6) - A(m-3, 6) * A(m+1,3) - A(m-4,6) * A(m+2,2) - A(m-5,6) * A(m+3,1)
    A(m-3, 7) = A(m-3, 7) - A(m-4, 7) * A(m+1,3) - A(m-5,7) * A(m+2,2) - A(m-6,7) * A(m+3,1)
    A(m-4, 8) = A(m-4, 8) - A(m-5, 8) * A(m+1,3) - A(m-6,8) * A(m+2,2) - A(m-7,8) * A(m+3,1)
    A(m-5, 9) = A(m-5, 9) - A(m-6, 9) * A(m+1,3) - A(m-7,9) * A(m+2,2)
    A(m-6,10) = A(m-6,10) - A(m-7,10) * A(m+1,3) 

    dummy = 1.d0 / A(m,4)
    A(m+1,4) = dummy * (A(m+1,4) - A(m-1,4) * A(m+2,3) - A(m-2,4) * A(m+3,2) - A(m-3,4) * A(m+4,1))
    A(m+2,4) = dummy * (A(m+2,4) - A(m-1,4) * A(m+3,3) - A(m-2,4) * A(m+4,2) - A(m-3,4) * A(m+5,1))
    A(m+3,4) = dummy * (A(m+3,4) - A(m-1,4) * A(m+4,3) - A(m-2,4) * A(m+5,2) - A(m-3,4) * A(m+6,1))
    A(m+4,4) = dummy * (A(m+4,4) - A(m-1,4) * A(m+5,3) - A(m-2,4) * A(m+6,2) - A(m-3,4) * A(m+7,1))
    A(m+5,4) = dummy * (A(m+5,4) - A(m-1,4) * A(m+6,3) - A(m-2,4) * A(m+7,2) - A(m-3,4) * A(m-1,1))
    A(m+6,4) = dummy * (A(m+6,4) - A(m-1,4) * A(m+7,3) - A(m-2,4) * A(m-2,2) - A(m-3,4) * A(m-2,1))
    A(m+7,4) = dummy * (A(m+7,4) - A(m-1,4) * A(m-3,3) - A(m-2,4) * A(m-3,2) - A(m-3,4) * A(m-3,1))
    A(m-4,4) = dummy * (A(m-4,4) - A(m-1,4) * A(m-4,3) - A(m-2,4) * A(m-4,2) - A(m-3,4) * A(m-4,1))
    A(m-5,4) = dummy * (A(m-5,4) - A(m-1,4) * A(m-5,3) - A(m-2,4) * A(m-5,2) - A(m-3,4) * A(m-5,1))
    A(m-6,4) = dummy * (A(m-6,4) - A(m-1,4) * A(m-6,3) - A(m-2,4) * A(m-6,2) - A(m-3,4) * A(m-6,1))
    A(m-7,4) = dummy * (A(m-7,4) - A(m-1,4) * A(m-7,3) - A(m-2,4) * A(m-7,2) - A(m-3,4) * A(m-7,1))

    A(m  , 5) = A(m  , 5) - A(m-1, 5) * A(m+1,4) - A(m-2, 5) * A(m+2,3) - A(m-3,5) * A(m+3,2) - A(m-4,5) * A(m+4,1)    
    A(m-1, 6) = A(m-1, 6) - A(m-2, 6) * A(m+1,4) - A(m-3, 6) * A(m+2,3) - A(m-4,6) * A(m+3,2) - A(m-5,6) * A(m+4,1)    
    A(m-2, 7) = A(m-2, 7) - A(m-3, 7) * A(m+1,4) - A(m-4, 7) * A(m+2,3) - A(m-5,7) * A(m+3,2) - A(m-6,7) * A(m+4,1) 
    A(m-3, 8) = A(m-3, 8) - A(m-4, 8) * A(m+1,4) - A(m-5, 8) * A(m+2,3) - A(m-6,8) * A(m+3,2) - A(m-7,8) * A(m+4,1) 
    A(m-4, 9) = A(m-4, 9) - A(m-5, 9) * A(m+1,4) - A(m-6, 9) * A(m+2,3) - A(m-7,9) * A(m+3,2)
    A(m-5,10) = A(m-5,10) - A(m-6,10) * A(m+1,4) - A(m-7,10) * A(m+2,3)
    A(m-6,11) = A(m-6,11) - A(m-7,11) * A(m+1,4) 

    dummy = 1.d0 / A(m,5)
    A(m+1,5) = dummy * (A(m+1,5) - A(m-1,5) * A(m+2,4) - A(m-2,5) * A(m+3,3) - A(m-3,5) * A(m+4,2) - A(m-4,5) * A(m+5,1))
    A(m+2,5) = dummy * (A(m+2,5) - A(m-1,5) * A(m+3,4) - A(m-2,5) * A(m+4,3) - A(m-3,5) * A(m+5,2) - A(m-4,5) * A(m+6,1))
    A(m+3,5) = dummy * (A(m+3,5) - A(m-1,5) * A(m+4,4) - A(m-2,5) * A(m+5,3) - A(m-3,5) * A(m+6,2) - A(m-4,5) * A(m+7,1))
    A(m+4,5) = dummy * (A(m+4,5) - A(m-1,5) * A(m+5,4) - A(m-2,5) * A(m+6,3) - A(m-3,5) * A(m+7,2) - A(m-4,5) * A(m-1,1))
    A(m+5,5) = dummy * (A(m+5,5) - A(m-1,5) * A(m+6,4) - A(m-2,5) * A(m+7,3) - A(m-3,5) * A(m-2,2) - A(m-4,5) * A(m-2,1))
    A(m+6,5) = dummy * (A(m+6,5) - A(m-1,5) * A(m+7,4) - A(m-2,5) * A(m-3,3) - A(m-3,5) * A(m-3,2) - A(m-4,5) * A(m-3,1))
    A(m+7,5) = dummy * (A(m+7,5) - A(m-1,5) * A(m-4,4) - A(m-2,5) * A(m-4,3) - A(m-3,5) * A(m-4,2) - A(m-4,5) * A(m-4,1))
    A(m-5,5) = dummy * (A(m-5,5) - A(m-1,5) * A(m-5,4) - A(m-2,5) * A(m-5,3) - A(m-3,5) * A(m-5,2) - A(m-4,5) * A(m-5,1))
    A(m-6,5) = dummy * (A(m-6,5) - A(m-1,5) * A(m-6,4) - A(m-2,5) * A(m-6,3) - A(m-3,5) * A(m-6,2) - A(m-4,5) * A(m-6,1))
    A(m-7,5) = dummy * (A(m-7,5) - A(m-1,5) * A(m-7,4) - A(m-2,5) * A(m-7,3) - A(m-3,5) * A(m-7,2) - A(m-4,5) * A(m-7,1))

    A(m  , 6) = A(m  , 6) - A(m-1, 6) * A(m+1,5) - A(m-2, 6) * A(m+2,4) - A(m-3, 6) * A(m+3,3) - A(m-4,6) * A(m+4,2)&
         - A(m-5,6) * A(m+5,1) 
    A(m-1, 7) = A(m-1, 7) - A(m-2, 7) * A(m+1,5) - A(m-3, 7) * A(m+2,4) - A(m-4, 7) * A(m+3,3) - A(m-5,7) * A(m+4,2)&
         - A(m-6,7) * A(m+5,1) 
    A(m-2, 8) = A(m-2, 8) - A(m-3, 8) * A(m+1,5) - A(m-4, 8) * A(m+2,4) - A(m-5, 8) * A(m+3,3) - A(m-6,8) * A(m+4,2)&
         - A(m-7,8) * A(m+5,1)
    A(m-3, 9) = A(m-3, 9) - A(m-4, 9) * A(m+1,5) - A(m-5, 9) * A(m+2,4) - A(m-6, 9) * A(m+3,3) - A(m-7,9) * A(m+4,2) 
    A(m-4,10) = A(m-4,10) - A(m-5,10) * A(m+1,5) - A(m-6,10) * A(m+2,4) - A(m-7,10) * A(m+3,3)
    A(m-5,11) = A(m-5,11) - A(m-6,11) * A(m+1,5) - A(m-7,11) * A(m+2,4)
    A(m-6,12) = A(m-6,12) - A(m-7,12) * A(m+1,5) 

    dummy = 1.d0 / A(m,6)
    A(m+1,6) = dummy * (A(m+1,6) - A(m-1,6) * A(m+2,5) - A(m-2,6) * A(m+3,4) - A(m-3,6) * A(m+4,3) - A(m-4,6) * A(m+5,2)&
         - A(m-5,6) * A(m+6,1))
    A(m+2,6) = dummy * (A(m+2,6) - A(m-1,6) * A(m+3,5) - A(m-2,6) * A(m+4,4) - A(m-3,6) * A(m+5,3) - A(m-4,6) * A(m+6,2)&
         - A(m-5,6) * A(m+7,1))
    A(m+3,6) = dummy * (A(m+3,6) - A(m-1,6) * A(m+4,5) - A(m-2,6) * A(m+5,4) - A(m-3,6) * A(m+6,3) - A(m-4,6) * A(m+7,2)&
         - A(m-5,6) * A(m-1,1))
    A(m+4,6) = dummy * (A(m+4,6) - A(m-1,6) * A(m+5,5) - A(m-2,6) * A(m+6,4) - A(m-3,6) * A(m+7,3) - A(m-4,6) * A(m-2,2)&
         - A(m-5,6) * A(m-2,1))
    A(m+5,6) = dummy * (A(m+5,6) - A(m-1,6) * A(m+6,5) - A(m-2,6) * A(m+7,4) - A(m-3,6) * A(m-3,3) - A(m-4,6) * A(m-3,2)&
         - A(m-5,6) * A(m-3,1))
    A(m+6,6) = dummy * (A(m+6,6) - A(m-1,6) * A(m+7,5) - A(m-2,6) * A(m-4,4) - A(m-3,6) * A(m-4,3) - A(m-4,6) * A(m-4,2)&
         - A(m-5,6) * A(m-4,1))
    A(m+7,6) = dummy * (A(m+7,6) - A(m-1,6) * A(m-5,5) - A(m-2,6) * A(m-5,4) - A(m-3,6) * A(m-5,3) - A(m-4,6) * A(m-5,2)&
         - A(m-5,6) * A(m-5,1))
    A(m-6,6) = dummy * (A(m-6,6) - A(m-1,6) * A(m-6,5) - A(m-2,6) * A(m-6,4) - A(m-3,6) * A(m-6,3) - A(m-4,6) * A(m-6,2)&
         - A(m-5,6) * A(m-6,1))
    A(m-7,6) = dummy * (A(m-7,6) - A(m-1,6) * A(m-7,5) - A(m-2,6) * A(m-7,4) - A(m-3,6) * A(m-7,3) - A(m-4,6) * A(m-7,2)&
         - A(m-5,6) * A(m-7,1))

    A(m  , 7) = A(m  , 7) - A(m-1, 7) * A(m+1,6) - A(m-2, 7) * A(m+2,5) - A(m-3, 7) * A(m+3,4) - A(m-4, 7) * A(m+4,3)&
         - A(m-5,7) * A(m+5,2) - A(m-6,7) * A(m+6,1) 
    A(m-1, 8) = A(m-1, 8) - A(m-2, 8) * A(m+1,6) - A(m-3, 8) * A(m+2,5) - A(m-4, 8) * A(m+3,4) - A(m-5, 8) * A(m+4,3)&
         - A(m-6,8) * A(m+5,2) - A(m-7,8) * A(m+6,1) 
    A(m-2, 9) = A(m-2, 9) - A(m-3, 9) * A(m+1,6) - A(m-4, 9) * A(m+2,5) - A(m-5, 9) * A(m+3,4) - A(m-6, 9) * A(m+4,3)&
         - A(m-7,9) * A(m+5,2)
    A(m-3,10) = A(m-3,10) - A(m-4,10) * A(m+1,6) - A(m-5,10) * A(m+2,5) - A(m-6,10) * A(m+3,4) - A(m-7,10) * A(m+4,3) 
    A(m-4,11) = A(m-4,11) - A(m-5,11) * A(m+1,6) - A(m-6,11) * A(m+2,5) - A(m-7,11) * A(m+3,4)
    A(m-5,12) = A(m-5,12) - A(m-6,12) * A(m+1,6) - A(m-7,12) * A(m+2,5)
    A(m-6,13) = A(m-6,13) - A(m-7,13) * A(m+1,6) 




    ! Subupper Region

    dummy = 1.d0 / A(m,7)
    A(m+1,7) = dummy * (A(m+1,7) - A(m-1,7) * A(m+2,6) - A(m-2,7) * A(m+3,5) - A(m-3,7) * A(m+4,4) - A(m-4,7) * A(m+5,3)&
         - A(m-5,7) * A(m+6,2) - A(m-6,7) * A(m+7,1))
    A(m+2,7) = dummy * (A(m+2,7) - A(m-1,7) * A(m+3,6) - A(m-2,7) * A(m+4,5) - A(m-3,7) * A(m+5,4) - A(m-4,7) * A(m+6,3)&
         - A(m-5,7) * A(m+7,2) - A(m-6,7) * A(m-1,1))
    A(m+3,7) = dummy * (A(m+3,7) - A(m-1,7) * A(m+4,6) - A(m-2,7) * A(m+5,5) - A(m-3,7) * A(m+6,4) - A(m-4,7) * A(m+7,3)&
         - A(m-5,7) * A(m-2,2) - A(m-6,7) * A(m-2,1))
    A(m+4,7) = dummy * (A(m+4,7) - A(m-1,7) * A(m+5,6) - A(m-2,7) * A(m+6,5) - A(m-3,7) * A(m+7,4) - A(m-4,7) * A(m-3,3)&
         - A(m-5,7) * A(m-3,2) - A(m-6,7) * A(m-3,1))
    A(m+5,7) = dummy * (A(m+5,7) - A(m-1,7) * A(m+6,6) - A(m-2,7) * A(m+7,5) - A(m-3,7) * A(m-4,4) - A(m-4,7) * A(m-4,3)&
         - A(m-5,7) * A(m-4,2) - A(m-6,7) * A(m-4,1))
    A(m+6,7) = dummy * (A(m+6,7) - A(m-1,7) * A(m+7,6) - A(m-2,7) * A(m-5,5) - A(m-3,7) * A(m-5,4) - A(m-4,7) * A(m-5,3)&
         - A(m-5,7) * A(m-5,2) - A(m-6,7) * A(m-5,1))
    A(m+7,7) = dummy * (A(m+7,7) - A(m-1,7) * A(m-6,6) - A(m-2,7) * A(m-6,5) - A(m-3,7) * A(m-6,4) - A(m-4,7) * A(m-6,3)&
         - A(m-5,7) * A(m-6,2) - A(m-6,7) * A(m-6,1))
    A(m-7,7) = dummy * (A(m-7,7) - A(m-1,7) * A(m-7,6) - A(m-2,7) * A(m-7,5) - A(m-3,7) * A(m-7,4) - A(m-4,7) * A(m-7,3)&
         - A(m-5,7) * A(m-7,2) - A(m-6,7) * A(m-7,1))

    A(m  , 8) = A(m  , 8) - A(m-1, 8) * A(m+1,7) - A(m-2, 8) * A(m+2,6) - A(m-3, 8) * A(m+3,5) - A(m-4, 8) * A(m+4,4)&
         - A(m-5, 8) * A(m+5,3) - A(m-6,8) * A(m+6,2) - A(m-7,8) * A(m+7,1)
    A(m-1, 9) = A(m-1, 9) - A(m-2, 9) * A(m+1,7) - A(m-3, 9) * A(m+2,6) - A(m-4, 9) * A(m+3,5) - A(m-5, 9) * A(m+4,4)&
         - A(m-6, 9) * A(m+5,3) - A(m-7,9) * A(m+6,2)
    A(m-2,10) = A(m-2,10) - A(m-3,10) * A(m+1,7) - A(m-4,10) * A(m+2,6) - A(m-5,10) * A(m+3,5) - A(m-6,10) * A(m+4,4)&
         - A(m-7,10) * A(m+5,3)
    A(m-3,11) = A(m-3,11) - A(m-4,11) * A(m+1,7) - A(m-5,11) * A(m+2,6) - A(m-6,11) * A(m+3,5) - A(m-7,11) * A(m+4,4)
    A(m-4,12) = A(m-4,12) - A(m-5,12) * A(m+1,7) - A(m-6,12) * A(m+2,6) - A(m-7,12) * A(m+3,5)
    A(m-5,13) = A(m-5,13) - A(m-6,13) * A(m+1,7) - A(m-7,13) * A(m+2,6) 
    A(m-6,14) = A(m-6,14) - A(m-7,14) * A(m+1,7) 


    dummy = 1.d0 / A(m,8)
    A(m+1,8) = dummy * (A(m+1,8) - A(m-1,8) * A(m+2,7) - A(m-2,8) * A(m+3,6) - A(m-3,8) * A(m+4,5) - A(m-4,8) * A(m+5,4)&
         - A(m-5,8) * A(m+6,3) - A(m-6,8) * A(m+7,2) - A(m-7,8) * A(m-1,1))
    A(m+2,8) = dummy * (A(m+2,8) - A(m-1,8) * A(m+3,7) - A(m-2,8) * A(m+4,6) - A(m-3,8) * A(m+5,5) - A(m-4,8) * A(m+6,4)&
         - A(m-5,8) * A(m+7,3) - A(m-6,8) * A(m-2,2) - A(m-7,8) * A(m-2,1))
    A(m+3,8) = dummy * (A(m+3,8) - A(m-1,8) * A(m+4,7) - A(m-2,8) * A(m+5,6) - A(m-3,8) * A(m+6,5) - A(m-4,8) * A(m+7,4)&
         - A(m-5,8) * A(m-3,3) - A(m-6,8) * A(m-3,2) - A(m-7,8) * A(m-3,1))
    A(m+4,8) = dummy * (A(m+4,8) - A(m-1,8) * A(m+5,7) - A(m-2,8) * A(m+6,6) - A(m-3,8) * A(m+7,5) - A(m-4,8) * A(m-4,4)&
         - A(m-5,8) * A(m-4,3) - A(m-6,8) * A(m-4,2) - A(m-7,8) * A(m-4,1))
    A(m+5,8) = dummy * (A(m+5,8) - A(m-1,8) * A(m+6,7) - A(m-2,8) * A(m+7,6) - A(m-3,8) * A(m-5,5) - A(m-4,8) * A(m-5,4)&
         - A(m-5,8) * A(m-5,3) - A(m-6,8) * A(m-5,2) - A(m-7,8) * A(m-5,1))
    A(m+6,8) = dummy * (A(m+6,8) - A(m-1,8) * A(m+7,7) - A(m-2,8) * A(m-6,6) - A(m-3,8) * A(m-6,5) - A(m-4,8) * A(m-6,4)&
         - A(m-5,8) * A(m-6,3) - A(m-6,8) * A(m-6,2) - A(m-7,8) * A(m-6,1))
    A(m+7,8) = dummy * (A(m+7,8) - A(m-1,8) * A(m-7,7) - A(m-2,8) * A(m-7,6) - A(m-3,8) * A(m-7,5) - A(m-4,8) * A(m-7,4)&
         - A(m-5,8) * A(m-7,3) - A(m-6,8) * A(m-7,2) - A(m-7,8) * A(m-7,1))

    A(m  , 9) = A(m  , 9) - A(m-1, 9) * A(m+1,8) - A(m-2, 9) * A(m+2,7) - A(m-3, 9) * A(m+3,6) - A(m-4, 9) * A(m+4,5)&
         - A(m-5, 9) * A(m+5,4) - A(m-6, 9) * A(m+6,3) - A(m-7,9) * A(m+7,2)
    A(m-1,10) = A(m-1,10) - A(m-2,10) * A(m+1,8) - A(m-3,10) * A(m+2,7) - A(m-4,10) * A(m+3,6) - A(m-5,10) * A(m+4,5)&
         - A(m-6,10) * A(m+5,4) - A(m-7,10) * A(m+6,3)
    A(m-2,11) = A(m-2,11) - A(m-3,11) * A(m+1,8) - A(m-4,11) * A(m+2,7) - A(m-5,11) * A(m+3,6) - A(m-6,11) * A(m+4,5)&
         - A(m-7,11) * A(m+5,4)
    A(m-3,12) = A(m-3,12) - A(m-4,12) * A(m+1,8) - A(m-5,12) * A(m+2,7) - A(m-6,12) * A(m+3,6) - A(m-7,12) * A(m+4,5)
    A(m-4,13) = A(m-4,13) - A(m-5,13) * A(m+1,8) - A(m-6,13) * A(m+2,7) - A(m-7,13) * A(m+3,6)
    A(m-5,14) = A(m-5,14) - A(m-6,14) * A(m+1,8) - A(m-7,14) * A(m+2,7) 
    A(m-6,15) = A(m-6,15) - A(m-7,15) * A(m+1,8) 

    dummy = 1.d0 / A(m,9)
    A(m+1,9) = dummy * (A(m+1,9) - A(m-1,9) * A(m+2,8) - A(m-2,9) * A(m+3,7) - A(m-3,9) * A(m+4,6) - A(m-4,9) * A(m+5,5)&
         - A(m-5,9) * A(m+6,4) - A(m-6,9) * A(m+7,3) - A(m-7,9) * A(m-2,2))
    A(m+2,9) = dummy * (A(m+2,9) - A(m-1,9) * A(m+3,8) - A(m-2,9) * A(m+4,7) - A(m-3,9) * A(m+5,6) - A(m-4,9) * A(m+6,5)&
         - A(m-5,9) * A(m+7,4) - A(m-6,9) * A(m-3,3) - A(m-7,9) * A(m-3,2))
    A(m+3,9) = dummy * (A(m+3,9) - A(m-1,9) * A(m+4,8) - A(m-2,9) * A(m+5,7) - A(m-3,9) * A(m+6,6) - A(m-4,9) * A(m+7,5)&
         - A(m-5,9) * A(m-4,4) - A(m-6,9) * A(m-4,3) - A(m-7,9) * A(m-4,2))
    A(m+4,9) = dummy * (A(m+4,9) - A(m-1,9) * A(m+5,8) - A(m-2,9) * A(m+6,7) - A(m-3,9) * A(m+7,6) - A(m-4,9) * A(m-5,5)&
         - A(m-5,9) * A(m-5,4) - A(m-6,9) * A(m-5,3) - A(m-7,9) * A(m-5,2))
    A(m+5,9) = dummy * (A(m+5,9) - A(m-1,9) * A(m+6,8) - A(m-2,9) * A(m+7,7) - A(m-3,9) * A(m-6,6) - A(m-4,9) * A(m-6,5)&
         - A(m-5,9) * A(m-6,4) - A(m-6,9) * A(m-6,3) - A(m-7,9) * A(m-6,2))
    A(m+6,9) = dummy * (A(m+6,9) - A(m-1,9) * A(m+7,8) - A(m-2,9) * A(m-7,7) - A(m-3,9) * A(m-7,6) - A(m-4,9) * A(m-7,5)&
         - A(m-5,9) * A(m-7,4) - A(m-6,9) * A(m-7,3) - A(m-7,9) * A(m-7,2))
    A(m+7,9) = dummy *  A(m+7,9) 

    A(m  ,10) = A(m  ,10) - A(m-1,10) * A(m+1,9) - A(m-2,10) * A(m+2,8) - A(m-3,10) * A(m+3,7) - A(m-4,10) * A(m+4,6)&
         - A(m-5,10) * A(m+5,5) - A(m-6,10) * A(m+6,4) - A(m-7,10) * A(m+7,3)
    A(m-1,11) = A(m-1,11) - A(m-2,11) * A(m+1,9) - A(m-3,11) * A(m+2,8) - A(m-4,11) * A(m+3,7) - A(m-5,11) * A(m+4,6)&
         - A(m-6,11) * A(m+5,5) - A(m-7,11) * A(m+6,4)
    A(m-2,12) = A(m-2,12) - A(m-3,12) * A(m+1,9) - A(m-4,12) * A(m+2,8) - A(m-5,12) * A(m+3,7) - A(m-6,12) * A(m+4,6)&
         - A(m-7,12) * A(m+5,5)
    A(m-3,13) = A(m-3,13) - A(m-4,13) * A(m+1,9) - A(m-5,13) * A(m+2,8) - A(m-6,13) * A(m+3,7) - A(m-7,13) * A(m+4,6)
    A(m-4,14) = A(m-4,14) - A(m-5,14) * A(m+1,9) - A(m-6,14) * A(m+2,8) - A(m-7,14) * A(m+3,7)
    A(m-5,15) = A(m-5,15) - A(m-6,15) * A(m+1,9) - A(m-7,15) * A(m+2,8) 
    A(m-6,16) = A(m-6,16) - A(m-7,16) * A(m+1,9) 

    dummy = 1.d0 / A(m,10)
    A(m+1,10) = dummy * (A(m+1,10) - A(m-1,10) * A(m+2,9) - A(m-2,10) * A(m+3,8) - A(m-3,10) * A(m+4,7) - A(m-4,10) * A(m+5,6)&
         - A(m-5,10) * A(m+6,5) - A(m-6,10) * A(m+7,4) - A(m-7,10) * A(m-3,3))
    A(m+2,10) = dummy * (A(m+2,10) - A(m-1,10) * A(m+3,9) - A(m-2,10) * A(m+4,8) - A(m-3,10) * A(m+5,7) - A(m-4,10) * A(m+6,6)&
         - A(m-5,10) * A(m+7,5) - A(m-6,10) * A(m-4,4) - A(m-7,10) * A(m-4,3))
    A(m+3,10) = dummy * (A(m+3,10) - A(m-1,10) * A(m+4,9) - A(m-2,10) * A(m+5,8) - A(m-3,10) * A(m+6,7) - A(m-4,10) * A(m+7,6)&
         - A(m-5,10) * A(m-5,5) - A(m-6,10) * A(m-5,4) - A(m-7,10) * A(m-5,3))
    A(m+4,10) = dummy * (A(m+4,10) - A(m-1,10) * A(m+5,9) - A(m-2,10) * A(m+6,8) - A(m-3,10) * A(m+7,7) - A(m-4,10) * A(m-6,6)&
         - A(m-5,10) * A(m-6,5) - A(m-6,10) * A(m-6,4) - A(m-7,10) * A(m-6,3))
    A(m+5,10) = dummy * (A(m+5,10) - A(m-1,10) * A(m+6,9) - A(m-2,10) * A(m+7,8) - A(m-3,10) * A(m-7,7) - A(m-4,10) * A(m-7,6)&
         - A(m-5,10) * A(m-7,5) - A(m-6,10) * A(m-7,4) - A(m-7,10) * A(m-7,3))
    A(m+6,10) = dummy * (A(m+6,10) - A(m-1,10) * A(m+7,9))
    A(m+7,10) = dummy *  A(m+7,10) 

    A(m  ,11) = A(m  ,11) - A(m-1,11) * A(m+1,10) - A(m-2,11) * A(m+2,9) - A(m-3,11) * A(m+3,8) - A(m-4,11) * A(m+4,7)&
         - A(m-5,11) * A(m+5,6) - A(m-6,11) * A(m+6,5) - A(m-7,11) * A(m+7,4)
    A(m-1,12) = A(m-1,12) - A(m-2,12) * A(m+1,10) - A(m-3,12) * A(m+2,9) - A(m-4,12) * A(m+3,8) - A(m-5,12) * A(m+4,7)&
         - A(m-6,12) * A(m+5,6) - A(m-7,12) * A(m+6,5)
    A(m-2,13) = A(m-2,13) - A(m-3,13) * A(m+1,10) - A(m-4,13) * A(m+2,9) - A(m-5,13) * A(m+3,8) - A(m-6,13) * A(m+4,7)&
         - A(m-7,13) * A(m+5,6)
    A(m-3,14) = A(m-3,14) - A(m-4,14) * A(m+1,10) - A(m-5,14) * A(m+2,9) - A(m-6,14) * A(m+3,8) - A(m-7,14) * A(m+4,7)
    A(m-4,15) = A(m-4,15) - A(m-5,15) * A(m+1,10) - A(m-6,15) * A(m+2,9) - A(m-7,15) * A(m+3,8)
    A(m-5,16) = A(m-5,16) - A(m-6,16) * A(m+1,10) - A(m-7,16) * A(m+2,9) 
    A(m-6,17) = A(m-6,17) - A(m-7,17) * A(m+1,10) 

    dummy = 1.d0 / A(m,11)
    A(m+1,11) = dummy * (A(m+1,11) - A(m-1,11) * A(m+2,10) - A(m-2,11) * A(m+3,9) - A(m-3,11) * A(m+4,8) - A(m-4,11) * A(m+5,7)&
         - A(m-5,11) * A(m+6,6) - A(m-6,11) * A(m+7,5) - A(m-7,11) * A(m-4,4))
    A(m+2,11) = dummy * (A(m+2,11) - A(m-1,11) * A(m+3,10) - A(m-2,11) * A(m+4,9) - A(m-3,11) * A(m+5,8) - A(m-4,11) * A(m+6,7)&
         - A(m-5,11) * A(m+7,6) - A(m-6,11) * A(m-5,5) - A(m-7,11) * A(m-5,4))
    A(m+3,11) = dummy * (A(m+3,11) - A(m-1,11) * A(m+4,10) - A(m-2,11) * A(m+5,9) - A(m-3,11) * A(m+6,8) - A(m-4,11) * A(m+7,7)&
         - A(m-5,11) * A(m-6,6) - A(m-6,11) * A(m-6,5) - A(m-7,11) * A(m-6,4))
    A(m+4,11) = dummy * (A(m+4,11) - A(m-1,11) * A(m+5,10) - A(m-2,11) * A(m+6,9) - A(m-3,11) * A(m+7,8) - A(m-4,11) * A(m-7,7)&
         - A(m-5,11) * A(m-7,6) - A(m-6,11) * A(m-7,5) - A(m-7,11) * A(m-7,4))
    A(m+5,11) = dummy * (A(m+5,11) - A(m-1,11) * A(m+6,10) - A(m-2,11) * A(m+7,9))
    A(m+6,11) = dummy * (A(m+6,11) - A(m-1,11) * A(m+7,10))
    A(m+7,11) = dummy *  A(m+7,11) 

    A(m  ,12) = A(m  ,12) - A(m-1,12) * A(m+1,11) - A(m-2,12) * A(m+2,10) - A(m-3,12) * A(m+3,9) - A(m-4,12) * A(m+4,8)&
         - A(m-5,12) * A(m+5,7) - A(m-6,12) * A(m+6,6) - A(m-7,12) * A(m+7,5)
    A(m-1,13) = A(m-1,13) - A(m-2,13) * A(m+1,11) - A(m-3,13) * A(m+2,10) - A(m-4,13) * A(m+3,9) - A(m-5,13) * A(m+4,8)&
         - A(m-6,13) * A(m+5,7) - A(m-7,13) * A(m+6,6)
    A(m-2,14) = A(m-2,14) - A(m-3,14) * A(m+1,11) - A(m-4,14) * A(m+2,10) - A(m-5,14) * A(m+3,9) - A(m-6,14) * A(m+4,8)&
         - A(m-7,14) * A(m+5,7)
    A(m-3,15) = A(m-3,15) - A(m-4,15) * A(m+1,11) - A(m-5,15) * A(m+2,10) - A(m-6,15) * A(m+3,9) - A(m-7,15) * A(m+4,8)
    A(m-4,16) = A(m-4,16) - A(m-5,16) * A(m+1,11) - A(m-6,16) * A(m+2,10) - A(m-7,16) * A(m+3,9)
    A(m-5,17) = A(m-5,17) - A(m-6,17) * A(m+1,11) - A(m-7,17) * A(m+2,10) 
    A(m-6,18) = A(m-6,18) - A(m-7,18) * A(m+1,11) 

    dummy = 1.d0 / A(m,12)
    A(m+1,12) = dummy * (A(m+1,12) - A(m-1,12) * A(m+2,11) - A(m-2,12) * A(m+3,10) - A(m-3,12) * A(m+4,9) - A(m-4,12) * A(m+5,8)&
         - A(m-5,12) * A(m+6,7) - A(m-6,12) * A(m+7,6) - A(m-7,12) * A(m-5,5))
    A(m+2,12) = dummy * (A(m+2,12) - A(m-1,12) * A(m+3,11) - A(m-2,12) * A(m+4,10) - A(m-3,12) * A(m+5,9) - A(m-4,12) * A(m+6,8)&
         - A(m-5,12) * A(m+7,7) - A(m-6,12) * A(m-6,6) - A(m-7,12) * A(m-6,5))
    A(m+3,12) = dummy * (A(m+3,12) - A(m-1,12) * A(m+4,11) - A(m-2,12) * A(m+5,10) - A(m-3,12) * A(m+6,9) - A(m-4,12) * A(m+7,8)&
         - A(m-5,12) * A(m-7,7) - A(m-6,12) * A(m-7,6) - A(m-7,12) * A(m-7,5))
    A(m+4,12) = dummy * (A(m+4,12) - A(m-1,12) * A(m+5,11) - A(m-2,12) * A(m+6,10) - A(m-3,12) * A(m+7,9))
    A(m+5,12) = dummy * (A(m+5,12) - A(m-1,12) * A(m+6,11) - A(m-2,12) * A(m+7,10))
    A(m+6,12) = dummy * (A(m+6,12) - A(m-1,12) * A(m+7,11))
    A(m+7,12) = dummy *  A(m+7,12) 

    A(m  ,13) = A(m  ,13) - A(m-1,13) * A(m+1,12) - A(m-2,13) * A(m+2,11) - A(m-3,13) * A(m+3,10) - A(m-4,13) * A(m+4,9)&
         - A(m-5,13) * A(m+5,8) - A(m-6,13) * A(m+6,7) - A(m-7,13) * A(m+7,6)
    A(m-1,14) = A(m-1,14) - A(m-2,14) * A(m+1,12) - A(m-3,14) * A(m+2,11) - A(m-4,14) * A(m+3,10) - A(m-5,14) * A(m+4,9)&
         - A(m-6,14) * A(m+5,8) - A(m-7,14) * A(m+6,7)
    A(m-2,15) = A(m-2,15) - A(m-3,15) * A(m+1,12) - A(m-4,15) * A(m+2,11) - A(m-5,15) * A(m+3,10) - A(m-6,15) * A(m+4,9)&
         - A(m-7,15) * A(m+5,8)
    A(m-3,16) = A(m-3,16) - A(m-4,16) * A(m+1,12) - A(m-5,16) * A(m+2,11) - A(m-6,16) * A(m+3,10) - A(m-7,16) * A(m+4,9)
    A(m-4,17) = A(m-4,17) - A(m-5,17) * A(m+1,12) - A(m-6,17) * A(m+2,11) - A(m-7,17) * A(m+3,10)
    A(m-5,18) = A(m-5,18) - A(m-6,18) * A(m+1,12) - A(m-7,18) * A(m+2,11) 
    A(m-6,19) = A(m-6,19) - A(m-7,19) * A(m+1,12) 

    dummy = 1.d0 / A(m,13)
    A(m+1,13) = dummy * (A(m+1,13) - A(m-1,13) * A(m+2,12) - A(m-2,13) * A(m+3,11) - A(m-3,13) * A(m+4,10) - A(m-4,13) * A(m+5,9)&
         - A(m-5,13) * A(m+6,8) - A(m-6,13) * A(m+7,7) - A(m-7,13) * A(m-6,6))
    A(m+2,13) = dummy * (A(m+2,13) - A(m-1,13) * A(m+3,12) - A(m-2,13) * A(m+4,11) - A(m-3,13) * A(m+5,10) - A(m-4,13) * A(m+6,9)&
         - A(m-5,13) * A(m+7,8) - A(m-6,13) * A(m-7,7) - A(m-7,13) * A(m-7,6))
    A(m+3,13) = dummy * (A(m+3,13) - A(m-1,13) * A(m+4,12) - A(m-2,13) * A(m+5,11) - A(m-3,13) * A(m+6,10) - A(m-4,13) * A(m+7,9))
    A(m+4,13) = dummy * (A(m+4,13) - A(m-1,13) * A(m+5,12) - A(m-2,13) * A(m+6,11) - A(m-3,13) * A(m+7,10))
    A(m+5,13) = dummy * (A(m+5,13) - A(m-1,13) * A(m+6,12) - A(m-2,13) * A(m+7,11))
    A(m+6,13) = dummy * (A(m+6,13) - A(m-1,13) * A(m+7,12))
    A(m+7,13) = dummy *  A(m+7,13)

    A(m  ,14) = A(m  ,14) - A(m-1,14) * A(m+1,13) - A(m-2,14) * A(m+2,12) - A(m-3,14) * A(m+3,11) - A(m-4,14) * A(m+4,10)&
         - A(m-5,14) * A(m+5,9) - A(m-6,14) * A(m+6,8) - A(m-7,14) * A(m+7,7)
    A(m-1,15) = A(m-1,15) - A(m-2,15) * A(m+1,13) - A(m-3,15) * A(m+2,12) - A(m-4,15) * A(m+3,11) - A(m-5,15) * A(m+4,10)&
         - A(m-6,15) * A(m+5,9) - A(m-7,15) * A(m+6,8)
    A(m-2,16) = A(m-2,16) - A(m-3,16) * A(m+1,13) - A(m-4,16) * A(m+2,12) - A(m-5,16) * A(m+3,11) - A(m-6,16) * A(m+4,10)&
         - A(m-7,16) * A(m+5,9)
    A(m-3,17) = A(m-3,17) - A(m-4,17) * A(m+1,13) - A(m-5,17) * A(m+2,12) - A(m-6,17) * A(m+3,11) - A(m-7,17) * A(m+4,10)
    A(m-4,18) = A(m-4,18) - A(m-5,18) * A(m+1,13) - A(m-6,18) * A(m+2,12) - A(m-7,18) * A(m+3,11)
    A(m-5,19) = A(m-5,19) - A(m-6,19) * A(m+1,13) - A(m-7,19) * A(m+2,12) 
    A(m-6,20) = A(m-6,20) - A(m-7,20) * A(m+1,13) 

    dummy = 1.d0 / A(m,14)
    A(m+1,14) = dummy * (A(m+1,14) - A(m-1,14) * A(m+2,13) - A(m-2,14) * A(m+3,12) - A(m-3,14) * A(m+4,11) - A(m-4,14) * A(m+5,10)&
         - A(m-5,14) * A(m+6,9) - A(m-6,14) * A(m+7,8) - A(m-7,14) * A(m-7,7))
    A(m+2,14) = dummy * (A(m+2,14) - A(m-1,14) * A(m+3,13) - A(m-2,14) * A(m+4,12) - A(m-3,14) * A(m+5,11) - A(m-4,14) * A(m+6,10)&
         - A(m-5,14) * A(m+7,9))
    A(m+3,14) = dummy * (A(m+3,14) - A(m-1,14) * A(m+4,13) - A(m-2,14) * A(m+5,12) - A(m-3,14) * A(m+6,11) - A(m-4,14) * A(m+7,10))
    A(m+4,14) = dummy * (A(m+4,14) - A(m-1,14) * A(m+5,13) - A(m-2,14) * A(m+6,12) - A(m-3,14) * A(m+7,11))
    A(m+5,14) = dummy * (A(m+5,14) - A(m-1,14) * A(m+6,13) - A(m-2,14) * A(m+7,12))
    A(m+6,14) = dummy * (A(m+6,14) - A(m-1,14) * A(m+7,13))
    A(m+7,14) = dummy *  A(m+7,14)

    !center region

    do i = 2*bw+1,n-2*bw-1
       A(m  ,i  ) = A(m  ,i  ) - A(m-1,i  ) * A(m+1,i-1) - A(m-2,i)   * A(m+2,i-2) - A(m-3,i  ) * A(m+3,i-3)&
            - A(m-4,i)   * A(m+4,i-4) - A(m-5,i)   * A(m+5,i-5) - A(m-6,i)   * A(m+6,i-6) - A(m-7,i) * A(m+7,i-7)
       A(m-1,i+1) = A(m-1,i+1) - A(m-2,i+1) * A(m+1,i-1) - A(m-3,i+1) * A(m+2,i-2) - A(m-4,i+1) * A(m+3,i-3) &
            - A(m-5,i+1) * A(m+4,i-4) - A(m-6,i+1) * A(m+5,i-5) - A(m-7,i+1) * A(m+6,i-6)
       A(m-2,i+2) = A(m-2,i+2) - A(m-3,i+2) * A(m+1,i-1) - A(m-4,i+2) * A(m+2,i-2) - A(m-5,i+2) * A(m+3,i-3) &
            - A(m-6,i+2) * A(m+4,i-4) - A(m-7,i+2) * A(m+5,i-5)
       A(m-3,i+3) = A(m-3,i+3) - A(m-4,i+3) * A(m+1,i-1) - A(m-5,i+3) * A(m+2,i-2) - A(m-6,i+3) * A(m+3,i-3)&
            - A(m-7,i+3) * A(m+4,i-4)
       A(m-4,i+4) = A(m-4,i+4) - A(m-5,i+4) * A(m+1,i-1) - A(m-6,i+4) * A(m+2,i-2) - A(m-7,i+4) * A(m+3,i-3)
       A(m-5,i+5) = A(m-5,i+5) - A(m-6,i+5) * A(m+1,i-1) - A(m-7,i+5) * A(m+2,i-2) 
       A(m-6,i+6) = A(m-6,i+6) - A(m-7,i+6) * A(m+1,i-1) 

       dummy = 1.d0 / A(m,i)
       A(m+1,i) = dummy * (A(m+1,i) - A(m-1,i) * A(m+2,i-1) - A(m-2,i) * A(m+3,i-2) - A(m-3,i) * A(m+4,i-3) &
            - A(m-4,i) * A(m+5,i-4) - A(m-5,i) * A(m+6,i-5) - A(m-6,i) * A(m+7,i-6))
       A(m+2,i) = dummy * (A(m+2,i) - A(m-1,i) * A(m+3,i-1) - A(m-2,i) * A(m+4,i-2) - A(m-3,i) * A(m+5,i-3)&
            - A(m-4,i) * A(m+6,i-4) - A(m-5,i) * A(m+7,i-5))
       A(m+3,i) = dummy * (A(m+3,i) - A(m-1,i) * A(m+4,i-1) - A(m-2,i) * A(m+5,i-2) - A(m-3,i) * A(m+6,i-3)&
            - A(m-4,i) * A(m+7,i-4))
       A(m+4,i) = dummy * (A(m+4,i) - A(m-1,i) * A(m+5,i-1) - A(m-2,i) * A(m+6,i-2) - A(m-3,i) * A(m+7,i-3))
       A(m+5,i) = dummy * (A(m+5,i) - A(m-1,i) * A(m+6,i-1) - A(m-2,i) * A(m+7,i-2))
       A(m+6,i) = dummy * (A(m+6,i) - A(m-1,i) * A(m+7,i-1))
       A(m+7,i) = dummy * (A(m+7,i))
    enddo

    ! sublower region
    A(m  ,n-2*bw)    = A(m  ,n-2*bw)    - A(m-1,n-2*bw)   * A(m+1,n-2*bw-1) - A(m-2,n-2*bw)   * A(m+2,n-2*bw-2) &
         - A(m-3,n-2*bw)   * A(m+3,n-2*bw-3) - A(m-4,n-2*bw)   * A(m+4,n-2*bw-4) - A(m-5,n-2*bw)   * A(m+5,n-2*bw-5) &
         - A(m-6,n-2*bw)   * A(m+6,n-2*bw-6) - A(m-7,n-2*bw) * A(m+7,n-2*bw-7)    
    A(m-1,n-2*bw+1)  = A(m-1,n-2*bw+1)  - A(m-2,n-2*bw+1) * A(m+1,n-2*bw-1) - A(m-3,n-2*bw+1) * A(m+2,n-2*bw-2)&
         - A(m-4,n-2*bw+1) * A(m+3,n-2*bw-3) - A(m-5,n-2*bw+1) * A(m+4,n-2*bw-4) - A(m-6,n-2*bw+1) * A(m+5,n-2*bw-5)&
         - A(m-7,n-2*bw+1) * A(m+6,n-2*bw-6)  
    A(m-2,n-2*bw+2)  = A(m-2,n-2*bw+2)  - A(m-3,n-2*bw+2) * A(m+1,n-2*bw-1) - A(m-4,n-2*bw+2) * A(m+2,n-2*bw-2)&
         - A(m-5,n-2*bw+2) * A(m+3,n-2*bw-3) - A(m-6,n-2*bw+2) * A(m+4,n-2*bw-4) - A(m-7,n-2*bw+2) * A(m+5,n-2*bw-5)    
    A(m-3,n-2*bw+3)  = A(m-3,n-2*bw+3)  - A(m-4,n-2*bw+3) * A(m+1,n-2*bw-1) - A(m-5,n-2*bw+3) * A(m+2,n-2*bw-2)&
         - A(m-6,n-2*bw+3) * A(m+3,n-2*bw-3) - A(m-7,n-2*bw+3) * A(m+4,n-2*bw-4) 
    A(m-4,n-2*bw+4)  = A(m-4,n-2*bw+4)  - A(m-5,n-2*bw+4) * A(m+1,n-2*bw-1) - A(m-6,n-2*bw+4) * A(m+2,n-2*bw-2)&
         - A(m-7,n-2*bw+4) * A(m+3,n-2*bw-3)
    A(m-5,n-2*bw+5)  = A(m-5,n-2*bw+5)  - A(m-6,n-2*bw+5) * A(m+1,n-2*bw-1) - A(m-7,n-2*bw+5) * A(m+2,n-2*bw-2)
    A(m-6,n-2*bw+6)  = A(m-6,n-2*bw+6)  - A(m-7,n-2*bw+6) * A(m+1,n-2*bw-1) 

    dummy = 1.d0 / A(m,n-2*bw)
    A(m+1,n-2*bw) = dummy * (A(m+1,n-2*bw) - A(m-1,n-2*bw) * A(m+2,n-2*bw-1) - A(m-2,n-2*bw) * A(m+3,n-2*bw-2) &
         - A(m-3,n-2*bw) * A(m+4,n-2*bw-3) - A(m-4,n-2*bw) * A(m+5,n-2*bw-4) - A(m-5,n-2*bw) * A(m+6,n-2*bw-5) &
         - A(m-6,n-2*bw) * A(m+7,n-2*bw-6))
    A(m+2,n-2*bw) = dummy * (A(m+2,n-2*bw) - A(m-1,n-2*bw) * A(m+3,n-2*bw-1) - A(m-2,n-2*bw) * A(m+4,n-2*bw-2) &
         - A(m-3,n-2*bw) * A(m+5,n-2*bw-3) - A(m-4,n-2*bw) * A(m+6,n-2*bw-4) - A(m-5,n-2*bw) * A(m+7,n-2*bw-5))
    A(m+3,n-2*bw) = dummy * (A(m+3,n-2*bw) - A(m-1,n-2*bw) * A(m+4,n-2*bw-1) - A(m-2,n-2*bw) * A(m+5,n-2*bw-2)&
         - A(m-3,n-2*bw) * A(m+6,n-2*bw-3) - A(m-4,n-2*bw) * A(m+7,n-2*bw-4))
    A(m+4,n-2*bw) = dummy * (A(m+4,n-2*bw) - A(m-1,n-2*bw) * A(m+5,n-2*bw-1) - A(m-2,n-2*bw) * A(m+6,n-2*bw-2)&
         - A(m-3,n-2*bw) * A(m+7,n-2*bw-3))
    A(m+5,n-2*bw) = dummy * (A(m+5,n-2*bw) - A(m-1,n-2*bw) * A(m+6,n-2*bw-1) - A(m-2,n-2*bw) * A(m+7,n-2*bw-2))
    A(m+6,n-2*bw) = dummy * (A(m+6,n-2*bw) - A(m-1,n-2*bw) * A(m+7,n-2*bw-1))
    A(m+7,n-2*bw) = dummy *  A(m+7,n-2*bw)

    A(m  ,n-2*bw+1)   = A(m  ,n-2*bw+1)  - A(m-1,n-2*bw+1) * A(m+1,n-2*bw) - A(m-2,n-2*bw+1) * A(m+2,n-2*bw-1) &
         - A(m-3,n-2*bw+1) * A(m+3,n-2*bw-2) - A(m-4,n-2*bw+1) * A(m+4,n-2*bw-3) - A(m-5,n-2*bw+1) * A(m+5,n-2*bw-4)&
         - A(m-6,n-2*bw+1) * A(m+6,n-2*bw-5) - A(m-7,n-2*bw+1) * A(m+7,n-2*bw-6)    
    A(m-1,n-2*bw+2)   = A(m-1,n-2*bw+2)  - A(m-2,n-2*bw+2) * A(m+1,n-2*bw) - A(m-3,n-2*bw+2) * A(m+2,n-2*bw-1) &
         - A(m-4,n-2*bw+2) * A(m+3,n-2*bw-2) - A(m-5,n-2*bw+2) * A(m+4,n-2*bw-3) - A(m-6,n-2*bw+2) * A(m+5,n-2*bw-4)&
         - A(m-7,n-2*bw+2) * A(m+6,n-2*bw-5)  
    A(m-2,n-2*bw+3)   = A(m-2,n-2*bw+3)  - A(m-3,n-2*bw+3) * A(m+1,n-2*bw) - A(m-4,n-2*bw+3) * A(m+2,n-2*bw-1) &
         - A(m-5,n-2*bw+3) * A(m+3,n-2*bw-2) - A(m-6,n-2*bw+3) * A(m+4,n-2*bw-3) - A(m-7,n-2*bw+3) * A(m+5,n-2*bw-4)    
    A(m-3,n-2*bw+4)   = A(m-3,n-2*bw+4)  - A(m-4,n-2*bw+4) * A(m+1,n-2*bw) - A(m-5,n-2*bw+4) * A(m+2,n-2*bw-1) &
         - A(m-6,n-2*bw+4) * A(m+3,n-2*bw-2) - A(m-7,n-2*bw+4) * A(m+4,n-2*bw-3) 
    A(m-4,n-2*bw+5)   = A(m-4,n-2*bw+5)  - A(m-5,n-2*bw+5) * A(m+1,n-2*bw) - A(m-6,n-2*bw+5) * A(m+2,n-2*bw-1)&
         - A(m-7,n-2*bw+5) * A(m+3,n-2*bw-2)
    A(m-5,n-2*bw+6)   = A(m-5,n-2*bw+6)  - A(m-6,n-2*bw+6) * A(m+1,n-2*bw) - A(m-7,n-2*bw+6) * A(m+2,n-2*bw-1)
    A(m-6,n-2*bw+7)   = A(m-6,n-2*bw+7)  - A(m-7,n-2*bw+7) * A(m+1,n-2*bw) 
    A(m-7,n-2*bw+8)   = A(m-7,n-2*bw+8)  - A(m+7,n-2*bw+8) * A(m+1,n-2*bw) 
    A(m+6,n-2*bw+9)   = A(m+6,n-2*bw+9)  - A(m+7,n-2*bw+9) * A(m+1,n-2*bw) 
    A(m+6,n-2*bw+10)  = A(m+6,n-2*bw+10)  - A(m+7,n-2*bw+10) * A(m+1,n-2*bw) 
    A(m+6,n-2*bw+11)  = A(m+6,n-2*bw+11)  - A(m+7,n-2*bw+11) * A(m+1,n-2*bw) 
    A(m+6,n-2*bw+12)  = A(m+6,n-2*bw+12)  - A(m+7,n-2*bw+12) * A(m+1,n-2*bw) 
    A(m+6,n-2*bw+13)  = A(m+6,n-2*bw+13)  - A(m+7,n-2*bw+13) * A(m+1,n-2*bw) 
    A(m+6,n-2*bw+14)  = A(m+6,n-2*bw+14)  - A(m+7,n-2*bw+14) * A(m+1,n-2*bw) 

    dummy = 1.d0 / A(m,n-2*bw+1)
    A(m+1,n-2*bw+1) = dummy * (A(m+1,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+2,n-2*bw) - A(m-2,n-2*bw+1) * A(m+3,n-2*bw-1)&
         - A(m-3,n-2*bw+1) * A(m+4,n-2*bw-2) - A(m-4,n-2*bw+1) * A(m+5,n-2*bw-3) - A(m-5,n-2*bw+1) * A(m+6,n-2*bw-4)&
         - A(m-6,n-2*bw+1) * A(m+7,n-2*bw-5))
    A(m+2,n-2*bw+1) = dummy * (A(m+2,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+3,n-2*bw) - A(m-2,n-2*bw+1) * A(m+4,n-2*bw-1)&
         - A(m-3,n-2*bw+1) * A(m+5,n-2*bw-2) - A(m-4,n-2*bw+1) * A(m+6,n-2*bw-3) - A(m-5,n-2*bw+1) * A(m+7,n-2*bw-4))
    A(m+3,n-2*bw+1) = dummy * (A(m+3,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+4,n-2*bw) - A(m-2,n-2*bw+1) * A(m+5,n-2*bw-1)&
         - A(m-3,n-2*bw+1) * A(m+6,n-2*bw-2) - A(m-4,n-2*bw+1) * A(m+7,n-2*bw-3))
    A(m+4,n-2*bw+1) = dummy * (A(m+4,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+5,n-2*bw) - A(m-2,n-2*bw+1) * A(m+6,n-2*bw-1)&
         - A(m-3,n-2*bw+1) * A(m+7,n-2*bw-2))
    A(m+5,n-2*bw+1) = dummy * (A(m+5,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+6,n-2*bw) - A(m-2,n-2*bw+1) * A(m+7,n-2*bw-1))
    A(m+6,n-2*bw+1) = dummy * (A(m+6,n-2*bw+1) - A(m-1,n-2*bw+1) * A(m+7,n-2*bw))
    A(m+7,n-2*bw+1) = dummy *  A(m+7,n-2*bw+1)

    A(m  ,n-2*bw+2)   = A(m  ,n-2*bw+2)  - A(m-1,n-2*bw+2) * A(m+1,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+2,n-2*bw) &
         - A(m-3,n-2*bw+2) * A(m+3,n-2*bw-1) - A(m-4,n-2*bw+2) * A(m+4,n-2*bw-2) - A(m-5,n-2*bw+2) * A(m+5,n-2*bw-3)&
         - A(m-6,n-2*bw+2) * A(m+6,n-2*bw-4) - A(m-7,n-2*bw+2) * A(m+7,n-2*bw-5)    
    A(m-1,n-2*bw+3)   = A(m-1,n-2*bw+3)  - A(m-2,n-2*bw+3) * A(m+1,n-2*bw+1) - A(m-3,n-2*bw+3) * A(m+2,n-2*bw) &
         - A(m-4,n-2*bw+3) * A(m+3,n-2*bw-1) - A(m-5,n-2*bw+3) * A(m+4,n-2*bw-2) - A(m-6,n-2*bw+3) * A(m+5,n-2*bw-3)&
         - A(m-7,n-2*bw+3) * A(m+6,n-2*bw-4)  
    A(m-2,n-2*bw+4)   = A(m-2,n-2*bw+4)  - A(m-3,n-2*bw+4) * A(m+1,n-2*bw+1) - A(m-4,n-2*bw+4) * A(m+2,n-2*bw)&
         - A(m-5,n-2*bw+4) * A(m+3,n-2*bw-1) - A(m-6,n-2*bw+4) * A(m+4,n-2*bw-2) - A(m-7,n-2*bw+4) * A(m+5,n-2*bw-3)    
    A(m-3,n-2*bw+5)   = A(m-3,n-2*bw+5)  - A(m-4,n-2*bw+5) * A(m+1,n-2*bw+1) - A(m-5,n-2*bw+5) * A(m+2,n-2*bw)&
         - A(m-6,n-2*bw+5) * A(m+3,n-2*bw-1) - A(m-7,n-2*bw+5) * A(m+4,n-2*bw-2) 
    A(m-4,n-2*bw+6)   = A(m-4,n-2*bw+6)  - A(m-5,n-2*bw+6) * A(m+1,n-2*bw+1) - A(m-6,n-2*bw+6) * A(m+2,n-2*bw) &
         - A(m-7,n-2*bw+6) * A(m+3,n-2*bw-1)
    A(m-5,n-2*bw+7)   = A(m-5,n-2*bw+7)  - A(m-6,n-2*bw+7) * A(m+1,n-2*bw+1) - A(m-7,n-2*bw+7) * A(m+2,n-2*bw)
    A(m-6,n-2*bw+8)   = A(m-6,n-2*bw+8)  - A(m-7,n-2*bw+8) * A(m+1,n-2*bw+1) - A(m+7,n-2*bw+8) * A(m+2,n-2*bw)
    A(m-7,n-2*bw+9)   = A(m-7,n-2*bw+9)  - A(m+6,n-2*bw+9) * A(m+1,n-2*bw+1) - A(m+7,n-2*bw+9) * A(m+2,n-2*bw)
    A(m+5,n-2*bw+10)   = A(m+5,n-2*bw+10)  - A(m+6,n-2*bw+10) * A(m+1,n-2*bw+1) - A(m+7,n-2*bw+10) * A(m+2,n-2*bw)
    A(m+5,n-2*bw+11)   = A(m+5,n-2*bw+11)  - A(m+6,n-2*bw+11) * A(m+1,n-2*bw+1) - A(m+7,n-2*bw+11) * A(m+2,n-2*bw)
    A(m+5,n-2*bw+12)   = A(m+5,n-2*bw+12)  - A(m+6,n-2*bw+12) * A(m+1,n-2*bw+1) - A(m+7,n-2*bw+12) * A(m+2,n-2*bw)
    A(m+5,n-2*bw+13)   = A(m+5,n-2*bw+13)  - A(m+6,n-2*bw+13) * A(m+1,n-2*bw+1) - A(m+7,n-2*bw+13) * A(m+2,n-2*bw)
    A(m+5,n-2*bw+14)   = A(m+5,n-2*bw+14)  - A(m+6,n-2*bw+14) * A(m+1,n-2*bw+1) - A(m+7,n-2*bw+14) * A(m+2,n-2*bw)

    dummy = 1.d0 / A(m,n-2*bw+2)
    A(m+1,n-2*bw+2) = dummy * (A(m+1,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+2,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+3,n-2*bw)&
         - A(m-3,n-2*bw+2) * A(m+4,n-2*bw-1) - A(m-4,n-2*bw+2) * A(m+5,n-2*bw-2) - A(m-5,n-2*bw+2) * A(m+6,n-2*bw-3) &
         - A(m-6,n-2*bw+2) * A(m+7,n-2*bw-4))
    A(m+2,n-2*bw+2) = dummy * (A(m+2,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+3,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+4,n-2*bw) &
         - A(m-3,n-2*bw+2) * A(m+5,n-2*bw-1) - A(m-4,n-2*bw+2) * A(m+6,n-2*bw-2) - A(m-5,n-2*bw+2) * A(m+7,n-2*bw-3))
    A(m+3,n-2*bw+2) = dummy * (A(m+3,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+4,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+5,n-2*bw)&
         - A(m-3,n-2*bw+2) * A(m+6,n-2*bw-1) - A(m-4,n-2*bw+2) * A(m+7,n-2*bw-2))
    A(m+4,n-2*bw+2) = dummy * (A(m+4,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+5,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+6,n-2*bw)&
         - A(m-3,n-2*bw+2) * A(m+7,n-2*bw-1))
    A(m+5,n-2*bw+2) = dummy * (A(m+5,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+6,n-2*bw+1) - A(m-2,n-2*bw+2) * A(m+7,n-2*bw))
    A(m+6,n-2*bw+2) = dummy * (A(m+6,n-2*bw+2) - A(m-1,n-2*bw+2) * A(m+7,n-2*bw+1))
    A(m+7,n-2*bw+2) = dummy *  A(m+7,n-2*bw+2)

    A(m  ,n-2*bw+3)   = A(m  ,n-2*bw+3)  - A(m-1,n-2*bw+3) * A(m+1,n-2*bw+2) - A(m-2,n-2*bw+3) * A(m+2,n-2*bw+1) &
         - A(m-3,n-2*bw+3) * A(m+3,n-2*bw) - A(m-4,n-2*bw+3) * A(m+4,n-2*bw-1) - A(m-5,n-2*bw+3) * A(m+5,n-2*bw-2) &
         - A(m-6,n-2*bw+3) * A(m+6,n-2*bw-3) - A(m-7,n-2*bw+3) * A(m+7,n-2*bw-4)    
    A(m-1,n-2*bw+4)   = A(m-1,n-2*bw+4)  - A(m-2,n-2*bw+4) * A(m+1,n-2*bw+2) - A(m-3,n-2*bw+4) * A(m+2,n-2*bw+1)&
         - A(m-4,n-2*bw+4) * A(m+3,n-2*bw) - A(m-5,n-2*bw+4) * A(m+4,n-2*bw-1) - A(m-6,n-2*bw+4) * A(m+5,n-2*bw-2)&
         - A(m-7,n-2*bw+4) * A(m+6,n-2*bw-3)  
    A(m-2,n-2*bw+5)   = A(m-2,n-2*bw+5)  - A(m-3,n-2*bw+5) * A(m+1,n-2*bw+2) - A(m-4,n-2*bw+5) * A(m+2,n-2*bw+1)&
         - A(m-5,n-2*bw+5) * A(m+3,n-2*bw) - A(m-6,n-2*bw+5) * A(m+4,n-2*bw-1) - A(m-7,n-2*bw+5) * A(m+5,n-2*bw-2)    
    A(m-3,n-2*bw+6)   = A(m-3,n-2*bw+6)  - A(m-4,n-2*bw+6) * A(m+1,n-2*bw+2) - A(m-5,n-2*bw+6) * A(m+2,n-2*bw+1)&
         - A(m-6,n-2*bw+6) * A(m+3,n-2*bw) - A(m-7,n-2*bw+6) * A(m+4,n-2*bw-1) 
    A(m-4,n-2*bw+7)   = A(m-4,n-2*bw+7)  - A(m-5,n-2*bw+7) * A(m+1,n-2*bw+2) - A(m-6,n-2*bw+7) * A(m+2,n-2*bw+1) &
         - A(m-7,n-2*bw+7) * A(m+3,n-2*bw)
    A(m-5,n-2*bw+8)    = A(m-5,n-2*bw+8)   - A(m-6,n-2*bw+8)  * A(m+1,n-2*bw+2) - A(m-7,n-2*bw+8)  * A(m+2,n-2*bw+1)&
         - A(m+7,n-2*bw+8) * A(m+3,n-2*bw)
    A(m-6,n-2*bw+9)    = A(m-6,n-2*bw+9)   - A(m-7,n-2*bw+9)  * A(m+1,n-2*bw+2) - A(m+6,n-2*bw+9)  * A(m+2,n-2*bw+1)&
         - A(m+7,n-2*bw+9) * A(m+3,n-2*bw)
    A(m-7,n-2*bw+10)   = A(m-7,n-2*bw+10)  - A(m+5,n-2*bw+10) * A(m+1,n-2*bw+2) - A(m+6,n-2*bw+10) * A(m+2,n-2*bw+1)&
         - A(m+7,n-2*bw+10) * A(m+3,n-2*bw)
    A(m+4,n-2*bw+11)   = A(m+4,n-2*bw+11)  - A(m+5,n-2*bw+11) * A(m+1,n-2*bw+2) - A(m+6,n-2*bw+11) * A(m+2,n-2*bw+1)&
         - A(m+7,n-2*bw+11) * A(m+3,n-2*bw)
    A(m+4,n-2*bw+12)   = A(m+4,n-2*bw+12)  - A(m+5,n-2*bw+12) * A(m+1,n-2*bw+2) - A(m+6,n-2*bw+12) * A(m+2,n-2*bw+1)&
         - A(m+7,n-2*bw+12) * A(m+3,n-2*bw)
    A(m+4,n-2*bw+13)   = A(m+4,n-2*bw+13)  - A(m+5,n-2*bw+13) * A(m+1,n-2*bw+2) - A(m+6,n-2*bw+13) * A(m+2,n-2*bw+1)&
         - A(m+7,n-2*bw+13) * A(m+3,n-2*bw)
    A(m+4,n-2*bw+14)   = A(m+4,n-2*bw+14)  - A(m+5,n-2*bw+14) * A(m+1,n-2*bw+2) - A(m+6,n-2*bw+14) * A(m+2,n-2*bw+1)&
         - A(m+7,n-2*bw+14) * A(m+3,n-2*bw)

    dummy = 1.d0 / A(m,n-2*bw+3)
    A(m+1,n-2*bw+3) = dummy * (A(m+1,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+2,n-2*bw+2) - A(m-2,n-2*bw+3) * A(m+3,n-2*bw+1)&
         - A(m-3,n-2*bw+3) * A(m+4,n-2*bw) - A(m-4,n-2*bw+3) * A(m+5,n-2*bw-1) - A(m-5,n-2*bw+3) * A(m+6,n-2*bw-2) &
         - A(m-6,n-2*bw+3) * A(m+7,n-2*bw-3))
    A(m+2,n-2*bw+3) = dummy * (A(m+2,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+3,n-2*bw+2) - A(m-2,n-2*bw+3) * A(m+4,n-2*bw+1) &
         - A(m-3,n-2*bw+3) * A(m+5,n-2*bw) - A(m-4,n-2*bw+3) * A(m+6,n-2*bw-1) - A(m-5,n-2*bw+3) * A(m+7,n-2*bw-2))
    A(m+3,n-2*bw+3) = dummy * (A(m+3,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+4,n-2*bw+2) - A(m-2,n-2*bw+3) * A(m+5,n-2*bw+1)&
         - A(m-3,n-2*bw+3) * A(m+6,n-2*bw) - A(m-4,n-2*bw+3) * A(m+7,n-2*bw-1))
    A(m+4,n-2*bw+3) = dummy * (A(m+4,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+5,n-2*bw+2) - A(m-2,n-2*bw+3) * A(m+6,n-2*bw+1) &
         - A(m-3,n-2*bw+3) * A(m+7,n-2*bw))
    A(m+5,n-2*bw+3) = dummy * (A(m+5,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+6,n-2*bw+2) - A(m-2,n-2*bw+3) * A(m+7,n-2*bw+1))
    A(m+6,n-2*bw+3) = dummy * (A(m+6,n-2*bw+3) - A(m-1,n-2*bw+3) * A(m+7,n-2*bw+2))
    A(m+7,n-2*bw+3) = dummy *  A(m+7,n-2*bw+3)

    A(m  ,n-2*bw+4)   = A(m  ,n-2*bw+4)  - A(m-1,n-2*bw+4) * A(m+1,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+2,n-2*bw+2) &
         - A(m-3,n-2*bw+4) * A(m+3,n-2*bw+1) - A(m-4,n-2*bw+4) * A(m+4,n-2*bw) - A(m-5,n-2*bw+4) * A(m+5,n-2*bw-1) &
         - A(m-6,n-2*bw+4) * A(m+6,n-2*bw-2) - A(m-7,n-2*bw+4) * A(m+7,n-2*bw-3)    
    A(m-1,n-2*bw+5)   = A(m-1,n-2*bw+5)  - A(m-2,n-2*bw+5) * A(m+1,n-2*bw+3) - A(m-3,n-2*bw+5) * A(m+2,n-2*bw+2) &
         - A(m-4,n-2*bw+5) * A(m+3,n-2*bw+1) - A(m-5,n-2*bw+5) * A(m+4,n-2*bw) - A(m-6,n-2*bw+5) * A(m+5,n-2*bw-1) &
         - A(m-7,n-2*bw+5) * A(m+6,n-2*bw-2)  
    A(m-2,n-2*bw+6)   = A(m-2,n-2*bw+6)  - A(m-3,n-2*bw+6) * A(m+1,n-2*bw+3) - A(m-4,n-2*bw+6) * A(m+2,n-2*bw+2) &
         - A(m-5,n-2*bw+6) * A(m+3,n-2*bw+1) - A(m-6,n-2*bw+6) * A(m+4,n-2*bw) - A(m-7,n-2*bw+6) * A(m+5,n-2*bw-1)    
    A(m-3,n-2*bw+7)   = A(m-3,n-2*bw+7)  - A(m-4,n-2*bw+7) * A(m+1,n-2*bw+3) - A(m-5,n-2*bw+7) * A(m+2,n-2*bw+2) &
         - A(m-6,n-2*bw+7) * A(m+3,n-2*bw+1) - A(m-7,n-2*bw+7) * A(m+4,n-2*bw) 
    A(m-4,n-2*bw+8)   = A(m-4,n-2*bw+8)  - A(m-5,n-2*bw+8) * A(m+1,n-2*bw+3) - A(m-6,n-2*bw+8) * A(m+2,n-2*bw+2) &
         - A(m-7,n-2*bw+8) * A(m+3,n-2*bw+1) - A(m+7,n-2*bw+8) * A(m+4,n-2*bw)
    A(m-5,n-2*bw+9)    = A(m-5,n-2*bw+9)   - A(m-6,n-2*bw+9)  * A(m+1,n-2*bw+3) - A(m-7,n-2*bw+9)  * A(m+2,n-2*bw+2)&
         - A(m+6,n-2*bw+9)  * A(m+3,n-2*bw+1) - A(m+7,n-2*bw+9) * A(m+4,n-2*bw)
    A(m-6,n-2*bw+10)   = A(m-6,n-2*bw+10)  - A(m-7,n-2*bw+10) * A(m+1,n-2*bw+3) - A(m+5,n-2*bw+10) * A(m+2,n-2*bw+2)&
         - A(m+6,n-2*bw+10) * A(m+3,n-2*bw+1)- A(m+7,n-2*bw+10) * A(m+4,n-2*bw)
    A(m-7,n-2*bw+11)   = A(m-7,n-2*bw+11)  - A(m+4,n-2*bw+11) * A(m+1,n-2*bw+3) - A(m+5,n-2*bw+11) * A(m+2,n-2*bw+2)&
         - A(m+6,n-2*bw+11) * A(m+3,n-2*bw+1)- A(m+7,n-2*bw+11) * A(m+4,n-2*bw)
    A(m+3,n-2*bw+12)   = A(m+3,n-2*bw+12)  - A(m+4,n-2*bw+12) * A(m+1,n-2*bw+3) - A(m+5,n-2*bw+12) * A(m+2,n-2*bw+2)&
         - A(m+6,n-2*bw+12) * A(m+3,n-2*bw+1)- A(m+7,n-2*bw+12) * A(m+4,n-2*bw)
    A(m+3,n-2*bw+13)   = A(m+3,n-2*bw+13)  - A(m+4,n-2*bw+13) * A(m+1,n-2*bw+3) - A(m+5,n-2*bw+13) * A(m+2,n-2*bw+2)&
         - A(m+6,n-2*bw+13) * A(m+3,n-2*bw+1)- A(m+7,n-2*bw+13) * A(m+4,n-2*bw)
    A(m+3,n-2*bw+14)   = A(m+3,n-2*bw+14)  - A(m+4,n-2*bw+14) * A(m+1,n-2*bw+3) - A(m+5,n-2*bw+14) * A(m+2,n-2*bw+2)&
         - A(m+6,n-2*bw+14) * A(m+3,n-2*bw+1)- A(m+7,n-2*bw+14) * A(m+4,n-2*bw)


    dummy = 1.d0 / A(m,n-2*bw+4)
    A(m+1,n-2*bw+4) = dummy * (A(m+1,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+2,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+3,n-2*bw+2)&
         - A(m-3,n-2*bw+4) * A(m+4,n-2*bw+1) - A(m-4,n-2*bw+4) * A(m+5,n-2*bw) - A(m-5,n-2*bw+4) * A(m+6,n-2*bw-1)&
         - A(m-6,n-2*bw+4) * A(m+7,n-2*bw-2))
    A(m+2,n-2*bw+4) = dummy * (A(m+2,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+3,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+4,n-2*bw+2)&
         - A(m-3,n-2*bw+4) * A(m+5,n-2*bw+1) - A(m-4,n-2*bw+4) * A(m+6,n-2*bw) - A(m-5,n-2*bw+4) * A(m+7,n-2*bw-1))
    A(m+3,n-2*bw+4) = dummy * (A(m+3,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+4,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+5,n-2*bw+2)&
         - A(m-3,n-2*bw+4) * A(m+6,n-2*bw+1) - A(m-4,n-2*bw+4) * A(m+7,n-2*bw))
    A(m+4,n-2*bw+4) = dummy * (A(m+4,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+5,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+6,n-2*bw+2) &
         - A(m-3,n-2*bw+4) * A(m+7,n-2*bw+1))
    A(m+5,n-2*bw+4) = dummy * (A(m+5,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+6,n-2*bw+3) - A(m-2,n-2*bw+4) * A(m+7,n-2*bw+2))
    A(m+6,n-2*bw+4) = dummy * (A(m+6,n-2*bw+4) - A(m-1,n-2*bw+4) * A(m+7,n-2*bw+3))
    A(m+7,n-2*bw+4) = dummy *  A(m+7,n-2*bw+4)

    A(m  ,n-2*bw+5)   = A(m  ,n-2*bw+5)  - A(m-1,n-2*bw+5)  * A(m+1,n-2*bw+4) - A(m-2,n-2*bw+5)  * A(m+2,n-2*bw+3) &
         - A(m-3,n-2*bw+5)  * A(m+3,n-2*bw+2) - A(m-4,n-2*bw+5)  * A(m+4,n-2*bw+1) - A(m-5,n-2*bw+5)  * A(m+5,n-2*bw)&
         - A(m-6,n-2*bw+5) * A(m+6,n-2*bw-1) - A(m-7,n-2*bw+5) * A(m+7,n-2*bw-2)    
    A(m-1,n-2*bw+6)   = A(m-1,n-2*bw+6)  - A(m-2,n-2*bw+6)  * A(m+1,n-2*bw+4) - A(m-3,n-2*bw+6)  * A(m+2,n-2*bw+3) &
         - A(m-4,n-2*bw+6)  * A(m+3,n-2*bw+2) - A(m-5,n-2*bw+6)  * A(m+4,n-2*bw+1) - A(m-6,n-2*bw+6)  * A(m+5,n-2*bw)&
         - A(m-7,n-2*bw+6) * A(m+6,n-2*bw-1)  
    A(m-2,n-2*bw+7)   = A(m-2,n-2*bw+7)  - A(m-3,n-2*bw+7)  * A(m+1,n-2*bw+4) - A(m-4,n-2*bw+7)  * A(m+2,n-2*bw+3)&
         - A(m-5,n-2*bw+7)  * A(m+3,n-2*bw+2) - A(m-6,n-2*bw+7)  * A(m+4,n-2*bw+1) - A(m-7,n-2*bw+7)  * A(m+5,n-2*bw)    
    A(m-3,n-2*bw+8)   = A(m-3,n-2*bw+8)  - A(m-4,n-2*bw+8)  * A(m+1,n-2*bw+4) - A(m-5,n-2*bw+8)  * A(m+2,n-2*bw+3) &
         - A(m-6,n-2*bw+8)  * A(m+3,n-2*bw+2) - A(m-7,n-2*bw+8)  * A(m+4,n-2*bw+1) - A(m+7,n-2*bw+8)  * A(m+5,n-2*bw)    
    A(m-4,n-2*bw+9)   = A(m-4,n-2*bw+9)  - A(m-5,n-2*bw+9)  * A(m+1,n-2*bw+4) - A(m-6,n-2*bw+9)  * A(m+2,n-2*bw+3) &
         - A(m-7,n-2*bw+9)  * A(m+3,n-2*bw+2) - A(m+6,n-2*bw+9)  * A(m+4,n-2*bw+1) - A(m+7,n-2*bw+9)  * A(m+5,n-2*bw)    
    A(m-5,n-2*bw+10)  = A(m-5,n-2*bw+10) - A(m-6,n-2*bw+10) * A(m+1,n-2*bw+4) - A(m-7,n-2*bw+10) * A(m+2,n-2*bw+3) &
         - A(m+5,n-2*bw+10) * A(m+3,n-2*bw+2) - A(m+6,n-2*bw+10) * A(m+4,n-2*bw+1) - A(m+7,n-2*bw+10) * A(m+5,n-2*bw)    
    A(m-6,n-2*bw+11)  = A(m-6,n-2*bw+11) - A(m-7,n-2*bw+11) * A(m+1,n-2*bw+4) - A(m+4,n-2*bw+11) * A(m+2,n-2*bw+3) &
         - A(m+5,n-2*bw+11) * A(m+3,n-2*bw+2) - A(m+6,n-2*bw+11) * A(m+4,n-2*bw+1) - A(m+7,n-2*bw+11) * A(m+5,n-2*bw)    
    A(m-7,n-2*bw+12)  = A(m-7,n-2*bw+12) - A(m+3,n-2*bw+12) * A(m+1,n-2*bw+4) - A(m+4,n-2*bw+12) * A(m+2,n-2*bw+3) &
         - A(m+5,n-2*bw+12) * A(m+3,n-2*bw+2) - A(m+6,n-2*bw+12) * A(m+4,n-2*bw+1) - A(m+7,n-2*bw+12) * A(m+5,n-2*bw)    
    A(m+2,n-2*bw+13)  = A(m+2,n-2*bw+13) - A(m+3,n-2*bw+13) * A(m+1,n-2*bw+4) - A(m+4,n-2*bw+13) * A(m+2,n-2*bw+3) &
         - A(m+5,n-2*bw+13) * A(m+3,n-2*bw+2) - A(m+6,n-2*bw+13) * A(m+4,n-2*bw+1) - A(m+7,n-2*bw+13) * A(m+5,n-2*bw)    
    A(m+2,n-2*bw+14)  = A(m+2,n-2*bw+14) - A(m+3,n-2*bw+14) * A(m+1,n-2*bw+4) - A(m+4,n-2*bw+14) * A(m+2,n-2*bw+3) &
         - A(m+5,n-2*bw+14) * A(m+3,n-2*bw+2) - A(m+6,n-2*bw+14) * A(m+4,n-2*bw+1) - A(m+7,n-2*bw+14) * A(m+5,n-2*bw)    

    dummy = 1.d0 / A(m,n-2*bw+5)
    A(m+1,n-2*bw+5) = dummy * (A(m+1,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+2,n-2*bw+4) - A(m-2,n-2*bw+5) * A(m+3,n-2*bw+3)&
         - A(m-3,n-2*bw+5) * A(m+4,n-2*bw+2) - A(m-4,n-2*bw+5) * A(m+5,n-2*bw+1) - A(m-5,n-2*bw+5) * A(m+6,n-2*bw) &
         - A(m-6,n-2*bw+5) * A(m+7,n-2*bw-1))
    A(m+2,n-2*bw+5) = dummy * (A(m+2,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+3,n-2*bw+4) - A(m-2,n-2*bw+5) * A(m+4,n-2*bw+3)&
         - A(m-3,n-2*bw+5) * A(m+5,n-2*bw+2) - A(m-4,n-2*bw+5) * A(m+6,n-2*bw+1) - A(m-5,n-2*bw+5) * A(m+7,n-2*bw))
    A(m+3,n-2*bw+5) = dummy * (A(m+3,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+4,n-2*bw+4) - A(m-2,n-2*bw+5) * A(m+5,n-2*bw+3)&
         - A(m-3,n-2*bw+5) * A(m+6,n-2*bw+2) - A(m-4,n-2*bw+5) * A(m+7,n-2*bw+1))
    A(m+4,n-2*bw+5) = dummy * (A(m+4,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+5,n-2*bw+4) - A(m-2,n-2*bw+5) * A(m+6,n-2*bw+3)&
         - A(m-3,n-2*bw+5) * A(m+7,n-2*bw+2))
    A(m+5,n-2*bw+5) = dummy * (A(m+5,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+6,n-2*bw+4) - A(m-2,n-2*bw+5) * A(m+7,n-2*bw+3))
    A(m+6,n-2*bw+5) = dummy * (A(m+6,n-2*bw+5) - A(m-1,n-2*bw+5) * A(m+7,n-2*bw+4))
    A(m+7,n-2*bw+5) = dummy *  A(m+7,n-2*bw+5)


    A(m  ,n-2*bw+6)   = A(m  ,n-2*bw+6)  - A(m-1,n-2*bw+6)  * A(m+1,n-2*bw+5) - A(m-2,n-2*bw+6)  * A(m+2,n-2*bw+4) &
         - A(m-3,n-2*bw+6)  * A(m+3,n-2*bw+3) - A(m-4,n-2*bw+6)  * A(m+4,n-2*bw+2) - A(m-5,n-2*bw+6)  * A(m+5,n-2*bw+1) &
         - A(m-6,n-2*bw+6)  * A(m+6,n-2*bw) - A(m-7,n-2*bw+6) * A(m+7,n-2*bw-1)    
    A(m-1,n-2*bw+7)   = A(m-1,n-2*bw+7)  - A(m-2,n-2*bw+7)  * A(m+1,n-2*bw+5) - A(m-3,n-2*bw+7)  * A(m+2,n-2*bw+4)&
         - A(m-4,n-2*bw+7)  * A(m+3,n-2*bw+3) - A(m-5,n-2*bw+7)  * A(m+4,n-2*bw+2) - A(m-6,n-2*bw+7)  * A(m+5,n-2*bw+1)&
         - A(m-7,n-2*bw+7)  * A(m+6,n-2*bw)  
    A(m-2,n-2*bw+8)   = A(m-2,n-2*bw+8)  - A(m-3,n-2*bw+8)  * A(m+1,n-2*bw+5) - A(m-4,n-2*bw+8)  * A(m+2,n-2*bw+4) &
         - A(m-5,n-2*bw+8)  * A(m+3,n-2*bw+3) - A(m-6,n-2*bw+8)  * A(m+4,n-2*bw+2) - A(m-7,n-2*bw+8)  * A(m+5,n-2*bw+1) &
         - A(m+7,n-2*bw+8)  * A(m+6,n-2*bw)    
    A(m-3,n-2*bw+9)   = A(m-3,n-2*bw+9)  - A(m-4,n-2*bw+9)  * A(m+1,n-2*bw+5) - A(m-5,n-2*bw+9)  * A(m+2,n-2*bw+4) &
         - A(m-6,n-2*bw+9)  * A(m+3,n-2*bw+3) - A(m-7,n-2*bw+9)  * A(m+4,n-2*bw+2) - A(m+6,n-2*bw+9)  * A(m+5,n-2*bw+1)&
         - A(m+7,n-2*bw+9)  * A(m+6,n-2*bw)    
    A(m-4,n-2*bw+10)  = A(m-4,n-2*bw+10) - A(m-5,n-2*bw+10) * A(m+1,n-2*bw+5) - A(m-6,n-2*bw+10) * A(m+2,n-2*bw+4) &
         - A(m-7,n-2*bw+10) * A(m+3,n-2*bw+3) - A(m+5,n-2*bw+10) * A(m+4,n-2*bw+2) - A(m+6,n-2*bw+10) * A(m+5,n-2*bw+1) &
         - A(m+7,n-2*bw+10) * A(m+6,n-2*bw)    
    A(m-5,n-2*bw+11)  = A(m-5,n-2*bw+11) - A(m-6,n-2*bw+11) * A(m+1,n-2*bw+5) - A(m-7,n-2*bw+11) * A(m+2,n-2*bw+4)&
         - A(m+4,n-2*bw+11) * A(m+3,n-2*bw+3) - A(m+5,n-2*bw+11) * A(m+4,n-2*bw+2) - A(m+6,n-2*bw+11) * A(m+5,n-2*bw+1) &
         - A(m+7,n-2*bw+11) * A(m+6,n-2*bw)    
    A(m-6,n-2*bw+12)  = A(m-6,n-2*bw+12) - A(m-7,n-2*bw+12) * A(m+1,n-2*bw+5) - A(m+3,n-2*bw+12) * A(m+2,n-2*bw+4) &
         - A(m+4,n-2*bw+12) * A(m+3,n-2*bw+3) - A(m+5,n-2*bw+12) * A(m+4,n-2*bw+2) - A(m+6,n-2*bw+12) * A(m+5,n-2*bw+1)&
         - A(m+7,n-2*bw+12) * A(m+6,n-2*bw)    
    A(m-7,n-2*bw+13)  = A(m-7,n-2*bw+13) - A(m+2,n-2*bw+13) * A(m+1,n-2*bw+5) - A(m+3,n-2*bw+13) * A(m+2,n-2*bw+4)&
         - A(m+4,n-2*bw+13) * A(m+3,n-2*bw+3) - A(m+5,n-2*bw+13) * A(m+4,n-2*bw+2) - A(m+6,n-2*bw+13) * A(m+5,n-2*bw+1)&
         - A(m+7,n-2*bw+13) * A(m+6,n-2*bw)    
    A(m+1,n-2*bw+14)  = A(m+1,n-2*bw+14) - A(m+2,n-2*bw+14) * A(m+1,n-2*bw+5) - A(m+3,n-2*bw+14) * A(m+2,n-2*bw+4) &
         - A(m+4,n-2*bw+14) * A(m+3,n-2*bw+3) - A(m+5,n-2*bw+14) * A(m+4,n-2*bw+2) - A(m+6,n-2*bw+14) * A(m+5,n-2*bw+1)&
         - A(m+7,n-2*bw+14) * A(m+6,n-2*bw)

    dummy = 1.d0 / A(m,n-2*bw+6)
    A(m+1,n-2*bw+6) = dummy * (A(m+1,n-2*bw+6) - A(m-1,n-2*bw+6) * A(m+2,n-2*bw+5) - A(m-2,n-2*bw+6) * A(m+3,n-2*bw+4)&
         - A(m-3,n-2*bw+6) * A(m+4,n-2*bw+3) - A(m-4,n-2*bw+6) * A(m+5,n-2*bw+2) - A(m-5,n-2*bw+6) * A(m+6,n-2*bw+1)&
         - A(m-6,n-2*bw+6) * A(m+7,n-2*bw))
    A(m+2,n-2*bw+6) = dummy * (A(m+2,n-2*bw+6) - A(m-1,n-2*bw+6) * A(m+3,n-2*bw+5) - A(m-2,n-2*bw+6) * A(m+4,n-2*bw+4)&
         - A(m-3,n-2*bw+6) * A(m+5,n-2*bw+3) - A(m-4,n-2*bw+6) * A(m+6,n-2*bw+2) - A(m-5,n-2*bw+6) * A(m+7,n-2*bw+1))
    A(m+3,n-2*bw+6) = dummy * (A(m+3,n-2*bw+6) - A(m-1,n-2*bw+6) * A(m+4,n-2*bw+5) - A(m-2,n-2*bw+6) * A(m+5,n-2*bw+4)&
         - A(m-3,n-2*bw+6) * A(m+6,n-2*bw+3) - A(m-4,n-2*bw+6) * A(m+7,n-2*bw+2))
    A(m+4,n-2*bw+6) = dummy * (A(m+4,n-2*bw+6) - A(m-1,n-2*bw+6) * A(m+5,n-2*bw+5) - A(m-2,n-2*bw+6) * A(m+6,n-2*bw+4) &
         - A(m-3,n-2*bw+6) * A(m+7,n-2*bw+3))
    A(m+5,n-2*bw+6) = dummy * (A(m+5,n-2*bw+6) - A(m-1,n-2*bw+6) * A(m+6,n-2*bw+5) - A(m-2,n-2*bw+6) * A(m+7,n-2*bw+4))
    A(m+6,n-2*bw+6) = dummy * (A(m+6,n-2*bw+6) - A(m-1,n-2*bw+6) * A(m+7,n-2*bw+5))
    A(m+7,n-2*bw+6) = dummy *  A(m+7,n-2*bw+6)

    A(m  ,n-2*bw+7)   = A(m  ,n-2*bw+7)  - A(m-1,n-2*bw+7)  * A(m+1,n-2*bw+6) - A(m-2,n-2*bw+7)  * A(m+2,n-2*bw+5) &
         - A(m-3,n-2*bw+7)  * A(m+3,n-2*bw+4) - A(m-4,n-2*bw+7)  * A(m+4,n-2*bw+3) - A(m-5,n-2*bw+7)  * A(m+5,n-2*bw+2)&
         - A(m-6,n-2*bw+7)  * A(m+6,n-2*bw+1) - A(m-7,n-2*bw+7) * A(m+7,n-2*bw)    
    A(m-1,n-2*bw+8)   = A(m-1,n-2*bw+8)  - A(m-2,n-2*bw+8)  * A(m+1,n-2*bw+6) - A(m-3,n-2*bw+8)  * A(m+2,n-2*bw+5) &
         - A(m-4,n-2*bw+8)  * A(m+3,n-2*bw+4) - A(m-5,n-2*bw+8)  * A(m+4,n-2*bw+3) - A(m-6,n-2*bw+8)  * A(m+5,n-2*bw+2) &
         - A(m-7,n-2*bw+8)  * A(m+6,n-2*bw+1) - A(m+7,n-2*bw+8) * A(m+7,n-2*bw)
    A(m-2,n-2*bw+9)   = A(m-2,n-2*bw+9)  - A(m-3,n-2*bw+9)  * A(m+1,n-2*bw+6) - A(m-4,n-2*bw+9)  * A(m+2,n-2*bw+5) &
         - A(m-5,n-2*bw+9)  * A(m+3,n-2*bw+4) - A(m-6,n-2*bw+9)  * A(m+4,n-2*bw+3) - A(m-7,n-2*bw+9)  * A(m+5,n-2*bw+2)&
         - A(m+6,n-2*bw+9)  * A(m+6,n-2*bw+1) - A(m+7,n-2*bw+9) * A(m+7,n-2*bw)
    A(m-3,n-2*bw+10)  = A(m-3,n-2*bw+10) - A(m-4,n-2*bw+10) * A(m+1,n-2*bw+6) - A(m-5,n-2*bw+10) * A(m+2,n-2*bw+5) &
         - A(m-6,n-2*bw+10) * A(m+3,n-2*bw+4) - A(m-7,n-2*bw+10) * A(m+4,n-2*bw+3) - A(m+5,n-2*bw+10) * A(m+5,n-2*bw+2) &
         - A(m+6,n-2*bw+10) * A(m+6,n-2*bw+1) - A(m+7,n-2*bw+10) * A(m+7,n-2*bw)
    A(m-4,n-2*bw+11)  = A(m-4,n-2*bw+11) - A(m-5,n-2*bw+11) * A(m+1,n-2*bw+6) - A(m-6,n-2*bw+11) * A(m+2,n-2*bw+5) &
         - A(m-7,n-2*bw+11) * A(m+3,n-2*bw+4) - A(m+4,n-2*bw+11) * A(m+4,n-2*bw+3) - A(m+5,n-2*bw+11) * A(m+5,n-2*bw+2) &
         - A(m+6,n-2*bw+11) * A(m+6,n-2*bw+1) - A(m+7,n-2*bw+11) * A(m+7,n-2*bw)
    A(m-5,n-2*bw+12)  = A(m-5,n-2*bw+12) - A(m-6,n-2*bw+12) * A(m+1,n-2*bw+6) - A(m-7,n-2*bw+12) * A(m+2,n-2*bw+5) &
         - A(m+3,n-2*bw+12) * A(m+3,n-2*bw+4) - A(m+4,n-2*bw+12) * A(m+4,n-2*bw+3) - A(m+5,n-2*bw+12) * A(m+5,n-2*bw+2)&
         - A(m+6,n-2*bw+12) * A(m+6,n-2*bw+1) - A(m+7,n-2*bw+12) * A(m+7,n-2*bw)
    A(m-6,n-2*bw+13)  = A(m-6,n-2*bw+13) - A(m-7,n-2*bw+13) * A(m+1,n-2*bw+6) - A(m+2,n-2*bw+13) * A(m+2,n-2*bw+5) &
         - A(m+3,n-2*bw+13) * A(m+3,n-2*bw+4) - A(m+4,n-2*bw+13) * A(m+4,n-2*bw+3) - A(m+5,n-2*bw+13) * A(m+5,n-2*bw+2)&
         - A(m+6,n-2*bw+13) * A(m+6,n-2*bw+1) - A(m+7,n-2*bw+13) * A(m+7,n-2*bw)
    A(m-7,n-2*bw+14)  = A(m-7,n-2*bw+14) - A(m+1,n-2*bw+14) * A(m+1,n-2*bw+6) - A(m+2,n-2*bw+14) * A(m+2,n-2*bw+5) &
         - A(m+3,n-2*bw+14) * A(m+3,n-2*bw+4) - A(m+4,n-2*bw+14) * A(m+4,n-2*bw+3) - A(m+5,n-2*bw+14) * A(m+5,n-2*bw+2) &
         - A(m+6,n-2*bw+14) * A(m+6,n-2*bw+1) - A(m+7,n-2*bw+14) * A(m+7,n-2*bw)

    dummy = 1.d0 / A(m,n-2*bw+7)
    A(m+1,n-2*bw+7) = dummy * (A(m+1,n-2*bw+7) - A(m-1,n-2*bw+7) * A(m+2,n-2*bw+6) - A(m-2,n-2*bw+7) * A(m+3,n-2*bw+5)&
         - A(m-3,n-2*bw+7) * A(m+4,n-2*bw+4) - A(m-4,n-2*bw+7) * A(m+5,n-2*bw+3) - A(m-5,n-2*bw+7) * A(m+6,n-2*bw+2) &
         - A(m-6,n-2*bw+7) * A(m+7,n-2*bw+1))
    A(m+2,n-2*bw+7) = dummy * (A(m+2,n-2*bw+7) - A(m-1,n-2*bw+7) * A(m+3,n-2*bw+6) - A(m-2,n-2*bw+7) * A(m+4,n-2*bw+5)&
         - A(m-3,n-2*bw+7) * A(m+5,n-2*bw+4) - A(m-4,n-2*bw+7) * A(m+6,n-2*bw+3) - A(m-5,n-2*bw+7) * A(m+7,n-2*bw+2))
    A(m+3,n-2*bw+7) = dummy * (A(m+3,n-2*bw+7) - A(m-1,n-2*bw+7) * A(m+4,n-2*bw+6) - A(m-2,n-2*bw+7) * A(m+5,n-2*bw+5)&
         - A(m-3,n-2*bw+7) * A(m+6,n-2*bw+4) - A(m-4,n-2*bw+7) * A(m+7,n-2*bw+3))
    A(m+4,n-2*bw+7) = dummy * (A(m+4,n-2*bw+7) - A(m-1,n-2*bw+7) * A(m+5,n-2*bw+6) - A(m-2,n-2*bw+7) * A(m+6,n-2*bw+5) &
         - A(m-3,n-2*bw+7) * A(m+7,n-2*bw+4))
    A(m+5,n-2*bw+7) = dummy * (A(m+5,n-2*bw+7) - A(m-1,n-2*bw+7) * A(m+6,n-2*bw+6) - A(m-2,n-2*bw+7) * A(m+7,n-2*bw+5))
    A(m+6,n-2*bw+7) = dummy * (A(m+6,n-2*bw+7) - A(m-1,n-2*bw+7) * A(m+7,n-2*bw+6))
    A(m+7,n-2*bw+7) = dummy *  A(m+7,n-2*bw+7)

    ! Lower region
    A(m  ,n-2*bw+8)   = A(m  ,n-2*bw+8)  - A(m-1,n-2*bw+8)  * A(m+1,n-2*bw+7) - A(m-2,n-2*bw+8)  * A(m+2,n-2*bw+6) &
         - A(m-3,n-2*bw+8)  * A(m+3,n-2*bw+5) - A(m-4,n-2*bw+8)  * A(m+4,n-2*bw+4) - A(m-5,n-2*bw+8)  * A(m+5,n-2*bw+3) &
         - A(m-6,n-2*bw+8)  * A(m+6,n-2*bw+2) - A(m-7,n-2*bw+8) * A(m+7,n-2*bw+1)    
    A(m-1,n-2*bw+9)   = A(m-1,n-2*bw+9)  - A(m-2,n-2*bw+9)  * A(m+1,n-2*bw+7) - A(m-3,n-2*bw+9)  * A(m+2,n-2*bw+6)&
         - A(m-4,n-2*bw+9)  * A(m+3,n-2*bw+5) - A(m-5,n-2*bw+9)  * A(m+4,n-2*bw+4) - A(m-6,n-2*bw+9)  * A(m+5,n-2*bw+3) &
         - A(m-7,n-2*bw+9)  * A(m+6,n-2*bw+2) - A(m+6,n-2*bw+9) * A(m+7,n-2*bw+1)
    A(m-2,n-2*bw+10)  = A(m-2,n-2*bw+10) - A(m-3,n-2*bw+10) * A(m+1,n-2*bw+7) - A(m-4,n-2*bw+10) * A(m+2,n-2*bw+6) &
         - A(m-5,n-2*bw+10) * A(m+3,n-2*bw+5) - A(m-6,n-2*bw+10) * A(m+4,n-2*bw+4) - A(m-7,n-2*bw+10) * A(m+5,n-2*bw+3) &
         - A(m+5,n-2*bw+10) * A(m+6,n-2*bw+2) - A(m+6,n-2*bw+10) * A(m+7,n-2*bw+1)
    A(m-3,n-2*bw+11)  = A(m-3,n-2*bw+11) - A(m-4,n-2*bw+11) * A(m+1,n-2*bw+7) - A(m-5,n-2*bw+11) * A(m+2,n-2*bw+6) &
         - A(m-6,n-2*bw+11) * A(m+3,n-2*bw+5) - A(m-7,n-2*bw+11) * A(m+4,n-2*bw+4) - A(m+4,n-2*bw+11) * A(m+5,n-2*bw+3) &
         - A(m+5,n-2*bw+11) * A(m+6,n-2*bw+2) - A(m+6,n-2*bw+11) * A(m+7,n-2*bw+1)
    A(m-4,n-2*bw+12)  = A(m-4,n-2*bw+12) - A(m-5,n-2*bw+12) * A(m+1,n-2*bw+7) - A(m-6,n-2*bw+12) * A(m+2,n-2*bw+6) &
         - A(m-7,n-2*bw+12) * A(m+3,n-2*bw+5) - A(m+3,n-2*bw+12) * A(m+4,n-2*bw+4) - A(m+4,n-2*bw+12) * A(m+5,n-2*bw+3) &
         - A(m+5,n-2*bw+12) * A(m+6,n-2*bw+2) - A(m+6,n-2*bw+12) * A(m+7,n-2*bw+1)
    A(m-5,n-2*bw+13)  = A(m-5,n-2*bw+13) - A(m-6,n-2*bw+13) * A(m+1,n-2*bw+7) - A(m-7,n-2*bw+13) * A(m+2,n-2*bw+6) &
         - A(m+2,n-2*bw+13) * A(m+3,n-2*bw+5) - A(m+3,n-2*bw+13) * A(m+4,n-2*bw+4) - A(m+4,n-2*bw+13) * A(m+5,n-2*bw+3) &
         - A(m+5,n-2*bw+13) * A(m+6,n-2*bw+2) - A(m+6,n-2*bw+13) * A(m+7,n-2*bw+1)
    A(m-6,n-2*bw+14)  = A(m-6,n-2*bw+14) - A(m-7,n-2*bw+14) * A(m+1,n-2*bw+7) - A(m+1,n-2*bw+14) * A(m+2,n-2*bw+6) &
         - A(m+2,n-2*bw+14) * A(m+3,n-2*bw+5) - A(m+3,n-2*bw+14) * A(m+4,n-2*bw+4) - A(m+4,n-2*bw+14) * A(m+5,n-2*bw+3) &
         - A(m+5,n-2*bw+14) * A(m+6,n-2*bw+2) - A(m+6,n-2*bw+14) * A(m+7,n-2*bw+1)

    dummy = 1.d0 / A(m,n-2*bw+8)
    A(m+1,n-2*bw+8) = dummy * (A(m+1,n-2*bw+8) - A(m-1,n-2*bw+8) * A(m+2,n-2*bw+7) - A(m-2,n-2*bw+8) * A(m+3,n-2*bw+6)&
         - A(m-3,n-2*bw+8) * A(m+4,n-2*bw+5) - A(m-4,n-2*bw+8) * A(m+5,n-2*bw+4) - A(m-5,n-2*bw+8) * A(m+6,n-2*bw+3) &
         - A(m-6,n-2*bw+8) * A(m+7,n-2*bw+2))
    A(m+2,n-2*bw+8) = dummy * (A(m+2,n-2*bw+8) - A(m-1,n-2*bw+8) * A(m+3,n-2*bw+7) - A(m-2,n-2*bw+8) * A(m+4,n-2*bw+6)&
         - A(m-3,n-2*bw+8) * A(m+5,n-2*bw+5) - A(m-4,n-2*bw+8) * A(m+6,n-2*bw+4) - A(m-5,n-2*bw+8) * A(m+7,n-2*bw+3))
    A(m+3,n-2*bw+8) = dummy * (A(m+3,n-2*bw+8) - A(m-1,n-2*bw+8) * A(m+4,n-2*bw+7) - A(m-2,n-2*bw+8) * A(m+5,n-2*bw+6)&
         - A(m-3,n-2*bw+8) * A(m+6,n-2*bw+5) - A(m-4,n-2*bw+8) * A(m+7,n-2*bw+4))
    A(m+4,n-2*bw+8) = dummy * (A(m+4,n-2*bw+8) - A(m-1,n-2*bw+8) * A(m+5,n-2*bw+7) - A(m-2,n-2*bw+8) * A(m+6,n-2*bw+6)&
         - A(m-3,n-2*bw+8) * A(m+7,n-2*bw+5))
    A(m+5,n-2*bw+8) = dummy * (A(m+5,n-2*bw+8) - A(m-1,n-2*bw+8) * A(m+6,n-2*bw+7) - A(m-2,n-2*bw+8) * A(m+7,n-2*bw+6))
    A(m+6,n-2*bw+8) = dummy * (A(m+6,n-2*bw+8) - A(m-1,n-2*bw+8) * A(m+7,n-2*bw+7))

    A(m  ,n-2*bw+9)   = A(m  ,n-2*bw+9)  - A(m-1,n-2*bw+9)  * A(m+1,n-2*bw+8) - A(m-2,n-2*bw+9)  * A(m+2,n-2*bw+7)&
         - A(m-3,n-2*bw+9)  * A(m+3,n-2*bw+6) - A(m-4,n-2*bw+9)  * A(m+4,n-2*bw+5) - A(m-5,n-2*bw+9)  * A(m+5,n-2*bw+4)&
         - A(m-6,n-2*bw+9)  * A(m+6,n-2*bw+3) - A(m-7,n-2*bw+9)  * A(m+7,n-2*bw+2)    
    A(m-1,n-2*bw+10)  = A(m-1,n-2*bw+10) - A(m-2,n-2*bw+10) * A(m+1,n-2*bw+8) - A(m-3,n-2*bw+10) * A(m+2,n-2*bw+7) &
         - A(m-4,n-2*bw+10) * A(m+3,n-2*bw+6) - A(m-5,n-2*bw+10) * A(m+4,n-2*bw+5) - A(m-6,n-2*bw+10) * A(m+5,n-2*bw+4)&
         - A(m-7,n-2*bw+10) * A(m+6,n-2*bw+3) - A(m+5,n-2*bw+10) * A(m+7,n-2*bw+2)
    A(m-2,n-2*bw+11)  = A(m-2,n-2*bw+11) - A(m-3,n-2*bw+11) * A(m+1,n-2*bw+8) - A(m-4,n-2*bw+11) * A(m+2,n-2*bw+7) &
         - A(m-5,n-2*bw+11) * A(m+3,n-2*bw+6) - A(m-6,n-2*bw+11) * A(m+4,n-2*bw+5) - A(m-7,n-2*bw+11) * A(m+5,n-2*bw+4) &
         - A(m+4,n-2*bw+11) * A(m+6,n-2*bw+3) - A(m+5,n-2*bw+11) * A(m+7,n-2*bw+2)
    A(m-3,n-2*bw+12)  = A(m-3,n-2*bw+12) - A(m-4,n-2*bw+12) * A(m+1,n-2*bw+8) - A(m-5,n-2*bw+12) * A(m+2,n-2*bw+7) &
         - A(m-6,n-2*bw+12) * A(m+3,n-2*bw+6) - A(m-7,n-2*bw+12) * A(m+4,n-2*bw+5) - A(m+3,n-2*bw+12) * A(m+5,n-2*bw+4) &
         - A(m+4,n-2*bw+12) * A(m+6,n-2*bw+3) - A(m+5,n-2*bw+12) * A(m+7,n-2*bw+2)
    A(m-4,n-2*bw+13)  = A(m-4,n-2*bw+13) - A(m-5,n-2*bw+13) * A(m+1,n-2*bw+8) - A(m-6,n-2*bw+13) * A(m+2,n-2*bw+7) &
         - A(m-7,n-2*bw+13) * A(m+3,n-2*bw+6) - A(m+2,n-2*bw+13) * A(m+4,n-2*bw+5) - A(m+3,n-2*bw+13) * A(m+5,n-2*bw+4) &
         - A(m+4,n-2*bw+13) * A(m+6,n-2*bw+3) - A(m+5,n-2*bw+13) * A(m+7,n-2*bw+2)
    A(m-5,n-2*bw+14)  = A(m-5,n-2*bw+14) - A(m-6,n-2*bw+14) * A(m+1,n-2*bw+8) - A(m-7,n-2*bw+14) * A(m+2,n-2*bw+7) &
         - A(m+1,n-2*bw+14) * A(m+3,n-2*bw+6) - A(m+2,n-2*bw+14) * A(m+4,n-2*bw+5) - A(m+3,n-2*bw+14) * A(m+5,n-2*bw+4)&
         - A(m+4,n-2*bw+14) * A(m+6,n-2*bw+3) - A(m+5,n-2*bw+14) * A(m+7,n-2*bw+2)

    dummy = 1.d0 / A(m,n-2*bw+9)
    A(m+1,n-2*bw+9) = dummy * (A(m+1,n-2*bw+9) - A(m-1,n-2*bw+9) * A(m+2,n-2*bw+8) - A(m-2,n-2*bw+9) * A(m+3,n-2*bw+7)&
         - A(m-3,n-2*bw+9) * A(m+4,n-2*bw+6) - A(m-4,n-2*bw+9) * A(m+5,n-2*bw+5) - A(m-5,n-2*bw+9) * A(m+6,n-2*bw+4) &
         - A(m-6,n-2*bw+9) * A(m+7,n-2*bw+3))
    A(m+2,n-2*bw+9) = dummy * (A(m+2,n-2*bw+9) - A(m-1,n-2*bw+9) * A(m+3,n-2*bw+8) - A(m-2,n-2*bw+9) * A(m+4,n-2*bw+7) &
         - A(m-3,n-2*bw+9) * A(m+5,n-2*bw+6) - A(m-4,n-2*bw+9) * A(m+6,n-2*bw+5) - A(m-5,n-2*bw+9) * A(m+7,n-2*bw+4))
    A(m+3,n-2*bw+9) = dummy * (A(m+3,n-2*bw+9) - A(m-1,n-2*bw+9) * A(m+4,n-2*bw+8) - A(m-2,n-2*bw+9) * A(m+5,n-2*bw+7)&
         - A(m-3,n-2*bw+9) * A(m+6,n-2*bw+6) - A(m-4,n-2*bw+9) * A(m+7,n-2*bw+5))
    A(m+4,n-2*bw+9) = dummy * (A(m+4,n-2*bw+9) - A(m-1,n-2*bw+9) * A(m+5,n-2*bw+8) - A(m-2,n-2*bw+9) * A(m+6,n-2*bw+7) &
         - A(m-3,n-2*bw+9) * A(m+7,n-2*bw+6))
    A(m+5,n-2*bw+9) = dummy * (A(m+5,n-2*bw+9) - A(m-1,n-2*bw+9) * A(m+6,n-2*bw+8) - A(m-2,n-2*bw+9) * A(m+7,n-2*bw+7))

    A(m  ,n-2*bw+10)  = A(m  ,n-2*bw+10) - A(m-1,n-2*bw+10) * A(m+1,n-2*bw+9) - A(m-2,n-2*bw+10) * A(m+2,n-2*bw+8) &
         - A(m-3,n-2*bw+10) * A(m+3,n-2*bw+7) - A(m-4,n-2*bw+10) * A(m+4,n-2*bw+6) - A(m-5,n-2*bw+10) * A(m+5,n-2*bw+5) &
         - A(m-6,n-2*bw+10) * A(m+6,n-2*bw+4) - A(m-7,n-2*bw+10) * A(m+7,n-2*bw+3)    
    A(m-1,n-2*bw+11)  = A(m-1,n-2*bw+11) - A(m-2,n-2*bw+11) * A(m+1,n-2*bw+9) - A(m-3,n-2*bw+11) * A(m+2,n-2*bw+8) &
         - A(m-4,n-2*bw+11) * A(m+3,n-2*bw+7) - A(m-5,n-2*bw+11) * A(m+4,n-2*bw+6) - A(m-6,n-2*bw+11) * A(m+5,n-2*bw+5) &
         - A(m-7,n-2*bw+11) * A(m+6,n-2*bw+4) - A(m+4,n-2*bw+11) * A(m+7,n-2*bw+3)
    A(m-2,n-2*bw+12)  = A(m-2,n-2*bw+12) - A(m-3,n-2*bw+12) * A(m+1,n-2*bw+9) - A(m-4,n-2*bw+12) * A(m+2,n-2*bw+8)&
         - A(m-5,n-2*bw+12) * A(m+3,n-2*bw+7) - A(m-6,n-2*bw+12) * A(m+4,n-2*bw+6) - A(m-7,n-2*bw+12) * A(m+5,n-2*bw+5)&
         - A(m+3,n-2*bw+12) * A(m+6,n-2*bw+4) - A(m+4,n-2*bw+12) * A(m+7,n-2*bw+3)
    A(m-3,n-2*bw+13)  = A(m-3,n-2*bw+13) - A(m-4,n-2*bw+13) * A(m+1,n-2*bw+9) - A(m-5,n-2*bw+13) * A(m+2,n-2*bw+8) &
         - A(m-6,n-2*bw+13) * A(m+3,n-2*bw+7) - A(m-7,n-2*bw+13) * A(m+4,n-2*bw+6) - A(m+2,n-2*bw+13) * A(m+5,n-2*bw+5) &
         - A(m+3,n-2*bw+13) * A(m+6,n-2*bw+4) - A(m+4,n-2*bw+13) * A(m+7,n-2*bw+3)
    A(m-4,n-2*bw+14)  = A(m-4,n-2*bw+14) - A(m-5,n-2*bw+14) * A(m+1,n-2*bw+9) - A(m-6,n-2*bw+14) * A(m+2,n-2*bw+8)&
         - A(m-7,n-2*bw+14) * A(m+3,n-2*bw+7) - A(m+1,n-2*bw+14) * A(m+4,n-2*bw+6) - A(m+2,n-2*bw+14) * A(m+5,n-2*bw+5)&
         - A(m+3,n-2*bw+14) * A(m+6,n-2*bw+4) - A(m+4,n-2*bw+14) * A(m+7,n-2*bw+3)

    dummy = 1.d0 / A(m,n-2*bw+10)
    A(m+1,n-2*bw+10) = dummy * (A(m+1,n-2*bw+10) - A(m-1,n-2*bw+10) * A(m+2,n-2*bw+9) - A(m-2,n-2*bw+10) * A(m+3,n-2*bw+8)&
         - A(m-3,n-2*bw+10) * A(m+4,n-2*bw+7) - A(m-4,n-2*bw+10) * A(m+5,n-2*bw+6) - A(m-5,n-2*bw+10) * A(m+6,n-2*bw+5) &
         - A(m-6,n-2*bw+10) * A(m+7,n-2*bw+4))
    A(m+2,n-2*bw+10) = dummy * (A(m+2,n-2*bw+10) - A(m-1,n-2*bw+10) * A(m+3,n-2*bw+9) - A(m-2,n-2*bw+10) * A(m+4,n-2*bw+8)&
         - A(m-3,n-2*bw+10) * A(m+5,n-2*bw+7) - A(m-4,n-2*bw+10) * A(m+6,n-2*bw+6) - A(m-5,n-2*bw+10) * A(m+7,n-2*bw+5))
    A(m+3,n-2*bw+10) = dummy * (A(m+3,n-2*bw+10) - A(m-1,n-2*bw+10) * A(m+4,n-2*bw+9) - A(m-2,n-2*bw+10) * A(m+5,n-2*bw+8)&
         - A(m-3,n-2*bw+10) * A(m+6,n-2*bw+7) - A(m-4,n-2*bw+10) * A(m+7,n-2*bw+6))
    A(m+4,n-2*bw+10) = dummy * (A(m+4,n-2*bw+10) - A(m-1,n-2*bw+10) * A(m+5,n-2*bw+9) - A(m-2,n-2*bw+10) * A(m+6,n-2*bw+8) &
         - A(m-3,n-2*bw+10) * A(m+7,n-2*bw+7))

    A(m  ,n-2*bw+11)  = A(m  ,n-2*bw+11) - A(m-1,n-2*bw+11) * A(m+1,n-2*bw+10) - A(m-2,n-2*bw+11) * A(m+2,n-2*bw+9)&
         - A(m-3,n-2*bw+11) * A(m+3,n-2*bw+8) - A(m-4,n-2*bw+11) * A(m+4,n-2*bw+7) - A(m-5,n-2*bw+11) * A(m+5,n-2*bw+6) &
         - A(m-6,n-2*bw+11) * A(m+6,n-2*bw+5) - A(m-7,n-2*bw+11) * A(m+7,n-2*bw+4)    
    A(m-1,n-2*bw+12)  = A(m-1,n-2*bw+12) - A(m-2,n-2*bw+12) * A(m+1,n-2*bw+10) - A(m-3,n-2*bw+12) * A(m+2,n-2*bw+9)&
         - A(m-4,n-2*bw+12) * A(m+3,n-2*bw+8) - A(m-5,n-2*bw+12) * A(m+4,n-2*bw+7) - A(m-6,n-2*bw+12) * A(m+5,n-2*bw+6)&
         - A(m-7,n-2*bw+12) * A(m+6,n-2*bw+5) - A(m+3,n-2*bw+12) * A(m+7,n-2*bw+4)
    A(m-2,n-2*bw+13)  = A(m-2,n-2*bw+13) - A(m-3,n-2*bw+13) * A(m+1,n-2*bw+10) - A(m-4,n-2*bw+13) * A(m+2,n-2*bw+9)&
         - A(m-5,n-2*bw+13) * A(m+3,n-2*bw+8) - A(m-6,n-2*bw+13) * A(m+4,n-2*bw+7) - A(m-7,n-2*bw+13) * A(m+5,n-2*bw+6)&
         - A(m+2,n-2*bw+13) * A(m+6,n-2*bw+5) - A(m+3,n-2*bw+13) * A(m+7,n-2*bw+4)
    A(m-3,n-2*bw+14)  = A(m-3,n-2*bw+14) - A(m-4,n-2*bw+14) * A(m+1,n-2*bw+10) - A(m-5,n-2*bw+14) * A(m+2,n-2*bw+9) &
         - A(m-6,n-2*bw+14) * A(m+3,n-2*bw+8) - A(m-7,n-2*bw+14) * A(m+4,n-2*bw+7) - A(m+1,n-2*bw+14) * A(m+5,n-2*bw+6) &
         - A(m+2,n-2*bw+14) * A(m+6,n-2*bw+5) - A(m+3,n-2*bw+14) * A(m+7,n-2*bw+4)

    dummy = 1.d0 / A(m,n-2*bw+11)
    A(m+1,n-2*bw+11) = dummy * (A(m+1,n-2*bw+11) - A(m-1,n-2*bw+11) * A(m+2,n-2*bw+10) - A(m-2,n-2*bw+11) * A(m+3,n-2*bw+9)&
         - A(m-3,n-2*bw+11) * A(m+4,n-2*bw+8) - A(m-4,n-2*bw+11) * A(m+5,n-2*bw+7) - A(m-5,n-2*bw+11) * A(m+6,n-2*bw+6) &
         - A(m-6,n-2*bw+11) * A(m+7,n-2*bw+5))
    A(m+2,n-2*bw+11) = dummy * (A(m+2,n-2*bw+11) - A(m-1,n-2*bw+11) * A(m+3,n-2*bw+10) - A(m-2,n-2*bw+11) * A(m+4,n-2*bw+9)&
         - A(m-3,n-2*bw+11) * A(m+5,n-2*bw+8) - A(m-4,n-2*bw+11) * A(m+6,n-2*bw+7) - A(m-5,n-2*bw+11) * A(m+7,n-2*bw+6))
    A(m+3,n-2*bw+11) = dummy * (A(m+3,n-2*bw+11) - A(m-1,n-2*bw+11) * A(m+4,n-2*bw+10) - A(m-2,n-2*bw+11) * A(m+5,n-2*bw+9) &
         - A(m-3,n-2*bw+11) * A(m+6,n-2*bw+8) - A(m-4,n-2*bw+11) * A(m+7,n-2*bw+7))

    A(m  ,n-2*bw+12)  = A(m  ,n-2*bw+12) - A(m-1,n-2*bw+12) * A(m+1,n-2*bw+11) - A(m-2,n-2*bw+12) * A(m+2,n-2*bw+10) &
         - A(m-3,n-2*bw+12) * A(m+3,n-2*bw+9) - A(m-4,n-2*bw+12) * A(m+4,n-2*bw+8) - A(m-5,n-2*bw+12) * A(m+5,n-2*bw+7) &
         - A(m-6,n-2*bw+12) * A(m+6,n-2*bw+6) - A(m-7,n-2*bw+12) * A(m+7,n-2*bw+5)    
    A(m-1,n-2*bw+13)  = A(m-1,n-2*bw+13) - A(m-2,n-2*bw+13) * A(m+1,n-2*bw+11) - A(m-3,n-2*bw+13) * A(m+2,n-2*bw+10) &
         - A(m-4,n-2*bw+13) * A(m+3,n-2*bw+9) - A(m-5,n-2*bw+13) * A(m+4,n-2*bw+8) - A(m-6,n-2*bw+13) * A(m+5,n-2*bw+7)&
         - A(m-7,n-2*bw+13) * A(m+6,n-2*bw+6) - A(m+2,n-2*bw+13) * A(m+7,n-2*bw+5)
    A(m-2,n-2*bw+14)  = A(m-2,n-2*bw+14) - A(m-3,n-2*bw+14) * A(m+1,n-2*bw+11) - A(m-4,n-2*bw+14) * A(m+2,n-2*bw+10) &
         - A(m-5,n-2*bw+14) * A(m+3,n-2*bw+9) - A(m-6,n-2*bw+14) * A(m+4,n-2*bw+8) - A(m-7,n-2*bw+14) * A(m+5,n-2*bw+7) &
         - A(m+1,n-2*bw+14) * A(m+6,n-2*bw+6) - A(m+2,n-2*bw+14) * A(m+7,n-2*bw+5)

    dummy = 1.d0 / A(m,n-2*bw+12)
    A(m+1,n-2*bw+12) = dummy * (A(m+1,n-2*bw+12) - A(m-1,n-2*bw+12) * A(m+2,n-2*bw+11) - A(m-2,n-2*bw+12) * A(m+3,n-2*bw+10)&
         - A(m-3,n-2*bw+12) * A(m+4,n-2*bw+9) - A(m-4,n-2*bw+12) * A(m+5,n-2*bw+8) - A(m-5,n-2*bw+12) * A(m+6,n-2*bw+7)&
         - A(m-6,n-2*bw+12) * A(m+7,n-2*bw+6))
    A(m+2,n-2*bw+12) = dummy * (A(m+2,n-2*bw+12) - A(m-1,n-2*bw+12) * A(m+3,n-2*bw+11) - A(m-2,n-2*bw+12) * A(m+4,n-2*bw+10)&
         - A(m-3,n-2*bw+12) * A(m+5,n-2*bw+9) - A(m-4,n-2*bw+12) * A(m+6,n-2*bw+8) - A(m-5,n-2*bw+12) * A(m+7,n-2*bw+7))

    A(m  ,n-2*bw+13)  = A(m  ,n-2*bw+13) - A(m-1,n-2*bw+13) * A(m+1,n-2*bw+12) - A(m-2,n-2*bw+13) * A(m+2,n-2*bw+11)&
         - A(m-3,n-2*bw+13) * A(m+3,n-2*bw+10) - A(m-4,n-2*bw+13) * A(m+4,n-2*bw+9) - A(m-5,n-2*bw+13) * A(m+5,n-2*bw+8)&
         - A(m-6,n-2*bw+13) * A(m+6,n-2*bw+7) - A(m-7,n-2*bw+13) * A(m+7,n-2*bw+6)    
    A(m-1,n-2*bw+14)  = A(m-1,n-2*bw+14) - A(m-2,n-2*bw+14) * A(m+1,n-2*bw+12) - A(m-3,n-2*bw+14) * A(m+2,n-2*bw+11)&
         - A(m-4,n-2*bw+14) * A(m+3,n-2*bw+10) - A(m-5,n-2*bw+14) * A(m+4,n-2*bw+9) - A(m-6,n-2*bw+14) * A(m+5,n-2*bw+8)&
         - A(m-7,n-2*bw+14) * A(m+6,n-2*bw+7) - A(m+1,n-2*bw+14) * A(m+7,n-2*bw+6)

    dummy = 1.d0 / A(m,n-2*bw+13)
    A(m+1,n-2*bw+13) = dummy * (A(m+1,n-2*bw+13) - A(m-1,n-2*bw+13) * A(m+2,n-2*bw+12) - A(m-2,n-2*bw+13) * A(m+3,n-2*bw+11)&
         - A(m-3,n-2*bw+13) * A(m+4,n-2*bw+10) - A(m-4,n-2*bw+13) * A(m+5,n-2*bw+9) - A(m-5,n-2*bw+13) * A(m+6,n-2*bw+8) &
         - A(m-6,n-2*bw+13) * A(m+7,n-2*bw+7))

    A(m  ,n-2*bw+14)  = A(m  ,n-2*bw+14) - A(m-1,n-2*bw+14) * A(m+1,n-2*bw+13) - A(m-2,n-2*bw+14) * A(m+2,n-2*bw+12)&
         - A(m-3,n-2*bw+14) * A(m+3,n-2*bw+11) - A(m-4,n-2*bw+14) * A(m+4,n-2*bw+10) - A(m-5,n-2*bw+14) * A(m+5,n-2*bw+9)&
         - A(m-6,n-2*bw+14) * A(m+6,n-2*bw+8) - A(m-7,n-2*bw+14) * A(m+7,n-2*bw+7)    
    return

  end subroutine fixed_bandec7


  subroutine fixed_bandec8(A,n)
    implicit none
    integer :: n
    real(c_double) :: A(17,n)
    integer :: i,j,k
    real(c_double) :: dummy
    integer, parameter :: m = 9
    integer, parameter :: bw = 8

    ! Upper region
    do i=1,bw
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j=1,bw+1-i
          do k=1,i-1
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)           
          enddo
          A(m+j,i) = dummy * A(m+j,i)
       enddo
       do j=bw+2-i,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          do k=bw-j+1,i-1
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+bw-i-j+1,i-k)
          enddo
          A(m+j,i) = dummy * A(m+j,i)
       enddo
       do j=i,bw
          do k=1,i-1
             A(m-j,i) = A(m-j,i) - A(m-k,i)*A(m-j,i-k)           
          enddo
          A(m-j,i) = dummy * A(m-j,i)
       enddo
       ! Downward direction
       do j=1,bw+1-i
          do k=1,i
             A(m+1-j,i+j) = A(m+1-j,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
       do j=bw+2-i,bw
          do k=1,bw+1-j
             A(m+1-j,i+j) = A(m+1-j,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo

    ! Subupper Region
    do i=m,2*bw
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j=1,2*bw+1-i
          do k=1, bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          do k=bw-j+1, bw
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(2*m-i-j,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       do j=2*bw+2-i,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,bw
          do k=1,bw+1-j
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m-j+1-k,i+j)*A(m+k,i-k+1)
          enddo
       enddo
    enddo



    !center region
    do i = 2*bw+1,n-2*bw-1
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j=1,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,bw
          do k=1,bw-j+1
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo



    ! sublower region
    do i = n-2*bw, n-bw-1
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j = 1,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,n-bw-i
          do k=1,bw-j+1
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
       do j=n-bw-i+1,bw+1
          do k=1,bw+1-j
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
          do k=bw+2-j,i-n+2*bw+1 
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+n-bw-1+k-i,i+j)*A(m+k,i+1-k)
          enddo
       enddo
       do j=bw+2,n-i
          do k=1,i-n+2*bw+1
             A(m+n-bw-1-i,i+j) = A(m+n-bw-1-i,i+j) - A(m+n-bw-1-i+k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo

    ! lower region
    do i = n-bw, n-1
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j = 1,n-i
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,n-i
          do k=1,1+bw-j
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
          do k=2+bw-j,bw
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+n-bw-1-i+k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo
    return

  end subroutine fixed_bandec8


  subroutine fixed_bandec9(A,n)
    implicit none
    integer :: n
    real(c_double) :: A(19,n)
    integer :: i,j,k
    real(c_double) :: dummy
    integer, parameter :: m = 10
    integer, parameter :: bw = 9

    ! Upper region
    do i=1,bw
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j=1,bw+1-i
          do k=1,i-1
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)           
          enddo
          A(m+j,i) = dummy * A(m+j,i)
       enddo
       do j=bw+2-i,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          do k=bw-j+1,i-1
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+bw-i-j+1,i-k)
          enddo
          A(m+j,i) = dummy * A(m+j,i)
       enddo
       do j=i,bw
          do k=1,i-1
             A(m-j,i) = A(m-j,i) - A(m-k,i)*A(m-j,i-k)           
          enddo
          A(m-j,i) = dummy * A(m-j,i)
       enddo
       ! Downward direction
       do j=1,bw+1-i
          do k=1,i
             A(m+1-j,i+j) = A(m+1-j,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
       do j=bw+2-i,bw
          do k=1,bw+1-j
             A(m+1-j,i+j) = A(m+1-j,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo

    ! Subupper Region
    do i=m,2*bw
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j=1,2*bw+1-i
          do k=1, bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          do k=bw-j+1, bw
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(2*m-i-j,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       do j=2*bw+2-i,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,bw
          do k=1,bw+1-j
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m-j+1-k,i+j)*A(m+k,i-k+1)
          enddo
       enddo
    enddo



    !center region
    do i = 2*bw+1,n-2*bw-1
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j=1,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,bw
          do k=1,bw-j+1
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo



    ! sublower region
    do i = n-2*bw, n-bw-1
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j = 1,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,n-bw-i
          do k=1,bw-j+1
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
       do j=n-bw-i+1,bw+1
          do k=1,bw+1-j
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
          do k=bw+2-j,i-n+2*bw+1 
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+n-bw-1+k-i,i+j)*A(m+k,i+1-k)
          enddo
       enddo
       do j=bw+2,n-i
          do k=1,i-n+2*bw+1
             A(m+n-bw-1-i,i+j) = A(m+n-bw-1-i,i+j) - A(m+n-bw-1-i+k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo

    ! lower region
    do i = n-bw, n-1
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j = 1,n-i
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,n-i
          do k=1,1+bw-j
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
          do k=2+bw-j,bw
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+n-bw-1-i+k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo
    return

  end subroutine fixed_bandec9


  subroutine fixed_bandec10(A,n)
    implicit none
    integer :: n
    real(c_double) :: A(21,n)
    integer :: i,j,k
    real(c_double) :: dummy
    integer, parameter :: m = 11
    integer, parameter :: bw = 10

    ! Upper region
    do i=1,bw
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j=1,bw+1-i
          do k=1,i-1
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)           
          enddo
          A(m+j,i) = dummy * A(m+j,i)
       enddo
       do j=bw+2-i,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          do k=bw-j+1,i-1
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+bw-i-j+1,i-k)
          enddo
          A(m+j,i) = dummy * A(m+j,i)
       enddo
       do j=i,bw
          do k=1,i-1
             A(m-j,i) = A(m-j,i) - A(m-k,i)*A(m-j,i-k)           
          enddo
          A(m-j,i) = dummy * A(m-j,i)
       enddo
       ! Downward direction
       do j=1,bw+1-i
          do k=1,i
             A(m+1-j,i+j) = A(m+1-j,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
       do j=bw+2-i,bw
          do k=1,bw+1-j
             A(m+1-j,i+j) = A(m+1-j,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo

    ! Subupper Region
    do i=m,2*bw
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j=1,2*bw+1-i
          do k=1, bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          do k=bw-j+1, bw
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(2*m-i-j,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       do j=2*bw+2-i,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,bw
          do k=1,bw+1-j
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m-j+1-k,i+j)*A(m+k,i-k+1)
          enddo
       enddo
    enddo



    !center region
    do i = 2*bw+1,n-2*bw-1
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j=1,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,bw
          do k=1,bw-j+1
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo



    ! sublower region
    do i = n-2*bw, n-bw-1
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j = 1,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,n-bw-i
          do k=1,bw-j+1
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
       do j=n-bw-i+1,bw+1
          do k=1,bw+1-j
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
          do k=bw+2-j,i-n+2*bw+1 
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+n-bw-1+k-i,i+j)*A(m+k,i+1-k)
          enddo
       enddo
       do j=bw+2,n-i
          do k=1,i-n+2*bw+1
             A(m+n-bw-1-i,i+j) = A(m+n-bw-1-i,i+j) - A(m+n-bw-1-i+k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo

    ! lower region
    do i = n-bw, n-1
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j = 1,n-i
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,n-i
          do k=1,1+bw-j
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
          do k=2+bw-j,bw
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+n-bw-1-i+k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo
    return

  end subroutine fixed_bandec10


  subroutine fixed_bandec(A,n,bw)
    implicit none
    integer :: n,bw
    real(c_double) :: A(2*bw+1,n)
    integer :: i,j,k,m
    real(c_double) :: dummy

    m = bw+1

    ! Upper region
    do i=1,bw
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j=1,bw+1-i
          do k=1,i-1
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)           
          enddo
          A(m+j,i) = dummy * A(m+j,i)
       enddo
       do j=bw+2-i,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          do k=bw-j+1,i-1
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+bw-i-j+1,i-k)
          enddo
          A(m+j,i) = dummy * A(m+j,i)
       enddo
       do j=i,bw
          do k=1,i-1
             A(m-j,i) = A(m-j,i) - A(m-k,i)*A(m-j,i-k)           
          enddo
          A(m-j,i) = dummy * A(m-j,i)
       enddo
       ! Downward direction
       do j=1,bw+1-i
          do k=1,i
             A(m+1-j,i+j) = A(m+1-j,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
       do j=bw+2-i,bw
          do k=1,bw+1-j
             A(m+1-j,i+j) = A(m+1-j,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo

    ! Subupper Region
    do i=m,2*bw
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j=1,2*bw+1-i
          do k=1, bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          do k=bw-j+1, bw
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(2*m-i-j,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       do j=2*bw+2-i,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,bw
          do k=1,bw+1-j
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m-j+1-k,i+j)*A(m+k,i-k+1)
          enddo
       enddo
    enddo



    !center region
    do i = 2*bw+1,n-2*bw-1
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j=1,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,bw
          do k=1,bw-j+1
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo



    ! sublower region
    do i = n-2*bw, n-bw-1
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j = 1,bw
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,n-bw-i
          do k=1,bw-j+1
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
       do j=n-bw-i+1,bw+1
          do k=1,bw+1-j
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
          do k=bw+2-j,i-n+2*bw+1 
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+n-bw-1+k-i,i+j)*A(m+k,i+1-k)
          enddo
       enddo
       do j=bw+2,n-i
          do k=1,i-n+2*bw+1
             A(m+n-bw-1-i,i+j) = A(m+n-bw-1-i,i+j) - A(m+n-bw-1-i+k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo

    ! lower region
    do i = n-bw, n-1
       ! right direction
       dummy = 1.d0 / A(m,i)
       do j = 1,n-i
          do k=1,bw-j
             A(m+j,i) = A(m+j,i) - A(m-k,i)*A(m+j+k,i-k)
          enddo
          A(m+j,i) = dummy*A(m+j,i)
       enddo
       ! downward direction
       do j=1,n-i
          do k=1,1+bw-j
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+1-j-k,i+j)*A(m+k,i+1-k)
          enddo
          do k=2+bw-j,bw
             A(m-j+1,i+j) = A(m-j+1,i+j) - A(m+n-bw-1-i+k,i+j)*A(m+k,i+1-k)
          enddo
       enddo
    enddo
    return

  end subroutine fixed_bandec



!------------------------------------------------------------------------------------------------------



  subroutine fixed_banbks1(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(3,n)
    complex(c_double_complex) :: b(n)
    integer :: m, bw
    integer :: i
    m = 2
    bw = 1

    ! dummy = L^-1 b  -- Upper part
    b(1) = b(1) / A(m,1) 

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       b(i) = (b(i) - b(i-bw)*A(m-1,i)) / A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    b(n  ) = (b(n) - b(n-2)*A(m+1,n  ) - b(n-1)*A(m-1,n  )) / A(m,n  )

!!$  ! now dummy is in b
!!$  ! f = U^-1 dummy

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       b(i) = b(i) - b(i+1)*A(m+1,i) 
    enddo

    !b(n) = b(n) -- upper part

    b(1) = b(1) - b(2)*A(m+1,1) - b(3)*A(m-1,1)  

    return
  end subroutine fixed_banbks1

  subroutine fixed_banbks2(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(5,n)
    complex(c_double_complex) :: b(n)
    integer :: m, bw
    integer :: i
    m = 3
    bw = 2

    ! dummy = L^-1 b  -- Upper part
    b(1) = b(1) / A(m,1) 
    b(2) = (b(2) - b(1)*A(m-1,2)) / A(m,2)

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       b(i) = (b(i) - b(i-bw)*A(m-2,i) - b(i-bw+1)*A(m-1,i)) / A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    b(n-1) = (b(n-1) - b(n-4)*A(m+2,n-1) - b(n-3)*A(m-2,n-1) - b(n-2)*A(m-1,n-1)) / A(m,n-1)
    b(n  ) = (b(n  ) - b(n-4)*A(m+2,n  ) - b(n-3)*A(m+1,n  ) - b(n-2)*A(m-2,n  ) - b(n-1)*A(m-1,n  )) / A(m,n  )

!!$  ! now dummy is in b
!!$  ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    b(n-1) = b(n-1) - b(n)  *A(m+1,n-1)

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       b(i) = b(i) - b(i+1)*A(m+1,i) - b(i+2)*A(m+2,i) 
    enddo

    !b(n) = b(n) -- upper part

    b(2) = b(2) - b(3)*A(m+1,2) - b(4)*A(m+2,2) -  b(5)*A(m-2,2) 
    b(1) = b(1) - b(2)*A(m+1,1) - b(3)*A(m+2,1) -  b(4)*A(m-1,1) - b(5)*A(m-2,1) 
    return
  end subroutine fixed_banbks2

  subroutine fixed_banbks3(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(7,n)
    complex(c_double_complex) :: b(n)
    integer :: m, bw
    integer :: i
    m = 4
    bw = 3

    ! dummy = L^-1 b  -- Upper part
    b(1) = b(1) / A(m,1) 
    b(2) = (b(2) - b(1)*A(m-1,2)) / A(m,2)
    b(3) = (b(3) - b(1)*A(m-2,3) - b(2)*A(m-1,3)) / A(m,3)

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       b(i) = (b(i) - b(i-bw)*A(m-3,i) - b(i-bw+1)*A(m-2,i) - b(i-bw+2)*A(m-1,i)) / A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    b(n-2) = (b(n-2) - b(n-6)*A(m+3,n-2) - b(n-5)*A(m-3,n-2) - b(n-4)*A(m-2,n-2) - b(n-3)*A(m-1,n-2)) / A(m,n-2)
    b(n-1) = (b(n-1) - b(n-6)*A(m+3,n-1) - b(n-5)*A(m+2,n-1) - b(n-4)*A(m-3,n-1) - b(n-3)*A(m-2,n-1) &
         - b(n-2)*A(m-1,n-1)) / A(m,n-1)
    b(n  ) = (b(n  ) - b(n-6)*A(m+3,n  ) - b(n-5)*A(m+2,n  ) - b(n-4)*A(m+1,n  ) - b(n-3)*A(m-3,n  ) &
         - b(n-2)*A(m-2,n  ) - b(n-1)*A(m-1,n  )) / A(m,n  )

!!$  ! now dummy is in b
!!$  ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    b(n-1) = b(n-1) - b(n)  *A(m+1,n-1)
    b(n-2) = b(n-2) - b(n-1)*A(m+1,n-2) - b(n)  *A(m+2,n-2)

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       b(i) = b(i) - b(i+1)*A(m+1,i) - b(i+2)*A(m+2,i) - b(i+3)*A(m+3,i)
    enddo

    !b(n) = b(n) -- upper part

    b(3) = b(3) - b(4)*A(m+1,3) - b(5)*A(m+2,3) -  b(6)*A(m+3,3) -  b(7)*A(m-3,3)
    b(2) = b(2) - b(3)*A(m+1,2) - b(4)*A(m+2,2) -  b(5)*A(m+3,2) -  b(6)*A(m-2,2) -  b(7)*A(m-3,2) 
    b(1) = b(1) - b(2)*A(m+1,1) - b(3)*A(m+2,1) -  b(4)*A(m+3,1) -  b(5)*A(m-1,1) -  b(6)*A(m-2,1) -  b(7)*A(m-3,1) 

    return
  end subroutine fixed_banbks3

  subroutine fixed_banbks4(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(9,n)
    complex(c_double_complex) :: b(n)
!!$  integer :: m, bw
    integer :: i
!!$  m = 5
!!$  bw = 4

    ! dummy = L^-1 b  -- Upper part
    b(1) = b(1) / A(5,1) 
    b(2) = (b(2) - b(1)*A(4,2)) / A(5,2)
    b(3) = (b(3) - b(1)*A(3,3) - b(2)*A(4,3)) / A(5,3)
    b(4) = (b(4) - b(1)*A(2,4) - b(2)*A(3,4) - b(3)*A(4,4)) / A(5,4)

    ! dummy =  L^-1 b  -- Center part
    do i=5,n-4    
       b(i) = (b(i) - b(i-4)*A(1,i) - b(i-3)*A(2,i) - b(i-2)*A(3,i) - b(i-1)*A(4,i)) / A(5,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    b(n-3) = (b(n-3) - b(n-8)*A(9,n-3) - b(n-7)*A(1,n-3) - b(n-6)*A(2,n-3) - b(n-5)*A(3,n-3) - b(n-4)*A(4,n-3))&
         / A(5,n-3)
    b(n-2) = (b(n-2) - b(n-8)*A(9,n-2) - b(n-7)*A(8,n-2) - b(n-6)*A(1,n-2) - b(n-5)*A(2,n-2) - b(n-4)*A(3,n-2)&
         - b(n-3)*A(4,n-2)) / A(5,n-2)
    b(n-1) = (b(n-1) - b(n-8)*A(9,n-1) - b(n-7)*A(8,n-1) - b(n-6)*A(7,n-1) - b(n-5)*A(1,n-1) - b(n-4)*A(2,n-1)&
         - b(n-3)*A(3,n-1) - b(n-2)*A(4,n-1)) / A(5,n-1)
    b(n  ) = (b(n  ) - b(n-8)*A(9,n  ) - b(n-7)*A(8,n  ) - b(n-6)*A(7,n  ) - b(n-5)*A(6,n  ) - b(n-4)*A(1,n  ) &
         - b(n-3)*A(2,n  ) - b(n-2)*A(3,n  ) - b(n-1)*A(4,n  )) / A(5,n  )

!!$  ! now dummy is in b
!!$  ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    b(n-1) = b(n-1) - b(n)  *A(6,n-1)
    b(n-2) = b(n-2) - b(n-1)*A(6,n-2) - b(n)  *A(7,n-2)
    b(n-3) = b(n-3) - b(n-2)*A(6,n-3) - b(n-1)*A(7,n-3) - b(n)  *A(8,n-3)

    !b(n) = b(n) -- center part
    do i=n-4,5,-1
       b(i) = b(i) - b(i+1)*A(6,i) - b(i+2)*A(7,i) - b(i+3)*A(8,i) - b(i+4)*A(9,i) 
    enddo

    !b(n) = b(n) -- upper part

    b(4) = b(4) - b(5)*A(6,4) - b(6)*A(7,4) -  b(7)*A(8,4) -  b(8)*A(9,4) -  b(9)*A(1,4)
    b(3) = b(3) - b(4)*A(6,3) - b(5)*A(7,3) -  b(6)*A(8,3) -  b(7)*A(9,3) -  b(8)*A(2,3) -  b(9)*A(1,3)
    b(2) = b(2) - b(3)*A(6,2) - b(4)*A(7,2) -  b(5)*A(8,2) -  b(6)*A(9,2) -  b(7)*A(3,2) -  b(8)*A(2,2)&
         -  b(9)*A(1,2)
    b(1) = b(1) - b(2)*A(6,1) - b(3)*A(7,1) -  b(4)*A(8,1) -  b(5)*A(9,1) -  b(6)*A(4,1) -  b(7)*A(3,1)&
         -  b(8)*A(2,1) - b(9)*A(1,1)

    return
  end subroutine fixed_banbks4

  subroutine fixed_banbks5(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(11,n)
    complex(c_double_complex) :: b(n)
    integer :: m, bw
    integer :: i
    m = 6
    bw = 5

    ! dummy = L^-1 b  -- Upper part
    b(1) = b(1) / A(m,1) 
    b(2) = (b(2) - b(1)*A(m-1,2)) / A(m,2)
    b(3) = (b(3) - b(1)*A(m-2,3) - b(2)*A(m-1,3)) / A(m,3)
    b(4) = (b(4) - b(1)*A(m-3,4) - b(2)*A(m-2,4) - b(3)*A(m-1,4)) / A(m,4)
    b(5) = (b(5) - b(1)*A(m-4,5) - b(2)*A(m-3,5) - b(3)*A(m-2,5) - b(4)*A(m-1,5)) / A(m,5)  

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       b(i) = (b(i) - b(i-bw)*A(m-5,i) - b(i-bw+1)*A(m-4,i) - b(i-bw+2)*A(m-3,i) - b(i-bw+3)*A(m-2,i) &
            - b(i-bw+4)*A(m-1,i)) / A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    b(n-4) = (b(n-4) - b(n-10)*A(m+5,n-4) - b(n-9)*A(m-5,n-4) - b(n-8)*A(m-4,n-4) - b(n-7)*A(m-3,n-4) &
         - b(n-6)*A(m-2,n-4) - b(n-5)*A(m-1,n-4)) / A(m,n-4)
    b(n-3) = (b(n-3) - b(n-10)*A(m+5,n-3) - b(n-9)*A(m+4,n-3) - b(n-8)*A(m-5,n-3) - b(n-7)*A(m-4,n-3) &
         - b(n-6)*A(m-3,n-3) - b(n-5)*A(m-2,n-3) - b(n-4)*A(m-1,n-3)) / A(m,n-3)
    b(n-2) = (b(n-2) - b(n-10)*A(m+5,n-2) - b(n-9)*A(m+4,n-2) - b(n-8)*A(m+3,n-2) - b(n-7)*A(m-5,n-2)&
         - b(n-6)*A(m-4,n-2) - b(n-5)*A(m-3,n-2) - b(n-4)*A(m-2,n-2) - b(n-3)*A(m-1,n-2)) / A(m,n-2)
    b(n-1) = (b(n-1) - b(n-10)*A(m+5,n-1) - b(n-9)*A(m+4,n-1) - b(n-8)*A(m+3,n-1) - b(n-7)*A(m+2,n-1) &
         - b(n-6)*A(m-5,n-1) - b(n-5)*A(m-4,n-1) - b(n-4)*A(m-3,n-1) - b(n-3)*A(m-2,n-1) - b(n-2)*A(m-1,n-1)) / A(m,n-1)
    b(n  ) = (b(n  ) - b(n-10)*A(m+5,n  ) - b(n-9)*A(m+4,n  ) - b(n-8)*A(m+3,n  ) - b(n-7)*A(m+2,n  ) &
         - b(n-6)*A(m+1,n  ) - b(n-5)*A(m-5,n  ) - b(n-4)*A(m-4,n  ) - b(n-3)*A(m-3,n  ) - b(n-2)*A(m-2,n  ) &
         - b(n-1)*A(m-1,n  )) / A(m,n  )

!!$  ! now dummy is in b
!!$  ! f = U^-1 dummy


    !b(n) = b(n) -- Lower part
    b(n-1) = b(n-1) - b(n)  *A(m+1,n-1)
    b(n-2) = b(n-2) - b(n-1)*A(m+1,n-2) - b(n)  *A(m+2,n-2)
    b(n-3) = b(n-3) - b(n-2)*A(m+1,n-3) - b(n-1)*A(m+2,n-3) - b(n)  *A(m+3,n-3)
    b(n-4) = b(n-4) - b(n-3)*A(m+1,n-4) - b(n-2)*A(m+2,n-4) - b(n-1)*A(m+3,n-4) - b(n)  *A(m+4,n-4)

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       b(i) = b(i) - b(i+1)*A(m+1,i) - b(i+2)*A(m+2,i) - b(i+3)*A(m+3,i) - b(i+4)*A(m+4,i) - b(i+5)*A(m+5,i) 
    enddo

    !b(n) = b(n) -- upper part

    b(5) = b(5) - b(6)*A(m+1,5) - b(7)*A(m+2,5) -  b(8)*A(m+3,5) -  b(9)*A(m+4,5) - b(10)*A(m+5,5)&
         - b(11)*A(m-5,5)
    b(4) = b(4) - b(5)*A(m+1,4) - b(6)*A(m+2,4) -  b(7)*A(m+3,4) -  b(8)*A(m+4,4) -  b(9)*A(m+5,4)&
         - b(10)*A(m-4,4) - b(11)*A(m-5,4) 
    b(3) = b(3) - b(4)*A(m+1,3) - b(5)*A(m+2,3) -  b(6)*A(m+3,3) -  b(7)*A(m+4,3) -  b(8)*A(m+5,3)&
         -  b(9)*A(m-3,3) - b(10)*A(m-4,3) -b(11)*A(m-5,3) 
    b(2) = b(2) - b(3)*A(m+1,2) - b(4)*A(m+2,2) -  b(5)*A(m+3,2) -  b(6)*A(m+4,2) -  b(7)*A(m+5,2)&
         -  b(8)*A(m-2,2) -  b(9)*A(m-3,2) -b(10)*A(m-4,2) -b(11)*A(m-5,2) 
    b(1) = b(1) - b(2)*A(m+1,1) - b(3)*A(m+2,1) -  b(4)*A(m+3,1) -  b(5)*A(m+4,1) -  b(6)*A(m+5,1) &
         -  b(7)*A(m-1,1) -  b(8)*A(m-2,1) - b(9)*A(m-3,1) -b(10)*A(m-4,1) -b(11)*A(m-5,1) 

    return
  end subroutine fixed_banbks5

  subroutine fixed_banbks6(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(13,n)
    complex(c_double_complex) :: b(n)
    integer :: m, bw
    integer :: i
    m = 7
    bw = 6

    ! dummy = L^-1 b  -- Upper part
    b(1) = b(1) / A(m,1) 
    b(2) = (b(2) - b(1)*A(m-1,2)) / A(m,2)
    b(3) = (b(3) - b(1)*A(m-2,3) - b(2)*A(m-1,3)) / A(m,3)
    b(4) = (b(4) - b(1)*A(m-3,4) - b(2)*A(m-2,4) - b(3)*A(m-1,4)) / A(m,4)
    b(5) = (b(5) - b(1)*A(m-4,5) - b(2)*A(m-3,5) - b(3)*A(m-2,5) - b(4)*A(m-1,5)) / A(m,5)  
    b(6) = (b(6) - b(1)*A(m-5,6) - b(2)*A(m-4,6) - b(3)*A(m-3,6) - b(4)*A(m-2,6) - b(5)*A(m-1,6)) / A(m,6)  

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       b(i) = (b(i) - b(i-bw)*A(m-6,i) - b(i-bw+1)*A(m-5,i) - b(i-bw+2)*A(m-4,i) - b(i-bw+3)*A(m-3,i) - b(i-bw+4)*A(m-2,i)&
            - b(i-bw+5)*A(m-1,i)) / A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    b(n-5) = (b(n-5) - b(n-12)*A(m+6,n-5) - b(n-11)*A(m-6,n-5) - b(n-10)*A(m-5,n-5) - b(n-9)*A(m-4,n-5)&
         - b(n-8)*A(m-3,n-5) - b(n-7)*A(m-2,n-5) - b(n-6)*A(m-1,n-5)) / A(m,n-5)
    b(n-4) = (b(n-4) - b(n-12)*A(m+6,n-4) - b(n-11)*A(m+5,n-4) - b(n-10)*A(m-6,n-4) - b(n-9)*A(m-5,n-4) &
         - b(n-8)*A(m-4,n-4) - b(n-7)*A(m-3,n-4) - b(n-6)*A(m-2,n-4) - b(n-5)*A(m-1,n-4)) / A(m,n-4)
    b(n-3) = (b(n-3) - b(n-12)*A(m+6,n-3) - b(n-11)*A(m+5,n-3) - b(n-10)*A(m+4,n-3) - b(n-9)*A(m-6,n-3)&
         - b(n-8)*A(m-5,n-3) - b(n-7)*A(m-4,n-3) - b(n-6)*A(m-3,n-3) - b(n-5)*A(m-2,n-3) - b(n-4)*A(m-1,n-3)) / A(m,n-3)
    b(n-2) = (b(n-2) - b(n-12)*A(m+6,n-2) - b(n-11)*A(m+5,n-2) - b(n-10)*A(m+4,n-2) - b(n-9)*A(m+3,n-2) &
         - b(n-8)*A(m-6,n-2) - b(n-7)*A(m-5,n-2) - b(n-6)*A(m-4,n-2) - b(n-5)*A(m-3,n-2) - b(n-4)*A(m-2,n-2) - b(n-3)*A(m-1,n-2)) &
         / A(m,n-2)
    b(n-1) = (b(n-1) - b(n-12)*A(m+6,n-1) - b(n-11)*A(m+5,n-1) - b(n-10)*A(m+4,n-1) - b(n-9)*A(m+3,n-1) &
         - b(n-8)*A(m+2,n-1) - b(n-7)*A(m-6,n-1) - b(n-6)*A(m-5,n-1) - b(n-5)*A(m-4,n-1) - b(n-4)*A(m-3,n-1) - b(n-3)*A(m-2,n-1)&
         - b(n-2)*A(m-1,n-1)) / A(m,n-1)
    b(n  ) = (b(n  ) - b(n-12)*A(m+6,n  ) - b(n-11)*A(m+5,n  ) - b(n-10)*A(m+4,n  ) - b(n-9)*A(m+3,n  ) &
         - b(n-8)*A(m+2,n  ) - b(n-7)*A(m+1,n  ) - b(n-6)*A(m-6,n  ) - b(n-5)*A(m-5,n  ) - b(n-4)*A(m-4,n  ) - b(n-3)*A(m-3,n  )&
         - b(n-2)*A(m-2,n  ) - b(n-1)*A(m-1,n  )) / A(m,n  )

!!$  ! now dummy is in b
!!$  ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    b(n-1) = b(n-1) - b(n)  *A(m+1,n-1)
    b(n-2) = b(n-2) - b(n-1)*A(m+1,n-2) - b(n)  *A(m+2,n-2)
    b(n-3) = b(n-3) - b(n-2)*A(m+1,n-3) - b(n-1)*A(m+2,n-3) - b(n)  *A(m+3,n-3)
    b(n-4) = b(n-4) - b(n-3)*A(m+1,n-4) - b(n-2)*A(m+2,n-4) - b(n-1)*A(m+3,n-4) - b(n)  *A(m+4,n-4)
    b(n-5) = b(n-5) - b(n-4)*A(m+1,n-5) - b(n-3)*A(m+2,n-5) - b(n-2)*A(m+3,n-5) - b(n-1)*A(m+4,n-5) - b(n)  *A(m+5,n-5)

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       b(i) = b(i) - b(i+1)*A(m+1,i) - b(i+2)*A(m+2,i) - b(i+3)*A(m+3,i) - b(i+4)*A(m+4,i) - b(i+5)*A(m+5,i) - b(i+6)*A(m+6,i)
    enddo

    !b(n) = b(n) -- upper part

    b(6) = b(6) - b(7)*A(m+1,6) - b(8)*A(m+2,6) -  b(9)*A(m+3,6) - b(10)*A(m+4,6) - b(11)*A(m+5,6) - b(12)*A(m+6,6)&
         - b(13)*A(m-6,6)
    b(5) = b(5) - b(6)*A(m+1,5) - b(7)*A(m+2,5) -  b(8)*A(m+3,5) -  b(9)*A(m+4,5) - b(10)*A(m+5,5) - b(11)*A(m+6,5)&
         - b(12)*A(m-5,5) -b(13)*A(m-6,5) 
    b(4) = b(4) - b(5)*A(m+1,4) - b(6)*A(m+2,4) -  b(7)*A(m+3,4) -  b(8)*A(m+4,4) -  b(9)*A(m+5,4) - b(10)*A(m+6,4)&
         - b(11)*A(m-4,4) -b(12)*A(m-5,4) -b(13)*A(m-6,4) 
    b(3) = b(3) - b(4)*A(m+1,3) - b(5)*A(m+2,3) -  b(6)*A(m+3,3) -  b(7)*A(m+4,3) -  b(8)*A(m+5,3) -  b(9)*A(m+6,3)&
         - b(10)*A(m-3,3) -b(11)*A(m-4,3) -b(12)*A(m-5,3) -b(13)*A(m-6,3) 
    b(2) = b(2) - b(3)*A(m+1,2) - b(4)*A(m+2,2) -  b(5)*A(m+3,2) -  b(6)*A(m+4,2) -  b(7)*A(m+5,2) -  b(8)*A(m+6,2)&
         -  b(9)*A(m-2,2) -b(10)*A(m-3,2) -b(11)*A(m-4,2) -b(12)*A(m-5,2) -b(13)*A(m-6,2)
    b(1) = b(1) - b(2)*A(m+1,1) - b(3)*A(m+2,1) -  b(4)*A(m+3,1) -  b(5)*A(m+4,1) -  b(6)*A(m+5,1) -  b(7)*A(m+6,1)&
         -  b(8)*A(m-1,1) - b(9)*A(m-2,1) -b(10)*A(m-3,1) -b(11)*A(m-4,1) -b(12)*A(m-5,1) -b(13)*A(m-6,1) 

    return
  end subroutine fixed_banbks6

  subroutine fixed_banbks7(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(15,n)
    complex(c_double_complex) :: b(n)
    integer :: m, bw
    integer :: i
    m = 8
    bw = 7

    ! dummy = L^-1 b  -- Upper part
    b(1) = b(1) / A(m,1) 
    b(2) = (b(2) - b(1)*A(m-1,2)) / A(m,2)
    b(3) = (b(3) - b(1)*A(m-2,3) - b(2)*A(m-1,3)) / A(m,3)
    b(4) = (b(4) - b(1)*A(m-3,4) - b(2)*A(m-2,4) - b(3)*A(m-1,4)) / A(m,4)
    b(5) = (b(5) - b(1)*A(m-4,5) - b(2)*A(m-3,5) - b(3)*A(m-2,5) - b(4)*A(m-1,5)) / A(m,5)  
    b(6) = (b(6) - b(1)*A(m-5,6) - b(2)*A(m-4,6) - b(3)*A(m-3,6) - b(4)*A(m-2,6) - b(5)*A(m-1,6)) / A(m,6)  
    b(7) = (b(7) - b(1)*A(m-6,7) - b(2)*A(m-5,7) - b(3)*A(m-4,7) - b(4)*A(m-3,7) - b(5)*A(m-2,7) - b(6)*A(m-1,7)) / A(m,7)  

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       b(i) = (b(i) - b(i-bw)*A(m-7,i) - b(i-bw+1)*A(m-6,i) - b(i-bw+2)*A(m-5,i) - b(i-bw+3)*A(m-4,i) - b(i-bw+4)*A(m-3,i) &
            - b(i-bw+5)*A(m-2,i) - b(i-bw+6)*A(m-1,i)) / A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    b(n-6) = (b(n-6) - b(n-14)*A(m+7,n-6) - b(n-13)*A(m-7,n-6) - b(n-12)*A(m-6,n-6) - b(n-11)*A(m-5,n-6) - b(n-10)*A(m-4,n-6)&
         - b(n-9)*A(m-3,n-6) - b(n-8)*A(m-2,n-6) - b(n-7)*A(m-1,n-6)) / A(m,n-6)
    b(n-5) = (b(n-5) - b(n-14)*A(m+7,n-5) - b(n-13)*A(m+6,n-5) - b(n-12)*A(m-7,n-5) - b(n-11)*A(m-6,n-5) - b(n-10)*A(m-5,n-5)&
         - b(n-9)*A(m-4,n-5) - b(n-8)*A(m-3,n-5) - b(n-7)*A(m-2,n-5) - b(n-6)*A(m-1,n-5)) / A(m,n-5)
    b(n-4) = (b(n-4) - b(n-14)*A(m+7,n-4) - b(n-13)*A(m+6,n-4) - b(n-12)*A(m+5,n-4) - b(n-11)*A(m-7,n-4) - b(n-10)*A(m-6,n-4)&
         - b(n-9)*A(m-5,n-4) - b(n-8)*A(m-4,n-4) - b(n-7)*A(m-3,n-4) - b(n-6)*A(m-2,n-4) - b(n-5)*A(m-1,n-4)) / A(m,n-4)
    b(n-3) = (b(n-3) - b(n-14)*A(m+7,n-3) - b(n-13)*A(m+6,n-3) - b(n-12)*A(m+5,n-3) - b(n-11)*A(m+4,n-3) - b(n-10)*A(m-7,n-3)&
         - b(n-9)*A(m-6,n-3) - b(n-8)*A(m-5,n-3) - b(n-7)*A(m-4,n-3) - b(n-6)*A(m-3,n-3) - b(n-5)*A(m-2,n-3) - b(n-4)*A(m-1,n-3))&
         / A(m,n-3)
    b(n-2) = (b(n-2) - b(n-14)*A(m+7,n-2) - b(n-13)*A(m+6,n-2) - b(n-12)*A(m+5,n-2) - b(n-11)*A(m+4,n-2) - b(n-10)*A(m+3,n-2)&
         - b(n-9)*A(m-7,n-2) - b(n-8)*A(m-6,n-2) - b(n-7)*A(m-5,n-2) - b(n-6)*A(m-4,n-2) - b(n-5)*A(m-3,n-2) - b(n-4)*A(m-2,n-2)&
         - b(n-3)*A(m-1,n-2)) / A(m,n-2)
    b(n-1) = (b(n-1) - b(n-14)*A(m+7,n-1) - b(n-13)*A(m+6,n-1) - b(n-12)*A(m+5,n-1) - b(n-11)*A(m+4,n-1) - b(n-10)*A(m+3,n-1)&
         - b(n-9)*A(m+2,n-1) - b(n-8)*A(m-7,n-1) - b(n-7)*A(m-6,n-1) - b(n-6)*A(m-5,n-1) - b(n-5)*A(m-4,n-1) - b(n-4)*A(m-3,n-1)&
         - b(n-3)*A(m-2,n-1) - b(n-2)*A(m-1,n-1)) / A(m,n-1)
    b(n  ) = (b(n  ) - b(n-14)*A(m+7,n  ) - b(n-13)*A(m+6,n  ) - b(n-12)*A(m+5,n  ) - b(n-11)*A(m+4,n  ) - b(n-10)*A(m+3,n  )&
         - b(n-9)*A(m+2,n  ) - b(n-8)*A(m+1,n  ) - b(n-7)*A(m-7,n  ) - b(n-6)*A(m-6,n  ) - b(n-5)*A(m-5,n  ) - b(n-4)*A(m-4,n  )&
         - b(n-3)*A(m-3,n  ) - b(n-2)*A(m-2,n  ) - b(n-1)*A(m-1,n  )) / A(m,n  )

!!$  ! now dummy is in b
!!$  ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    b(n-1) = b(n-1) - b(n)  *A(m+1,n-1)
    b(n-2) = b(n-2) - b(n-1)*A(m+1,n-2) - b(n)  *A(m+2,n-2)
    b(n-3) = b(n-3) - b(n-2)*A(m+1,n-3) - b(n-1)*A(m+2,n-3) - b(n)  *A(m+3,n-3)
    b(n-4) = b(n-4) - b(n-3)*A(m+1,n-4) - b(n-2)*A(m+2,n-4) - b(n-1)*A(m+3,n-4) - b(n)  *A(m+4,n-4)
    b(n-5) = b(n-5) - b(n-4)*A(m+1,n-5) - b(n-3)*A(m+2,n-5) - b(n-2)*A(m+3,n-5) - b(n-1)*A(m+4,n-5) - b(n)  *A(m+5,n-5)
    b(n-6) = b(n-6) - b(n-5)*A(m+1,n-6) - b(n-4)*A(m+2,n-6) - b(n-3)*A(m+3,n-6) - b(n-2)*A(m+4,n-6) - b(n-1)*A(m+5,n-6)&
         - b(n)*A(m+6,n-6)

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       b(i) = b(i) - b(i+1)*A(m+1,i) - b(i+2)*A(m+2,i) - b(i+3)*A(m+3,i) - b(i+4)*A(m+4,i) - b(i+5)*A(m+5,i) - b(i+6)*A(m+6,i) &
            - b(i+7)*A(m+7,i)
    enddo

    !b(n) = b(n) -- upper part

    b(7) = b(7) - b(8)*A(m+1,7) - b(9)*A(m+2,7) - b(10)*A(m+3,7) - b(11)*A(m+4,7) - b(12)*A(m+5,7) - b(13)*A(m+6,7) &
         - b(14)*A(m+7,7) -b(15)*A(m-7,7)
    b(6) = b(6) - b(7)*A(m+1,6) - b(8)*A(m+2,6) -  b(9)*A(m+3,6) - b(10)*A(m+4,6) - b(11)*A(m+5,6) - b(12)*A(m+6,6) &
         - b(13)*A(m+7,6) -b(14)*A(m-6,6) -b(15)*A(m-7,6)
    b(5) = b(5) - b(6)*A(m+1,5) - b(7)*A(m+2,5) -  b(8)*A(m+3,5) - b(9)*A(m+4,5) - b(10)*A(m+5,5) - b(11)*A(m+6,5) &
         - b(12)*A(m+7,5) -b(13)*A(m-5,5) -b(14)*A(m-6,5) -b(15)*A(m-7,5)
    b(4) = b(4) - b(5)*A(m+1,4) - b(6)*A(m+2,4) -  b(7)*A(m+3,4) -  b(8)*A(m+4,4) -  b(9)*A(m+5,4) - b(10)*A(m+6,4) &
         - b(11)*A(m+7,4) -b(12)*A(m-4,4) -b(13)*A(m-5,4) -b(14)*A(m-6,4) -b(15)*A(m-7,4)
    b(3) = b(3) - b(4)*A(m+1,3) - b(5)*A(m+2,3) -  b(6)*A(m+3,3) -  b(7)*A(m+4,3) -  b(8)*A(m+5,3) -  b(9)*A(m+6,3) &
         - b(10)*A(m+7,3) -b(11)*A(m-3,3) -b(12)*A(m-4,3) -b(13)*A(m-5,3) -b(14)*A(m-6,3) -b(15)*A(m-7,3)
    b(2) = b(2) - b(3)*A(m+1,2) - b(4)*A(m+2,2) -  b(5)*A(m+3,2) -  b(6)*A(m+4,2) -  b(7)*A(m+5,2) -  b(8)*A(m+6,2) &
         -  b(9)*A(m+7,2) -b(10)*A(m-2,2) -b(11)*A(m-3,2) -b(12)*A(m-4,2) -b(13)*A(m-5,2) -b(14)*A(m-6,2) -b(15)*A(m-7,2)
    b(1) = b(1) - b(2)*A(m+1,1) - b(3)*A(m+2,1) -  b(4)*A(m+3,1) -  b(5)*A(m+4,1) -  b(6)*A(m+5,1) -  b(7)*A(m+6,1)&
         -  b(8)*A(m+7,1) - b(9)*A(m-1,1) -b(10)*A(m-2,1) -b(11)*A(m-3,1) -b(12)*A(m-4,1) -b(13)*A(m-5,1) -b(14)*A(m-6,1) &
         -b(15)*A(m-7,1)


    return
  end subroutine fixed_banbks7

  subroutine fixed_banbks8(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(17,n)
    complex(c_double_complex) :: b(n)
    integer :: i,j,k
    integer, parameter :: m = 9
    integer, parameter :: bw = 8

    ! dummy = L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i-1
          b(i) = b(i) - b(j)*A(m-i+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       do j=1,bw
          b(i) = b(i) - b(i-bw-1+j)*A(m-bw-1+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m+bw+1-j,n-bw+i)
       enddo
       do j=i+1,i+bw
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m-bw+j-i-1,n-bw+i)
       enddo
       b(n-bw+i) = b(n-bw+i)/A(m,n-bw+i)
    enddo


    ! now dummy is in b
    ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    do i=1,bw-1
       do j=1,i
          b(n-i) = b(n-i) - b(n-i+j)*A(m+j,n-i)
       enddo
    enddo

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       do j=1,bw
          b(i) = b(i) - b(i+j)*A(m+j,i)
       enddo
    enddo

    !b(n) = b(n) -- upper part
    do i=1,bw
       do j=1,bw
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+j,bw+1-i)
       enddo
       do j=bw+1,bw+i
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+i-j,bw+1-i)
       enddo
    enddo

    return
  end subroutine fixed_banbks8

  subroutine fixed_banbks9(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(19,n)
    complex(c_double_complex) :: b(n)
    integer :: i,j,k
    integer, parameter :: m = 10
    integer, parameter :: bw = 9

    ! dummy = L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i-1
          b(i) = b(i) - b(j)*A(m-i+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       do j=1,bw
          b(i) = b(i) - b(i-bw-1+j)*A(m-bw-1+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m+bw+1-j,n-bw+i)
       enddo
       do j=i+1,i+bw
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m-bw+j-i-1,n-bw+i)
       enddo
       b(n-bw+i) = b(n-bw+i)/A(m,n-bw+i)
    enddo


    ! now dummy is in b
    ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    do i=1,bw-1
       do j=1,i
          b(n-i) = b(n-i) - b(n-i+j)*A(m+j,n-i)
       enddo
    enddo

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       do j=1,bw
          b(i) = b(i) - b(i+j)*A(m+j,i)
       enddo
    enddo

    !b(n) = b(n) -- upper part
    do i=1,bw
       do j=1,bw
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+j,bw+1-i)
       enddo
       do j=bw+1,bw+i
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+i-j,bw+1-i)
       enddo
    enddo

    return
  end subroutine fixed_banbks9

  subroutine fixed_banbks10(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(21,n)
    complex(c_double_complex) :: b(n)
    integer :: i,j,k
    integer, parameter :: m = 11
    integer, parameter :: bw = 10

    ! dummy = L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i-1
          b(i) = b(i) - b(j)*A(m-i+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       do j=1,bw
          b(i) = b(i) - b(i-bw-1+j)*A(m-bw-1+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m+bw+1-j,n-bw+i)
       enddo
       do j=i+1,i+bw
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m-bw+j-i-1,n-bw+i)
       enddo
       b(n-bw+i) = b(n-bw+i)/A(m,n-bw+i)
    enddo


    ! now dummy is in b
    ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    do i=1,bw-1
       do j=1,i
          b(n-i) = b(n-i) - b(n-i+j)*A(m+j,n-i)
       enddo
    enddo

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       do j=1,bw
          b(i) = b(i) - b(i+j)*A(m+j,i)
       enddo
    enddo

    !b(n) = b(n) -- upper part
    do i=1,bw
       do j=1,bw
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+j,bw+1-i)
       enddo
       do j=bw+1,bw+i
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+i-j,bw+1-i)
       enddo
    enddo

    return
  end subroutine fixed_banbks10

  subroutine fixed_banbks(A,n,b,bw)
    implicit none
    integer :: n,bw
    real(c_double) :: A(2*bw+1,n)
    integer :: i,j,k,m
    complex(c_double_complex) :: b(n)
    real(c_double) :: dummy

    m = bw+1

    ! dummy = L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i-1
          b(i) = b(i) - b(j)*A(m-i+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       do j=1,bw
          b(i) = b(i) - b(i-bw-1+j)*A(m-bw-1+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m+bw+1-j,n-bw+i)
       enddo
       do j=i+1,i+bw
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m-bw+j-i-1,n-bw+i)
       enddo
       b(n-bw+i) = b(n-bw+i)/A(m,n-bw+i)
    enddo


    ! now dummy is in b
    ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    do i=1,bw-1
       do j=1,i
          b(n-i) = b(n-i) - b(n-i+j)*A(m+j,n-i)
       enddo
    enddo

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       do j=1,bw
          b(i) = b(i) - b(i+j)*A(m+j,i)
       enddo
    enddo

    !b(n) = b(n) -- upper part
    do i=1,bw
       do j=1,bw
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+j,bw+1-i)
       enddo
       do j=bw+1,bw+i
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+i-j,bw+1-i)
       enddo
    enddo

    return
  end subroutine fixed_banbks






  subroutine LABS_mat_mul_1(A,x,b,n)
    implicit none
    real(c_double) :: A(3,n)
    complex(c_double_complex) :: x(n), b(n)
    integer :: i,j,k,n

    ! Upper part
    b(1) = A(2,1)*x(1) + A(3,1)*x(2) + A(1,1)*x(3)     

    ! center part
    do j=2,n-1
       b(j) = + A(1,j)*x(j-1) + A(2,j)*x(j) + A(3,j)*x(j+1) 
    enddo

    ! Lower part
    b(n  ) = A(3,n  )*x(n-2) + A(1,n-1)*x(n-1) + A(2,n  )*x(n)

    return
  end subroutine LABS_mat_mul_1


  subroutine LABS_mat_mul_2(A,x,b,n)
    implicit none
    real(c_double) :: A(5,n)
    complex(c_double_complex) :: x(n), b(n)
    integer :: i,j,k,n

    ! Upper part
    b(1) = A(3,1)*x(1)+A(4,1)*x(2)+A(5,1)*x(3)+A(2,1)*x(4)+A(1,1)*x(5)    
    b(2) = A(2,2)*x(1)+A(3,2)*x(2)+A(4,2)*x(3)+A(5,2)*x(4)+A(1,2)*x(5)    

    ! center part
    do j=3,n-2
       b(j) = A(1,j)*x(j-2)+A(2,j)*x(j-1)+A(3,j)*x(j)+A(4,j)*x(j+1)+A(5,j)*x(j+2) 
    enddo

    ! Lower part
    b(n-1) = A(5,n-1)*x(n-4)+A(1,n-1)*x(n-3)+A(2,n-1)*x(n-2)+A(3,n-1)*x(n-1)+A(4,n-1)*x(n)
    b(n  ) = A(5,n  )*x(n-4)+A(4,n  )*x(n-3)+A(1,n  )*x(n-2)+A(2,n  )*x(n-1)+A(3,n  )*x(n)

    return
  end subroutine LABS_mat_mul_2



  subroutine LABS_mat_mul_3(A,x,b,n)
    implicit none
    real(c_double) :: A(7,n)
    complex(c_double_complex) :: x(n), b(n)
    integer :: i,j,k,n

    ! Upper part
    b(1) = A(4,1)*x(1)+A(5,1)*x(2)+A(6,1)*x(3)+A(7,1)*x(4)+A(3,1)*x(5)+A(2,1)*x(6)+A(1,1)*x(7)    
    b(2) = A(3,2)*x(1)+A(4,2)*x(2)+A(5,2)*x(3)+A(6,2)*x(4)+A(7,2)*x(5)+A(2,2)*x(6)+A(1,2)*x(7)    
    b(3) = A(2,3)*x(1)+A(3,3)*x(2)+A(4,3)*x(3)+A(5,3)*x(4)+A(6,3)*x(5)+A(7,3)*x(6)+A(1,3)*x(7)   

    ! center part
    do j=4,n-3
       b(j) = A(1,j)*x(j-3)+A(2,j)*x(j-2)+A(3,j)*x(j-1)+A(4,j)*x(j)+A(5,j)*x(j+1)+A(6,j)*x(j+2)+A(7,j)*x(j+3)
    enddo

    ! Lower part
    b(n-2) = A(7,n-2)*x(n-6)+A(1,n-2)*x(n-5)+A(2,n-2)*x(n-4)+A(3,n-2)*x(n-3)+A(4,n-2)*x(n-2)+A(5,n-2)*x(n-1)+A(6,n-2)*x(n)
    b(n-1) = A(7,n-1)*x(n-6)+A(6,n-1)*x(n-5)+A(1,n-1)*x(n-4)+A(2,n-1)*x(n-3)+A(3,n-1)*x(n-2)+A(4,n-1)*x(n-1)+A(5,n-1)*x(n)
    b(n  ) = A(7,n  )*x(n-6)+A(6,n  )*x(n-5)+A(5,n  )*x(n-4)+A(1,n  )*x(n-3)+A(2,n  )*x(n-2)+A(3,n  )*x(n-1)+A(4,n  )*x(n)

    return
  end subroutine LABS_mat_mul_3


  subroutine LABS_mat_mul_4(A,x,b,n)
    implicit none
    real(c_double) :: A(9,n)
    complex(c_double_complex) :: x(n), b(n)
    integer :: i,j,k,n
    integer :: m, bw

    m = 5
    bw = 4


    ! Upper part
    b(1) = A(5,1)*x(1)+A(6,1)*x(2)+A(7,1)*x(3)+A(8,1)*x(4)+A(9,1)*x(5)+A(4,1)*x(6)+A(3,1)*x(7)+A(2,1)*x(8)+A(1,1)*x(9)   
    b(2) = A(4,2)*x(1)+A(5,2)*x(2)+A(6,2)*x(3)+A(7,2)*x(4)+A(8,2)*x(5)+A(9,2)*x(6)+A(3,2)*x(7)+A(2,2)*x(8)+A(1,2)*x(9)   
    b(3) = A(3,3)*x(1)+A(4,3)*x(2)+A(5,3)*x(3)+A(6,3)*x(4)+A(7,3)*x(5)+A(8,3)*x(6)+A(9,3)*x(7)+A(2,3)*x(8)+A(1,3)*x(9)  
    b(4) = A(2,4)*x(1)+A(3,4)*x(2)+A(4,4)*x(3)+A(5,4)*x(4)+A(6,4)*x(5)+A(7,4)*x(6)+A(8,4)*x(7)+A(9,4)*x(8)+A(1,4)*x(9)   

    ! center part
    do j=5,n-4
       b(j) = A(1,j)*x(j-4)+A(2,j)*x(j-3)+A(3,j)*x(j-2)+A(4,j)*x(j-1)+A(5,j)*x(j)+A(6,j)*x(j+1)+A(7,j)*x(j+2)+A(8,j)*x(j+3) &
            +A(9,j)*x(j+4)
    enddo


    ! Lower part
    b(n-3) = A(9,n-3)*x(n-8)+A(1,n-3)*x(n-7)+A(2,n-3)*x(n-6)+A(3,n-3)*x(n-5)+A(4,n-3)*x(n-4)+A(5,n-3)*x(n-3)+A(6,n-3)*x(n-2)&
         +A(7,n-3)*x(n-1)+A(8,n-3)*x(n)
    b(n-2) = A(9,n-2)*x(n-8)+A(8,n-2)*x(n-7)+A(1,n-2)*x(n-6)+A(2,n-2)*x(n-5)+A(3,n-2)*x(n-4)+A(4,n-2)*x(n-3)+A(5,n-2)*x(n-2)&
         +A(6,n-2)*x(n-1)+A(7,n-2)*x(n)
    b(n-1) = A(9,n-1)*x(n-8)+A(8,n-1)*x(n-7)+A(7,n-1)*x(n-6)+A(1,n-1)*x(n-5)+A(2,n-1)*x(n-4)+A(3,n-1)*x(n-3)+A(4,n-1)*x(n-2)&
         +A(5,n-1)*x(n-1)+A(6,n-1)*x(n)
    b(n  ) = A(9,n  )*x(n-8)+A(8,n  )*x(n-7)+A(7,n  )*x(n-6)+A(6,n  )*x(n-5)+A(1,n  )*x(n-4)+A(2,n  )*x(n-3)+A(3,n  )*x(n-2)&
         +A(4,n  )*x(n-1)+A(5,n  )*x(n)

    return
  end subroutine LABS_mat_mul_4



  subroutine LABS_mat_mul_5(A,x,b,n)
    implicit none
    real(c_double) :: A(11,n)
    complex(c_double_complex) :: x(n), b(n)
    integer :: i,j,k,n

    ! Upper part
    b(1) = A(6,1)*x(1)+A(7,1)*x(2)+A(8,1)*x(3)+A(9,1)*x(4)+A(10,1)*x(5)+A(11,1)*x(6)+A(5,1)*x(7)+A(4,1)*x(8)+A(3,1)*x(9)&
         +A(2,1)*x(10)+A(1,1)*x(11)  
    b(2) = A(5,2)*x(1)+A(6,2)*x(2)+A(7,2)*x(3)+A(8,2)*x(4)+A(9,2)*x(5)+A(10,2)*x(6)+A(11,2)*x(7)+A(4,2)*x(8)+A(3,2)*x(9)&
         +A(2,2)*x(10)+A(1,2)*x(11)  
    b(3) = A(4,3)*x(1)+A(5,3)*x(2)+A(6,3)*x(3)+A(7,3)*x(4)+A(8,3)*x(5)+A(9,3)*x(6)+A(10,3)*x(7)+A(11,3)*x(8)+A(3,3)*x(9)&
         +A(2,3)*x(10)+A(1,3)*x(11)  
    b(4) = A(3,4)*x(1)+A(4,4)*x(2)+A(5,4)*x(3)+A(6,4)*x(4)+A(7,4)*x(5)+A(8,4)*x(6)+A(9,4)*x(7)+A(10,4)*x(8)+A(11,4)*x(9)&
         +A(2,4)*x(10)+A(1,4)*x(11)  
    b(5) = A(2,5)*x(1)+A(3,5)*x(2)+A(4,5)*x(3)+A(5,5)*x(4)+A(6,5)*x(5)+A(7,5)*x(6)+A(8,5)*x(7)+A(9,5)*x(8)+A(10,5)*x(9)&
         +A(11,5)*x(10)+A(1,5)*x(11)  

    ! center part
    do j=6,n-5
       b(j) = A(1,j)*x(j-5)+A(2,j)*x(j-4)+A(3,j)*x(j-3)+A(4,j)*x(j-2)+A(5,j)*x(j-1)+A(6,j)*x(j)+A(7,j)*x(j+1)+A(8,j)*x(j+2)&
            +A(9,j)*x(j+3)+A(10,j)*x(j+4)+A(11,j)*x(j+5)
    enddo


    ! Lower part
    b(n-4) = A(11,n-4)*x(n-10)+A(1,n-4)*x(n-9)+A(2,n-4)*x(n-8)+A(3,n-4)*x(n-7)+A(4,n-4)*x(n-6)+A(5,n-4)*x(n-5)+A(6,n-4)*x(n-4)&
         +A(7,n-4)*x(n-3)+A(8,n-4)*x(n-2)+A(9,n-4)*x(n-1)+A(10,n-4)*x(n) 
    b(n-3) = A(11,n-3)*x(n-10)+A(10,n-3)*x(n-9)+A(1,n-3)*x(n-8)+A(2,n-3)*x(n-7)+A(3,n-3)*x(n-6)+A(4,n-3)*x(n-5)+A(5,n-3)*x(n-4)&
         +A(6,n-3)*x(n-3)+A(7,n-3)*x(n-2)+A(8,n-3)*x(n-1)+A(9,n-3)*x(n)
    b(n-2) = A(11,n-2)*x(n-10)+A(10,n-2)*x(n-9)+A(9,n-2)*x(n-8)+A(1,n-2)*x(n-7)+A(2,n-2)*x(n-6)+A(3,n-2)*x(n-5)+A(4,n-2)*x(n-4)&
         +A(5,n-2)*x(n-3)+A(6,n-2)*x(n-2)+A(7,n-2)*x(n-1)+A(8,n-2)*x(n)
    b(n-1) = A(11,n-1)*x(n-10)+A(10,n-1)*x(n-9)+A(9,n-1)*x(n-8)+A(8,n-1)*x(n-7)+A(1,n-1)*x(n-6)+A(2,n-1)*x(n-5)+A(3,n-1)*x(n-4)&
         +A(4,n-1)*x(n-3)+A(5,n-1)*x(n-2)+A(6,n-1)*x(n-1)+A(7,n-1)*x(n)
    b(n  ) = A(11,n  )*x(n-10)+A(10,n  )*x(n-9)+A(9,n  )*x(n-8)+A(8,n  )*x(n-7)+A(7,n  )*x(n-6)+A(1,n  )*x(n-5)+A(2,n  )*x(n-4)&
         +A(3,n  )*x(n-3)+A(4,n  )*x(n-2)+A(5,n  )*x(n-1)+A(6,n  )*x(n)

    return
  end subroutine LABS_mat_mul_5


  subroutine LABS_mat_mul_6(A,x,b,n)
    implicit none
    real(c_double) :: A(13,n)
    complex(c_double_complex) :: x(n), b(n)
    integer :: i,j,k,n


    ! Upper part
    b(1) = A(7,1)*x(1)+A(8,1)*x(2)+A(9,1)*x(3)+A(10,1)*x(4)+A(11,1)*x(5)+A(12,1)*x(6)+A(13,1)*x(7)+A(6,1)*x(8)+A(5,1)*x(9)&
         +A(4,1)*x(10)+A(3,1)*x(11)+A(2,1)*x(12)+A(1,1)*x(13)  
    b(2) = A(6,2)*x(1)+A(7,2)*x(2)+A(8,2)*x(3)+A(9,2)*x(4)+A(10,2)*x(5)+A(11,2)*x(6)+A(12,2)*x(7)+A(13,2)*x(8)+A(5,2)*x(9)&
         +A(4,2)*x(10)+A(3,2)*x(11)+A(2,2)*x(12)+A(1,2)*x(13)  
    b(3) = A(5,3)*x(1)+A(6,3)*x(2)+A(7,3)*x(3)+A(8,3)*x(4)+A(9,3)*x(5)+A(10,3)*x(6)+A(11,3)*x(7)+A(12,3)*x(8)+A(13,3)*x(9)&
         +A(4,3)*x(10)+A(3,3)*x(11)+A(2,3)*x(12)+A(1,3)*x(13)  
    b(4) = A(4,4)*x(1)+A(5,4)*x(2)+A(6,4)*x(3)+A(7,4)*x(4)+A(8,4)*x(5)+A(9,4)*x(6)+A(10,4)*x(7)+A(11,4)*x(8)+A(12,4)*x(9)&
         +A(13,4)*x(10)+A(3,4)*x(11)+A(2,4)*x(12)+A(1,4)*x(13)  
    b(5) = A(3,5)*x(1)+A(4,5)*x(2)+A(5,5)*x(3)+A(6,5)*x(4)+A(7,5)*x(5)+A(8,5)*x(6)+A(9,5)*x(7)+A(10,5)*x(8)+A(11,5)*x(9)&
         +A(12,5)*x(10)+A(13,5)*x(11)+A(2,5)*x(12)+A(1,5)*x(13)  
    b(6) = A(2,6)*x(1)+A(3,6)*x(2)+A(4,6)*x(3)+A(5,6)*x(4)+A(6,6)*x(5)+A(7,6)*x(6)+A(8,6)*x(7)+A(9,6)*x(8)+A(10,6)*x(9)&
         +A(11,6)*x(10)+A(12,6)*x(11)+A(13,6)*x(12)+A(1,6)*x(13)  

    ! center part
    do j=7,n-6
       b(j) = A(1,j)*x(j-6)+A(2,j)*x(j-5)+A(3,j)*x(j-4)+A(4,j)*x(j-3)+A(5,j)*x(j-2)+A(6,j)*x(j-1)+A(7,j)*x(j)+A(8,j)*x(j+1)&
            +A(9,j)*x(j+2)+A(10,j)*x(j+3)+A(11,j)*x(j+4)+A(12,j)*x(j+5)+A(13,j)*x(j+6)  
    enddo


    ! Lower part
    b(n-5) = A(13,n-5)*x(n-12)+A(1,n-5)*x(n-11)+A(2,n-5)*x(n-10)+A(3,n-5)*x(n-9)+A(4,n-5)*x(n-8)+A(5,n-5)*x(n-7)+A(6,n-5)*x(n-6)&
         +A(7,n-5)*x(n-5)+A(8,n-5)*x(n-4)+A(9,n-5)*x(n-3)+A(10,n-5)*x(n-2)+A(11,n-5)*x(n-1)+A(12,n-5)*x(n) 
    b(n-4) = A(13,n-4)*x(n-12)+A(12,n-4)*x(n-11)+A(1,n-4)*x(n-10)+A(2,n-4)*x(n-9)+A(3,n-4)*x(n-8)+A(4,n-4)*x(n-7)+A(5,n-4)*x(n-6)&
         +A(6,n-4)*x(n-5)+A(7,n-4)*x(n-4)+A(8,n-4)*x(n-3)+A(9,n-4)*x(n-2)+A(10,n-4)*x(n-1)+A(11,n-4)*x(n) 
    b(n-3) = A(13,n-3)*x(n-12)+A(12,n-3)*x(n-11)+A(11,n-3)*x(n-10)+A(1,n-3)*x(n-9)+A(2,n-3)*x(n-8)+A(3,n-3)*x(n-7)+A(4,n-3)*x(n-6)&
         +A(5,n-3)*x(n-5)+A(6,n-3)*x(n-4)+A(7,n-3)*x(n-3)+A(8,n-3)*x(n-2)+A(9,n-3)*x(n-1)+A(10,n-3)*x(n) 
    b(n-2) = A(13,n-2)*x(n-12)+A(12,n-2)*x(n-11)+A(11,n-2)*x(n-10)+A(10,n-2)*x(n-9)+A(1,n-2)*x(n-8)+A(2,n-2)*x(n-7)+A(3,n-2)*x(n-6)&
         +A(4,n-2)*x(n-5)+A(5,n-2)*x(n-4)+A(6,n-2)*x(n-3)+A(7,n-2)*x(n-2)+A(8,n-2)*x(n-1)+A(9,n-2)*x(n) 
    b(n-1) = A(13,n-1)*x(n-12)+A(12,n-1)*x(n-11)+A(11,n-1)*x(n-10)+A(10,n-1)*x(n-9)+A(9,n-1)*x(n-8)+A(1,n-1)*x(n-7)+A(2,n-1)*x(n-6)&
         +A(3,n-1)*x(n-5)+A(4,n-1)*x(n-4)+A(5,n-1)*x(n-3)+A(6,n-1)*x(n-2)+A(7,n-1)*x(n-1)+A(8,n-1)*x(n) 
    b(n  ) = A(13,n  )*x(n-12)+A(12,n  )*x(n-11)+A(11,n  )*x(n-10)+A(10,n  )*x(n-9)+A(9,n  )*x(n-8)+A(8,n  )*x(n-7)+A(1,n  )*x(n-6)&
         +A(2,n  )*x(n-5)+A(3,n  )*x(n-4)+A(4,n  )*x(n-3)+A(5,n  )*x(n-2)+A(6,n  )*x(n-1)+A(7,n  )*x(n) 

    return
  end subroutine LABS_mat_mul_6

  subroutine LABS_mat_mul_7(A,x,b,n)
    implicit none
    real(c_double) :: A(15,n)
    complex(c_double_complex) :: x(n), b(n)
    integer :: i,j,k,n

    ! Upper part
    b(1) = A(8,1)*x(1)+A(9,1)*x(2)+A(10,1)*x(3)+A(11,1)*x(4)+A(12,1)*x(5)+A(13,1)*x(6)+A(14,1)*x(7)+A(15,1)*x(8)+A(7,1)*x(9)&
         +A(6,1)*x(10)+A(5,1)*x(11)+A(4,1)*x(12)+A(3,1)*x(13)+A(2,1)*x(14)+A(1,1)*x(15) 
    b(2) = A(7,2)*x(1)+A(8,2)*x(2)+A(9,2)*x(3)+A(10,2)*x(4)+A(11,2)*x(5)+A(12,2)*x(6)+A(13,2)*x(7)+A(14,2)*x(8)+A(15,2)*x(9)&
         +A(6,2)*x(10)+A(5,2)*x(11)+A(4,2)*x(12)+A(3,2)*x(13)+A(2,2)*x(14)+A(1,2)*x(15) 
    b(3) = A(6,3)*x(1)+A(7,3)*x(2)+A(8,3)*x(3)+A(9,3)*x(4)+A(10,3)*x(5)+A(11,3)*x(6)+A(12,3)*x(7)+A(13,3)*x(8)+A(14,3)*x(9)&
         +A(15,3)*x(10)+A(5,3)*x(11)+A(4,3)*x(12)+A(3,3)*x(13)+A(2,3)*x(14)+A(1,3)*x(15) 
    b(4) = A(5,4)*x(1)+A(6,4)*x(2)+A(7,4)*x(3)+A(8,4)*x(4)+A(9,4)*x(5)+A(10,4)*x(6)+A(11,4)*x(7)+A(12,4)*x(8)+A(13,4)*x(9)&
         +A(14,4)*x(10)+A(15,4)*x(11)+A(4,4)*x(12)+A(3,4)*x(13)+A(2,4)*x(14)+A(1,4)*x(15) 
    b(5) = A(4,5)*x(1)+A(5,5)*x(2)+A(6,5)*x(3)+A(7,5)*x(4)+A(8,5)*x(5)+A(9,5)*x(6)+A(10,5)*x(7)+A(11,5)*x(8)+A(12,5)*x(9)&
         +A(13,5)*x(10)+A(14,5)*x(11)+A(15,5)*x(12)+A(3,5)*x(13)+A(2,5)*x(14)+A(1,5)*x(15) 
    b(6) = A(3,6)*x(1)+A(4,6)*x(2)+A(5,6)*x(3)+A(6,6)*x(4)+A(7,6)*x(5)+A(8,6)*x(6)+A(9,6)*x(7)+A(10,6)*x(8)+A(11,6)*x(9)&
         +A(12,6)*x(10)+A(13,6)*x(11)+A(14,6)*x(12)+A(15,6)*x(13)+A(2,6)*x(14)+A(1,6)*x(15) 
    b(7) = A(2,7)*x(1)+A(3,7)*x(2)+A(4,7)*x(3)+A(5,7)*x(4)+A(6,7)*x(5)+A(7,7)*x(6)+A(8,7)*x(7)+A(9,7)*x(8)+A(10,7)*x(9)&
         +A(11,7)*x(10)+A(12,7)*x(11)+A(13,7)*x(12)+A(14,7)*x(13)+A(15,7)*x(14)+A(1,7)*x(15) 

    ! center part
    do j=8,n-7
       b(j) = A(1,j)*x(j-7)+A(2,j)*x(j-6)+A(3,j)*x(j-5)+A(4,j)*x(j-4)+A(5,j)*x(j-3)+A(6,j)*x(j-2)+A(7,j)*x(j-1)+A(8,j)*x(j)&
            +A(9,j)*x(j+1)+A(10,j)*x(j+2)+A(11,j)*x(j+3)+A(12,j)*x(j+4)+A(13,j)*x(j+5)+A(14,j)*x(j+6)+A(15,j)*x(j+7) 
    enddo


    ! Lower part
    b(n-6) = A(15,n-6)*x(n-14)+A(1,n-6)*x(n-13)+A(2,n-6)*x(n-12)+A(3,n-6)*x(n-11)+A(4,n-6)*x(n-10)+A(5,n-6)*x(n-9)+A(6,n-6)*x(n-8)&
         +A(7,n-6)*x(n-7)+A(8,n-6)*x(n-6)+A(9,n-6)*x(n-5)+A(10,n-6)*x(n-4)+A(11,n-6)*x(n-3)+A(12,n-6)*x(n-2)+A(13,n-6)*x(n-1)&
         +A(14,n-6)*x(n) 
    b(n-5) = A(15,n-5)*x(n-14)+A(14,n-5)*x(n-13)+A(1,n-5)*x(n-12)+A(2,n-5)*x(n-11)+A(3,n-5)*x(n-10)+A(4,n-5)*x(n-9)+A(5,n-5)*x(n-8)&
         +A(6,n-5)*x(n-7)+A(7,n-5)*x(n-6)+A(8,n-5)*x(n-5)+A(9,n-5)*x(n-4)+A(10,n-5)*x(n-3)+A(11,n-5)*x(n-2)+A(12,n-5)*x(n-1)&
         +A(13,n-5)*x(n) 
    b(n-4) = A(15,n-4)*x(n-14)+A(14,n-4)*x(n-13)+A(13,n-4)*x(n-12)+A(1,n-4)*x(n-11)+A(2,n-4)*x(n-10)+A(3,n-4)*x(n-9) &
         +A(4,n-4)*x(n-8)+ A(5,n-4)*x(n-7)+A(6,n-4)*x(n-6)+A(7,n-4)*x(n-5)+A(8,n-4)*x(n-4)+A(9,n-4)*x(n-3)+A(10,n-4)*x(n-2)&
         +A(11,n-4)*x(n-1)+ A(12,n-4)*x(n) 
    b(n-3) = A(15,n-3)*x(n-14)+A(14,n-3)*x(n-13)+A(13,n-3)*x(n-12)+A(12,n-3)*x(n-11)+A(1,n-3)*x(n-10)+A(2,n-3)*x(n-9)&
         +A(3,n-3)*x(n-8)+A(4,n-3)*x(n-7)+A(5,n-3)*x(n-6)+A(6,n-3)*x(n-5)+A(7,n-3)*x(n-4)+A(8,n-3)*x(n-3)+A(9,n-3)*x(n-2)&
         +A(10,n-3)*x(n-1) +A(11,n-3)*x(n) 
    b(n-2) = A(15,n-2)*x(n-14)+A(14,n-2)*x(n-13)+A(13,n-2)*x(n-12)+A(12,n-2)*x(n-11)+A(11,n-2)*x(n-10)+A(1,n-2)*x(n-9)&
         +A(2,n-2)*x(n-8)+A(3,n-2)*x(n-7)+A(4,n-2)*x(n-6)+A(5,n-2)*x(n-5)+A(6,n-2)*x(n-4)+A(7,n-2)*x(n-3)+A(8,n-2)*x(n-2)&
         +A(9,n-2)*x(n-1)+A(10,n-2)*x(n) 
    b(n-1) = A(15,n-1)*x(n-14)+A(14,n-1)*x(n-13)+A(13,n-1)*x(n-12)+A(12,n-1)*x(n-11)+A(11,n-1)*x(n-10)+A(10,n-1)*x(n-9)&
         +A(1,n-1)*x(n-8)+A(2,n-1)*x(n-7)+A(3,n-1)*x(n-6)+A(4,n-1)*x(n-5)+A(5,n-1)*x(n-4)+A(6,n-1)*x(n-3)+A(7,n-1)*x(n-2)&
         +A(8,n-1)*x(n-1)+A(9,n-1)*x(n) 
    b(n  ) = A(15,n  )*x(n-14)+A(14,n  )*x(n-13)+A(13,n  )*x(n-12)+A(12,n  )*x(n-11)+A(11,n)*x(n-10)+A(10,n)*x(n-9)+A(9,n)*x(n-8)&
         +A(1,n  )*x(n-7)+A(2,n  )*x(n-6)+A(3,n  )*x(n-5)+A(4,n  )*x(n-4)+A(5,n  )*x(n-3)+A(6,n  )*x(n-2)+A(7,n  )*x(n-1)&
         +A(8,n  )*x(n) 

    return
  end subroutine LABS_mat_mul_7


  subroutine LABS_mat_mul_8(A,x,b,n)
    implicit none
    integer, parameter :: bw = 8
    integer :: n
    real(c_double) :: A(2*bw+1,n)
    complex(c_double_complex) :: x(n), b(n)
    integer :: i,j,k

    do j=1,bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=1-j,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,bw+1-j
          b(j) = b(j) + A(bw+2-j-i,j)*x(j+bw+i)
       enddo
    enddo

    do j=bw+1,n-bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
    enddo

    do j=n-bw+1,n
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,n-j
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,j-n+bw
          b(j) = b(j) + A(n+bw+1-j+i,j)*x(j-bw-i)
       enddo
    enddo
    return
  end subroutine LABS_mat_mul_8

  subroutine LABS_mat_mul_9(A,x,b,n)
    implicit none
    integer, parameter :: bw = 9
    integer :: n
    real(c_double) :: A(2*bw+1,n)
    complex(c_double_complex) :: x(n), b(n)
    integer :: i,j,k

    do j=1,bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=1-j,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,bw+1-j
          b(j) = b(j) + A(bw+2-j-i,j)*x(j+bw+i)
       enddo
    enddo

    do j=bw+1,n-bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
    enddo

    do j=n-bw+1,n
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,n-j
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,j-n+bw
          b(j) = b(j) + A(n+bw+1-j+i,j)*x(j-bw-i)
       enddo
    enddo
    return
  end subroutine LABS_mat_mul_9

  subroutine LABS_mat_mul_10(A,x,b,n)
    implicit none
    integer, parameter :: bw = 10
    integer :: n
    real(c_double) :: A(2*bw+1,n)
    complex(c_double_complex) :: x(n), b(n)
    integer :: i,j,k

    do j=1,bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=1-j,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,bw+1-j
          b(j) = b(j) + A(bw+2-j-i,j)*x(j+bw+i)
       enddo
    enddo

    do j=bw+1,n-bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
    enddo

    do j=n-bw+1,n
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,n-j
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,j-n+bw
          b(j) = b(j) + A(n+bw+1-j+i,j)*x(j-bw-i)
       enddo
    enddo
    return
  end subroutine LABS_mat_mul_10

  subroutine LABS_mat_mul_bw(A,x,b,n,bw)
    implicit none
    integer :: bw
    integer :: n
    real(c_double) :: A(2*bw+1,n)
    complex(c_double_complex) :: x(n), b(n)
    integer :: i,j,k

    do j=1,bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=1-j,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,bw+1-j
          b(j) = b(j) + A(bw+2-j-i,j)*x(j+bw+i)
       enddo
    enddo

    do j=bw+1,n-bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
    enddo

    do j=n-bw+1,n
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,n-j
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,j-n+bw
          b(j) = b(j) + A(n+bw+1-j+i,j)*x(j-bw-i)
       enddo
    enddo
    return
  end subroutine LABS_mat_mul_bw









  subroutine fixed_banbks1_r(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(3,n)
    real(c_double) :: b(n)
    integer :: m, bw
    integer :: i
    m = 2
    bw = 1

    ! dummy = L^-1 b  -- Upper part
    b(1) = b(1) / A(m,1) 

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       b(i) = (b(i) - b(i-bw)*A(m-1,i)) / A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    b(n  ) = (b(n) - b(n-2)*A(m+1,n  ) - b(n-1)*A(m-1,n  )) / A(m,n  )

!!$  ! now dummy is in b
!!$  ! f = U^-1 dummy

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       b(i) = b(i) - b(i+1)*A(m+1,i) 
    enddo

    !b(n) = b(n) -- upper part

    b(1) = b(1) - b(2)*A(m+1,1) - b(3)*A(m-1,1)  

    return
  end subroutine fixed_banbks1_r

  subroutine fixed_banbks2_r(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(5,n)
    real(c_double) :: b(n)
    integer :: m, bw
    integer :: i
    m = 3
    bw = 2

    ! dummy = L^-1 b  -- Upper part
    b(1) = b(1) / A(m,1) 
    b(2) = (b(2) - b(1)*A(m-1,2)) / A(m,2)

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       b(i) = (b(i) - b(i-bw)*A(m-2,i) - b(i-bw+1)*A(m-1,i)) / A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    b(n-1) = (b(n-1) - b(n-4)*A(m+2,n-1) - b(n-3)*A(m-2,n-1) - b(n-2)*A(m-1,n-1)) / A(m,n-1)
    b(n  ) = (b(n  ) - b(n-4)*A(m+2,n  ) - b(n-3)*A(m+1,n  ) - b(n-2)*A(m-2,n  ) - b(n-1)*A(m-1,n  )) / A(m,n  )

!!$  ! now dummy is in b
!!$  ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    b(n-1) = b(n-1) - b(n)  *A(m+1,n-1)

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       b(i) = b(i) - b(i+1)*A(m+1,i) - b(i+2)*A(m+2,i) 
    enddo

    !b(n) = b(n) -- upper part

    b(2) = b(2) - b(3)*A(m+1,2) - b(4)*A(m+2,2) -  b(5)*A(m-2,2) 
    b(1) = b(1) - b(2)*A(m+1,1) - b(3)*A(m+2,1) -  b(4)*A(m-1,1) - b(5)*A(m-2,1) 
    return
  end subroutine fixed_banbks2_r

  subroutine fixed_banbks3_r(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(7,n)
    real(c_double) :: b(n)
    integer :: m, bw
    integer :: i
    m = 4
    bw = 3

    ! dummy = L^-1 b  -- Upper part
    b(1) = b(1) / A(m,1) 
    b(2) = (b(2) - b(1)*A(m-1,2)) / A(m,2)
    b(3) = (b(3) - b(1)*A(m-2,3) - b(2)*A(m-1,3)) / A(m,3)

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       b(i) = (b(i) - b(i-bw)*A(m-3,i) - b(i-bw+1)*A(m-2,i) - b(i-bw+2)*A(m-1,i)) / A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    b(n-2) = (b(n-2) - b(n-6)*A(m+3,n-2) - b(n-5)*A(m-3,n-2) - b(n-4)*A(m-2,n-2) - b(n-3)*A(m-1,n-2)) / A(m,n-2)
    b(n-1) = (b(n-1) - b(n-6)*A(m+3,n-1) - b(n-5)*A(m+2,n-1) - b(n-4)*A(m-3,n-1) - b(n-3)*A(m-2,n-1) &
         - b(n-2)*A(m-1,n-1)) / A(m,n-1)
    b(n  ) = (b(n  ) - b(n-6)*A(m+3,n  ) - b(n-5)*A(m+2,n  ) - b(n-4)*A(m+1,n  ) - b(n-3)*A(m-3,n  ) &
         - b(n-2)*A(m-2,n  ) - b(n-1)*A(m-1,n  )) / A(m,n  )

!!$  ! now dummy is in b
!!$  ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    b(n-1) = b(n-1) - b(n)  *A(m+1,n-1)
    b(n-2) = b(n-2) - b(n-1)*A(m+1,n-2) - b(n)  *A(m+2,n-2)

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       b(i) = b(i) - b(i+1)*A(m+1,i) - b(i+2)*A(m+2,i) - b(i+3)*A(m+3,i)
    enddo

    !b(n) = b(n) -- upper part

    b(3) = b(3) - b(4)*A(m+1,3) - b(5)*A(m+2,3) -  b(6)*A(m+3,3) -  b(7)*A(m-3,3)
    b(2) = b(2) - b(3)*A(m+1,2) - b(4)*A(m+2,2) -  b(5)*A(m+3,2) -  b(6)*A(m-2,2) -  b(7)*A(m-3,2) 
    b(1) = b(1) - b(2)*A(m+1,1) - b(3)*A(m+2,1) -  b(4)*A(m+3,1) -  b(5)*A(m-1,1) -  b(6)*A(m-2,1) -  b(7)*A(m-3,1) 

    return
  end subroutine fixed_banbks3_r

  subroutine fixed_banbks4_r(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(9,n)
    real(c_double) :: b(n)
!!$  integer :: m, bw
    integer :: i
!!$  m = 5
!!$  bw = 4

    ! dummy = L^-1 b  -- Upper part
    b(1) = b(1) / A(5,1) 
    b(2) = (b(2) - b(1)*A(4,2)) / A(5,2)
    b(3) = (b(3) - b(1)*A(3,3) - b(2)*A(4,3)) / A(5,3)
    b(4) = (b(4) - b(1)*A(2,4) - b(2)*A(3,4) - b(3)*A(4,4)) / A(5,4)

    ! dummy =  L^-1 b  -- Center part
    do i=5,n-4    
       b(i) = (b(i) - b(i-4)*A(1,i) - b(i-3)*A(2,i) - b(i-2)*A(3,i) - b(i-1)*A(4,i)) / A(5,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    b(n-3) = (b(n-3) - b(n-8)*A(9,n-3) - b(n-7)*A(1,n-3) - b(n-6)*A(2,n-3) - b(n-5)*A(3,n-3) - b(n-4)*A(4,n-3))&
         / A(5,n-3)
    b(n-2) = (b(n-2) - b(n-8)*A(9,n-2) - b(n-7)*A(8,n-2) - b(n-6)*A(1,n-2) - b(n-5)*A(2,n-2) - b(n-4)*A(3,n-2)&
         - b(n-3)*A(4,n-2)) / A(5,n-2)
    b(n-1) = (b(n-1) - b(n-8)*A(9,n-1) - b(n-7)*A(8,n-1) - b(n-6)*A(7,n-1) - b(n-5)*A(1,n-1) - b(n-4)*A(2,n-1)&
         - b(n-3)*A(3,n-1) - b(n-2)*A(4,n-1)) / A(5,n-1)
    b(n  ) = (b(n  ) - b(n-8)*A(9,n  ) - b(n-7)*A(8,n  ) - b(n-6)*A(7,n  ) - b(n-5)*A(6,n  ) - b(n-4)*A(1,n  ) &
         - b(n-3)*A(2,n  ) - b(n-2)*A(3,n  ) - b(n-1)*A(4,n  )) / A(5,n  )

!!$  ! now dummy is in b
!!$  ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    b(n-1) = b(n-1) - b(n)  *A(6,n-1)
    b(n-2) = b(n-2) - b(n-1)*A(6,n-2) - b(n)  *A(7,n-2)
    b(n-3) = b(n-3) - b(n-2)*A(6,n-3) - b(n-1)*A(7,n-3) - b(n)  *A(8,n-3)

    !b(n) = b(n) -- center part
    do i=n-4,5,-1
       b(i) = b(i) - b(i+1)*A(6,i) - b(i+2)*A(7,i) - b(i+3)*A(8,i) - b(i+4)*A(9,i) 
    enddo

    !b(n) = b(n) -- upper part

    b(4) = b(4) - b(5)*A(6,4) - b(6)*A(7,4) -  b(7)*A(8,4) -  b(8)*A(9,4) -  b(9)*A(1,4)
    b(3) = b(3) - b(4)*A(6,3) - b(5)*A(7,3) -  b(6)*A(8,3) -  b(7)*A(9,3) -  b(8)*A(2,3) -  b(9)*A(1,3)
    b(2) = b(2) - b(3)*A(6,2) - b(4)*A(7,2) -  b(5)*A(8,2) -  b(6)*A(9,2) -  b(7)*A(3,2) -  b(8)*A(2,2)&
         -  b(9)*A(1,2)
    b(1) = b(1) - b(2)*A(6,1) - b(3)*A(7,1) -  b(4)*A(8,1) -  b(5)*A(9,1) -  b(6)*A(4,1) -  b(7)*A(3,1)&
         -  b(8)*A(2,1) - b(9)*A(1,1)

    return
  end subroutine fixed_banbks4_r

  subroutine fixed_banbks5_r(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(11,n)
    real(c_double) :: b(n)
    integer :: m, bw
    integer :: i
    m = 6
    bw = 5

    ! dummy = L^-1 b  -- Upper part
    b(1) = b(1) / A(m,1) 
    b(2) = (b(2) - b(1)*A(m-1,2)) / A(m,2)
    b(3) = (b(3) - b(1)*A(m-2,3) - b(2)*A(m-1,3)) / A(m,3)
    b(4) = (b(4) - b(1)*A(m-3,4) - b(2)*A(m-2,4) - b(3)*A(m-1,4)) / A(m,4)
    b(5) = (b(5) - b(1)*A(m-4,5) - b(2)*A(m-3,5) - b(3)*A(m-2,5) - b(4)*A(m-1,5)) / A(m,5)  

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       b(i) = (b(i) - b(i-bw)*A(m-5,i) - b(i-bw+1)*A(m-4,i) - b(i-bw+2)*A(m-3,i) - b(i-bw+3)*A(m-2,i) &
            - b(i-bw+4)*A(m-1,i)) / A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    b(n-4) = (b(n-4) - b(n-10)*A(m+5,n-4) - b(n-9)*A(m-5,n-4) - b(n-8)*A(m-4,n-4) - b(n-7)*A(m-3,n-4) &
         - b(n-6)*A(m-2,n-4) - b(n-5)*A(m-1,n-4)) / A(m,n-4)
    b(n-3) = (b(n-3) - b(n-10)*A(m+5,n-3) - b(n-9)*A(m+4,n-3) - b(n-8)*A(m-5,n-3) - b(n-7)*A(m-4,n-3) &
         - b(n-6)*A(m-3,n-3) - b(n-5)*A(m-2,n-3) - b(n-4)*A(m-1,n-3)) / A(m,n-3)
    b(n-2) = (b(n-2) - b(n-10)*A(m+5,n-2) - b(n-9)*A(m+4,n-2) - b(n-8)*A(m+3,n-2) - b(n-7)*A(m-5,n-2)&
         - b(n-6)*A(m-4,n-2) - b(n-5)*A(m-3,n-2) - b(n-4)*A(m-2,n-2) - b(n-3)*A(m-1,n-2)) / A(m,n-2)
    b(n-1) = (b(n-1) - b(n-10)*A(m+5,n-1) - b(n-9)*A(m+4,n-1) - b(n-8)*A(m+3,n-1) - b(n-7)*A(m+2,n-1) &
         - b(n-6)*A(m-5,n-1) - b(n-5)*A(m-4,n-1) - b(n-4)*A(m-3,n-1) - b(n-3)*A(m-2,n-1) - b(n-2)*A(m-1,n-1)) / A(m,n-1)
    b(n  ) = (b(n  ) - b(n-10)*A(m+5,n  ) - b(n-9)*A(m+4,n  ) - b(n-8)*A(m+3,n  ) - b(n-7)*A(m+2,n  ) &
         - b(n-6)*A(m+1,n  ) - b(n-5)*A(m-5,n  ) - b(n-4)*A(m-4,n  ) - b(n-3)*A(m-3,n  ) - b(n-2)*A(m-2,n  ) &
         - b(n-1)*A(m-1,n  )) / A(m,n  )

!!$  ! now dummy is in b
!!$  ! f = U^-1 dummy


    !b(n) = b(n) -- Lower part
    b(n-1) = b(n-1) - b(n)  *A(m+1,n-1)
    b(n-2) = b(n-2) - b(n-1)*A(m+1,n-2) - b(n)  *A(m+2,n-2)
    b(n-3) = b(n-3) - b(n-2)*A(m+1,n-3) - b(n-1)*A(m+2,n-3) - b(n)  *A(m+3,n-3)
    b(n-4) = b(n-4) - b(n-3)*A(m+1,n-4) - b(n-2)*A(m+2,n-4) - b(n-1)*A(m+3,n-4) - b(n)  *A(m+4,n-4)

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       b(i) = b(i) - b(i+1)*A(m+1,i) - b(i+2)*A(m+2,i) - b(i+3)*A(m+3,i) - b(i+4)*A(m+4,i) - b(i+5)*A(m+5,i) 
    enddo

    !b(n) = b(n) -- upper part

    b(5) = b(5) - b(6)*A(m+1,5) - b(7)*A(m+2,5) -  b(8)*A(m+3,5) -  b(9)*A(m+4,5) - b(10)*A(m+5,5)&
         - b(11)*A(m-5,5)
    b(4) = b(4) - b(5)*A(m+1,4) - b(6)*A(m+2,4) -  b(7)*A(m+3,4) -  b(8)*A(m+4,4) -  b(9)*A(m+5,4)&
         - b(10)*A(m-4,4) - b(11)*A(m-5,4) 
    b(3) = b(3) - b(4)*A(m+1,3) - b(5)*A(m+2,3) -  b(6)*A(m+3,3) -  b(7)*A(m+4,3) -  b(8)*A(m+5,3)&
         -  b(9)*A(m-3,3) - b(10)*A(m-4,3) -b(11)*A(m-5,3) 
    b(2) = b(2) - b(3)*A(m+1,2) - b(4)*A(m+2,2) -  b(5)*A(m+3,2) -  b(6)*A(m+4,2) -  b(7)*A(m+5,2)&
         -  b(8)*A(m-2,2) -  b(9)*A(m-3,2) -b(10)*A(m-4,2) -b(11)*A(m-5,2) 
    b(1) = b(1) - b(2)*A(m+1,1) - b(3)*A(m+2,1) -  b(4)*A(m+3,1) -  b(5)*A(m+4,1) -  b(6)*A(m+5,1) &
         -  b(7)*A(m-1,1) -  b(8)*A(m-2,1) - b(9)*A(m-3,1) -b(10)*A(m-4,1) -b(11)*A(m-5,1) 

    return
  end subroutine fixed_banbks5_r

  subroutine fixed_banbks6_r(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(13,n)
    real(c_double) :: b(n)
    integer :: m, bw
    integer :: i
    m = 7
    bw = 6

    ! dummy = L^-1 b  -- Upper part
    b(1) = b(1) / A(m,1) 
    b(2) = (b(2) - b(1)*A(m-1,2)) / A(m,2)
    b(3) = (b(3) - b(1)*A(m-2,3) - b(2)*A(m-1,3)) / A(m,3)
    b(4) = (b(4) - b(1)*A(m-3,4) - b(2)*A(m-2,4) - b(3)*A(m-1,4)) / A(m,4)
    b(5) = (b(5) - b(1)*A(m-4,5) - b(2)*A(m-3,5) - b(3)*A(m-2,5) - b(4)*A(m-1,5)) / A(m,5)  
    b(6) = (b(6) - b(1)*A(m-5,6) - b(2)*A(m-4,6) - b(3)*A(m-3,6) - b(4)*A(m-2,6) - b(5)*A(m-1,6)) / A(m,6)  

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       b(i) = (b(i) - b(i-bw)*A(m-6,i) - b(i-bw+1)*A(m-5,i) - b(i-bw+2)*A(m-4,i) - b(i-bw+3)*A(m-3,i) - b(i-bw+4)*A(m-2,i)&
            - b(i-bw+5)*A(m-1,i)) / A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    b(n-5) = (b(n-5) - b(n-12)*A(m+6,n-5) - b(n-11)*A(m-6,n-5) - b(n-10)*A(m-5,n-5) - b(n-9)*A(m-4,n-5)&
         - b(n-8)*A(m-3,n-5) - b(n-7)*A(m-2,n-5) - b(n-6)*A(m-1,n-5)) / A(m,n-5)
    b(n-4) = (b(n-4) - b(n-12)*A(m+6,n-4) - b(n-11)*A(m+5,n-4) - b(n-10)*A(m-6,n-4) - b(n-9)*A(m-5,n-4) &
         - b(n-8)*A(m-4,n-4) - b(n-7)*A(m-3,n-4) - b(n-6)*A(m-2,n-4) - b(n-5)*A(m-1,n-4)) / A(m,n-4)
    b(n-3) = (b(n-3) - b(n-12)*A(m+6,n-3) - b(n-11)*A(m+5,n-3) - b(n-10)*A(m+4,n-3) - b(n-9)*A(m-6,n-3)&
         - b(n-8)*A(m-5,n-3) - b(n-7)*A(m-4,n-3) - b(n-6)*A(m-3,n-3) - b(n-5)*A(m-2,n-3) - b(n-4)*A(m-1,n-3)) / A(m,n-3)
    b(n-2) = (b(n-2) - b(n-12)*A(m+6,n-2) - b(n-11)*A(m+5,n-2) - b(n-10)*A(m+4,n-2) - b(n-9)*A(m+3,n-2) &
         - b(n-8)*A(m-6,n-2) - b(n-7)*A(m-5,n-2) - b(n-6)*A(m-4,n-2) - b(n-5)*A(m-3,n-2) - b(n-4)*A(m-2,n-2) - b(n-3)*A(m-1,n-2)) &
         / A(m,n-2)
    b(n-1) = (b(n-1) - b(n-12)*A(m+6,n-1) - b(n-11)*A(m+5,n-1) - b(n-10)*A(m+4,n-1) - b(n-9)*A(m+3,n-1) &
         - b(n-8)*A(m+2,n-1) - b(n-7)*A(m-6,n-1) - b(n-6)*A(m-5,n-1) - b(n-5)*A(m-4,n-1) - b(n-4)*A(m-3,n-1) - b(n-3)*A(m-2,n-1)&
         - b(n-2)*A(m-1,n-1)) / A(m,n-1)
    b(n  ) = (b(n  ) - b(n-12)*A(m+6,n  ) - b(n-11)*A(m+5,n  ) - b(n-10)*A(m+4,n  ) - b(n-9)*A(m+3,n  ) &
         - b(n-8)*A(m+2,n  ) - b(n-7)*A(m+1,n  ) - b(n-6)*A(m-6,n  ) - b(n-5)*A(m-5,n  ) - b(n-4)*A(m-4,n  ) - b(n-3)*A(m-3,n  )&
         - b(n-2)*A(m-2,n  ) - b(n-1)*A(m-1,n  )) / A(m,n  )

!!$  ! now dummy is in b
!!$  ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    b(n-1) = b(n-1) - b(n)  *A(m+1,n-1)
    b(n-2) = b(n-2) - b(n-1)*A(m+1,n-2) - b(n)  *A(m+2,n-2)
    b(n-3) = b(n-3) - b(n-2)*A(m+1,n-3) - b(n-1)*A(m+2,n-3) - b(n)  *A(m+3,n-3)
    b(n-4) = b(n-4) - b(n-3)*A(m+1,n-4) - b(n-2)*A(m+2,n-4) - b(n-1)*A(m+3,n-4) - b(n)  *A(m+4,n-4)
    b(n-5) = b(n-5) - b(n-4)*A(m+1,n-5) - b(n-3)*A(m+2,n-5) - b(n-2)*A(m+3,n-5) - b(n-1)*A(m+4,n-5) - b(n)  *A(m+5,n-5)

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       b(i) = b(i) - b(i+1)*A(m+1,i) - b(i+2)*A(m+2,i) - b(i+3)*A(m+3,i) - b(i+4)*A(m+4,i) - b(i+5)*A(m+5,i) - b(i+6)*A(m+6,i)
    enddo

    !b(n) = b(n) -- upper part

    b(6) = b(6) - b(7)*A(m+1,6) - b(8)*A(m+2,6) -  b(9)*A(m+3,6) - b(10)*A(m+4,6) - b(11)*A(m+5,6) - b(12)*A(m+6,6)&
         - b(13)*A(m-6,6)
    b(5) = b(5) - b(6)*A(m+1,5) - b(7)*A(m+2,5) -  b(8)*A(m+3,5) -  b(9)*A(m+4,5) - b(10)*A(m+5,5) - b(11)*A(m+6,5)&
         - b(12)*A(m-5,5) -b(13)*A(m-6,5) 
    b(4) = b(4) - b(5)*A(m+1,4) - b(6)*A(m+2,4) -  b(7)*A(m+3,4) -  b(8)*A(m+4,4) -  b(9)*A(m+5,4) - b(10)*A(m+6,4)&
         - b(11)*A(m-4,4) -b(12)*A(m-5,4) -b(13)*A(m-6,4) 
    b(3) = b(3) - b(4)*A(m+1,3) - b(5)*A(m+2,3) -  b(6)*A(m+3,3) -  b(7)*A(m+4,3) -  b(8)*A(m+5,3) -  b(9)*A(m+6,3)&
         - b(10)*A(m-3,3) -b(11)*A(m-4,3) -b(12)*A(m-5,3) -b(13)*A(m-6,3) 
    b(2) = b(2) - b(3)*A(m+1,2) - b(4)*A(m+2,2) -  b(5)*A(m+3,2) -  b(6)*A(m+4,2) -  b(7)*A(m+5,2) -  b(8)*A(m+6,2)&
         -  b(9)*A(m-2,2) -b(10)*A(m-3,2) -b(11)*A(m-4,2) -b(12)*A(m-5,2) -b(13)*A(m-6,2)
    b(1) = b(1) - b(2)*A(m+1,1) - b(3)*A(m+2,1) -  b(4)*A(m+3,1) -  b(5)*A(m+4,1) -  b(6)*A(m+5,1) -  b(7)*A(m+6,1)&
         -  b(8)*A(m-1,1) - b(9)*A(m-2,1) -b(10)*A(m-3,1) -b(11)*A(m-4,1) -b(12)*A(m-5,1) -b(13)*A(m-6,1) 

    return
  end subroutine fixed_banbks6_r

  subroutine fixed_banbks7_r(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(15,n)
    real(c_double) :: b(n)
    integer :: m, bw
    integer :: i
    m = 8
    bw = 7

    ! dummy = L^-1 b  -- Upper part
    b(1) = b(1) / A(m,1) 
    b(2) = (b(2) - b(1)*A(m-1,2)) / A(m,2)
    b(3) = (b(3) - b(1)*A(m-2,3) - b(2)*A(m-1,3)) / A(m,3)
    b(4) = (b(4) - b(1)*A(m-3,4) - b(2)*A(m-2,4) - b(3)*A(m-1,4)) / A(m,4)
    b(5) = (b(5) - b(1)*A(m-4,5) - b(2)*A(m-3,5) - b(3)*A(m-2,5) - b(4)*A(m-1,5)) / A(m,5)  
    b(6) = (b(6) - b(1)*A(m-5,6) - b(2)*A(m-4,6) - b(3)*A(m-3,6) - b(4)*A(m-2,6) - b(5)*A(m-1,6)) / A(m,6)  
    b(7) = (b(7) - b(1)*A(m-6,7) - b(2)*A(m-5,7) - b(3)*A(m-4,7) - b(4)*A(m-3,7) - b(5)*A(m-2,7) - b(6)*A(m-1,7)) / A(m,7)  

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       b(i) = (b(i) - b(i-bw)*A(m-7,i) - b(i-bw+1)*A(m-6,i) - b(i-bw+2)*A(m-5,i) - b(i-bw+3)*A(m-4,i) - b(i-bw+4)*A(m-3,i) &
            - b(i-bw+5)*A(m-2,i) - b(i-bw+6)*A(m-1,i)) / A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    b(n-6) = (b(n-6) - b(n-14)*A(m+7,n-6) - b(n-13)*A(m-7,n-6) - b(n-12)*A(m-6,n-6) - b(n-11)*A(m-5,n-6) - b(n-10)*A(m-4,n-6)&
         - b(n-9)*A(m-3,n-6) - b(n-8)*A(m-2,n-6) - b(n-7)*A(m-1,n-6)) / A(m,n-6)
    b(n-5) = (b(n-5) - b(n-14)*A(m+7,n-5) - b(n-13)*A(m+6,n-5) - b(n-12)*A(m-7,n-5) - b(n-11)*A(m-6,n-5) - b(n-10)*A(m-5,n-5)&
         - b(n-9)*A(m-4,n-5) - b(n-8)*A(m-3,n-5) - b(n-7)*A(m-2,n-5) - b(n-6)*A(m-1,n-5)) / A(m,n-5)
    b(n-4) = (b(n-4) - b(n-14)*A(m+7,n-4) - b(n-13)*A(m+6,n-4) - b(n-12)*A(m+5,n-4) - b(n-11)*A(m-7,n-4) - b(n-10)*A(m-6,n-4)&
         - b(n-9)*A(m-5,n-4) - b(n-8)*A(m-4,n-4) - b(n-7)*A(m-3,n-4) - b(n-6)*A(m-2,n-4) - b(n-5)*A(m-1,n-4)) / A(m,n-4)
    b(n-3) = (b(n-3) - b(n-14)*A(m+7,n-3) - b(n-13)*A(m+6,n-3) - b(n-12)*A(m+5,n-3) - b(n-11)*A(m+4,n-3) - b(n-10)*A(m-7,n-3)&
         - b(n-9)*A(m-6,n-3) - b(n-8)*A(m-5,n-3) - b(n-7)*A(m-4,n-3) - b(n-6)*A(m-3,n-3) - b(n-5)*A(m-2,n-3) - b(n-4)*A(m-1,n-3))&
         / A(m,n-3)
    b(n-2) = (b(n-2) - b(n-14)*A(m+7,n-2) - b(n-13)*A(m+6,n-2) - b(n-12)*A(m+5,n-2) - b(n-11)*A(m+4,n-2) - b(n-10)*A(m+3,n-2)&
         - b(n-9)*A(m-7,n-2) - b(n-8)*A(m-6,n-2) - b(n-7)*A(m-5,n-2) - b(n-6)*A(m-4,n-2) - b(n-5)*A(m-3,n-2) - b(n-4)*A(m-2,n-2)&
         - b(n-3)*A(m-1,n-2)) / A(m,n-2)
    b(n-1) = (b(n-1) - b(n-14)*A(m+7,n-1) - b(n-13)*A(m+6,n-1) - b(n-12)*A(m+5,n-1) - b(n-11)*A(m+4,n-1) - b(n-10)*A(m+3,n-1)&
         - b(n-9)*A(m+2,n-1) - b(n-8)*A(m-7,n-1) - b(n-7)*A(m-6,n-1) - b(n-6)*A(m-5,n-1) - b(n-5)*A(m-4,n-1) - b(n-4)*A(m-3,n-1)&
         - b(n-3)*A(m-2,n-1) - b(n-2)*A(m-1,n-1)) / A(m,n-1)
    b(n  ) = (b(n  ) - b(n-14)*A(m+7,n  ) - b(n-13)*A(m+6,n  ) - b(n-12)*A(m+5,n  ) - b(n-11)*A(m+4,n  ) - b(n-10)*A(m+3,n  )&
         - b(n-9)*A(m+2,n  ) - b(n-8)*A(m+1,n  ) - b(n-7)*A(m-7,n  ) - b(n-6)*A(m-6,n  ) - b(n-5)*A(m-5,n  ) - b(n-4)*A(m-4,n  )&
         - b(n-3)*A(m-3,n  ) - b(n-2)*A(m-2,n  ) - b(n-1)*A(m-1,n  )) / A(m,n  )

!!$  ! now dummy is in b
!!$  ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    b(n-1) = b(n-1) - b(n)  *A(m+1,n-1)
    b(n-2) = b(n-2) - b(n-1)*A(m+1,n-2) - b(n)  *A(m+2,n-2)
    b(n-3) = b(n-3) - b(n-2)*A(m+1,n-3) - b(n-1)*A(m+2,n-3) - b(n)  *A(m+3,n-3)
    b(n-4) = b(n-4) - b(n-3)*A(m+1,n-4) - b(n-2)*A(m+2,n-4) - b(n-1)*A(m+3,n-4) - b(n)  *A(m+4,n-4)
    b(n-5) = b(n-5) - b(n-4)*A(m+1,n-5) - b(n-3)*A(m+2,n-5) - b(n-2)*A(m+3,n-5) - b(n-1)*A(m+4,n-5) - b(n)  *A(m+5,n-5)
    b(n-6) = b(n-6) - b(n-5)*A(m+1,n-6) - b(n-4)*A(m+2,n-6) - b(n-3)*A(m+3,n-6) - b(n-2)*A(m+4,n-6) - b(n-1)*A(m+5,n-6)&
         - b(n)*A(m+6,n-6)

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       b(i) = b(i) - b(i+1)*A(m+1,i) - b(i+2)*A(m+2,i) - b(i+3)*A(m+3,i) - b(i+4)*A(m+4,i) - b(i+5)*A(m+5,i) - b(i+6)*A(m+6,i) &
            - b(i+7)*A(m+7,i)
    enddo

    !b(n) = b(n) -- upper part

    b(7) = b(7) - b(8)*A(m+1,7) - b(9)*A(m+2,7) - b(10)*A(m+3,7) - b(11)*A(m+4,7) - b(12)*A(m+5,7) - b(13)*A(m+6,7) &
         - b(14)*A(m+7,7) -b(15)*A(m-7,7)
    b(6) = b(6) - b(7)*A(m+1,6) - b(8)*A(m+2,6) -  b(9)*A(m+3,6) - b(10)*A(m+4,6) - b(11)*A(m+5,6) - b(12)*A(m+6,6) &
         - b(13)*A(m+7,6) -b(14)*A(m-6,6) -b(15)*A(m-7,6)
    b(5) = b(5) - b(6)*A(m+1,5) - b(7)*A(m+2,5) -  b(8)*A(m+3,5) - b(9)*A(m+4,5) - b(10)*A(m+5,5) - b(11)*A(m+6,5) &
         - b(12)*A(m+7,5) -b(13)*A(m-5,5) -b(14)*A(m-6,5) -b(15)*A(m-7,5)
    b(4) = b(4) - b(5)*A(m+1,4) - b(6)*A(m+2,4) -  b(7)*A(m+3,4) -  b(8)*A(m+4,4) -  b(9)*A(m+5,4) - b(10)*A(m+6,4) &
         - b(11)*A(m+7,4) -b(12)*A(m-4,4) -b(13)*A(m-5,4) -b(14)*A(m-6,4) -b(15)*A(m-7,4)
    b(3) = b(3) - b(4)*A(m+1,3) - b(5)*A(m+2,3) -  b(6)*A(m+3,3) -  b(7)*A(m+4,3) -  b(8)*A(m+5,3) -  b(9)*A(m+6,3) &
         - b(10)*A(m+7,3) -b(11)*A(m-3,3) -b(12)*A(m-4,3) -b(13)*A(m-5,3) -b(14)*A(m-6,3) -b(15)*A(m-7,3)
    b(2) = b(2) - b(3)*A(m+1,2) - b(4)*A(m+2,2) -  b(5)*A(m+3,2) -  b(6)*A(m+4,2) -  b(7)*A(m+5,2) -  b(8)*A(m+6,2) &
         -  b(9)*A(m+7,2) -b(10)*A(m-2,2) -b(11)*A(m-3,2) -b(12)*A(m-4,2) -b(13)*A(m-5,2) -b(14)*A(m-6,2) -b(15)*A(m-7,2)
    b(1) = b(1) - b(2)*A(m+1,1) - b(3)*A(m+2,1) -  b(4)*A(m+3,1) -  b(5)*A(m+4,1) -  b(6)*A(m+5,1) -  b(7)*A(m+6,1)&
         -  b(8)*A(m+7,1) - b(9)*A(m-1,1) -b(10)*A(m-2,1) -b(11)*A(m-3,1) -b(12)*A(m-4,1) -b(13)*A(m-5,1) -b(14)*A(m-6,1) &
         -b(15)*A(m-7,1)


    return
  end subroutine fixed_banbks7_r

  subroutine fixed_banbks8_r(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(17,n)
    real(c_double) :: b(n)
    integer :: i,j,k
    integer, parameter :: m = 9
    integer, parameter :: bw = 8

    ! dummy = L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i-1
          b(i) = b(i) - b(j)*A(m-i+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       do j=1,bw
          b(i) = b(i) - b(i-bw-1+j)*A(m-bw-1+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m+bw+1-j,n-bw+i)
       enddo
       do j=i+1,i+bw
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m-bw+j-i-1,n-bw+i)
       enddo
       b(n-bw+i) = b(n-bw+i)/A(m,n-bw+i)
    enddo


    ! now dummy is in b
    ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    do i=1,bw-1
       do j=1,i
          b(n-i) = b(n-i) - b(n-i+j)*A(m+j,n-i)
       enddo
    enddo

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       do j=1,bw
          b(i) = b(i) - b(i+j)*A(m+j,i)
       enddo
    enddo

    !b(n) = b(n) -- upper part
    do i=1,bw
       do j=1,bw
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+j,bw+1-i)
       enddo
       do j=bw+1,bw+i
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+i-j,bw+1-i)
       enddo
    enddo

    return
  end subroutine fixed_banbks8_r

  subroutine fixed_banbks9_r(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(19,n)
    real(c_double) :: b(n)
    integer :: i,j,k
    integer, parameter :: m = 10
    integer, parameter :: bw = 9

    ! dummy = L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i-1
          b(i) = b(i) - b(j)*A(m-i+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       do j=1,bw
          b(i) = b(i) - b(i-bw-1+j)*A(m-bw-1+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m+bw+1-j,n-bw+i)
       enddo
       do j=i+1,i+bw
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m-bw+j-i-1,n-bw+i)
       enddo
       b(n-bw+i) = b(n-bw+i)/A(m,n-bw+i)
    enddo


    ! now dummy is in b
    ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    do i=1,bw-1
       do j=1,i
          b(n-i) = b(n-i) - b(n-i+j)*A(m+j,n-i)
       enddo
    enddo

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       do j=1,bw
          b(i) = b(i) - b(i+j)*A(m+j,i)
       enddo
    enddo

    !b(n) = b(n) -- upper part
    do i=1,bw
       do j=1,bw
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+j,bw+1-i)
       enddo
       do j=bw+1,bw+i
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+i-j,bw+1-i)
       enddo
    enddo

    return
  end subroutine fixed_banbks9_r

  subroutine fixed_banbks10_r(A,n,b)
    implicit none
    integer :: n
    real(c_double) :: A(21,n)
    real(c_double) :: b(n)
    integer :: i,j,k
    integer, parameter :: m = 11
    integer, parameter :: bw = 10

    ! dummy = L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i-1
          b(i) = b(i) - b(j)*A(m-i+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       do j=1,bw
          b(i) = b(i) - b(i-bw-1+j)*A(m-bw-1+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m+bw+1-j,n-bw+i)
       enddo
       do j=i+1,i+bw
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m-bw+j-i-1,n-bw+i)
       enddo
       b(n-bw+i) = b(n-bw+i)/A(m,n-bw+i)
    enddo


    ! now dummy is in b
    ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    do i=1,bw-1
       do j=1,i
          b(n-i) = b(n-i) - b(n-i+j)*A(m+j,n-i)
       enddo
    enddo

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       do j=1,bw
          b(i) = b(i) - b(i+j)*A(m+j,i)
       enddo
    enddo

    !b(n) = b(n) -- upper part
    do i=1,bw
       do j=1,bw
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+j,bw+1-i)
       enddo
       do j=bw+1,bw+i
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+i-j,bw+1-i)
       enddo
    enddo

    return
  end subroutine fixed_banbks10_r

  subroutine fixed_banbks_r(A,n,b,bw)
    implicit none
    integer :: n,bw
    real(c_double) :: A(2*bw+1,n)
    integer :: i,j,k,m
    real(c_double) :: b(n)
    real(c_double) :: dummy

    m = bw+1

    ! dummy = L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i-1
          b(i) = b(i) - b(j)*A(m-i+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    ! dummy =  L^-1 b  -- Center part
    do i=bw+1,n-BW    
       do j=1,bw
          b(i) = b(i) - b(i-bw-1+j)*A(m-bw-1+j,i)
       enddo
       b(i) = b(i)/A(m,i)
    enddo

    !  dummy =  L^-1 b  -- Upper part
    do i=1,bw
       do j=1,i
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m+bw+1-j,n-bw+i)
       enddo
       do j=i+1,i+bw
          b(n-bw+i) = b(n-bw+i) - b(n-2*bw-1+j)*A(m-bw+j-i-1,n-bw+i)
       enddo
       b(n-bw+i) = b(n-bw+i)/A(m,n-bw+i)
    enddo


    ! now dummy is in b
    ! f = U^-1 dummy

    !b(n) = b(n) -- Lower part
    do i=1,bw-1
       do j=1,i
          b(n-i) = b(n-i) - b(n-i+j)*A(m+j,n-i)
       enddo
    enddo

    !b(n) = b(n) -- center part
    do i=n-bw,bw+1,-1
       do j=1,bw
          b(i) = b(i) - b(i+j)*A(m+j,i)
       enddo
    enddo

    !b(n) = b(n) -- upper part
    do i=1,bw
       do j=1,bw
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+j,bw+1-i)
       enddo
       do j=bw+1,bw+i
          b(bw+1-i) = b(bw+1-i) - b(bw+1-i+j)*A(m+i-j,bw+1-i)
       enddo
    enddo

    return
  end subroutine fixed_banbks_r






  subroutine LABS_mat_mul_1_r(A,x,b,n)
    implicit none
    real(c_double) :: A(3,n)
    real(c_double) :: x(n), b(n)
    integer :: i,j,k,n

    ! Upper part
    b(1) = A(2,1)*x(1) + A(3,1)*x(2) + A(1,1)*x(3)     

    ! center part
    do j=2,n-1
       b(j) = + A(1,j)*x(j-1) + A(2,j)*x(j) + A(3,j)*x(j+1) 
    enddo

    ! Lower part
    b(n  ) = A(3,n  )*x(n-2) + A(1,n-1)*x(n-1) + A(2,n  )*x(n)

    return
  end subroutine LABS_mat_mul_1_r


  subroutine LABS_mat_mul_2_r(A,x,b,n)
    implicit none
    real(c_double) :: A(5,n)
    real(c_double) :: x(n), b(n)
    integer :: i,j,k,n

    ! Upper part
    b(1) = A(3,1)*x(1)+A(4,1)*x(2)+A(5,1)*x(3)+A(2,1)*x(4)+A(1,1)*x(5)    
    b(2) = A(2,2)*x(1)+A(3,2)*x(2)+A(4,2)*x(3)+A(5,2)*x(4)+A(1,2)*x(5)    

    ! center part
    do j=3,n-2
       b(j) = A(1,j)*x(j-2)+A(2,j)*x(j-1)+A(3,j)*x(j)+A(4,j)*x(j+1)+A(5,j)*x(j+2) 
    enddo

    ! Lower part
    b(n-1) = A(5,n-1)*x(n-4)+A(1,n-1)*x(n-3)+A(2,n-1)*x(n-2)+A(3,n-1)*x(n-1)+A(4,n-1)*x(n)
    b(n  ) = A(5,n  )*x(n-4)+A(4,n  )*x(n-3)+A(1,n  )*x(n-2)+A(2,n  )*x(n-1)+A(3,n  )*x(n)

    return
  end subroutine LABS_mat_mul_2_r



  subroutine LABS_mat_mul_3_r(A,x,b,n)
    implicit none
    real(c_double) :: A(7,n)
    real(c_double) :: x(n), b(n)
    integer :: i,j,k,n

    ! Upper part
    b(1) = A(4,1)*x(1)+A(5,1)*x(2)+A(6,1)*x(3)+A(7,1)*x(4)+A(3,1)*x(5)+A(2,1)*x(6)+A(1,1)*x(7)    
    b(2) = A(3,2)*x(1)+A(4,2)*x(2)+A(5,2)*x(3)+A(6,2)*x(4)+A(7,2)*x(5)+A(2,2)*x(6)+A(1,2)*x(7)    
    b(3) = A(2,3)*x(1)+A(3,3)*x(2)+A(4,3)*x(3)+A(5,3)*x(4)+A(6,3)*x(5)+A(7,3)*x(6)+A(1,3)*x(7)   

    ! center part
    do j=4,n-3
       b(j) = A(1,j)*x(j-3)+A(2,j)*x(j-2)+A(3,j)*x(j-1)+A(4,j)*x(j)+A(5,j)*x(j+1)+A(6,j)*x(j+2)+A(7,j)*x(j+3)
    enddo

    ! Lower part
    b(n-2) = A(7,n-2)*x(n-6)+A(1,n-2)*x(n-5)+A(2,n-2)*x(n-4)+A(3,n-2)*x(n-3)+A(4,n-2)*x(n-2)+A(5,n-2)*x(n-1)+A(6,n-2)*x(n)
    b(n-1) = A(7,n-1)*x(n-6)+A(6,n-1)*x(n-5)+A(1,n-1)*x(n-4)+A(2,n-1)*x(n-3)+A(3,n-1)*x(n-2)+A(4,n-1)*x(n-1)+A(5,n-1)*x(n)
    b(n  ) = A(7,n  )*x(n-6)+A(6,n  )*x(n-5)+A(5,n  )*x(n-4)+A(1,n  )*x(n-3)+A(2,n  )*x(n-2)+A(3,n  )*x(n-1)+A(4,n  )*x(n)

    return
  end subroutine LABS_mat_mul_3_r


  subroutine LABS_mat_mul_4_r(A,x,b,n)
    implicit none
    real(c_double) :: A(9,n)
    real(c_double) :: x(n), b(n)
    integer :: i,j,k,n
    integer :: m, bw

    m = 5
    bw = 4


    ! Upper part
    b(1) = A(5,1)*x(1)+A(6,1)*x(2)+A(7,1)*x(3)+A(8,1)*x(4)+A(9,1)*x(5)+A(4,1)*x(6)+A(3,1)*x(7)+A(2,1)*x(8)+A(1,1)*x(9)   
    b(2) = A(4,2)*x(1)+A(5,2)*x(2)+A(6,2)*x(3)+A(7,2)*x(4)+A(8,2)*x(5)+A(9,2)*x(6)+A(3,2)*x(7)+A(2,2)*x(8)+A(1,2)*x(9)   
    b(3) = A(3,3)*x(1)+A(4,3)*x(2)+A(5,3)*x(3)+A(6,3)*x(4)+A(7,3)*x(5)+A(8,3)*x(6)+A(9,3)*x(7)+A(2,3)*x(8)+A(1,3)*x(9)  
    b(4) = A(2,4)*x(1)+A(3,4)*x(2)+A(4,4)*x(3)+A(5,4)*x(4)+A(6,4)*x(5)+A(7,4)*x(6)+A(8,4)*x(7)+A(9,4)*x(8)+A(1,4)*x(9)   

    ! center part
    do j=5,n-4
       b(j) = A(1,j)*x(j-4)+A(2,j)*x(j-3)+A(3,j)*x(j-2)+A(4,j)*x(j-1)+A(5,j)*x(j)+A(6,j)*x(j+1)+A(7,j)*x(j+2)+A(8,j)*x(j+3) &
            +A(9,j)*x(j+4)
    enddo


    ! Lower part
    b(n-3) = A(9,n-3)*x(n-8)+A(1,n-3)*x(n-7)+A(2,n-3)*x(n-6)+A(3,n-3)*x(n-5)+A(4,n-3)*x(n-4)+A(5,n-3)*x(n-3)+A(6,n-3)*x(n-2)&
         +A(7,n-3)*x(n-1)+A(8,n-3)*x(n)
    b(n-2) = A(9,n-2)*x(n-8)+A(8,n-2)*x(n-7)+A(1,n-2)*x(n-6)+A(2,n-2)*x(n-5)+A(3,n-2)*x(n-4)+A(4,n-2)*x(n-3)+A(5,n-2)*x(n-2)&
         +A(6,n-2)*x(n-1)+A(7,n-2)*x(n)
    b(n-1) = A(9,n-1)*x(n-8)+A(8,n-1)*x(n-7)+A(7,n-1)*x(n-6)+A(1,n-1)*x(n-5)+A(2,n-1)*x(n-4)+A(3,n-1)*x(n-3)+A(4,n-1)*x(n-2)&
         +A(5,n-1)*x(n-1)+A(6,n-1)*x(n)
    b(n  ) = A(9,n  )*x(n-8)+A(8,n  )*x(n-7)+A(7,n  )*x(n-6)+A(6,n  )*x(n-5)+A(1,n  )*x(n-4)+A(2,n  )*x(n-3)+A(3,n  )*x(n-2)&
         +A(4,n  )*x(n-1)+A(5,n  )*x(n)

    return
  end subroutine LABS_mat_mul_4_r



  subroutine LABS_mat_mul_5_r(A,x,b,n)
    implicit none
    real(c_double) :: A(11,n)
    real(c_double) :: x(n), b(n)
    integer :: i,j,k,n

    ! Upper part
    b(1) = A(6,1)*x(1)+A(7,1)*x(2)+A(8,1)*x(3)+A(9,1)*x(4)+A(10,1)*x(5)+A(11,1)*x(6)+A(5,1)*x(7)+A(4,1)*x(8)+A(3,1)*x(9)&
         +A(2,1)*x(10)+A(1,1)*x(11)  
    b(2) = A(5,2)*x(1)+A(6,2)*x(2)+A(7,2)*x(3)+A(8,2)*x(4)+A(9,2)*x(5)+A(10,2)*x(6)+A(11,2)*x(7)+A(4,2)*x(8)+A(3,2)*x(9)&
         +A(2,2)*x(10)+A(1,2)*x(11)  
    b(3) = A(4,3)*x(1)+A(5,3)*x(2)+A(6,3)*x(3)+A(7,3)*x(4)+A(8,3)*x(5)+A(9,3)*x(6)+A(10,3)*x(7)+A(11,3)*x(8)+A(3,3)*x(9)&
         +A(2,3)*x(10)+A(1,3)*x(11)  
    b(4) = A(3,4)*x(1)+A(4,4)*x(2)+A(5,4)*x(3)+A(6,4)*x(4)+A(7,4)*x(5)+A(8,4)*x(6)+A(9,4)*x(7)+A(10,4)*x(8)+A(11,4)*x(9)&
         +A(2,4)*x(10)+A(1,4)*x(11)  
    b(5) = A(2,5)*x(1)+A(3,5)*x(2)+A(4,5)*x(3)+A(5,5)*x(4)+A(6,5)*x(5)+A(7,5)*x(6)+A(8,5)*x(7)+A(9,5)*x(8)+A(10,5)*x(9)&
         +A(11,5)*x(10)+A(1,5)*x(11)  

    ! center part
    do j=6,n-5
       b(j) = A(1,j)*x(j-5)+A(2,j)*x(j-4)+A(3,j)*x(j-3)+A(4,j)*x(j-2)+A(5,j)*x(j-1)+A(6,j)*x(j)+A(7,j)*x(j+1)+A(8,j)*x(j+2)&
            +A(9,j)*x(j+3)+A(10,j)*x(j+4)+A(11,j)*x(j+5)
    enddo


    ! Lower part
    b(n-4) = A(11,n-4)*x(n-10)+A(1,n-4)*x(n-9)+A(2,n-4)*x(n-8)+A(3,n-4)*x(n-7)+A(4,n-4)*x(n-6)+A(5,n-4)*x(n-5)+A(6,n-4)*x(n-4)&
         +A(7,n-4)*x(n-3)+A(8,n-4)*x(n-2)+A(9,n-4)*x(n-1)+A(10,n-4)*x(n) 
    b(n-3) = A(11,n-3)*x(n-10)+A(10,n-3)*x(n-9)+A(1,n-3)*x(n-8)+A(2,n-3)*x(n-7)+A(3,n-3)*x(n-6)+A(4,n-3)*x(n-5)+A(5,n-3)*x(n-4)&
         +A(6,n-3)*x(n-3)+A(7,n-3)*x(n-2)+A(8,n-3)*x(n-1)+A(9,n-3)*x(n)
    b(n-2) = A(11,n-2)*x(n-10)+A(10,n-2)*x(n-9)+A(9,n-2)*x(n-8)+A(1,n-2)*x(n-7)+A(2,n-2)*x(n-6)+A(3,n-2)*x(n-5)+A(4,n-2)*x(n-4)&
         +A(5,n-2)*x(n-3)+A(6,n-2)*x(n-2)+A(7,n-2)*x(n-1)+A(8,n-2)*x(n)
    b(n-1) = A(11,n-1)*x(n-10)+A(10,n-1)*x(n-9)+A(9,n-1)*x(n-8)+A(8,n-1)*x(n-7)+A(1,n-1)*x(n-6)+A(2,n-1)*x(n-5)+A(3,n-1)*x(n-4)&
         +A(4,n-1)*x(n-3)+A(5,n-1)*x(n-2)+A(6,n-1)*x(n-1)+A(7,n-1)*x(n)
    b(n  ) = A(11,n  )*x(n-10)+A(10,n  )*x(n-9)+A(9,n  )*x(n-8)+A(8,n  )*x(n-7)+A(7,n  )*x(n-6)+A(1,n  )*x(n-5)+A(2,n  )*x(n-4)&
         +A(3,n  )*x(n-3)+A(4,n  )*x(n-2)+A(5,n  )*x(n-1)+A(6,n  )*x(n)

    return
  end subroutine LABS_mat_mul_5_r


  subroutine LABS_mat_mul_6_r(A,x,b,n)
    implicit none
    real(c_double) :: A(13,n)
    real(c_double) :: x(n), b(n)
    integer :: i,j,k,n


    ! Upper part
    b(1) = A(7,1)*x(1)+A(8,1)*x(2)+A(9,1)*x(3)+A(10,1)*x(4)+A(11,1)*x(5)+A(12,1)*x(6)+A(13,1)*x(7)+A(6,1)*x(8)+A(5,1)*x(9)&
         +A(4,1)*x(10)+A(3,1)*x(11)+A(2,1)*x(12)+A(1,1)*x(13)  
    b(2) = A(6,2)*x(1)+A(7,2)*x(2)+A(8,2)*x(3)+A(9,2)*x(4)+A(10,2)*x(5)+A(11,2)*x(6)+A(12,2)*x(7)+A(13,2)*x(8)+A(5,2)*x(9)&
         +A(4,2)*x(10)+A(3,2)*x(11)+A(2,2)*x(12)+A(1,2)*x(13)  
    b(3) = A(5,3)*x(1)+A(6,3)*x(2)+A(7,3)*x(3)+A(8,3)*x(4)+A(9,3)*x(5)+A(10,3)*x(6)+A(11,3)*x(7)+A(12,3)*x(8)+A(13,3)*x(9)&
         +A(4,3)*x(10)+A(3,3)*x(11)+A(2,3)*x(12)+A(1,3)*x(13)  
    b(4) = A(4,4)*x(1)+A(5,4)*x(2)+A(6,4)*x(3)+A(7,4)*x(4)+A(8,4)*x(5)+A(9,4)*x(6)+A(10,4)*x(7)+A(11,4)*x(8)+A(12,4)*x(9)&
         +A(13,4)*x(10)+A(3,4)*x(11)+A(2,4)*x(12)+A(1,4)*x(13)  
    b(5) = A(3,5)*x(1)+A(4,5)*x(2)+A(5,5)*x(3)+A(6,5)*x(4)+A(7,5)*x(5)+A(8,5)*x(6)+A(9,5)*x(7)+A(10,5)*x(8)+A(11,5)*x(9)&
         +A(12,5)*x(10)+A(13,5)*x(11)+A(2,5)*x(12)+A(1,5)*x(13)  
    b(6) = A(2,6)*x(1)+A(3,6)*x(2)+A(4,6)*x(3)+A(5,6)*x(4)+A(6,6)*x(5)+A(7,6)*x(6)+A(8,6)*x(7)+A(9,6)*x(8)+A(10,6)*x(9)&
         +A(11,6)*x(10)+A(12,6)*x(11)+A(13,6)*x(12)+A(1,6)*x(13)  

    ! center part
    do j=7,n-6
       b(j) = A(1,j)*x(j-6)+A(2,j)*x(j-5)+A(3,j)*x(j-4)+A(4,j)*x(j-3)+A(5,j)*x(j-2)+A(6,j)*x(j-1)+A(7,j)*x(j)+A(8,j)*x(j+1)&
            +A(9,j)*x(j+2)+A(10,j)*x(j+3)+A(11,j)*x(j+4)+A(12,j)*x(j+5)+A(13,j)*x(j+6)  
    enddo


    ! Lower part
    b(n-5) = A(13,n-5)*x(n-12)+A(1,n-5)*x(n-11)+A(2,n-5)*x(n-10)+A(3,n-5)*x(n-9)+A(4,n-5)*x(n-8)+A(5,n-5)*x(n-7)+A(6,n-5)*x(n-6)&
         +A(7,n-5)*x(n-5)+A(8,n-5)*x(n-4)+A(9,n-5)*x(n-3)+A(10,n-5)*x(n-2)+A(11,n-5)*x(n-1)+A(12,n-5)*x(n) 
    b(n-4) = A(13,n-4)*x(n-12)+A(12,n-4)*x(n-11)+A(1,n-4)*x(n-10)+A(2,n-4)*x(n-9)+A(3,n-4)*x(n-8)+A(4,n-4)*x(n-7)+A(5,n-4)*x(n-6)&
         +A(6,n-4)*x(n-5)+A(7,n-4)*x(n-4)+A(8,n-4)*x(n-3)+A(9,n-4)*x(n-2)+A(10,n-4)*x(n-1)+A(11,n-4)*x(n) 
    b(n-3) = A(13,n-3)*x(n-12)+A(12,n-3)*x(n-11)+A(11,n-3)*x(n-10)+A(1,n-3)*x(n-9)+A(2,n-3)*x(n-8)+A(3,n-3)*x(n-7)+A(4,n-3)*x(n-6)&
         +A(5,n-3)*x(n-5)+A(6,n-3)*x(n-4)+A(7,n-3)*x(n-3)+A(8,n-3)*x(n-2)+A(9,n-3)*x(n-1)+A(10,n-3)*x(n) 
    b(n-2) = A(13,n-2)*x(n-12)+A(12,n-2)*x(n-11)+A(11,n-2)*x(n-10)+A(10,n-2)*x(n-9)+A(1,n-2)*x(n-8)+A(2,n-2)*x(n-7)+A(3,n-2)*x(n-6)&
         +A(4,n-2)*x(n-5)+A(5,n-2)*x(n-4)+A(6,n-2)*x(n-3)+A(7,n-2)*x(n-2)+A(8,n-2)*x(n-1)+A(9,n-2)*x(n) 
    b(n-1) = A(13,n-1)*x(n-12)+A(12,n-1)*x(n-11)+A(11,n-1)*x(n-10)+A(10,n-1)*x(n-9)+A(9,n-1)*x(n-8)+A(1,n-1)*x(n-7)+A(2,n-1)*x(n-6)&
         +A(3,n-1)*x(n-5)+A(4,n-1)*x(n-4)+A(5,n-1)*x(n-3)+A(6,n-1)*x(n-2)+A(7,n-1)*x(n-1)+A(8,n-1)*x(n) 
    b(n  ) = A(13,n  )*x(n-12)+A(12,n  )*x(n-11)+A(11,n  )*x(n-10)+A(10,n  )*x(n-9)+A(9,n  )*x(n-8)+A(8,n  )*x(n-7)+A(1,n  )*x(n-6)&
         +A(2,n  )*x(n-5)+A(3,n  )*x(n-4)+A(4,n  )*x(n-3)+A(5,n  )*x(n-2)+A(6,n  )*x(n-1)+A(7,n  )*x(n) 

    return
  end subroutine LABS_mat_mul_6_r
  
  subroutine LABS_mat_mul_7_r(A,x,b,n)
    implicit none
    real(c_double) :: A(15,n)
    real(c_double) :: x(n), b(n)
    integer :: i,j,k,n

    ! Upper part
    b(1) = A(8,1)*x(1)+A(9,1)*x(2)+A(10,1)*x(3)+A(11,1)*x(4)+A(12,1)*x(5)+A(13,1)*x(6)+A(14,1)*x(7)+A(15,1)*x(8)+A(7,1)*x(9)&
         +A(6,1)*x(10)+A(5,1)*x(11)+A(4,1)*x(12)+A(3,1)*x(13)+A(2,1)*x(14)+A(1,1)*x(15) 
    b(2) = A(7,2)*x(1)+A(8,2)*x(2)+A(9,2)*x(3)+A(10,2)*x(4)+A(11,2)*x(5)+A(12,2)*x(6)+A(13,2)*x(7)+A(14,2)*x(8)+A(15,2)*x(9)&
         +A(6,2)*x(10)+A(5,2)*x(11)+A(4,2)*x(12)+A(3,2)*x(13)+A(2,2)*x(14)+A(1,2)*x(15) 
    b(3) = A(6,3)*x(1)+A(7,3)*x(2)+A(8,3)*x(3)+A(9,3)*x(4)+A(10,3)*x(5)+A(11,3)*x(6)+A(12,3)*x(7)+A(13,3)*x(8)+A(14,3)*x(9)&
         +A(15,3)*x(10)+A(5,3)*x(11)+A(4,3)*x(12)+A(3,3)*x(13)+A(2,3)*x(14)+A(1,3)*x(15) 
    b(4) = A(5,4)*x(1)+A(6,4)*x(2)+A(7,4)*x(3)+A(8,4)*x(4)+A(9,4)*x(5)+A(10,4)*x(6)+A(11,4)*x(7)+A(12,4)*x(8)+A(13,4)*x(9)&
         +A(14,4)*x(10)+A(15,4)*x(11)+A(4,4)*x(12)+A(3,4)*x(13)+A(2,4)*x(14)+A(1,4)*x(15) 
    b(5) = A(4,5)*x(1)+A(5,5)*x(2)+A(6,5)*x(3)+A(7,5)*x(4)+A(8,5)*x(5)+A(9,5)*x(6)+A(10,5)*x(7)+A(11,5)*x(8)+A(12,5)*x(9)&
         +A(13,5)*x(10)+A(14,5)*x(11)+A(15,5)*x(12)+A(3,5)*x(13)+A(2,5)*x(14)+A(1,5)*x(15) 
    b(6) = A(3,6)*x(1)+A(4,6)*x(2)+A(5,6)*x(3)+A(6,6)*x(4)+A(7,6)*x(5)+A(8,6)*x(6)+A(9,6)*x(7)+A(10,6)*x(8)+A(11,6)*x(9)&
         +A(12,6)*x(10)+A(13,6)*x(11)+A(14,6)*x(12)+A(15,6)*x(13)+A(2,6)*x(14)+A(1,6)*x(15) 
    b(7) = A(2,7)*x(1)+A(3,7)*x(2)+A(4,7)*x(3)+A(5,7)*x(4)+A(6,7)*x(5)+A(7,7)*x(6)+A(8,7)*x(7)+A(9,7)*x(8)+A(10,7)*x(9)&
         +A(11,7)*x(10)+A(12,7)*x(11)+A(13,7)*x(12)+A(14,7)*x(13)+A(15,7)*x(14)+A(1,7)*x(15) 

    ! center part
    do j=8,n-7
       b(j) = A(1,j)*x(j-7)+A(2,j)*x(j-6)+A(3,j)*x(j-5)+A(4,j)*x(j-4)+A(5,j)*x(j-3)+A(6,j)*x(j-2)+A(7,j)*x(j-1)+A(8,j)*x(j)&
            +A(9,j)*x(j+1)+A(10,j)*x(j+2)+A(11,j)*x(j+3)+A(12,j)*x(j+4)+A(13,j)*x(j+5)+A(14,j)*x(j+6)+A(15,j)*x(j+7) 
    enddo


    ! Lower part
    b(n-6) = A(15,n-6)*x(n-14)+A(1,n-6)*x(n-13)+A(2,n-6)*x(n-12)+A(3,n-6)*x(n-11)+A(4,n-6)*x(n-10)+A(5,n-6)*x(n-9)+A(6,n-6)*x(n-8)&
         +A(7,n-6)*x(n-7)+A(8,n-6)*x(n-6)+A(9,n-6)*x(n-5)+A(10,n-6)*x(n-4)+A(11,n-6)*x(n-3)+A(12,n-6)*x(n-2)+A(13,n-6)*x(n-1)&
         +A(14,n-6)*x(n) 
    b(n-5) = A(15,n-5)*x(n-14)+A(14,n-5)*x(n-13)+A(1,n-5)*x(n-12)+A(2,n-5)*x(n-11)+A(3,n-5)*x(n-10)+A(4,n-5)*x(n-9)+A(5,n-5)*x(n-8)&
         +A(6,n-5)*x(n-7)+A(7,n-5)*x(n-6)+A(8,n-5)*x(n-5)+A(9,n-5)*x(n-4)+A(10,n-5)*x(n-3)+A(11,n-5)*x(n-2)+A(12,n-5)*x(n-1)&
         +A(13,n-5)*x(n) 
    b(n-4) = A(15,n-4)*x(n-14)+A(14,n-4)*x(n-13)+A(13,n-4)*x(n-12)+A(1,n-4)*x(n-11)+A(2,n-4)*x(n-10)+A(3,n-4)*x(n-9) &
         +A(4,n-4)*x(n-8)+ A(5,n-4)*x(n-7)+A(6,n-4)*x(n-6)+A(7,n-4)*x(n-5)+A(8,n-4)*x(n-4)+A(9,n-4)*x(n-3)+A(10,n-4)*x(n-2)&
         +A(11,n-4)*x(n-1)+ A(12,n-4)*x(n) 
    b(n-3) = A(15,n-3)*x(n-14)+A(14,n-3)*x(n-13)+A(13,n-3)*x(n-12)+A(12,n-3)*x(n-11)+A(1,n-3)*x(n-10)+A(2,n-3)*x(n-9)&
         +A(3,n-3)*x(n-8)+A(4,n-3)*x(n-7)+A(5,n-3)*x(n-6)+A(6,n-3)*x(n-5)+A(7,n-3)*x(n-4)+A(8,n-3)*x(n-3)+A(9,n-3)*x(n-2)&
         +A(10,n-3)*x(n-1) +A(11,n-3)*x(n) 
    b(n-2) = A(15,n-2)*x(n-14)+A(14,n-2)*x(n-13)+A(13,n-2)*x(n-12)+A(12,n-2)*x(n-11)+A(11,n-2)*x(n-10)+A(1,n-2)*x(n-9)&
         +A(2,n-2)*x(n-8)+A(3,n-2)*x(n-7)+A(4,n-2)*x(n-6)+A(5,n-2)*x(n-5)+A(6,n-2)*x(n-4)+A(7,n-2)*x(n-3)+A(8,n-2)*x(n-2)&
         +A(9,n-2)*x(n-1)+A(10,n-2)*x(n) 
    b(n-1) = A(15,n-1)*x(n-14)+A(14,n-1)*x(n-13)+A(13,n-1)*x(n-12)+A(12,n-1)*x(n-11)+A(11,n-1)*x(n-10)+A(10,n-1)*x(n-9)&
         +A(1,n-1)*x(n-8)+A(2,n-1)*x(n-7)+A(3,n-1)*x(n-6)+A(4,n-1)*x(n-5)+A(5,n-1)*x(n-4)+A(6,n-1)*x(n-3)+A(7,n-1)*x(n-2)&
         +A(8,n-1)*x(n-1)+A(9,n-1)*x(n) 
    b(n  ) = A(15,n  )*x(n-14)+A(14,n  )*x(n-13)+A(13,n  )*x(n-12)+A(12,n  )*x(n-11)+A(11,n)*x(n-10)+A(10,n)*x(n-9)+A(9,n)*x(n-8)&
         +A(1,n  )*x(n-7)+A(2,n  )*x(n-6)+A(3,n  )*x(n-5)+A(4,n  )*x(n-4)+A(5,n  )*x(n-3)+A(6,n  )*x(n-2)+A(7,n  )*x(n-1)&
         +A(8,n  )*x(n) 

    return
  end subroutine LABS_mat_mul_7_r


  subroutine LABS_mat_mul_8_r(A,x,b,n)
    implicit none
    integer, parameter :: bw = 8
    integer :: n
    real(c_double) :: A(2*bw+1,n)
    real(c_double) :: x(n), b(n)
    integer :: i,j,k

    do j=1,bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=1-j,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,bw+1-j
          b(j) = b(j) + A(bw+2-j-i,j)*x(j+bw+i)
       enddo
    enddo

    do j=bw+1,n-bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
    enddo

    do j=n-bw+1,n
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,n-j
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,j-n+bw
          b(j) = b(j) + A(n+bw+1-j+i,j)*x(j-bw-i)
       enddo
    enddo
    return
  end subroutine LABS_mat_mul_8_r

  subroutine LABS_mat_mul_9_r(A,x,b,n)
    implicit none
    integer, parameter :: bw = 9
    integer :: n
    real(c_double) :: A(2*bw+1,n)
    real(c_double) :: x(n), b(n)
    integer :: i,j,k

    do j=1,bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=1-j,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,bw+1-j
          b(j) = b(j) + A(bw+2-j-i,j)*x(j+bw+i)
       enddo
    enddo

    do j=bw+1,n-bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
    enddo

    do j=n-bw+1,n
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,n-j
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,j-n+bw
          b(j) = b(j) + A(n+bw+1-j+i,j)*x(j-bw-i)
       enddo
    enddo
    return
  end subroutine LABS_mat_mul_9_r

  subroutine LABS_mat_mul_10_r(A,x,b,n)
    implicit none
    integer, parameter :: bw = 10
    integer :: n
    real(c_double) :: A(2*bw+1,n)
    real(c_double) :: x(n), b(n)
    integer :: i,j,k

    do j=1,bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=1-j,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,bw+1-j
          b(j) = b(j) + A(bw+2-j-i,j)*x(j+bw+i)
       enddo
    enddo

    do j=bw+1,n-bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
    enddo

    do j=n-bw+1,n
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,n-j
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,j-n+bw
          b(j) = b(j) + A(n+bw+1-j+i,j)*x(j-bw-i)
       enddo
    enddo
    return
  end subroutine LABS_mat_mul_10_r

  subroutine LABS_mat_mul_bw_r(A,x,b,n,bw)
    implicit none
    integer :: bw
    integer :: n
    real(c_double) :: A(2*bw+1,n)
    real(c_double) :: x(n), b(n)
    integer :: i,j,k

    do j=1,bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=1-j,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,bw+1-j
          b(j) = b(j) + A(bw+2-j-i,j)*x(j+bw+i)
       enddo
    enddo

    do j=bw+1,n-bw
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,bw
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
    enddo

    do j=n-bw+1,n
       b(j) = dcmplx(0.d0,0.d0)
       do i=-bw,n-j
          b(j) = b(j) + A(bw+1+i,j)*x(j+i)
       enddo
       do i=1,j-n+bw
          b(j) = b(j) + A(n+bw+1-j+i,j)*x(j-bw-i)
       enddo
    enddo
    return
  end subroutine LABS_mat_mul_bw_r

end module LABS








