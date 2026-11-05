!--------------------------------------------------------------------
! This module contains variables to be used 
! for time advancing with Runge Kutta method.
! Now only support low-storage implicit/explicit 3rd order scheme.
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------
module RUNGE_KUTTA
  use,intrinsic :: iso_c_binding
  implicit none

  integer :: RK_SCHEME
  integer :: RK_STEP
  integer :: RK_SUBSTEP
  real(c_double) :: dt, time

  real(c_double), allocatable, dimension(:,:) :: RK_3_COEF
  real(c_double), allocatable, dimension(:) :: phi_dpdx
  real(c_double), allocatable, dimension(:) :: phi_dtdx
  

contains


  subroutine RK3_initialize
    use, intrinsic :: iso_c_binding
    implicit none

    allocate(RK_3_COEF(4,3))
  
    RK_3_COEF(1,1) =  29d0 /  96d0 ! Alpha 1
    RK_3_COEF(1,2) = - 3d0 /  40d0 ! Alpha 2
    RK_3_COEF(1,3) =   1d0 /   6d0 ! Alpha 3

    RK_3_COEF(2,1) =  37d0 / 160d0 ! Beta  1
    RK_3_COEF(2,2) =   5d0 /  24d0 ! Beta  2
    RK_3_COEF(2,3) =   1d0 /   6d0 ! Beta  3

    RK_3_COEF(3,1) =   8d0 /  15d0 ! Gamma 1
    RK_3_COEF(3,2) =   5d0 /  12d0 ! Gamma 2
    RK_3_COEF(3,3) =   3d0 /   4d0 ! Gamma 3

    RK_3_COEF(4,1) =   0d0         ! Zeta  1
    RK_3_COEF(4,2) = -17d0 /  60d0 ! Zeta  2
    RK_3_COEF(4,3) = - 5d0 /  12d0 ! Zeta  3
    
    allocate(phi_dpdx(4))
    allocate(phi_dtdx(4))

    return
  end subroutine RK3_initialize


  subroutine RK3_finalize
    use, intrinsic :: iso_c_binding
    implicit none

    deallocate(RK_3_COEF)
    deallocate(phi_dpdx)
    deallocate(phi_dtdx)

    return
  end subroutine RK3_finalize



end module RUNGE_KUTTA
