!--------------------------------------------------------------------
! This subroutine rescale the data after DFT
!
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------
subroutine rescale_data(U,ysz,zsz,xsz)
  use, intrinsic :: iso_c_binding
  use runge_kutta
  use global_variable
  implicit none
  integer(C_INTPTR_T) :: ysz,zsz,xsz
  real(C_DOUBLE) :: U(2*ysz*xsz*zsz*10)
!!$  complex(C_DOUBLE_COMPLEX) :: U(ysz*xsz*zsz*6)
  integer :: n_total,i



  if(stat_store_or_not .eqv. .true.) then
     n_total = 2*zsz*xsz*ysz*10
  else
     n_total = 2*zsz*xsz*ysz*10
  endif
  call omp_set_num_threads(omp_compute)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i) 
  !$OMP DO SCHEDULE(STATIC)
  do i=1,n_total
     U(i) = rescale*U(i)
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL

  return

end subroutine rescale_data
