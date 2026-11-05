!--------------------------------------------------------------------
! This subroutine perform inverse DFT in z direction 
! and data transpose from y pencil to z pencil 
! and from z pencil to x pencil
! 
! It is necessary to use this subroutine because data type conversion

! Detailed Procedure
! 
! - Global data transpose from y pencil to z pencil
! - Inverse DFT in Z direction with data padding for 3/2 dealiasing
! - Global data transpose from z pencil to x pencil
! 
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------

subroutine y_to_x(U_r,U_c,Upad_r,Upad_c)
  use ppdfftw
  implicit none

  complex(C_DOUBLE_COMPLEX) :: Upad_c(nz32,Zxsz*Zysz)
  complex(C_DOUBLE_COMPLEX) :: U_c(nz,Zxsz*Zysz)

  real(C_DOUBLE) :: Upad_r(alloc_size1)
  real(C_DOUBLE) :: U_r(alloc_size2)


 
  ! Execute Transpose from y pencil to z pencil
  call ppdfftw_y_to_z(u_r)
  
  ! Execute FFT complex to complex in z pencil (inverse)
!!$  if(myrank.eq.0) time_transform_start = MPI_WTIME()
  call ppdfftw_dft_z_inverse(u_c,Upad_c,buf_line_z)
!!$  if(myrank.eq.0) time_transform_end = MPI_WTIME()
!!$  if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "DFT in Z (inverse)", &
!!$       timestep,dt,time,time_transform_end-time_transform_start
  ! Now, data is in Upad.

  ! Execute Transpose from z pencil to x pencil
  call ppdfftw_z_to_x(Upad_r)



  return
end subroutine y_to_x
  
  
  
