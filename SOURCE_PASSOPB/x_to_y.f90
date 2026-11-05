!--------------------------------------------------------------------
! This subroutine perform forward DFT in z direction 
! and data transpose from x pencil to z pencil 
! and from z pencil to y pencil
! 
! It is necessary to use this subroutine because data type conversion

! Detailed Procedure
! 
! - Global data transpose from x pencil to z pencil
! - Forward DFT in Z direction with data packing for 3/2 dealiasing
! - Global data transpose from z pencil to y pencil
! 
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------

subroutine x_to_y(Upad_r,Upad_c,U_r,U_c)
  use ppdfftw
  implicit none

  complex(C_DOUBLE_COMPLEX) :: Upad_c(nz32,Zxsz*Zysz)
  complex(C_DOUBLE_COMPLEX) :: U_c(nz,Zxsz*Zysz)

  real(C_DOUBLE) :: Upad_r(alloc_size1)
  real(C_DOUBLE) :: U_r(alloc_size2)

  ! Execute Transpose from x pencil to z pencil
  call ppdfftw_x_to_z(Upad_r)

  ! Execute FFT complex to complex in z pencil (forward)
!!$  if(myrank.eq.0) time_transform_start = MPI_WTIME()
  call ppdfftw_dft_z_forward(Upad_c,u_c,buf_line_z)
!!$  if(myrank.eq.0) time_transform_end = MPI_WTIME()
!!$  if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "DFT in Z (forward)", &
!!$       timestep,dt,time,time_transform_end-time_transform_start
  ! Now, data is in U.

  ! Execute Transpose from z pencil to y pencil
  call ppdfftw_z_to_y(u_r)

  return
end subroutine x_to_y
  
  
  
