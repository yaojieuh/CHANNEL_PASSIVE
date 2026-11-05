!-----------------------------------------------------------
! This routine calculate the U and W
!
! input 
! u1 = g(=du/dz - dw/dx)
! u3 = f(=dv/dy)
! kx, kz : wave numbers
! 
! output 
! u1 = u
! u3 = w
! Developer : M.K. Lee (mk@ices.utexas.edu) 
! Last update : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!-----------------------------------------------------------

subroutine solve_u_w(g,f,kx,kz)
  use labs
  implicit none
  complex(c_double_complex), dimension(LABS_NY) :: g,f
  real(c_double) :: kx,kz
  complex(c_double_complex), parameter :: ci = (0d0,1d0)
  integer :: i
  complex(c_double_complex) :: kxx, kzz, dxg, dxf, dzg, dzf

  kxx = ci*kx/(kx*kx+kz*kz)
  kzz = ci*kz/(kx*kx+kz*kz)
 

  do i=1,LABS_NY
     dxg = kxx*g(i)
     dzg = kzz*g(i)
     dxf = kxx*f(i)
     dzf = kzz*f(i)
     g(i) = dxf-dzg
     f(i) = dxg+dzf
  enddo


  return
end subroutine solve_u_w
