!-----------------------------------------------------------
! This routine calculate the Hg and Hv

! input 
! u1 = u^2 - v^2
! u2 = uv
! u3 = w^2 - v^2
! u4 = vw
! u5 = uw
! D0_lu : Factorized banded matrix for D0
! D1 : B-Spline matrix for d/dy operation
! D2 : B-Spline matrix for d^2/dy^2 operation
! kx, kz : wave numbers
! 
! output 
! u1 = hg
! u2 = garbage
! u3 = garbage
! u4 = garbage
! u5 = hv
!
! Developer : M.K. Lee (mk@ices.utexas.edu) 
! Last update : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!-----------------------------------------------------------

subroutine hvhg(u1,u2,u3,u4,u5,u7,u8,u9,D0_lu,D1,D2,kx,kz,nx32,nz32)
  use labs
  implicit none
  complex(c_double_complex), dimension(LABS_NY) :: u1,u2,u3,u4,u5,u7,u8,u9
  real(c_double), dimension(2*LABS_BW+1,LABS_NY) :: D0_lu,D1,D2
  real(c_double) :: kx,kz,kz2,kx2,kxkz,scaling
  integer(c_intptr_t) :: nx32, nz32
  complex(c_double_complex), parameter :: ci = dcmplx(0.d0,1.d0)
  integer :: i
  complex(c_double_complex) :: dummy,ikx,ikz,ikx3,ikxkz2,ikx2kz,ikz3


  kz2 = kz**2
  kx2 = kx**2
  kxkz = kx*kz
  ikx = ci*kx
  ikz = ci*kz
  ikx3 = ci*kx**3
  ikxkz2 = ci*kx*(kz**2)
  ikx2kz = ci*(kx**2)*kz
  ikz3 = ci*(kz**3)
  
!  u1 = u^2 - v^2
!  u2 = uv
!  u3 = w^2 - v^2
!  u4 = vw
!  u5 = uw

  do i=1,LABS_NY
     dummy = u1(i)
     u1(i) = kxkz*(u1(i) - u3(i)) + (kz2-kx2)*u5(i)
     u3(i) = -kx2*dummy - 2d0*kxkz*u5(i) -kz2*u3(i)
  enddo

!  u1 = kx*kz*(u^2-w^2) + (kz^2 - kx^2)*uw
!  u2 = uv
!  u3 = -kx^2*(u^2 - v^2) -2*kx*kz*uw - kz^2*(w^2 - v^2)
!  u4 = vw
!  u5 = uw <- Don't care

  call LABS_bs(D0_lu,u3)
  call LABS_mat_mul(D1,u3,u5)

!  u1 = kx*kz*(u^2-w^2) + (kz^2 - kx^2)*uw
!  u2 = uv
!  u3 = Don't care
!  u4 = vw
!  u5 = d/dy (-kx^2*(u^2 - v^2) -2*kx*kz*uw - kz^2*(w^2 - v^2))

  do i=1,LABS_NY
     u5(i) = u5(i) + (ikx3+ikxkz2)*u2(i) + (ikx2kz + ikz3)*u4(i)
  enddo

!  u1 = kx*kz*(u^2-w^2) + (kz^2 - kx^2)*uw
!  u2 = uv
!  u3 = Don't care
!  u4 = vw
!  u5 = d/dy (-kx^2*(u^2 - v^2) -2*kx*kz*uw - kz^2*(w^2 - v^2))
!        + (i*kx^3 + i*kx*kz^2)*uv + (i*kx^2*kz + i*kz^3)*vw

  do i=1,LABS_NY
     u3(i) = ikx*u4(i) - ikz*u2(i)
     u4(i) = ikz*u4(i) + ikx*u2(i)
  enddo

!  u1 = kx*kz*(u^2-w^2) + (kz^2 - kx^2)*uw
!  u2 = Don't care
!  u3 = i*kx*vw - i*kz*uv
!  u4 = i*kz*vw + i*kx*uv
!  u5 = d/dy (-kx^2*(u^2 - v^2) -2*kx*kz*uw - kz^2*(w^2 - v^2))
!        + (i*kx^3 + i*kx*kz^2)*uv + (i*kx^2*kz + i*kz^3)*vw

  call LABS_bs(D0_lu,u3)
  call LABS_mat_mul(D1,u3,u2)
  
!  u1 = kx*kz*(u^2-w^2) + (kz^2 - kx^2)*uw
!  u2 = d/dy (i*kx*vw - i*kz*uv)
!  u3 = Don't care
!  u4 = i*kz*vw + i*kx*uv
!  u5 = d/dy (-kx^2*(u^2 - v^2) -2*kx*kz*uw - kz^2*(w^2 - v^2))
!        + (i*kx^3 + i*kx*kz^2)*uv + (i*kx^2*kz + i*kz^3)*vw

  call LABS_bs(D0_lu,u4)
  call LABS_mat_mul(D2,u4,u3)
  
!  u1 = kx*kz*(u^2-w^2) + (kz^2 - kx^2)*uw
!  u2 = d/dy (i*kx*vw - i*kz*uv)
!  u3 = d^2/dy^2 (i*kz*vw + i*kx*uv) 
!  u4 = Don't care
!  u5 = d/dy (-kx^2*(u^2 - v^2) -2*kx*kz*uw - kz^2*(w^2 - v^2))
!        + (i*kx^3 + i*kx*kz^2)*uv + (i*kx^2*kz + i*kz^3)*vw

  do i=1,LABS_NY
     u5(i) = u5(i) + u3(i)
     u1(i) = u1(i) + u2(i)
  enddo

!  u1 = kx*kz*(u^2-w^2) + (kz^2 - kx^2)*uw
!        + d/dy (i*kx*vw - i*kz*uv)
!     = hg  
!
!  u2 = Don't care
!  u3 = Don't care
!  u4 = Don't care
!
!  u5 = d/dy (-kx^2*(u^2 - v^2) -2*kx*kz*uw - kz^2*(w^2 - v^2))
!        + (i*kx^3 + i*kx*kz^2)*uv + (i*kx^2*kz + i*kz^3)*vw
!        + d^2/dy^2 (i*kz*vw + i*kx*uv) 
!     = hv


!  u7 = ut
!  u8 = vt
!  u9 = wt

  do i=1,LABS_NY
     u7(i) = ikx*u7(i) + ikz*u9(i)
  enddo
!  u7 = i*kx*ut+i*kz*wt
!  u8 = vt
!  u9 = Don't care

  call LABS_bs(D0_lu,u8)
  call LABS_mat_mul(D1,u8,u9)

!  u7 = i*kx*ut+i*kz*wt
!  u8 = vt
!  u9 = dvt/dy

  do i=1,LABS_NY
     u7(i) = -(u7(i) + u9(i))
  enddo
!  u7 = ht

  return
end subroutine hvhg
