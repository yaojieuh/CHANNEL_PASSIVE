!-----------------------------------------------------------
! This routine calculate the g (=du/dz - dw/dx)
!
! input 
! u1 = hg
! store_g = Non linear storage term for g from previous step (except for rk_step=0)
! buf1 = Don't care
! buf2 = Don't care
! D0_lu : Factorized banded matrix for D0
! D0 : B-Spline matrix 
! D2 : B-Spline matrix for d^2/dy^2 operation
! D2_lu : Factorized banded matrix for D2
! D_buf : Scratch space to solve linear equation for time step
! nu : viscosity
! kx, kz : wave numbers
! ny : number of grid points in y direction
!
! output 
! u1 = g_new (=du/dz - dw/dx))
! store_g = Low storage terms of g for next step 
! D_buf : LU-decomposed matrix (will be reused to calculate phi)
!
! Developer : M.K. Lee (mk@ices.utexas.edu) 
! Last update : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!-----------------------------------------------------------
subroutine solve_g(u1,store_g,buf1,buf2, D0_LU, D0, D2, D_buf, kx, kz, nu, ny)

  use labs
  use runge_kutta
  implicit none
  complex(c_double_complex), dimension(ny) :: u1,store_g,buf1,buf2
  real(c_double), dimension(2*LABS_BW+1,ny) :: D0_LU,D0,D2,D_buf
  real(c_double) :: kx,kz,k2, alpha, beta, gamma, zeta, dummy,nu
  integer :: i,j,ny
  k2 = kx**2 + kz**2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part1 .Calculate the LHS matrix and get LU decomposed matrix of it
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  beta  = dt*RK_3_COEF(2,RK_SUBSTEP)*nu
  dummy = 1d0 + beta*k2

  ! To enforce boundary condition of first element
  do i=1,LABS_BW
     D_buf(i,1) = 0d0
  enddo
  D_buf(LABS_BW+1,1) = dummy*D0(LABS_BW+1,1) - beta*D2(LABS_BW+1,1)
  do i=LABS_BW+2,2*LABS_BW+1
     D_buf(i,1) = 0d0
  enddo

  do j = 2,ny-1
     do i=1,2*LABS_BW+1
        D_buf(i,j) = (dummy*D0(i,j) - beta*D2(i,j))
     enddo
  enddo

  ! To enforce boundary condition of last element
  do i=1,LABS_BW
     D_buf(i,ny) = 0d0
  enddo
  D_buf(LABS_BW+1,ny) = dummy*D0(LABS_BW+1,ny) - beta*D2(LABS_BW+1,ny)
  do i=LABS_BW+2,2*LABS_BW+1
     D_buf(i,ny) = 0d0
  enddo
  
  
  ! LU decomposition of LHS mass matrix
  call LABS_lu(D_buf)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part1 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=





  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part2. Calculate the RHS to calculate new G
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (RK_SUBSTEP.eq.1) then
     alpha = dt*RK_3_COEF(1,RK_SUBSTEP)*nu
     gamma = dt*RK_3_COEF(3,RK_SUBSTEP)
     !---------------------
     ! u1 = hg
     ! store_g = g_old
     ! buf1 = Don't care
     ! buf2 = Don't care
     !---------------------
     dummy = 1.d0-alpha*k2
     do i=1,ny
        buf1(i) = dummy*store_g(i) + gamma*u1(i)
     enddo
     call LABS_bs(D0_lu,store_g)
     call LABS_mat_mul(D2,store_g,buf2)
     do i=1,ny
        buf1(i) = buf1(i) + alpha*buf2(i)
     enddo
  else
     !---------------------
     ! u1 = hg
     ! store_g = Non linear storage term for g from previous step (except for rk_step=1)
     ! buf1 = Don't care
     ! buf2 = Don't care
     !---------------------
     gamma = dt*RK_3_COEF(3,RK_SUBSTEP)
     do i=1,ny
        buf1(i) = store_g(i) + gamma*u1(i)
     enddo
     !---------------------
     ! u1 = hg
     ! store_g = Don't care
     ! buf1 = RHS for g
     ! buf2 = Don't care
     !---------------------
  endif

  ! boundary condition of g(y=+- 1) = 0
  buf1(1) = (0d0,0d0)
  buf1(ny) = (0d0,0d0)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  
  !---------------------
  ! u1 = hg
  ! store_g = Don't care
  ! buf1 = RHS for g
  ! buf2 = Don't care
  !---------------------
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part3. save the hg to RK-store (only in the case RK-substep == RK-step
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (RK_SUBSTEP.ne.RK_STEP) then
     zeta = dt*RK_3_COEF(4,RK_SUBSTEP+1)
     do i=1,ny
        store_g(i) = zeta*u1(i)
     enddo
  endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part3
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !---------------------
  ! u1 = Don't care
  ! store_g = dt*zeta*hg
  ! buf1 = RHS for g
  ! buf2 = Don't care
  !---------------------
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part4. Solve for g(new) and d2g (new)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  call labs_bs(D_buf,buf1)
  call labs_mat_mul(D0,buf1,u1)
  call labs_mat_mul(D2,buf1,buf2)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part4
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !---------------------
  ! u1 = g(new)
  ! store_g = dt*zeta*hg
  ! buf1 = RHS for g
  ! buf2 = d2g(new)
  !---------------------

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part5. Update the store_g for next step
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (RK_SUBSTEP.ne.RK_STEP) then
     alpha = dt*RK_3_COEF(1,RK_SUBSTEP+1)*nu
     dummy = 1.d0-alpha*k2
     do i=1,ny
        store_g(i) = store_g(i) + dummy*u1(i) + alpha*buf2(i)
     enddo
  else
     do i=1,ny
        store_g(i) = u1(i)
     enddo
  endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part5
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !---------------------
  ! u1 = g(new)
  ! buf1 = Don't care
  ! buf2 = Don't care
  ! store_g = dt*zeta*hg + dt*alpha*nu*(1-k^2)*g + dt*alpha*nu*d2g (when RK_substep .ne. RK_step)
  !         = g      (when RK_substep .eq. RK_step)
  !---------------------



  return
end subroutine solve_g
