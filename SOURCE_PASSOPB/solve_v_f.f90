!-----------------------------------------------------------
! This routine calculate the v and dv/dy
!
! input 
! v = Don't care
! dv = Don't care
! hv = hv_new
! store_phi = Non linear storage term for phi from previous step (except for rk_step=0)
! buf1 = Don't care
! buf2 = Don't care
! buf3 = Don't care
! D0_lu : Factorized banded matrix for D0
! D0 : B-Spline matrix 
! D1 : B-Spline matrix for d/dy operation
! D2 : B-Spline matrix for d^2/dy^2 operation
! D2_lu : Factorized banded matrix for D2
! D_buf : Scratch space to solve linear equation for time step (precalculated)
! kx, kz : wave numbers
! nu : viscosity
! ny : number of grids in y-direction
!
! output 
! v = v_new 
! dv= f (=dv/dy)
! hv = Don't care
! u5 = Don't care
! store_phi= Low storage terms of phi for next step 
! buf1 = Don't care
! buf2 = Don't care
! buf3 = Don't care
!
! Developer : M.K. Lee (mk@ices.utexas.edu) 
! Last update : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!-----------------------------------------------------------

subroutine solve_v_f(v,dv,hv,store_phi,buf1,buf2,buf3, &
     vold,D0_LU, D0, D1, D2, D_buf, kx, kz, nu, ny)

  use labs
  use runge_kutta
  implicit none

  complex(c_double_complex), dimension(ny) :: v,dv,hv,store_phi,buf1,buf2,buf3,vold
  real(c_double), dimension(2*LABS_BW+1,ny) :: D0_LU,D0,D1,D2,D_buf
  real(c_double) :: kx,kz,k2, alpha, beta, gamma, zeta, dummy, nu
  complex(c_double_complex) :: c1,c2
  integer :: i,j,ny

  k2 = kx**2 + kz**2

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part1. Calculate the RHS to calculate new Phi
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (RK_SUBSTEP.eq.1) then
     !---------------------
     ! v = Don't care
     ! dv = Don't care
     ! hv = hv_new
     ! store_phi = phi_old
     ! buf1 = Don't care
     ! buf2 = Don't care
     ! buf3 = Don't care
     !---------------------
     alpha = dt*RK_3_COEF(1,RK_SUBSTEP)*nu
     gamma = dt*RK_3_COEF(3,RK_SUBSTEP)
     dummy = 1.d0-alpha*k2
     do i=1,ny
        buf1(i) = dummy*store_phi(i) + gamma*hv(i)
     enddo
     call LABS_bs(D0_lu,store_phi)
     call LABS_mat_mul(D2,store_phi,buf2)
     do i=1,ny
        buf1(i) = buf1(i) + alpha*buf2(i)
     enddo
  else
     !---------------------
     ! v = Don't care
     ! dv = Don't care
     ! hv = hv_new
     ! store_phi = Non linear storage term for phi from previous step (except for rk_step=0)
     ! buf1 = Don't care
     ! buf2 = Don't care
     ! buf3 = Don't care
     !---------------------
     gamma = dt*RK_3_COEF(3,RK_SUBSTEP)
     do i=1,ny
        buf1(i) = store_phi(i) + gamma*hv(i)
     enddo
  endif

  ! boundary condition of phi(y=+- 1) = 0
  buf1(1) = dcmplx(0.d0,0.d0)
  buf1(ny) = dcmplx(0.d0,0.d0)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part1
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  
  !---------------------
  ! v = Don't care
  ! dv = Don't care
  ! hv = hv_new
  ! store_phi = Don't care
  ! buf1 = RHS to calculate new phi
  ! buf2 = Don't care
  ! buf3 = Don't care
  !---------------------


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part2. save the hv to RK-store (only in the case RK-substep == RK-step
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (RK_SUBSTEP.ne.RK_STEP) then
     zeta = dt*RK_3_COEF(4,RK_SUBSTEP+1)
     do i=1,ny
        store_phi(i) = zeta*hv(i)
     enddo
  endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !---------------------
  ! v = Don't care
  ! dv = Don't care
  ! hv = Don't care
  ! buf1 = RHS to calculate new phi
  ! buf2 = Don't care
  ! buf3 = Don't care
  ! store_phi = (RK-0) 0 
  ! store_phi = (RK-1,2) dt*zeta*hv_new 
  !---------------------



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part3. Calculate Coeff of Phi_particular / homogeneous 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  buf2(1) = dcmplx(D_buf(LABS_BW+1,1),D_buf(LABS_BW+1,1))
  do i=2,ny
     buf2(i) = dcmplx(0.d0,0.d0)
  enddo

  call labs_bs(D_buf,buf1)
  call labs_bs(D_buf,buf2)

  call labs_mat_mul(D0,buf1,hv)
  call labs_mat_mul(D0,buf2,buf3)
  call labs_mat_mul(D2,buf1,dv)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part 3
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !---------------------
  ! v = Don't care
  ! dv = d2 phi_particular
  ! hv = phi_particular
  ! buf1 = C_phi_particular
  ! buf2 = C_phi_homogeneous_1
  ! buf3 = phi_homogeneous_1
  ! store_phi = (RK-0) 0 
  ! store_phi = (RK-1,2) dt*zeta*hv_new 
  !---------------------


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part4. Update store_phi with phi_particular
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if(RK_SUBSTEP.ne.RK_STEP) then
     alpha = dt*nu*RK_3_COEF(1,RK_SUBSTEP+1)
     dummy = 1d0 - alpha*k2
     do i=1,ny
        store_phi(i) = store_phi(i) + dummy*hv(i) + alpha*dv(i)
     enddo
  else
     store_phi(:) = hv(:)
  endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part 4
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



  !---------------------
  ! v = Don't care
  ! dv = Don't care
  ! hv = phi_particular
  ! buf1 = Don't care
  ! buf2 = C_phi_homogeneous_1
  ! buf3 = phi_homogeneous_1
  ! store_phi = (RK-0) phi_particular
  ! store_phi = (RK-1,2) dt*zeta*hv_new + phi_particular + dt*alpha*nu*(Del)phi_particular 
  !---------------------



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part5. Make a Laplacian Operator / LU-Decompose
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! To enforce boundary condition of first element
  do i=1,LABS_BW
     D_buf(i,1) = 0d0
  enddo
  !D_buf(LABS_BW+1,1) = -k2*D0(LABS_BW+1,1) + D2(LABS_BW+1,1)
  D_buf(LABS_BW+1,1) = D0(LABS_BW+1,1) 
  do i=LABS_BW+2,2*LABS_BW+1
     D_buf(i,1) = 0d0
  enddo
  
  do j = 2,ny-1
     do i=1,2*LABS_BW+1
        D_buf(i,j) = -k2*D0(i,j) + D2(i,j)
     enddo
  enddo

  ! To enforce boundary condition of last element
  do i=1,LABS_BW
     D_buf(i,ny) = 0d0
  enddo
  !D_buf(LABS_BW+1,ny) = -k2*D0(LABS_BW+1,ny) + D2(LABS_BW+1,ny)
  D_buf(LABS_BW+1,ny) = D0(LABS_BW+1,ny)
  do i=LABS_BW+2,2*LABS_BW+1
     D_buf(i,ny) = 0d0
  enddo

  call labs_lu(D_buf)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part 5
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !---------------------
  ! v = Don't care
  ! dv = Don't care
  ! hv = phi_particular
  ! buf1 = Don't care
  ! buf2 = C_phi_homogeneous_1
  ! buf3 = phi_homogeneous_1
  ! store_phi = (RK-0) phi_particular
  ! store_phi = (RK-1,2) dt*zeta*hv_new + phi_particular + dt*alpha*nu*(Del)phi_particular 
  !---------------------


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part6. Calculate dyV and calculate the coefficent
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  hv(1) = (0d0,0d0)
  hv(ny) = (0d0,0d0)
  !hv(1)=-vold(23)
  !hv(ny)=-vold(ny-22)
  !hv(1)=-0.1d0*vold(22)+0.9*vold(1)
  !hv(ny)=-0.1d0*vold(ny-21)+0.9*vold(ny)
  !write(*,*) 'vld',vold(15),'hv', hv(1)
  call labs_bs(D_buf,hv)
  buf1(:) = buf3(:)
  buf1(1) = (0d0,0d0)
  buf1(ny) = (0d0,0d0)
  call labs_bs(D_buf,buf1)
  !hv(1) = (0d0,0d0)
  !hv(ny) = (0d0,0d0)
  !do i=1,2*LABS_BW+1
  !      hv(1)= hv(1)+ D0(i,15)*hv(11+i)
       ! hv(ny)=hv(ny)+D0(i,ny-14)*hv(ny-17+i)
  !enddo 
  call labs_mat_mul(D0,hv,v)
  call labs_mat_mul(D1,hv,dv)
  call labs_mat_mul(D1,buf1,hv)
  !---------------------
  ! v = v_particular
  ! dv = dv_particular
  ! hv = dv_homogeneous_1
  ! buf1 = C_v_homogeneous_1
  ! buf2 = C_phi_homogeneous_1
  ! buf3 = phi_homogeneous_1
  ! store_phi = (RK-0) phi_particular
  ! store_phi = (RK-1,2) dt*zeta*hv_new + phi_particular + dt*alpha*nu*(Del)phi_particular 
  !---------------------

  ! calculate C1, C2
  ! Vp'(1) + C1*Vh1'(1) + C2*vh2'(1) = 0
  ! Vp'(n) + C1*Vh1'(n) + C2*vh2'(n) = 0
  ! Using the antisymmetry property of homogeneous v'
  ! Vp'(1) + C1*vh1'(1) - C2*vh1'(n) = 0
  ! Vp'(n) + C1*vh1'(n) - C2*vh1'(1) = 0
  ! Hence 
  ! C1 = [Vh1'(n)*Vp'(n) - Vh1'(1)*Vp'(1)] /  [Vh1'(1)^2 - Vh1'(n)^2] 
  ! C2 = [Vh1'(1)*Vp'(n) - Vh1'(n)*Vp'(1)] /  [Vh1'(1)^2 - Vh1'(n)^2]

  C1 = (hv(ny)*dv(ny) - hv(1)*dv(1)) /  (hv(1)**2 - hv(ny)**2) 
  C2 = (hv(1)*dv(ny) - hv(ny)*dv(1)) /  (hv(1)**2 - hv(ny)**2)

  do i=1,ny
     dv(i) = dv(i) + C1*hv(i) - C2*hv(ny+1-i)
  enddo
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part 6
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !---------------------
  ! v = v_particular
  ! dv = dv  (DONE)
  ! hv = Don't care
  ! buf1 = C_v_homogeneous_1
  ! buf2 = C_phi_homogeneous_1
  ! buf3 = phi_homogeneous_1
  ! store_phi = (RK-0) phi_particular
  ! store_phi = (RK-1,2) dt*zeta*hv_new + phi_particular + dt*alpha*nu*(Del)phi_particular 
  !---------------------


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part7. Update store_phi 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  call labs_mat_mul(D2,buf2,hv)

  if(RK_SUBSTEP.ne.RK_STEP) then
     alpha = dt*nu*RK_3_COEF(1,RK_SUBSTEP+1)
     dummy = 1d0 - alpha*k2
     do i=1,ny
        store_phi(i) = store_phi(i) + C1*dummy*buf3(i) + C2*dummy*buf3(ny+1-i) + C1*alpha*hv(i) + C2*alpha*hv(ny+1-i)
     enddo
  else
     do i=1,ny
        store_phi(i) = store_phi(i) + C1*buf3(i) + C2*buf3(ny+1-i)
     enddo
  endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of part 7
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !---------------------
  ! v = v_particular
  ! dv = dv  (DONE)
  ! hv = Don't care
  ! buf1 = C_v_homogeneous_1
  ! buf2 = Don't care
  ! buf3 = Don't care
  ! store_phi = (RK-0) phi (DONE)
  ! store_phi = (RK-1,2) dt*zeta*hv_new + phi + dt*alpha*nu*(Del)phi (DONE) 
  !---------------------


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part8. Calculate v with homogeneous term
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  call labs_mat_mul(D0,buf1,hv)
  do i=1,ny
     v(i) = v(i) + C1*hv(i) + C2*hv(ny+1-i)
  enddo
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part8
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  
  !---------------------
  ! v = v (Done)
  ! dv = dv  (DONE)
  ! hv = Don't care
  ! buf1 = Don't care
  ! buf2 = Don't care
  ! buf3 = Don't care
  ! store_phi = (RK-0) phi (DONE)
  ! store_phi = (RK-1,2) dt*zeta*hv_new + phi + dt*alpha*nu*(Del)phi (DONE) 
  !---------------------
  

  

  return
end subroutine solve_v_f
