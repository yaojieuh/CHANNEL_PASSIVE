!--------------------------------------------------------------------
! This subroutine perform time advance calculation in wavespace
! at zero wave mode(kx = 0, kz = 0)
!
!      Solve simplified N-S eq
!      Forcing with constant mass flux
!      perform Galilean transform if u_gal is not zero
! 
! input 
! u1 = Don't care
! u2 = uv              
! u3 = Don't care
! u4 = vw
! u5 = Don't care
! store1 = Non linear storage term for u from previous step 
! store2 = Non linear storage term for w from previous step
! buf1 = Don't care
! buf2 = Don't care
! D0_lu : Factorized banded matrix for D0
! D0 : B-Spline matrix for U = Bc operation
! D1 : B-Spline matrix for d/dy operation
! D2 : B-Spline matrix for d^2/dy^2 operation
! D_buf : Scratch space to solve linear equation for time step
! 
! output 
! u1 = u_mean 
! u2 = v_mean (zero)
! u3 = w_mean 
! u4 = Don't care
! u5 = Don't care
! store1 = Low storage terms of u for next step 
! store2 = Low storage terms of w for next step 
! buf1 = Don't care
! buf2 = Don't care
!
! Developer : M.K. Lee (mk@ices.utexas.edu) 
! Last updated : 2013.01.09
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------
subroutine wave00(u1,u2,u3,u4,u5,u7,u8,u9,store1,store2,store3,buf1,buf2, &
     D0_LU, D0, D1, D2, D_buf, nu, pr, bulk_u, bulk_w, bulk_t, ny, nx32, nz32)
  use runge_kutta
  use global_variable, only : u_gal, wall_shear_stress, wall_theta_stress
  use labs
  implicit none
  complex(c_double_complex), dimension(ny) :: u1,u2,u3,u4,u5,u7,u8,u9,store1,store2,store3,buf1,buf2
  real(c_double), dimension(2*LABS_BW+1,ny) :: D0_LU,D0,D1,D2,D_buf
  real(c_double) :: alpha, beta, gamma, zeta, phi, nu,pr,bulk_u,bulk_w,bulk_t
  complex(c_double_complex) :: dummy_integration
  integer(c_intptr_t) :: nx32, nz32
  integer :: i,j,ny


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part1 .Calculate the LHS matrix and get LU decomposed matrix of it
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  beta = dt*RK_3_COEF(2,RK_SUBSTEP)*nu
  ! To enforce boundary condition of first element
  do i=1,LABS_BW
     D_buf(i,1) = 0.d0
  enddo
  D_buf(LABS_BW+1,1) = D0(LABS_BW+1,1) - beta*D2(LABS_BW+1,1)
  do i=LABS_BW+2,2*LABS_BW+1
     D_buf(i,1) = 0.d0
  enddo
  do j = 2,ny-1
     do i=1,2*LABS_BW+1
        D_buf(i,j) = (D0(i,j) - beta*D2(i,j))
     enddo
  enddo

  ! To enforce boundary condition of last element
  do i=1,LABS_BW
     D_buf(i,ny) = 0.d0
  enddo
  D_buf(LABS_BW+1,ny) = D0(LABS_BW+1,ny) - beta*D2(LABS_BW+1,ny)
  do i=LABS_BW+2,2*LABS_BW+1
     D_buf(i,ny) = 0.d0
  enddo

  call LABS_lu(D_buf)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part1 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



  ! u1 = Don't care
  ! u2 = uv        
  ! u3 = Don't care
  ! u4 = vw 
  ! u5 = Don't care
  ! store1 = Non linear storage term for u from previous step 
  ! store2 = Non linear storage term for w from previous step
  ! buf1 = Don't care
  ! buf2 = Don't care




  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part2. Calculate the RHS to calculate new U
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (RK_SUBSTEP.eq.1) then
     ! u1 = Don't care
     ! u2 = uv        
     ! u3 = Don't care
     ! u4 = vw 
     ! u5 = Don't care
     ! store1 = U_old
     ! store2 = W_old
     ! buf1 = Don't care
     ! buf2 = Don't care

     u1(:) = store1(:)
     call labs_bs(D0_lu,u1)
     call labs_bs(D0_lu,u2)
     call labs_mat_mul(D2,u1,buf1)
     call labs_mat_mul(D1,u2,u3)
     alpha = dt*RK_3_COEF(1,RK_SUBSTEP)*nu
     gamma = dt*RK_3_COEF(3,RK_SUBSTEP)
     do i=2,ny
        buf1(i) = store1(i) + alpha*buf1(i) - gamma*u3(i)
     enddo

     ! u1 = Don't care
     ! u2 = Don't care    
     ! u3 = d/dy uv 
     ! u4 = vw 
     ! u5 = Don't care
     ! store1 = Don't care
     ! store2 = W_old
     ! buf1 = RHS to solve U new
     ! buf2 = Don't care

  else
     
     ! u1 = Don't care
     ! u2 = uv        
     ! u3 = Don't care
     ! u4 = vw 
     ! u5 = Don't care
     ! store1 = U_old + dt*alpha*nu* d2/dy2 U_old - zeta*dt* d/dy uv
     ! store2 = Non linear storage term for w from previous step
     ! buf1 = Don't care
     ! buf2 = Don't care
     
     call labs_bs(D0_lu,u2)
     call labs_mat_mul(D1,u2,u3)
     gamma = dt*RK_3_COEF(3,RK_SUBSTEP)
     
     do i=2,ny
        buf1(i) = store1(i) - gamma*u3(i)
     enddo

     ! u1 = Don't care
     ! u2 = Don't care       
     ! u3 = d/dy uv 
     ! u4 = vw 
     ! u5 = Don't care
     ! store1 = Don't care
     ! store2 = Non linear storage term for w from previous step
     ! buf1 = RHS to solve U new
     ! buf2 = Don't care

  endif

  buf1(1) = (0.d0,0.d0)
  do i=2,ny-1
     buf1(i) = dcmplx(real(buf1(i))+u_gal,1.d0)
  enddo
  buf1(ny) = (0.d0,0.d0)

  ! u1 = Don't care
  ! u2 = Don't care
  ! u3 = d/dy uv 
  ! u4 = vw 
  ! u5 = Don't Care
  ! store1 = Don't care
  ! store2 = Non linear storage term for w from previous step
  ! buf1(real) = RHS for new bspline coeff of u
  ! buf1(imag) = [0,1,1,1,1,......,1,1,0]
  ! buf2 = Don't care

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of part 2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=






  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part 3. Caculate new U
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  call LABS_bs(D_buf,buf1)
  call labs_integration(buf1,dummy_integration)

  ! Calculate the forcing coeffcient 
  call calculate_phi(phi,bulk_u,dummy_integration)
  ! calculate new u
  do i=1,ny
     buf1(i) = dcmplx(real(buf1(i))+phi*aimag(buf1(i)),0.d0)
  enddo
!!$  phi_dpdx(RK_SUBSTEP) = (RK_3_COEF(3,RK_SUBSTEP)+RK_3_COEF(4,RK_SUBSTEP))*phi
  phi_dpdx(RK_SUBSTEP) = phi
  if(RK_SUBSTEP.eq.RK_STEP) then
     wall_shear_stress = 0.5d0*nu*(LABS_D1_B(LABS_BW+1,1,1)*real(buf1(1)) + LABS_D1_B(LABS_BW+2,1,1)*real(buf1(2)) &
          - LABS_D1_B(LABS_BW,ny,1)*real(buf1(Ny-1)) - LABS_D1_B(LABS_BW+1,ny,1)*real(buf1(Ny)))
  endif
  
  ! u1 = Don't care
  ! u2 = Don't care
  ! u3 = d/dy uv 
  ! u4 = vw 
  ! u5 = Don't Care
  ! store1 = Don't care
  ! store2 = Non linear storage term for w from previous step
  ! buf1(real) = bspline coeff of u after forcing
  ! buf1(imag) = 0
  ! buf2 = Don't care

  call LABS_mat_mul(D0,buf1,u1)
  do i=1,ny
     u1(i) = u1(i) - u_gal
  enddo
  ! u1 = u (new with forcing)
  ! u2 = Don't care
  ! u3 = d/dy uv 
  ! u4 = vw 
  ! u5 = Don't Care
  ! store1 = Don't care
  ! store2 = Non linear storage term for w from previous step
  ! buf1(real) = bspline coeff of u after forcing
  ! buf1(imag) = 0
  ! buf2 = Don't care
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part 3
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part 4. Store for next step ----- U
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if(RK_SUBSTEP.ne.RK_STEP) then
     call LABS_mat_mul(D2,buf1,store1)
     alpha = dt*RK_3_COEF(1,RK_SUBSTEP+1)*nu
     zeta = dt*RK_3_COEF(4,RK_SUBSTEP+1)
     do i=1,ny
        store1(i) = alpha*store1(i) + u1(i) - zeta*u3(i)
     enddo
  else
     store1(:) = u1(:)
  endif

  ! u1 = u (done)
  ! u2 = Don't care
  ! u3 = Don't care
  ! u4 = vw 
  ! u5 = Don't Care
  ! store1 = Non linear storage term for U next step (done)
  ! store2 = Non linear storage term for w from previous step
  ! buf1 = Don't care
  ! buf2 = Don't care

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part 4
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=








  ! u1 = u (done)
  ! u2 = Don't care
  ! u3 = Don't care
  ! u4 = vw 
  ! u5 = Don't Care
  ! store1 = Non linear storage term for U next step (done)
  ! store2 = Non linear storage term for w from previous step
  ! buf1 = Don't care
  ! buf2 = Don't care






  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part5. Calculate the RHS to calculate new W
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (RK_SUBSTEP.eq.1) then
     ! u1 = U_new (done)
     ! u2 = Don't care
     ! u3 = Don't care
     ! u4 = vw 
     ! u5 = Don't care
     ! store1 = Non linear storage term for U next step (done)
     ! store2 = W_old
     ! buf1 = Don't care
     ! buf2 = Don't care

     u2(:) = store2(:)
     call labs_bs(D0_lu,u2)
     call labs_mat_mul(D2,u2,buf1)
     call labs_bs(D0_lu,u4)
     call labs_mat_mul(D1,u4,u5)
     alpha = dt*RK_3_COEF(1,RK_SUBSTEP)*nu
     gamma = dt*RK_3_COEF(3,RK_SUBSTEP)
     do i=2,ny
        buf1(i) = store2(i) + alpha*buf1(i) - gamma*u5(i)
     enddo

     ! u1 = U_new (done)
     ! u2 = Don't care
     ! u3 = Don't care
     ! u4 = Don't care
     ! u5 = d/dy vw 
     ! store1 = Non linear storage term for U next step (done)
     ! store2 = W_old
     ! buf1 = RHS to calculate new W
     ! buf2 = Don't care

  else
     
     ! u1 = U_new (done)
     ! u2 = Don't care       
     ! u3 = Don't care
     ! u4 = vw 
     ! u5 = Don't care
     ! store1 = Non linear storage term for U next step (done)
     ! store2 = W_old + dt*alpha*nu* d2/dy2 W_old - zeta*dt* d/dy vw
     ! buf1 = Don't care
     ! buf2 = Don't care
     
     call labs_bs(D0_lu,u4)
     call labs_mat_mul(D1,u4,u5)
     gamma = dt*RK_3_COEF(3,RK_SUBSTEP)
     
     do i=2,ny
        buf1(i) = store2(i) - gamma*u5(i)
     enddo

     ! u1 = U_new (done)
     ! u2 = Don't care       
     ! u3 = Don't care
     ! u4 = Don't care
     ! u5 = d/dy vw 
     ! store1 = Non linear storage term for U next step (done)
     ! store2 = Don't care
     ! buf1 = RHS to solve W new
     ! buf2 = Don't care

  endif

  buf1(1) = (0.d0,0.d0)
  do i=2,ny-1
     buf1(i) = dcmplx(real(buf1(i)),1.d0)
  enddo
  buf1(ny) = (0.d0,0.d0)

  ! u1 = U_new (done)
  ! u2 = Don't care
  ! u3 = Don't care
  ! u4 = Don't care
  ! u5 = d/dy vw 
  ! store1 = Non linear storage term for U next step (done)
  ! store2 = Don't care
  ! buf1(real) = RHS for new bspline coeff of u
  ! buf1(imag) = [0,1,1,1,1,......,1,1,0]
  ! buf2 = Don't care

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of part 5
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=






  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part 6. Caculate new W
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  call LABS_bs(D_buf,buf1)
  call labs_integration(buf1,dummy_integration)

  ! Calculate the forcing coeffcient 
  call calculate_phi(phi,bulk_w,dummy_integration)
  phi = 0d0
  ! calculate new u
  do i=1,ny
     buf1(i) = dcmplx(real(buf1(i))+phi*aimag(buf1(i)),0.d0)
  enddo
  
  ! u1 = U_new (done)
  ! u2 = Don't care
  ! u3 = Don't care
  ! u4 = Don't care
  ! u5 = d/dy vw 
  ! store1 = Non linear storage term for U next step (done)
  ! store2 = Don't care
  ! buf1(real) = bspline coeff of u after forcing
  ! buf1(imag) = 0
  ! buf2 = Don't care

  call LABS_mat_mul(D0,buf1,u3)

  ! u1 = u (new with forcing)
  ! u2 = Don't care
  ! u3 = w (new with forcing)
  ! u4 = Don't care
  ! u5 = d/dy 
  ! store1 = Non linear storage term for U next step (done)
  ! store2 = Don't care
  ! buf1(real) = bspline coeff of u after forcing
  ! buf1(imag) = 0
  ! buf2 = Don't care
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part 6
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part 7. Store for next step ----- W
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if(RK_SUBSTEP.ne.RK_STEP) then
     call LABS_mat_mul(D2,buf1,store2)
     alpha = dt*RK_3_COEF(1,RK_SUBSTEP+1)*nu
     zeta = dt*RK_3_COEF(4,RK_SUBSTEP+1)
     do i=1,ny
        store2(i) = alpha*store2(i) + u3(i) - zeta*u5(i)
     enddo
  else
     store2(:) = u3(:)
  endif

  ! u1 = u (done)
  ! u2 = Don't care
  ! u3 = W (done)
  ! u4 = Don't care
  ! u5 = Don't Care
  ! store1 = Non linear storage term for U next step (done)
  ! store2 = Non linear storage term for W next step (done)
  ! buf1 = Don't care
  ! buf2 = Don't care

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part 7
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part 8. set V=0
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  u2(:) = (0d0,0d0)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of part 8
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Passive scalar
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part1 .Calculate the LHS matrix and get LU decomposed matrix of it
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  beta = dt*RK_3_COEF(2,RK_SUBSTEP)*nu/Pr
  ! To enforce boundary condition of first element
  do i=1,LABS_BW
     D_buf(i,1) = 0.d0
  enddo
  D_buf(LABS_BW+1,1) = D0(LABS_BW+1,1) - beta*D2(LABS_BW+1,1)
  do i=LABS_BW+2,2*LABS_BW+1
     D_buf(i,1) = 0.d0
  enddo
  do j = 2,ny-1
     do i=1,2*LABS_BW+1
        D_buf(i,j) = (D0(i,j) - beta*D2(i,j))
     enddo
  enddo

  ! To enforce boundary condition of last element
  do i=1,LABS_BW
     D_buf(i,ny) = 0.d0
  enddo
  D_buf(LABS_BW+1,ny) = D0(LABS_BW+1,ny) - beta*D2(LABS_BW+1,ny)
  do i=LABS_BW+2,2*LABS_BW+1
     D_buf(i,ny) = 0.d0
  enddo

  call LABS_lu(D_buf)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part1 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part2. Calculate the RHS to calculate new Passive
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (RK_SUBSTEP.eq.1) then
     ! u1 = U_new (done)
     ! u2 = Don't care
     ! u3 = Don't care
     ! u4 = vw 
     ! u5 = Don't care
     ! u8 = vt
     ! store1 = Non linear storage term for U next step (done)
     ! store2 = W_old
     ! store3 = T_old
     ! buf1 = Don't care
     ! buf2 = Don't care

     u7(:) = store3(:)
     call labs_bs(D0_lu,u7)
     call labs_mat_mul(D2,u7,buf1)
     call labs_bs(D0_lu,u8)
     call labs_mat_mul(D1,u8,u9)
     alpha = dt*RK_3_COEF(1,RK_SUBSTEP)*nu/Pr
     gamma = dt*RK_3_COEF(3,RK_SUBSTEP)
     do i=2,ny
        buf1(i) = store3(i) + alpha*buf1(i) - gamma*u9(i)
     enddo

     ! u1 = U_new (done)
     ! u2 = Don't care
     ! u3 = Don't care
     ! u4 = Don't care
     ! u5 = d/dy vw 
     ! u8 = vt
     ! u9 = d/dy vt
     ! store1 = Non linear storage term for U next step (done)
     ! store2 = W_old
     ! store3 = T_old
     ! buf1 = RHS to calculate new passive
     ! buf2 = Don't care

  else
     
     ! u1 = U_new (done)
     ! u2 = Don't care       
     ! u3 = Don't care
     ! u4 = vw 
     ! u5 = Don't care
     ! u8 = vt
     ! store1 = Non linear storage term for U next step (done)
     ! store2 = W_old + dt*alpha*nu* d2/dy2 W_old - zeta*dt* d/dy vw
     ! store3 = P_old + dt*alpha*nu/pr* d2/dy2 T_old - zeta*dt* d/dy vt
     ! buf1 = Don't care
     ! buf2 = Don't care
     
     call labs_bs(D0_lu,u8)
     call labs_mat_mul(D1,u8,u9)
     gamma = dt*RK_3_COEF(3,RK_SUBSTEP)
     
     do i=2,ny
        buf1(i) = store3(i) - gamma*u9(i)
     enddo

     ! u1 = U_new (done)
     ! u2 = Don't care       
     ! u3 = Don't care
     ! u4 = Don't care
     ! u5 = d/dy vw 
     ! u8 = vt
     ! u9 = d/dy vt
     ! store1 = Non linear storage term for U next step (done)
     ! store2 = Don't care
     ! buf1 = RHS to solve W new
     ! buf2 = Don't care

  endif

  buf1(1) = (0.d0,0.d0)
  do i=2,ny-1
     buf1(i) = dcmplx(real(buf1(i)),1.d0)
  enddo
  buf1(ny) = (0.d0,0.d0)

  ! u1 = U_new (done)
  ! u2 = Don't care
  ! u3 = Don't care
  ! u4 = Don't care
  ! u5 = d/dy vw 
  ! store1 = Non linear storage term for U next step (done)
  ! store2 = Don't care
  ! buf1(real) = RHS for new bspline coeff of u
  ! buf1(imag) = [0,1,1,1,1,......,1,1,0]
  ! buf2 = Don't care

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of part 2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=






  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part 3. Caculate new passive
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  call LABS_bs(D_buf,buf1)
  call labs_integration(buf1,dummy_integration)

  ! Calculate the forcing coeffcient 
  call calculate_phi(phi,bulk_t,dummy_integration)
  !phi = 0d0
  ! calculate new u
  do i=1,ny
     buf1(i) = dcmplx(real(buf1(i))+phi*aimag(buf1(i)),0.d0)
  enddo

  phi_dtdx(RK_SUBSTEP) = phi
  if(RK_SUBSTEP.eq.RK_STEP) then
     wall_theta_stress = 0.5d0*nu/pr*(LABS_D1_B(LABS_BW+1,1,1)*real(buf1(1)) + LABS_D1_B(LABS_BW+2,1,1)*real(buf1(2)) &
          - LABS_D1_B(LABS_BW,ny,1)*real(buf1(Ny-1)) - LABS_D1_B(LABS_BW+1,ny,1)*real(buf1(Ny)))
  endif
  
  ! u1 = U_new (done)
  ! u2 = Don't care
  ! u3 = Don't care
  ! u4 = Don't care
  ! u5 = d/dy vw 
  ! store1 = Non linear storage term for U next step (done)
  ! store2 = Don't care
  ! buf1(real) = bspline coeff of u after forcing
  ! buf1(imag) = 0
  ! buf2 = Don't care

  call LABS_mat_mul(D0,buf1,u7)

  ! u1 = u (new with forcing)
  ! u2 = Don't care
  ! u3 = w (new with forcing)
  ! u4 = Don't care
  ! u5 = d/dy 
  ! u7 = T
  ! store1 = Non linear storage term for U next step (done)
  ! store2 = Don't care
  ! buf1(real) = bspline coeff of u after forcing
  ! buf1(imag) = 0
  ! buf2 = Don't care
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part 3
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Part 4. Store for next step ----- T
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if(RK_SUBSTEP.ne.RK_STEP) then
     call LABS_mat_mul(D2,buf1,store3)
     alpha = dt*RK_3_COEF(1,RK_SUBSTEP+1)*nu/Pr
     zeta = dt*RK_3_COEF(4,RK_SUBSTEP+1)
     do i=1,ny
        store3(i) = alpha*store3(i) + u7(i) - zeta*u9(i)
     enddo
  else
     store3(:) = u7(:)
  endif

  ! u1 = u (done)
  ! u2 = Don't care
  ! u3 = W (done)
  ! u4 = Don't care
  ! u5 = Don't Care
  ! store1 = Non linear storage term for U next step (done)
  ! store2 = Non linear storage term for W next step (done)
  ! buf1 = Don't care
  ! buf2 = Don't care

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! End of Part 4
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

 
  return
end subroutine wave00








!--------------------------------------------------------------------
! This subroutine calculate the coefficient 
! to ensure the constant mass flux
! 
!
! Developer : M.K. Lee
! Last updated : 2013.01.04
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------
subroutine calculate_phi(phi,bulk_vel,dummy_integration)
  use iso_c_binding
  implicit none
  real(c_double) :: bulk_vel,phi
  complex(c_double_complex) :: dummy_integration
  
  phi = (bulk_vel - 0.5d0*real(dummy_integration))/(0.5d0*aimag(dummy_integration))
!!$  phi=0.d0

  return
end subroutine calculate_phi
