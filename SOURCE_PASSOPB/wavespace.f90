!--------------------------------------------------------------------
! This subroutine perform time advance calculation in wavespace
!
! When (kx = 0, kz = 0)
!      Solve simplified N-S eq
!      Forcing with constant mass flux
!      perform Galilean transform if u_gal is not zero
! 
! when (kx != 0 or kz != 0)
! 
!      Solve vorticity equation in y direction given wave numbers
!              g = W_y
!      Solve phi equation in y direction given wave numbers
!              phi = Del^2 v
!      Solve for u and w with g and dv/dy
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------
subroutine wavespace(U,RK_STORE,buf1,buf2,D0_LU,D0,D1,D2,D_buf,kx,kz,labs_ny,ny,omp_y,zsz,xsz,bw,nx32,nz32)
  use, intrinsic :: iso_c_binding
  use runge_kutta
  use global_variable
  implicit none
  integer(C_INTPTR_T) :: ny,zsz,xsz, nx32, nz32, omp_y
  integer :: kx(xsz), kz(zsz)
  complex(C_DOUBLE_COMPLEX) :: U(ny,xsz,zsz,10)
  complex(C_DOUBLE_COMPLEX) :: RK_STORE(ny,xsz,zsz,2)
  complex(C_DOUBLE_COMPLEX) :: buf1(ny,omp_y)
  complex(C_DOUBLE_COMPLEX) :: buf2(ny,omp_y)
  integer :: bw
  real(C_DOUBLE),dimension(2*bw+1,labs_ny,omp_compute) :: D0_LU,D0,D1,D2,D_buf
  integer :: n_total,xsz_dummy,labs_ny
  integer :: i,j,k,timer1, timer2,l
  ! This is needed to ensure compatibility of ISO_C_BINDING with some compilers (e.g. XL compiler)
  n_total = zsz*xsz
  xsz_dummy = xsz

  call omp_set_num_threads(omp_compute)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l) 
  l = omp_get_thread_num() + 1
  !$OMP DO SCHEDULE(STATIC)
  do i=1,n_total

     k=mod(int(i-1),xsz_dummy) + 1
     j = (i-1)/xsz_dummy + 1
     
     if((kx(k).eq.0) .and. (kz(j).eq.0)) then
        ! Different operation for [0,0] wave mode
        ! In this case v = 0 , Vorticity_y = 0
        call wave00(U(:,k,j,1),U(:,k,j,2),U(:,k,j,3),U(:,k,j,4),U(:,k,j,5), &
	     U(:,k,j,7),U(:,k,j,8),U(:,k,j,9),RK_store(:,k,j,1),&
             RK_store(:,k,j,2),RK_store(:,k,j,3),buf1(:,l),buf2(:,l), &
             D0_LU(:,:,l), D0(:,:,l), D1(:,:,l), D2(:,:,l), D_buf(:,:,l), &
             nu,pr, bulk_vel_u,bulk_vel_w,bulk_t, labs_ny, nx32, nz32)
     else

        ! u1 = u^2 - v^2
        ! u2 = uv
        ! u3 = w^2 - v^2
        ! u4 = vw
        ! u5 = uw

        call hvhg(U(:,k,j,1),U(:,k,j,2),U(:,k,j,3),U(:,k,j,4),U(:,k,j,5), &
		U(:,k,j,7),U(:,k,j,8),U(:,k,j,9),&
             D0_LU(:,:,l), D1(:,:,l), D2(:,:,l), kx(k)*LMFxr, kz(j)*LMFzr,nx32,nz32)

        ! u1 = hg
        ! u2 = Don't care
        ! u3 = Don't care
        ! u4 = Don't care
        ! u5 = hv

        call solve_g(U(:,k,j,1),RK_store(:,k,j,1),buf1(:,l),buf2(:,l),D0_LU(:,:,l),D0(:,:,l),D2(:,:,l),D_buf(:,:,l), &
             kx(k)*LMFxr, kz(j)*LMFzr,nu,labs_ny)

        ! u1 = g(new)
        ! u2 = Don't care
        ! u3 = Don't care
        ! u4 = Don't care
        ! u5 = hv

        call solve_v_f(U(:,k,j,2),U(:,k,j,3),U(:,k,j,5),RK_store(:,k,j,2),U(:,k,j,4),buf1(:,l),buf2(:,l), &
             U(:,k,j,10),D0_LU(:,:,l),D0(:,:,l),D1(:,:,l),D2(:,:,l),D_buf(:,:,l),kx(k)*LMFxr,kz(j)*LMFzr,nu,labs_ny)

        ! u1 = g(new)
        ! u2 = v(new)
        ! u3 = dvdy(new)
        ! u4 = Don't care
        ! u5 = Don't care

        call solve_u_w(U(:,k,j,1),U(:,k,j,3),kx(k)*LMFxr,kz(j)*LMFzr)

        call solve_t(U(:,k,j,7),RK_store(:,k,j,3),buf1(:,l),buf2(:,l),D0_LU(:,:,l),D0(:,:,l),D2(:,:,l),D_buf(:,:,l), &
             kx(k)*LMFxr, kz(j)*LMFzr,nu,pr,labs_ny)

        ! u1 = u(new)
        ! u2 = v(new)
        ! u3 = w(new)
        ! u4 = Don't care
        ! u5 = Don't care
        ! u7 = theta
     endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  return

end subroutine wavespace
