subroutine get_pressure

  use global_variable
  use ppdfftw
  use labs
  implicit none

  call cal_pressure(kx,kz,labs_ny, ny, omp_y, Yzsz, Yxsz, LABS_BW, U,stat_PPE_boundary,stat_vv, &
       buf_line_y_c1,buf_line_y_c2,LABS_D0_LU,LABS_D0_B,LABS_D1_B,LABS_D2_B,LABS_BUF)       

  return
end subroutine get_pressure




subroutine cal_pressure(kx,kz,labs_ny,ny,omp_y,zsz,xsz,bw,U,BC,vv,buf1,buf2,D0_LU,D0,D1,D2,D_buf)
  use, intrinsic :: iso_C_binding
  use global_variable
  implicit none
  integer(C_INTPTR_T) :: ny,zsz,xsz,omp_y
  integer :: kx(xsz), kz(zsz)
  integer :: n_total,xsz_dummy,labs_ny
  complex(C_DOUBLE_COMPLEX) :: U(ny,xsz,zsz,6)
  complex(C_DOUBLE_COMPLEX) :: BC(2,xsz,zsz)
  real(C_DOUBLE) :: vv(ny)
  complex(C_DOUBLE_COMPLEX) :: buf1(ny,omp_y)
  complex(C_DOUBLE_COMPLEX) :: buf2(ny,omp_y)
  integer :: bw
  real(C_DOUBLE),dimension(2*bw+1,labs_ny,omp_compute) :: D0_LU,D0,D1,D2,D_buf
  integer :: i,j,k,l,m

  xsz_dummy = xsz

  n_total = xsz*zsz
  
  call omp_set_num_threads(omp_compute)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l,m) 
  l = omp_get_thread_num() + 1
  !$OMP DO SCHEDULE(STATIC)
  do i=1,n_total

     k=mod(int(i-1),xsz_dummy) + 1
     j = (i-1)/xsz_dummy + 1

     if((kx(k).eq.0) .and. (kz(j).eq.0)) then
        do m=1,ny
           U(m,k,j,6) = dcmplx(-vv(m),0d0)
        enddo
     else

        call solve_ppe(U(:,k,j,1),U(:,k,j,2),U(:,k,j,3),U(:,k,j,4),U(:,k,j,5),U(:,k,j,6), &
             BC(1,k,j), BC(2,k,j), buf1(:,l), buf2(:,l), D0_LU(:,:,l), D0(:,:,l),D1(:,:,l),D2(:,:,l),D_buf(:,:,l), &
             kx(k)*LMFxr, kz(j)*LMFzr,labs_ny)
     endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  return
end subroutine cal_pressure



subroutine solve_ppe(u1,u2,u3,u4,u5,u6,BC1,BC2,buf1,buf2,D0_lu,D0,D1,D2,D_buf,kx,kz,ny)
  use labs
  use, intrinsic :: iso_c_binding

  implicit none

  complex(c_double_complex), dimension(ny) :: u1,u2,u3,u4,u5,u6,buf1,buf2
  real(c_double), dimension(2*LABS_BW+1,ny) :: D0_LU,D0,D1,D2,D_buf
  real(c_double) :: kx,kz,k2,kxkz
  integer :: i,j,ny
  complex(c_double_complex) :: ikx, ikz, BC1,BC2

  ikx = dcmplx(0d0,kx)
  ikz = dcmplx(0d0,kz)

  k2 = kx*kx + kz*kz
  kxkz = kx*kz
  ! u1 = u^2 - v^2
  ! u2 = uv
  ! u3 = w^2 - v^2
  ! u4 = vw
  ! u5 = uw
  ! u6 = v^2
  
  do i=2,ny-1
     buf1(i) = kx*kx*u1(i) + kz*kz*u3(i) + k2*u6(i) + 2d0*kxkz*u5(i)
  enddo
  
  call labs_bs(d0_lu(:,:),u6)
  call labs_mat_mul(d2(:,:),u6,buf2)
  do i=2,ny-1
     buf1(i) = buf1(i) - buf2(i)
     u6(i) = ikx*u2(i) + ikz*u4(i)
  enddo


  call labs_bs(d0_lu(:,:),u6)
  call labs_mat_mul(d1(:,:),u6,buf2)

  buf1(1) = BC1
  do i=2,ny-1
     buf1(i) = buf1(i) - 2d0*buf2(i)
  enddo
  buf1(ny) = BC2


  do i=1,2*LABS_BW+1
     D_buf(i,1) = D1(i,1)
  enddo
  do j=2,ny-1
     do i=1,2*LABS_BW+1
        D_buf(i,j) = -k2*D0(i,j) + D2(i,j)
     enddo
  enddo
  do i=1,2*LABS_BW+1
     D_buf(i,ny) = D1(i,ny)
  enddo


  call labs_lu(D_buf)
  call labs_bs(D_buf,buf1)
  call labs_mat_mul(D0,buf1,u6)

  ! u1 = u^2 - v^2
  ! u2 = uv
  ! u3 = w^2 - v^2
  ! u4 = vw
  ! u5 = uw
  ! u6 = p


  return
end subroutine solve_ppe

