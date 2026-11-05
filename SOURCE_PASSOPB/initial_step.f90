!--------------------------------------------------------------------
! This subroutine feed initial field
!
! - Fill U/V/W 
!     With restart file : read data file from 'restart_vel.h5'
!     Without restart file : fill non zero numbers in low wave number components
! - Fill RK_store field which is needed to solve time advance equation
! 
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------
subroutine initial_step
  use global_variable
  use ppdfftw
  use runge_kutta
  use labs
  implicit none

  if (restart_file_exist.eq.0) then
     call set_initial_value_in_y_pencil(U,u_gal,ny,Yxsz,Yxst,nz,Yzsz,Yzst,labs_y)
     call set_initial_data_uvw(U,RK_store,kx,kz,LMFxr,LMFzr,ny,Yxsz,Yxst,Yzsz,Yzst)
     timestep = 0
     time = 0d0
     dt = 0d0
     target_stat_time = stat_time_interval
     target_data_time = data_time_interval
  else 
     if (restart_option.eq.1) then
        call data_read_uvw(U)
        call set_initial_data_uvw(U,RK_store,kx,kz,LMFxr,LMFzr,ny,Yxsz,Yxst,Yzsz,Yzst)
     else if (restart_option.eq.2) then
        call data_read_gphi(U,RK_store)
        call set_initial_data_gphi(U,RK_store,kx,kz,LMFxr,LMFzr,ny,Yxsz,Yxst,Yzsz,Yzst,&
             LABS_D0_B(1,1,1),LABS_D1_B(1,1,1),LABS_D2_B(1,1,1))
     else
        call data_read_vg(U,RK_store)
        call set_initial_data_vg(U,RK_store,kx,kz,LMFxr,LMFzr,ny,Yxsz,Yxst,Yzsz,Yzst)
     endif
  endif
!!$  timestep = 0
!!$  time = 0d0
!  timestep = 0
!  time = 0d0
!!$  dt = 0d0
  timestep_code_start = timestep
  target_stat_time = (floor(time/stat_time_interval)+1)*stat_time_interval
  target_data_time = (floor(time/data_time_interval)+1)*data_time_interval
  if (time.le.3d-16) then
     stat_store_or_not = .true.
     target_stat_time = 0d0
  endif
  return
end subroutine initial_step






subroutine set_initial_value_in_y_pencil(U,u_gal,ny,Yxsz,Yxst,nz,Yzsz,Yzst,y)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_intptr_t) :: ny,Yxsz,Yxst,nz,Yzsz,Yzst
  complex(c_double_complex) :: U(ny,Yxsz,Yzsz,7)
  integer :: i,j,k,l,ny_dummy,xsz,size
  real(c_double) :: y(ny)
  real(c_double) :: u_gal
  ny_dummy = ny
  xsz = Yxsz
  size = Yxsz*Yzsz


  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l)
  !$OMP DO SCHEDULE(STATIC)
  do l=1,size
     k = mod(l-1,xsz) + 1
     j = (l-1)/xsz + 1
     if (((j+Yzst).eq.1) .and.((k+Yxst).eq.1)) then
        do i=1,ny_dummy
           U(i,k,j,1) = dcmplx( (1.d0-y(i)**2)*1.5d0-u_gal, 0.d0 )
        enddo
        U(:,k,j,2) = (0d0,0d0)
        U(:,k,j,3) = (0d0,0d0)
	U(:,k,j,7) = U(:,k,j,1)
     else if (((j+Yzst).eq.1) .and.((k+Yxst).eq.2)) then
        do i=1,ny_dummy
           U(i,k,j,2) = dcmplx(((1.d0-y(i))*(1d0+y(i)))**2 *2d-2,((1.d0-y(i))*(1d0+y(i)))**2 *1d-2)
           U(i,k,j,1) = (0d0,0d0)
           U(i,k,j,3) = (0d0,0d0)
	   U(:,k,j,7) = U(:,k,j,1)
        enddo
     else if (((j+Yzst).eq.2) .and.((k+Yxst).eq.1)) then
        do i=1,ny_dummy
           U(i,k,j,1) = dcmplx((1.d0-y(i)**2)*2d-2,(1.d0-y(i)**2)*1d-2)
           U(i,k,j,2) = dcmplx(-(1.d0-y(i)**2)*1d-2,(1.d0-y(i)**2)*1d-2)
           U(i,k,j,3) = dcmplx(-(1.d0-y(i)**2)*1d-2,-(1.d0-y(i)**2)*2d-2)
	   U(:,k,j,7) = U(:,k,j,1)
        enddo
     else if (((j+Yzst).eq.nz) .and.((k+Yxst).eq.1)) then
        do i=1,ny_dummy
           U(i,k,j,1) = dcmplx((1.d0-y(i)**2)*2d-2,-(1.d0-y(i)**2)*1d-2)
           U(i,k,j,2) = dcmplx(-(1.d0-y(i)**2)*1d-2,-(1.d0-y(i)**2)*1d-2)
           U(i,k,j,3) = dcmplx(-(1.d0-y(i)**2)*1d-2,(1.d0-y(i)**2)*2d-2)
	   U(:,k,j,7) = U(:,k,j,1)
        enddo
     else
        U(:,k,j,1) = (0d0,0d0)
        U(:,k,j,2) = (0d0,0d0)
        U(:,k,j,3) = (0d0,0d0)
	U(:,k,j,7) = U(:,k,j,1)
     endif

  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL

  return
end subroutine set_initial_value_in_y_pencil







!--------------------------------------------------------------------
! This subroutine calculate RK_store field which is needed to solve time advnace equation
! 
! u/v/w need to be provided.
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------
subroutine set_initial_data_uvw(U,RK_store,kx,kz,LMFxr,LMFzr,ny,Yxsz,Yxst,Yzsz,Yzst)
  use labs
  implicit none
  integer(c_intptr_t) :: ny,Yxsz,Yxst,Yzsz,Yzst
  complex(c_double_complex) :: U(ny,Yxsz,Yzsz,7)
  complex(c_double_complex) :: RK_store(ny,Yxsz,Yzsz,3)
  complex(c_double_complex) :: dummy_c(ny)
  integer :: kx(Yxsz), kz(Yzsz)
  real(c_double) :: LMFxr,LMFzr
  integer :: i,j,k,l,xsz,size


  size = Yzsz*Yxsz
  xsz = Yxsz

!!$  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l) FIRSTPRIVATE(LABS_D0_LU,LABS_D2_B,dummy_c,kx,LMFxr,kz,LMFzr)
!!$  !$OMP DO SCHEDULE(STATIC)
  do l=1,size
     k = mod(l-1,xsz) + 1
     j = (l-1)/xsz + 1

     if ((kx(k).eq.0).and.(kz(j).eq.0)) then
        RK_store(:,k,j,1) = U(:,k,j,1)
        RK_store(:,k,j,2) = U(:,k,j,3)
	RK_store(:,k,j,3) = U(:,k,j,7)
     else
        dummy_c(:) = U(:,k,j,2)
        call labs_bs(LABS_D0_LU(:,:,1),dummy_c)
        call labs_mat_mul(LABS_D2_B(:,:,1),dummy_c,RK_store(:,k,j,2))
        do i=1,ny
           RK_store(i,k,j,2) = RK_store(i,k,j,2) - ( (kx(k)*LMFxr)**2 + (kz(j)*LMFzr)**2 )*U(i,k,j,2)
           RK_store(i,k,j,1) = dcmplx(0.d0,kz(j)*LMFzr)*U(i,k,j,1) - dcmplx(0.d0,kx(k)*LMFxr)*U(i,k,j,3)
        enddo
	RK_store(:,k,j,3) = U(:,k,j,7)
     endif
  enddo
!!$  !$OMP END DO
!!$  !$OMP END PARALLEL

  return
end subroutine set_initial_data_uvw




subroutine data_read_uvw(Udata)

  use esio
  use global_variable
  use labs
  use ppdfftw
  use runge_kutta

  implicit none

  real(C_DOUBLE) :: Udata(2,ny,Yxsz,Yzsz,no_of_field)
  real(C_DOUBLE) :: dummy_real(1)
  integer :: dummy_integer(1)
  logical :: read_write
  type(esio_handle) :: h

  call esio_handle_initialize(h,MPI_COMM_WORLD)
  call esio_file_open(h,restart_file_name,.true.,ierr)

  call esio_line_establish(h,1,1,1,ierr)

  call esio_line_read_integer(h,"timestep",dummy_integer,0)
  timestep = dummy_integer(1)
  target_timestep = target_timestep + timestep
  call esio_line_read_double(h,"time",dummy_real,0)
  time = dummy_real(1)
  call esio_line_read_double(h,"dt",dummy_real,0)
  dt = dummy_real(1)

  call esio_field_establish (h, INT(ny), 1, INT(ny), INT(nx)/2, INT(Yxst)+1, INT(Yxsz), INT(nz), INT(Yzst)+1, INT(Yzsz), ierr)
  call esio_field_readv_double (h, "u",  Udata(1,1,1,1,1),2)
  call esio_field_readv_double (h, "v",  Udata(1,1,1,1,2),2)
  call esio_field_readv_double (h, "w",  Udata(1,1,1,1,3),2)
  call esio_field_readv_double (h, "theta",  Udata(1,1,1,1,7),2)
  !call esio_field_readv_double (h, "u",  Udata(1,1,1,1,7),2)

  call esio_file_close(h)
  call esio_handle_finalize(h)

  

 ! close off esio and datafile set
  if(myrank .eq. 0) then
     Udata(1,:,1,1,1) = Udata(1,:,1,1,1) - U_gal
  endif

  return
end subroutine data_read_uvw





subroutine data_read_vg(Udata,RKdata)

  use esio
  use global_variable
  use labs
  use ppdfftw
  use runge_kutta

  implicit none

  real(C_DOUBLE) :: Udata(2,ny,Yxsz,Yzsz,no_of_field)
  real(C_DOUBLE) :: RKdata(2,ny,Yxsz,Yzsz,3)
  real(C_DOUBLE) :: dummy_real(1)
  real(C_DOUBLE) :: dummy_u(ny), dummy_w(ny), dummy_t(ny)
  integer :: dummy_integer(1)
  logical :: read_write
  type(esio_handle) :: h


  call esio_handle_initialize(h,MPI_COMM_WORLD)
  call esio_file_open(h,restart_file_name,.true.,ierr)

  call esio_line_establish(h,1,1,1,ierr)

  call esio_line_read_integer(h,"timestep",dummy_integer,0)
  timestep = dummy_integer(1)
  target_timestep = target_timestep + timestep
  call esio_line_read_double(h,"time",dummy_real,0)
  time = dummy_real(1)
  call esio_line_read_double(h,"dt",dummy_real,0)
  dt = dummy_real(1)
!!$  dt = 0d0
!!$  time = 0d0
  call esio_field_establish (h, INT(ny), 1, INT(ny), INT(nx)/2, INT(Yxst)+1, INT(Yxsz), INT(nz), INT(Yzst)+1, INT(Yzsz), ierr)
  call esio_field_readv_double (h, "v",  Udata(1,1,1,1,2),2)
  call esio_field_readv_double (h, "g",  RKdata(1,1,1,1,1),2)
  call esio_field_readv_double (h, "theta",  RKdata(1,1,1,1,3),2)
  !call esio_field_readv_double (h, "u",  RKdata(1,1,1,1,3),2)

  if(myrank.eq.0) then
     call esio_line_establish (h, INT(ny), 1, INT(ny), ierr)
  else
     call esio_line_establish (h, INT(ny), 1, 0, ierr)
  endif
  call esio_line_read_double (h,"U_mean", dummy_u(1),1)
  call esio_line_read_double (h,"W_mean", dummy_w(1),1)
  call esio_line_read_double (h,"T_mean", dummy_t(1),1)
  !call esio_line_read_double (h,"U_mean", dummy_t(1),1)
  if(myrank.eq.0) then
     RKdata(1,:,1,1,1) = dummy_u(:)-U_gal
     RKdata(2,:,1,1,1) = 0d0
     RKdata(1,:,1,1,2) = dummy_w(:)
     RKdata(2,:,1,1,2) = 0d0
     RKdata(1,:,1,1,3) = dummy_t(:)
     RKdata(2,:,1,1,3) = 0d0
  endif

  call esio_file_close(h)
  call esio_handle_finalize(h)


  return
end subroutine data_read_vg

subroutine set_initial_data_vg(U,RK_store,kx,kz,LMFxr,LMFzr,ny,Yxsz,Yxst,Yzsz,Yzst)
  use labs
  implicit none
  integer(c_intptr_t) :: ny,Yxsz,Yxst,Yzsz,Yzst
  complex(c_double_complex) :: U(ny,Yxsz,Yzsz,7)
  complex(c_double_complex) :: RK_store(ny,Yxsz,Yzsz,3)
  integer :: kx(Yxsz), kz(Yzsz)
  real(c_double) :: LMFxr,LMFzr
  integer :: i,j,k,l,size,xsz
  complex(c_double_complex), parameter :: ci = (0d0,1d0)
  complex(c_double_complex) :: kxx, kzz, dxg, dxf, dzg, dzf

  size = Yxsz*Yzsz
  xsz = Yxsz
  

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l,kxx,kzz,dxg,dxf,dzg,dzf) &
  !$OMP FIRSTPRIVATE(LABS_D0_LU,LABS_D1_B,LABS_D2_B,kx,LMFxr,kz,LMFzr)
  !$OMP DO SCHEDULE(STATIC)
  do l=1,size
     k = mod(l-1,xsz) + 1
     j = (l-1)/xsz + 1

        if((kx(k).eq.0).and.(kz(j).eq.0)) then
           U(:,k,j,1) = RK_store(:,k,j,1)
           U(:,k,j,2) = (0d0,0d0)
           U(:,k,j,3) = RK_store(:,k,j,2)
	   U(:,k,j,7) = RK_store(:,k,j,3)
        else
           U(:,k,j,4) = U(:,k,j,2)
           call LABS_bs(LABS_D0_LU(:,:,1),U(:,k,j,4))
           call LABS_mat_mul(LABS_D1_B(:,:,1),U(:,k,j,4),U(:,k,j,5))
           call LABS_mat_mul(LABS_D2_B(:,:,1),U(:,k,j,4),RK_store(:,k,j,2))
           do i=1,ny
              RK_store(i,k,j,2) = RK_store(i,k,j,2) - ( (kx(k)*LMFxr)**2 + (kz(j)*LMFzr)**2 )*U(i,k,j,2)
           enddo
           kxx = ci*kx(k)*LMFxr/((kx(k)*LMFxr)**2 + (kz(j)*LMFzr)**2)
           kzz = ci*kz(j)*LMFzr/((kx(k)*LMFxr)**2 + (kz(j)*LMFzr)**2)
           do i=1,ny
              dxg = kxx*RK_STORE(i,k,j,1)
              dzg = kzz*RK_STORE(i,k,j,1)
              dxf = kxx*U(i,k,j,5)
              dzf = kzz*U(i,k,j,5)
              U(i,k,j,1) = dxf-dzg
              U(i,k,j,3) = dxg+dzf	      
           enddo
	   U(:,k,j,7) =RK_store(:,k,j,3)
        endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  return
end subroutine set_initial_data_vg









subroutine data_read_gphi(Udata,RKdata)

  use esio
  use global_variable
  use labs
  use ppdfftw
  use runge_kutta

  implicit none

  real(C_DOUBLE) :: Udata(2,ny,Yxsz,Yzsz,no_of_field)
  real(C_DOUBLE) :: RKdata(2,ny,Yxsz,Yzsz,3)
  real(C_DOUBLE) :: dummy_real(1)
  real(C_DOUBLE) :: dummy_u(ny), dummy_w(ny), dummy_t(ny)
  integer :: dummy_integer(1)
  logical :: read_write
  type(esio_handle) :: h


  call esio_handle_initialize(h,MPI_COMM_WORLD)
  call esio_file_open(h,restart_file_name,.true.,ierr)

  call esio_line_establish(h,1,1,1,ierr)

  call esio_line_read_integer(h,"timestep",dummy_integer,0)
  timestep = dummy_integer(1)
  target_timestep = target_timestep + timestep
  call esio_line_read_double(h,"time",dummy_real,0)
  time = dummy_real(1)
  call esio_line_read_double(h,"dt",dummy_real,0)
  dt = dummy_real(1)

  
  call esio_field_establish (h, INT(ny), 1, INT(ny), INT(nx)/2, INT(Yxst)+1, INT(Yxsz), INT(nz), INT(Yzst)+1, INT(Yzsz), ierr)
  call esio_field_readv_double (h, "g",  RKdata(1,1,1,1,1),2)
  call esio_field_readv_double (h, "phi",  RKdata(1,1,1,1,2),2)
  call esio_field_readv_double (h, "theta",  RKdata(1,1,1,1,3),2)

  if(myrank.eq.0) then
     call esio_line_establish (h, INT(ny), 1, INT(ny), ierr)
  else
     call esio_line_establish (h, INT(ny), 1, 0, ierr)
  endif
  call esio_line_read_double (h,"U_mean", dummy_u(:),1)
  call esio_line_read_double (h,"W_mean", dummy_w(:),1)
  call esio_line_read_double (h,"T_mean", dummy_t(:),1)
  if(myrank.eq.0) then
     RKdata(1,:,1,1,1) = dummy_u(:) - U_gal
     RKdata(2,:,1,1,1) = 0d0
     RKdata(1,:,1,1,2) = dummy_w(:)
     RKdata(2,:,1,1,2) = 0d0
     RKdata(1,:,1,1,3) = dummy_t(:) 
     RKdata(2,:,1,1,3) = 0d0
  endif

  call esio_file_close(h)
  call esio_handle_finalize(h)


  return
end subroutine data_read_gphi






subroutine set_initial_data_gphi(U,RK_store,kx,kz,LMFxr,LMFzr,ny,Yxsz,Yxst,Yzsz,Yzst,D0,D1,D2)
  use labs
  implicit none
  integer(c_intptr_t) :: ny,Yxsz,Yxst,Yzsz,Yzst
  complex(c_double_complex) :: U(ny,Yxsz,Yzsz,7)
  complex(c_double_complex) :: RK_store(ny,Yxsz,Yzsz,3)
  real(c_double), dimension(2*LABS_BW+1,LABS_NY) :: D0,D1,D2,D_buf
  integer :: kx(Yxsz), kz(Yzsz)
  real(c_double) :: LMFxr,LMFzr,k2
  integer :: i,j,k,l,m,size,xsz
  complex(c_double_complex), parameter :: ci = (0d0,1d0)
  complex(c_double_complex) :: kxx, kzz, dxg, dxf, dzg, dzf

  size = Yxsz*Yzsz
  xsz = Yxsz


  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l,m,kxx,kzz,dxg,dxf,dzg,dzf,k2) FIRSTPRIVATE(D0,D1,D2,D_buf,kx,LMFxr,kz,LMFzr)
  !$OMP DO SCHEDULE(STATIC)
  do l=1,size
     k = mod(l-1,xsz) + 1
     j = (l-1)/xsz + 1

     if((kx(k).eq.0).and.(kz(j).eq.0)) then
        U(:,k,j,1) = RK_store(:,k,j,1)
        U(:,k,j,3) = RK_store(:,k,j,2)
	U(:,k,j,7) = RK_store(:,k,j,3)
        U(:,k,j,2) = (0d0,0d0)	
     else
        k2 = (kz(j)*LMFzr)**2 + (kx(k)*LMFxr)**2
        kxx = ci*kx(k)*LMFxr/k2
        kzz = ci*kz(j)*LMFzr/k2

        do i=1,LABS_BW
           D_buf(i,1) = 0d0
        enddo
        D_buf(LABS_BW+1,1) = -k2*D0(LABS_BW+1,1) + D2(LABS_BW+1,1)
        do i=LABS_BW+2,2*LABS_BW+1
           D_buf(i,1) = 0d0
        enddo

        do m = 2,LABS_NY-1
           do i=1,2*LABS_BW+1
              D_buf(i,m) = -k2*D0(i,m) + D2(i,m)
           enddo
        enddo

        do i=1,LABS_BW
           D_buf(i,ny) = 0d0
        enddo
        D_buf(LABS_BW+1,labs_ny) = -k2*D0(LABS_BW+1,labs_ny) + D2(LABS_BW+1,labs_ny)
        do i=LABS_BW+2,2*LABS_BW+1
           D_buf(i,LABS_ny) = 0d0
        enddo

        call labs_lu(D_buf)

        U(:,k,j,5) = RK_store(:,k,j,2)
        U(1,k,j,5) = (0d0,0d0)
        U(LABS_NY,k,j,5) = (0d0,0d0)


        call labs_bs(D_buf,U(:,k,j,5))

        call labs_mat_mul(D0,U(:,k,j,5),U(:,k,j,2))
        call labs_mat_mul(D1,U(:,k,j,5),U(:,k,j,3))

        
        do i=1,LABS_NY
           dxg = kxx*RK_store(i,k,j,1)
           dzg = kzz*RK_store(i,k,j,1)
           dxf = kxx*U(i,k,j,3)
           dzf = kzz*U(i,k,j,3)
           U(i,k,j,1) = dxf-dzg
           U(i,k,j,3) = dxg+dzf
        enddo
	U(:,k,j,7) = RK_store(:,k,j,3)

     endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  return
end subroutine set_initial_data_gphi
