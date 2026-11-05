!--------------------------------------------------------------------
! This subroutine perform DFT in X direction 
! and nonlinear product
! When RK_substep = 1, it also calculate the dt (time inteval)
! with given CFL number. 
! 
! Detailed Procedure
! 
! - Padding from Nx to 3/2 Nx
! - Inverse DFT from wavespace to realspace in X direction
! - Calculate nonlinear product
!     ** if RK_substep = 1 calculate max(u/dx,v/dy,w/dz)
! - Forward DFT from realspace to wavespace in X direction
! - Packing from 3/2 Nx to Nx
! - if RK_substep = 1, calcualte dt and broadcast the data 
! 
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------

subroutine dft_x_realspace(Up,buf_c,buf_r)
  use ppdfftw
  use runge_kutta
  use labs
  use global_variable
  use, intrinsic :: iso_c_binding
  implicit none
  complex(C_DOUBLE_COMPLEX) :: Up(nxh,Xzsz*Xysz,no_of_field)
  complex(C_DOUBLE_COMPLEX) :: buf_c(nx32hpad,no_of_field,omp_x)
  real(C_DOUBLE) :: buf_r(nx32,no_of_field,omp_x)
  real(C_DOUBLE) :: temp1,temp2,temp3,temp4,bv2,velmax_dummy,vel, temp_u, temp_v, temp_w
  integer :: i,j,k,ntotal,yy,l

  velmax_dummy=0d0
  ! This is needed to ensure compatibility of ISO_C_BINDING with some compilers (e.g. XL compiler)
  ntotal = Xzsz*Xysz

  if (stat_store_or_not .eqv. .true.) then
     if (RK_SUBSTEP .eq. 1) then
        stat_skewness(:,:) = 0d0
        stat_flatness(:,:) = 0d0
        stat_uuv(:) = 0d0
        stat_vvv(:) = 0d0
        stat_uvv(:) = 0d0
        stat_vww(:) = 0d0
        stat_real_omp_buf(:,:,:) = 0d0 
        do i=1,omp_compute
           do j=1,ny
              stat_real_omp_buf(1,j,i) = stat_U(j)-u_gal
              stat_real_omp_buf(2,j,i) = stat_W(j)
           enddo
        enddo
     endif
  endif

  



  call omp_set_num_threads(omp_compute)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,l,yy,temp1,temp2,temp3,bv2,vel,temp_u,temp_v,temp_w)
  l = omp_get_thread_num()+1
  !$OMP DO SCHEDULE(STATIC) REDUCTION(max:velmax_dummy)
  do i=1,ntotal

     yy = (i-1)/Xzsz + 1

     ! padding / transform complex to real
     do j=1,nx/2
        buf_c(j,1,l) = Up(j,i,1)
     enddo
     do j=((nx/2)+1),((nx32/2)+1)
        buf_c(j,1,l) = (0d0,0d0)
     enddo
     call fftw_execute_dft_c2r(Plan_transform_inverse_x,buf_c(:,1,l),buf_r(:,1,l))


     do j=1,nx/2
        buf_c(j,2,l) = Up(j,i,2)
     enddo
     do j=((nx/2)+1),((nx32/2)+1)
        buf_c(j,2,l) = (0d0,0d0)
     enddo
     call fftw_execute_dft_c2r(Plan_transform_inverse_x,buf_c(:,2,l),buf_r(:,2,l))


     do j=1,nx/2
        buf_c(j,3,l) = Up(j,i,3)
     enddo
     do j=((nx/2)+1),((nx32/2)+1)
        buf_c(j,3,l) = (0d0,0d0)
     enddo
     call fftw_execute_dft_c2r(Plan_transform_inverse_x,buf_c(:,3,l),buf_r(:,3,l))

     do j=1,nx/2
        buf_c(j,7,l) = Up(j,i,7)
     enddo
     do j=((nx/2)+1),((nx32/2)+1)
        buf_c(j,7,l) = (0d0,0d0)
     enddo
     call fftw_execute_dft_c2r(Plan_transform_inverse_x,buf_c(:,7,l),buf_r(:,7,l))


     ! Calculating nonlinear product
     do j=1,nx32
        temp1 = buf_r(j,1,l)
        temp2 = buf_r(j,2,l)
        temp3 = buf_r(j,3,l)
        temp4 = buf_r(j,7,l)
        bv2=temp2*temp2
	buf_r(j,9,l)=temp3*temp4
	buf_r(j,8,l)=temp2*temp4
	buf_r(j,7,l)=temp1*temp4
        buf_r(j,5,l)=temp1*temp3
        buf_r(j,4,l)=temp2*temp3
        buf_r(j,3,l)=temp3*temp3-bv2
        buf_r(j,2,l)=temp1*temp2
        buf_r(j,1,l)=temp1*temp1-bv2
        buf_r(j,10,l)=temp2

        if (RK_SUBSTEP.eq.1) then
           vel = (dabs(temp1)*dxr)+(dabs(temp2)*dyr(yy))+(dabs(temp3)*dzr)
           velmax_dummy = max(vel,velmax_dummy)
           if (stat_store_or_not .eqv. .true.) then
              buf_r(j,6,l) = bv2
              temp_u = temp1 !- stat_real_omp_buf(1,yy+xyst,l)
              temp_v = temp2
              temp_w = temp3 !- stat_real_omp_buf(2,yy+xyst,l)
              
              stat_real_omp_buf(3,yy+xyst,l)  = stat_real_omp_buf(3,yy+xyst,l)  + temp_u*temp_u*temp_u
              stat_real_omp_buf(4,yy+xyst,l)  = stat_real_omp_buf(4,yy+xyst,l)  + temp_v*temp_v*temp_v
              stat_real_omp_buf(5,yy+xyst,l)  = stat_real_omp_buf(5,yy+xyst,l)  + temp_w*temp_w*temp_W
              stat_real_omp_buf(6,yy+xyst,l)  = stat_real_omp_buf(6,yy+xyst,l)  + temp_u*temp_u*temp_u*temp_u
              stat_real_omp_buf(7,yy+xyst,l)  = stat_real_omp_buf(7,yy+xyst,l)  + temp_v*temp_v*temp_v*temp_v
              stat_real_omp_buf(8,yy+xyst,l)  = stat_real_omp_buf(8,yy+xyst,l)  + temp_w*temp_w*temp_w*temp_w
              stat_real_omp_buf(9,yy+xyst,l)  = stat_real_omp_buf(9,yy+xyst,l)  - temp_u*temp_u*temp_v
              stat_real_omp_buf(10,yy+xyst,l) = stat_real_omp_buf(10,yy+xyst,l) - temp_v*temp_w*temp_w
              stat_real_omp_buf(11,yy+xyst,l) = stat_real_omp_buf(11,yy+xyst,l) - temp_v*temp_v*temp_v
              stat_real_omp_buf(12,yy+xyst,l) = stat_real_omp_buf(12,yy+xyst,l) - temp_u*temp_v*temp_v

              ! high order statistics
	      stat_real_omp_buf(13,yy+xyst,l)  = stat_real_omp_buf(13,yy+xyst,l)  + temp_u**6
	      stat_real_omp_buf(14,yy+xyst,l)  = stat_real_omp_buf(14,yy+xyst,l)  + temp_v**6
	      stat_real_omp_buf(15,yy+xyst,l)  = stat_real_omp_buf(15,yy+xyst,l)  + temp_w**6
	      stat_real_omp_buf(16,yy+xyst,l)  = stat_real_omp_buf(16,yy+xyst,l)  + temp_u**8
	      stat_real_omp_buf(17,yy+xyst,l)  = stat_real_omp_buf(17,yy+xyst,l)  + temp_v**8
	      stat_real_omp_buf(18,yy+xyst,l)  = stat_real_omp_buf(18,yy+xyst,l)  + temp_w**8
	      stat_real_omp_buf(19,yy+xyst,l)  = stat_real_omp_buf(19,yy+xyst,l)  + temp_u**10
	      stat_real_omp_buf(20,yy+xyst,l)  = stat_real_omp_buf(20,yy+xyst,l)  + temp_v**10
	      stat_real_omp_buf(21,yy+xyst,l)  = stat_real_omp_buf(21,yy+xyst,l)  + temp_w**10
           endif
        endif

     enddo
 

     ! packing / transform real to complex
     call fftw_execute_dft_r2c(Plan_transform_forward_x,buf_r(:,1,l),buf_c(:,1,l))
     do j=1,nx/2
        Up(j,i,1) = buf_c(j,1,l)
     enddo

     call fftw_execute_dft_r2c(Plan_transform_forward_x,buf_r(:,2,l),buf_c(:,2,l))
     do j=1,nx/2
        Up(j,i,2) = buf_c(j,2,l)
     enddo

     call fftw_execute_dft_r2c(Plan_transform_forward_x,buf_r(:,3,l),buf_c(:,3,l))
     do j=1,nx/2
        Up(j,i,3) = buf_c(j,3,l)
     enddo

     call fftw_execute_dft_r2c(Plan_transform_forward_x,buf_r(:,4,l),buf_c(:,4,l))
     do j=1,nx/2
        Up(j,i,4) = buf_c(j,4,l)
     enddo

     call fftw_execute_dft_r2c(Plan_transform_forward_x,buf_r(:,5,l),buf_c(:,5,l))
     do j=1,nx/2
        Up(j,i,5) = buf_c(j,5,l)
     enddo

     call fftw_execute_dft_r2c(Plan_transform_forward_x,buf_r(:,10,l),buf_c(:,10,l))
     do j=1,nx/2
        Up(j,i,10) = buf_c(j,10,l)
     enddo
     
     if ((RK_SUBSTEP.eq.1) .and. (stat_store_or_not .eqv. .true.)) then
        call fftw_execute_dft_r2c(Plan_transform_forward_x,buf_r(:,6,l),buf_c(:,6,l))
        do j=1,nx/2
           Up(j,i,6) = buf_c(j,6,l)
        enddo
     endif

     call fftw_execute_dft_r2c(Plan_transform_forward_x,buf_r(:,7,l),buf_c(:,7,l))
     do j=1,nx/2
        Up(j,i,7) = buf_c(j,7,l)
     enddo

     call fftw_execute_dft_r2c(Plan_transform_forward_x,buf_r(:,8,l),buf_c(:,8,l))
     do j=1,nx/2
        Up(j,i,8) = buf_c(j,8,l)
     enddo

     call fftw_execute_dft_r2c(Plan_transform_forward_x,buf_r(:,9,l),buf_c(:,9,l))
     do j=1,nx/2
        Up(j,i,9) = buf_c(j,9,l)
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL 


!!$  stat_skewness(:,1) = stat_real_omp_buf(:,3,1)
!!$  stat_skewness(:,2) = stat_real_omp_buf(:,4,1)
!!$  stat_skewness(:,3) = stat_real_omp_buf(:,5,1)
!!$  stat_flatness(:,1) = stat_real_omp_buf(:,6,1)
!!$  stat_flatness(:,2) = stat_real_omp_buf(:,7,1)
!!$  stat_flatness(:,3) = stat_real_omp_buf(:,8,1)
!!$  stat_h6(:,1) = stat_real_omp_buf(:,13,1)
!!$  stat_h6(:,2) = stat_real_omp_buf(:,14,1)
!!$  stat_h6(:,3) = stat_real_omp_buf(:,15,1)
!!$  stat_h8(:,1) = stat_real_omp_buf(:,16,1)
!!$  stat_h8(:,2) = stat_real_omp_buf(:,17,1)
!!$  stat_h8(:,3) = stat_real_omp_buf(:,18,1)
!!$  stat_h10(:,1) = stat_real_omp_buf(:,19,1)
!!$  stat_h10(:,2) = stat_real_omp_buf(:,20,1)
!!$  stat_h10(:,3) = stat_real_omp_buf(:,21,1)
!!$  stat_uuv(:) = stat_real_omp_buf(:,9,1)
!!$  stat_vww(:) = stat_real_omp_buf(:,10,1)
!!$  stat_vvv(:) = stat_real_omp_buf(:,11,1)
!!$  stat_uvv(:) = stat_real_omp_buf(:,12,1)

  if (stat_store_or_not .eqv. .true.) then
     if (RK_SUBSTEP .eq. 1) then
        ! make a local sum from omp threads data
        do l=1,omp_compute
           do i=1,ny
           stat_skewness(i,1) = stat_skewness(i,1) + stat_real_omp_buf(3,i,l)
           stat_skewness(i,2) = stat_skewness(i,2) + stat_real_omp_buf(4,i,l)
           stat_skewness(i,3) = stat_skewness(i,3) + stat_real_omp_buf(5,i,l)
           stat_flatness(i,1) = stat_flatness(i,1) + stat_real_omp_buf(6,i,l)
           stat_flatness(i,2) = stat_flatness(i,2) + stat_real_omp_buf(7,i,l)
           stat_flatness(i,3) = stat_flatness(i,3) + stat_real_omp_buf(8,i,l)

	   stat_hs6(i,1) = stat_hs6(i,1) + stat_real_omp_buf(13,i,l)
           stat_hs6(i,2) = stat_hs6(i,2) + stat_real_omp_buf(14,i,l)
           stat_hs6(i,3) = stat_hs6(i,3) + stat_real_omp_buf(15,i,l)
	   stat_hs8(i,1) = stat_hs8(i,1) + stat_real_omp_buf(16,i,l)
           stat_hs8(i,2) = stat_hs8(i,2) + stat_real_omp_buf(17,i,l)
           stat_hs8(i,3) = stat_hs8(i,3) + stat_real_omp_buf(18,i,l)
	   stat_hs10(i,1) = stat_hs10(i,1) + stat_real_omp_buf(19,i,l)
           stat_hs10(i,2) = stat_hs10(i,2) + stat_real_omp_buf(20,i,l)
           stat_hs10(i,3) = stat_hs10(i,3) + stat_real_omp_buf(21,i,l)

           stat_uuv(i) = stat_uuv(i) + stat_real_omp_buf(9,i,l)
           stat_vww(i) = stat_vww(i) + stat_real_omp_buf(10,i,l)
           stat_vvv(i) = stat_vvv(i) + stat_real_omp_buf(11,i,l)
           stat_uvv(i) = stat_uvv(i) + stat_real_omp_buf(12,i,l)
           enddo
        enddo

        ! calculate_skewness
        call mpi_reduce(stat_skewness(:,1),stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           stat_skewness(1,1) = 0d0
           do i=2,ny-1
              stat_skewness(i,1)  = rescale*stat_buf_line(i)
           enddo
        endif
           stat_skewness(ny,1) = 0d0
        call mpi_barrier(MPI_COMM_WORLD,ierr)

        call mpi_reduce(stat_skewness(:,2),stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           stat_skewness(1,2) = 0d0
           do i=2,ny-1
              stat_skewness(i,2)  = rescale*stat_buf_line(i)
           enddo
           stat_skewness(ny,2) = 0d0
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)

        call mpi_reduce(stat_skewness(:,3),stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           stat_skewness(1,3) = 0d0
           do i=2,ny-1
              stat_skewness(i,3)  = rescale*stat_buf_line(i)
           enddo
           stat_skewness(ny,3) = 0d0
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)

        ! calculate_flatness
        call mpi_reduce(stat_flatness(:,1),stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           stat_flatness(1,1) = 0d0
           do i=2,ny-1
              stat_flatness(i,1)  = rescale*stat_buf_line(i)
           enddo
           stat_flatness(ny,1) = 0d0
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)

        call mpi_reduce(stat_flatness(:,2),stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           stat_flatness(1,2) = 0d0
           do i=2,ny-1
              stat_flatness(i,2)  = rescale*stat_buf_line(i)
           enddo
           stat_flatness(ny,2) = 0d0
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)

        call mpi_reduce(stat_flatness(:,3),stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           stat_flatness(1,3) = 0d0
           do i=2,ny-1
              stat_flatness(i,3)  = rescale*stat_buf_line(i)
           enddo
           stat_flatness(ny,3) = 0d0
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)


	! calculate_highorder_stat_6th
        call mpi_reduce(stat_hs6(:,1),stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           stat_hs6(1,1) = 0d0
           do i=2,ny-1
              stat_hs6(i,1)  = rescale*stat_buf_line(i)
           enddo
        endif
           stat_hs6(ny,1) = 0d0
        call mpi_barrier(MPI_COMM_WORLD,ierr)

        call mpi_reduce(stat_hs6(:,2),stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           stat_hs6(1,2) = 0d0
           do i=2,ny-1
              stat_hs6(i,2)  = rescale*stat_buf_line(i)
           enddo
           stat_hs6(ny,2) = 0d0
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)

        call mpi_reduce(stat_hs6(:,3),stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           stat_hs6(1,3) = 0d0
           do i=2,ny-1
              stat_hs6(i,3)  = rescale*stat_buf_line(i)
           enddo
           stat_hs6(ny,3) = 0d0
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)

	! calculate_highorder_stat_8th
        call mpi_reduce(stat_hs8(:,1),stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           stat_hs8(1,1) = 0d0
           do i=2,ny-1
              stat_hs8(i,1)  = rescale*stat_buf_line(i)
           enddo
        endif
           stat_hs8(ny,1) = 0d0
        call mpi_barrier(MPI_COMM_WORLD,ierr)

        call mpi_reduce(stat_hs8(:,2),stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           stat_hs8(1,2) = 0d0
           do i=2,ny-1
              stat_hs8(i,2)  = rescale*stat_buf_line(i)
           enddo
           stat_hs8(ny,2) = 0d0
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)

        call mpi_reduce(stat_hs8(:,3),stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           stat_hs8(1,3) = 0d0
           do i=2,ny-1
              stat_hs8(i,3)  = rescale*stat_buf_line(i)
           enddo
           stat_hs8(ny,3) = 0d0
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)

	! calculate_highorder_stat_10th
        call mpi_reduce(stat_hs10(:,1),stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           stat_hs10(1,1) = 0d0
           do i=2,ny-1
              stat_hs10(i,1)  = rescale*stat_buf_line(i)
           enddo
        endif
           stat_hs10(ny,1) = 0d0
        call mpi_barrier(MPI_COMM_WORLD,ierr)

        call mpi_reduce(stat_hs10(:,2),stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           stat_hs10(1,2) = 0d0
           do i=2,ny-1
              stat_hs10(i,2)  = rescale*stat_buf_line(i)
           enddo
           stat_hs10(ny,2) = 0d0
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)

        call mpi_reduce(stat_hs10(:,3),stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           stat_hs10(1,3) = 0d0
           do i=2,ny-1
              stat_hs10(i,3)  = rescale*stat_buf_line(i)
           enddo
           stat_hs10(ny,3) = 0d0
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)



        ! Calculate turbulent transport terms for turb transport eqaution
        ! This will be stored in get_statistics2 subroutine
        ! Calculate -duuv/dy
        call mpi_reduce(stat_uuv,stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           do i=1,ny
              stat_buf_line_dummy1(i) = dcmplx(stat_buf_line(i)*rescale,0d0)
           enddo
           call LABS_bs(LABS_D0_LU(:,:,1),stat_buf_line_dummy1)
           call LABS_MAT_MUL(LABS_D1_B(:,:,1), stat_buf_line_dummy1, stat_buf_line_dummy2)
           do i=1,ny
              stat_uuv(i) = real(stat_buf_line_dummy2(i))
           enddo
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)

        ! Calculate -dvww/dy
        call mpi_reduce(stat_vww,stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           do i=1,ny
              stat_buf_line_dummy1(i) = dcmplx(stat_buf_line(i)*rescale,0d0)
           enddo
           call LABS_bs(LABS_D0_LU(:,:,1),stat_buf_line_dummy1)
           call LABS_MAT_MUL(LABS_D1_B(:,:,1), stat_buf_line_dummy1, stat_buf_line_dummy2)
           do i=1,ny
              stat_vww(i) = real(stat_buf_line_dummy2(i))
           enddo
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)

        ! Calculate -dvvv/dy
        call mpi_reduce(stat_vvv,stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           do i=1,ny
              stat_buf_line_dummy1(i) = dcmplx(stat_buf_line(i)*rescale,0d0)
           enddo
           call LABS_bs(LABS_D0_LU(:,:,1),stat_buf_line_dummy1)
           call LABS_MAT_MUL(LABS_D1_B(:,:,1), stat_buf_line_dummy1, stat_buf_line_dummy2)
           do i=1,ny
              stat_vvv(i) = real(stat_buf_line_dummy2(i))
           enddo
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)

        ! Calculate -duvv/dy
        call mpi_reduce(stat_uvv,stat_buf_line,int(ny),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           do i=1,ny
              stat_buf_line_dummy1(i) = dcmplx(stat_buf_line(i)*rescale,0d0)
           enddo
           call LABS_bs(LABS_D0_LU(:,:,1),stat_buf_line_dummy1)
           call LABS_MAT_MUL(LABS_D1_B(:,:,1), stat_buf_line_dummy1, stat_buf_line_dummy2)
           do i=1,ny
              stat_uvv(i) = real(stat_buf_line_dummy2(i))
           enddo
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)

     endif
  endif




  ! Calculate velocity maximum and dt
  if (RK_SUBSTEP.eq.1) then
     call mpi_allreduce(velmax_dummy,velmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr) 
     dt = dt_factor*dt + (1d0-dt_factor)*CFLx/velmax
     !dt = 0.01d0
     if (dt.le.1d-14 .or. dt.gt.1d5) then 
        if(myrank.eq.0) print*, "ABORT due to timestep dt =", dt
        call mpi_abort(MPI_COMM_WORLD,1,ierr)
     endif
  endif

  return
end subroutine dft_x_realspace
