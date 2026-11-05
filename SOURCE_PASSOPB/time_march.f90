!--------------------------------------------------------------------
! This subroutine finish 1 full RK substep
! 
! It is necessary to use this subroutine because data type conversion
!
! Detailed Procedure
! 
! - move data from y pencil to x pencil
!      global data transpose from y pencil to z pencil
!      inverse DFT in Z direction with padding
!      global data transpose from z pencil to x pencil
! - realspace calculation
!      inverse DFT in X direction with padding
!      calculate nonlinear product
!      calculate time interval with CFL
!      forward DFT in X direction with packing
!
! - move data from x pencil to y pencil
!      global data transpose from x pencil to z pencil
!      forward DFT in Z direction with packing
!      global data transpose from z pencil to y pencil
!
! - wavespace calculation
!     [0,0] mode : solve simplified N-S equation 
!                  Enforce constant mass flux (Forcing)
!     non [0,0] mode : calculate Hg / Hv
!                      solve g (=du/dz - dw/dx)
!                      solve v and dvdy
!                      solve u and w
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------
subroutine time_march


  use global_variable
  use ppdfftw
  use runge_kutta
  use labs
  use timer
  implicit none
  integer(C_INTPTR_T) :: i

  if (stat_store_or_not .eqv. .true.) then
     if (RK_SUBSTEP .eq. 1) then
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call get_statistics_1(U)
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call write_spectra(U)
        call mpi_barrier(MPI_COMM_WORLD,ierr)
     endif
  endif






!!$  if (myrank.eq.0) time_y_to_x_start = MPI_WTIME()
  do i=1,3
     call y_to_x(U(:,i),U(:,i),Upad(:,i),Upad(:,i))
     call mpi_barrier(MPI_COMM_WORLD,ierr)
  enddo
  call y_to_x(U(:,7),U(:,7),Upad(:,7),Upad(:,7))
  call mpi_barrier(MPI_COMM_WORLD,ierr)
!!$  if (myrank.eq.0) time_y_to_x_end = MPI_WTIME()
!!$  if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "y to x (with DFT in Z)", &
!!$       timestep,dt,time,time_y_to_x_end-time_y_to_x_start
!!$  call mpi_barrier(MPI_COMM_WORLD,ierr)

!!$  if (myrank.eq.0) time_realspace_start = MPI_WTIME()
  call dft_x_realspace(Upad,buf_line_x_c,buf_line_x_r)
!!$  if (myrank.eq.0) time_realspace_end = MPI_WTIME()
!!$  if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "Realspace (with DFT in X)", &
!!$       timestep,dt,time,time_realspace_end-time_realspace_start
!!$  call mpi_barrier(MPI_COMM_WORLD,ierr)
  
!!$  if (myrank.eq.0) time_x_to_y_start = MPI_WTIME()
  do i=1,5
     call x_to_y(Upad(:,i),Upad(:,i),U(:,i),U(:,i))
     call mpi_barrier(MPI_COMM_WORLD,ierr)
  enddo

  if (stat_store_or_not .eqv. .true.) then
     if (RK_SUBSTEP .eq. 1) then
        call x_to_y(Upad(:,6),Upad(:,6),U(:,6),U(:,6))
        call mpi_barrier(MPI_COMM_WORLD,ierr)
     endif
  endif

  do i=7,10
     call x_to_y(Upad(:,i),Upad(:,i),U(:,i),U(:,i))
     call mpi_barrier(MPI_COMM_WORLD,ierr)
  enddo
  
  call rescale_data(U,ny,Yzsz,Yxsz)

  if (stat_store_or_not .eqv. .true.) then
     if (RK_SUBSTEP .eq. 1) then
        call get_statistics_2(U)
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        target_stat_time = target_stat_time + stat_time_interval
        stat_store_or_not = .false.
        if(myrank.eq.0) write(*,'(a40,i12,4ES26.15)') "Statistics calculation", & 
             timestep,dt,time,U_center
        call mpi_barrier(MPI_COMM_WORLD,ierr)
     endif
  endif

!!$  if (myrank.eq.0) time_x_to_y_end = MPI_WTIME()
!!$  if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "x to z (with DFT in Z)", &
!!$       timestep,dt,time,time_x_to_y_end-time_x_to_y_start
!!$  call mpi_barrier(MPI_COMM_WORLD,ierr)
  
!!$  if (myrank.eq.0) time_wavespace_start = MPI_WTIME()
  call wavespace(U,RK_store,buf_line_y_c1,buf_line_y_c2,LABS_D0_LU, LABS_D0_B, LABS_D1_B, LABS_D2_B, LABS_BUF, &
       kx,kz,labs_ny,ny,omp_y,Yzsz,Yxsz,labs_bw,nx32,nz32)
!!$  if (myrank.eq.0) time_wavespace_end = MPI_WTIME()
!!$  if (myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "Wavespace", &
!!$       timestep,dt,time,time_wavespace_end-time_wavespace_start
!!$  call mpi_barrier(MPI_COMM_WORLD,ierr)


  return 
end subroutine time_march
