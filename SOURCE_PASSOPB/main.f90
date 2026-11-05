!---------------------------------------------------------------------------------
! 
! This program performs direct numerical simulation for turbulence flow 
! between two flat plates with passive scalar and oppostion control. 
! The wall-normal direction to the plate is Y-Dir. 
! 
! X-Dir : Streamwise Direction (Periodic Boundary Condition)
! Y-Dir : Wall Normal Direction (No Slip / No Mass Flux Boundary Condition)
! Z-Dir : Spanwise Direction (Periodic Boundary Condition)
! 
! This code is based on the paper 
! Kim Moin Moser J. Fluid Mech. (1987), vol. 177 pp. 133-166
! 
! B-spline is used for the operator of first and second derivatives estimation
! 
! 3/2 Dealliasing Technique is used for dealliasing for nonlinear products
! 
! FFTW library is used to transpose the global data and DFT of each data line
!
! Developer        :  Myoungkyu Lee   (mk@ices.utexas.edu)
! Software Support :  Nicholas Malaya (nick@ices.utexas.edu)
! PI               :  Robert Moser    (rmoser@ices.utexas.edu)
!
! Version : 0.0
! Date : 2013.01.04
!
!---------------------------------------------------------------------------------

program POONGBACK

  use ppdfftw
  use labs
  use runge_kutta
  use timer
  use global_variable
  implicit none
  integer :: i


  !--------------------------------------------------------
  ! initialize MPI envornment
  !--------------------------------------------------------
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,nprocs,ierr)
  call fftw_mpi_init 



  !--------------------------------------------------------
  ! Write program start time
  !--------------------------------------------------------
  if (myrank.eq.0)   call date_and_time(values=program_start_time)
  if (myrank.eq.0)   call print_program_info
  if (myrank.eq.0)   write(*,'(a9,i4,5(a1,i2))')  'Start at ',&
       program_start_time(1),'/', program_start_time(2),'/', program_start_time(3),' ', &
       program_start_time(5),':', program_start_time(6),':', program_start_time(7) 
  call mpi_barrier(MPI_COMM_WORLD,ierr)


  !--------------------------------------------------------
  ! Read Input File
  !--------------------------------------------------------
  call read_input_file
  call mpi_barrier(MPI_COMM_WORLD,ierr)


  !--------------------------------------------------------
  ! Initialize environment to use B-Spline
  !--------------------------------------------------------
  call labs_initialize(ny,stretch,Ly,omp_compute)
  call mpi_barrier(MPI_COMM_WORLD,ierr)


  !--------------------------------------------------------
  ! Initialize environment to use RK-method
  !--------------------------------------------------------
  call RK3_initialize
  call mpi_barrier(MPI_COMM_WORLD,ierr)


  !--------------------------------------------------------
  ! Set the Environment to use ppdfftw including memory allocation
  !--------------------------------------------------------
  call ppdfftw_initialize
  call mpi_barrier(MPI_COMM_WORLD,ierr)


  !--------------------------------------------------------
  ! If compiled with MASA, State this in the log
  !--------------------------------------------------------
!!$#ifdef MASA
!!$  if (myrank.eq.0)   write(*,*) "*------------------------*"
!!$  if (myrank.eq.0)   write(*,*) "  Initializing with MASA"
!!$  if (myrank.eq.0)   write(*,*) "*------------------------*"
!!$  if (myrank.eq.0)   write(*,*) " "
!!$#endif


  !--------------------------------------------------------
  ! Initializing the velocity field and proceed the first rk-substep
 !--------------------------------------------------------
  if(myrank.eq.0) write(*,*) "Program running log"
  if(myrank.eq.0) write(*,'(a40,a12,7a26)') "Event", "Timestep", "dt", "Simulation Time", &
       "Elapsed Time / Center vel","dP/dx", "Tau_wall = nu*dU/dy","dT/dx", "T_wall = nu/pr*dT/dy"
  if(myrank.eq.0) time_initial_step_start = MPI_WTIME()
  call initial_step
  if(myrank.eq.0) time_initial_step_end = MPI_WTIME()
  if(myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "Initial step", &
       timestep,dt,time,time_initial_step_end - time_initial_step_start
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  
  RK_SUBSTEP = 1

  !--------------------------------------------------------
  ! Write Initial Field and statistics
  !--------------------------------------------------------
  call print_labs_data
  if(myrank.eq.0) time_writing_field = 0d0
  if (restart_file_exist.eq.0) then
     if(myrank.eq.0) time_data_store_start = MPI_WTIME()
     call data_store(U,RK_STORE)
     if(myrank.eq.0) time_data_store_end = MPI_WTIME()
     if(myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "Data store in HDF5", &
          timestep,dt,time,time_data_store_end - time_data_store_start
     call mpi_barrier(MPI_COMM_WORLD,ierr)
  endif


  !--------------------------------------------------------
  ! Run Program
  !--------------------------------------------------------
  if (myrank.eq.0)   call date_and_time(values=simulation_start_time)
  do 
     if(myrank.eq.0 .and. RK_SUBSTEP.eq.1) time_one_full_step_start = MPI_WTIME()
     call time_march

     if (rk_substep.ne.rk_step) then
        
        rk_substep = rk_substep + 1

     else

        !---------------------------------------------------------
        ! Update status
        !---------------------------------------------------------
        if(myrank.eq.0) phi_dpdx(4) = (phi_dpdx(1) + phi_dpdx(2) + phi_dpdx(3))/dt
	if(myrank.eq.0) phi_dtdx(4) = (phi_dtdx(1) + phi_dtdx(2) + phi_dtdx(3))/dt
        if(myrank.eq.0) time_one_full_step_end = MPI_WTIME()
        if(myrank.eq.0) write(*,'(a40,i12,7ES26.15)') "One Full Timestep Computation", &
             timestep,dt,time,time_one_full_step_end-time_one_full_step_start,phi_dpdx(4), wall_shear_stress,&
		phi_dtdx(4), wall_theta_stress

        timestep = timestep + 1
        time = time + dt
        rk_substep = 1
        
        if(myrank.eq.0) call check_remaining_time
        call mpi_bcast(stop_program_or_not,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(data_store_or_not,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(stat_store_or_not,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call mpi_barrier(MPI_COMM_WORLD,ierr)

        !---------------------------------------------------------
        ! Write data if writing signal is on
        !---------------------------------------------------------
        if (data_store_or_not .eqv. .true.) then
           if(myrank.eq.0) time_data_store_start = MPI_WTIME()
           call data_store(U,RK_STORE)
           if(myrank.eq.0) time_data_store_end = MPI_WTIME()
           if(myrank.eq.0) write(*,'(a40,i12,3ES26.15)') "Data store in HDF5", &
                timestep,dt,time,time_data_store_end - time_data_store_start
           if(myrank.eq.0) time_writing_field = time_writing_field+time_data_store_end-time_data_store_start
           call mpi_barrier(MPI_COMM_WORLD,ierr)
           target_data_time = target_data_time + data_time_interval
        endif

	!--------------------------------------------
	! changed by jie to save final data with same name
	! comment: easy for restart 
	! data: 05/11/18
	!--------------------------------------------
        if (stop_program_or_not .eqv. .true.) then 
                !call data_store(U,RK_STORE)
		call data_store_f(U,RK_STORE)
		EXIT
	endif

     endif


     call mpi_barrier(MPI_COMM_WORLD,ierr)
  enddo
  


  !---------------------------------------------------------
  ! Finalize the program 
  ! (mostly deallocation and destroy FFTW plans)
  !---------------------------------------------------------
  call labs_finalize
  call mpi_barrier(MPI_COMM_WORLD,ierr) 
  call RK3_finalize
  call mpi_barrier(MPI_COMM_WORLD,ierr) 
  call ppdfftw_finalize
  call mpi_barrier(MPI_COMM_WORLD,ierr) 

  !---------------------------------------------------------
  ! Write Program Ending Time
  !---------------------------------------------------------
  if (myrank.eq.0)   call date_and_time(values=program_end_time)
  if (myrank.eq.0)   write(*,'(a9,i4,5(a,i2))') 'End at ',&
       program_end_time(1),'/', program_end_time(2),'/', program_end_time(3),' ', &
       program_end_time(5),':', program_end_time(6),':', program_end_time(7) 
  call mpi_finalize(ierr)

end program POONGBACK
