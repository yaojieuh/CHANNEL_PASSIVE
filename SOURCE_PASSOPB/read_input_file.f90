!--------------------------------------------------------------------
! This subroutine read input variable to run the program 
! from the file 'input'
! 
! It should be able to check the input variable is valid
! If not, it should be able to indicate which one is not valid.
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------

subroutine read_input_file
  use ppdfftw
  use runge_kutta
  use global_variable
  use labs
  implicit none
  integer :: nxtest,nz32test,nztest,nytest
  logical :: file_exists

  ! DO NOT USE MPI_BCAST to broadcast input data
  ! It may not work on some machines, e.g. BG/Q
  ! Because some values are pointer variables

  if (nprocs.eq.4096) then
     open(28,file='input_4096')
  else if (nprocs.eq.8192) then
     open(28,file='input_8192')
  else if (nprocs.eq.16384) then
     open(28,file='input_16384')
  else if (nprocs.eq.32768) then
     open(28,file='input_32768') 
  else
     open(28,file='input')
  endif

  read (28,fmt='(1x)')
  read (28,*) nx,ny,nz,no_of_field
  read (28,fmt='(1x)')
  read (28,*) dimsize(1),dimsize(2)
  read (28,fmt='(1x)')
  read (28,*) show_mpi_topology
  read (28,fmt='(1x)')
  read (28,*) rk_scheme
  read (28,fmt='(1x)')
  read (28,*) labs_k
  read (28,fmt='(1x)')
  read (28,*) LMFx,LMFy,LMFz
  read (28,fmt='(1x)')
  read (28,*) stretch
  read (28,fmt='(1x)')
  read (28,*) WRITE_LABS_INFO
  read (28,fmt='(1x)')
  read (28,*) nu,pr
  read (28,fmt='(1x)')
  read (28,*) bulk_vel_u,bulk_vel_w,bulk_t
  read (28,fmt='(1x)')
  read (28,*) u_gal
  read (28,fmt='(1x)')
  read (28,*) omp_reorder
  read (28,fmt='(1x)')
  read (28,*) omp_compute
  read (28,fmt='(1x)')
  read (28,*) CFL
  read (28,fmt='(1x)')
  read (28,*) dt_factor
  read (28,fmt='(1x)')
  read (28,*) dy_factor
  read (28,fmt='(1x)')
  read (28,*) program_run_hours
  read (28,fmt='(1x)')
  read (28,*) target_timestep
  read (28,fmt='(1x)')
  read (28,*) data_time_interval
  read (28,fmt='(1x)')
  read (28,*) write_u,write_v,write_w,write_g,write_phi,write_t
  read (28,fmt='(1x)')
  read (28,*) stat_time_interval
  read (28,fmt='(1x)')
  read (28,'(A)') output_file_path
  read (28,fmt='(1x)')
  read (28,*) restart_file_exist
  read (28,fmt='(1x)')
  read (28,*) restart_option
  read (28,fmt='(1x)')
  read (28,'(A)') restart_file_name
  close(28)

  restart_file_name=trim(adjustl(restart_file_name))
  inquire(FILE=restart_file_name, EXIST=file_exists)
  
  nxh = nx/2
  nxhpad = nxh + 1
  nz32 = nz*3/2
  nx32 = nx*3/2
  nx32hpad = nx*3/4+1
  abort_signal=0
  Lx = LMFx * 2 * pi
  Ly = LMFy * 2
  Lz = LMFz * 2 * pi
  Re = real(1.0,kind(Re))/nu ! To ensure we don't get 0 if nu > 1 (though nu really should never be this)
  rescale = 1.d0/(nx32*nz32)
  dxr = dble(nx)/Lx
  dzr = dble(nz)/Lz
  CFLx = CFL / pi

  ! Checking the load balance for transpose
  nxtest = nx
  nztest =nz
  nz32test = nz32
  nytest = ny

  LMFxr = 1d0/LMFx
  LMFzr = 1d0/LMFz

  
  if(rk_scheme.eq.3) rk_step=3

  if (myrank.eq.0) then
     print*, "---------------------------------------------------------------------"
     write(*,'(a6,i5)') " Nx = ",nx
     write(*,'(a6,i5)') " Ny = ",ny
     write(*,'(a6,i5)') " Nz = ",nz
     write(*,'(a10,i5)') " Nx_pad = ",nx32
     write(*,'(a10,i5)') " Nz_pad = ",nz32
     write(*,'(a20,i5)') " Number of fields = ",no_of_field
     write(*,'(a11,f4.1)') " Lx = 2*pi*",LMFx
     write(*,'(a11,f4.1)') " Ly =  2.0*",LMFy
     write(*,'(a11,f4.1)') " Lz = 2*pi*",LMFz
     write(*,'(a7,f4.1)') " CFL = ", CFL
     write(*,'(a28,f4.1)') " dy multiplication factor = ", dy_factor
     write(*,'(a6,es14.6)') " nu = ", nu
     write(*,'(a11,f10.2)') " Re_bulk = ", Re
     write(*,'(a22,f4.1,a18)') " Program will run for ", program_run_hours, " wall clock hours."
     write(*,'(a31)',advance='no') " Fields that will be written : " 
     if(write_u .ne. 0) write(*,'(a2)', advance='no') "u "
     if(write_v .ne. 0) write(*,'(a2)', advance='no') "v "
     if(write_w .ne. 0) write(*,'(a2)', advance='no') "w "
     if(write_g .ne. 0) write(*,'(a13)', advance='no') "g(= Omega y) "
     if(write_phi .ne. 0) write(*,'(a19)', advance='no') "phi(= Laplacian v) "
     write(*,'(a1)') " " 
     print*, "---------------------------------------------------------------------"

     if(target_timestep.le.0) then
        print*, "target_time_step is wrong!!!"
        print*, "---------------------------------------------------------------------"
        abort_signal=1
     else 
        write(*,'(a19,i10)') " Total Time Step = ", target_timestep
        print*, "---------------------------------------------------------------------"
     endif

     if(no_of_field.lt.6) then
        print*, "Minimum number of field is 6!!!"
        print*, "---------------------------------------------------------------------"
        abort_signal=1
     endif


!!$     if((iostep.le.0).or.(iostep.gt.target_timestep)) then
!!$        print*, "iostep is wrong or greater than target_timestep!!!"
!!$        print*, "---------------------------------------------------------------------"
!!$        abort_signal=1
!!$     else
!!$        write(*,'(a40,i0,a7)') " Velocity Field will be stored at every ",iostep," steps."
!!$        print*, "---------------------------------------------------------------------"
!!$     endif

        
     if(((mod(nxtest,dimsize(1)).eq.0).and.(mod(nz32test,dimsize(1)).eq.0)) .and. &
          ((mod(nztest,dimsize(2)).eq.0).and.(mod(nytest,dimsize(2)).eq.0))) then
        print*, "Good load balace for transpose!" 
        print*, "---------------------------------------------------------------------"
     else
        print*, "WARNING : BAD LOAD BALANCE!!"
        print*, "Transpose routine will be very slow."
        print*, "Stop running the code and check MPI topology." 
        print*, " "
        print*, "--Good Topology Condition--"
        print*, "mod(NX,dimsize1)     = 0"
        print*, "mod(NZ_3/2,dimsize1) = 0"
        print*, "mod(NZ,dimsize2)     = 0"
        print*, "mod(NY,dimsize2)     = 0"
        print*, "---------------------------------------------------------------------"
        abort_signal = 1
     endif

     if(rk_scheme.ne.3) then
        print*, "Currently only 3rd order Runge-Kutta methos is supported."
        print*, "Sorry....... :( "
        print*, "---------------------------------------------------------------------"
        abort_signal = 1
     endif


     if(restart_file_exist.eq.0) then
        print*, "Start without Restart File." 
        print*, "---------------------------------------------------------------------"
     else
        if(file_exists .eqv. .true.) then
           print*, "Restart file is provided. ", trim(adjustl(restart_file_name))
                write(*,'(a28)',advance='no') " Fields that will be read : " 
                if(restart_option .eq. 1) then 
                   write(*,'(a5)') "u v w"
                else if(restart_option .eq. 2) then
                   write(*,'(a5)') "g phi"
                else
                   write(*,'(a3)') "v g"
                endif
        else
           print*, "Restart file does not exist. ",  trim(adjustl(restart_file_name))
        endif
        print*, "---------------------------------------------------------------------"
     endif

     if(nprocs.ne.dimsize(1)*dimsize(2)) then
        print*, "Dimension for MPI topology is not correct!!!" 
        print*, "---------------------------------------------------------------------"
        abort_signal = 1
     endif

     if(ny.le.2*labs_k) then
        print*, "K, order of B-spline, is too large !!!" 
        print*, "---------------------------------------------------------------------"
        abort_signal = 1
     endif
     
     if (abort_signal.ne.0) then
        print*, "Check your input file"
        print*, " "
        call system('cat input')
        ierr = 2
        print*, " "
        print*, " "
        call mpi_abort(MPI_COMM_WORLD,ierr)
     endif
     open (99, file='RUNNING')
     write(99,*) 'This file indicates a job running in this directory.'
     write(99,*) 'Delete this file to cleanly terminate the process.'
     close(99)
  endif



  return
end subroutine read_input_file
