!--------------------------------------------------------------------
! This subroutine calculate the running time / remaining time 
! to make dicisions for
! 1. program running (or not)
! 2. Data store (or not)
! 3. Calculate Statistics (or not)
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------
subroutine check_remaining_time
  use global_variable
  use timer
  use runge_kutta
  implicit none
  logical :: file_exist
  integer :: timestep_code_run
  integer :: timestep_code_remain
  real(C_DOUBLE) :: avg_time_per_step
  real(C_DOUBLE) :: elapsed_time1, elapsed_time2, remained_time
  call date_and_time(values=program_current_time)

  timestep_code_run = timestep - timestep_code_start

  if(timestep_code_run .eq. 1) then
     if(program_start_time(3).eq.1) then
        if((program_start_time(2) .eq. 1 ) .or. &
             (program_start_time(2) .eq. 2) .or. &
             (program_start_time(2) .eq. 4) .or. &
             (program_start_time(2) .eq. 6) .or. &
             (program_start_time(2) .eq. 8) .or. &
             (program_start_time(2) .eq. 9) .or. &
             (program_start_time(2) .eq. 11) ) then 
           program_start_time(3) = program_start_time(3) + 31
        else if((program_start_time(2) .eq. 5) .or. &
             (program_start_time(2) .eq. 7) .or. &
             (program_start_time(2) .eq. 10) .or. &
             (program_start_time(2) .eq. 12)) then
           program_start_time(3) = program_start_time(3) + 30
        else 
           program_start_time(3) = program_start_time(3) + 28
        endif
     endif
     if(simulation_start_time(3).eq.1) then
        if((simulation_start_time(2) .eq. 1 ) .or. &
             (simulation_start_time(2) .eq. 2) .or. &
             (simulation_start_time(2) .eq. 4) .or. &
             (simulation_start_time(2) .eq. 6) .or. &
             (simulation_start_time(2) .eq. 8) .or. &
             (simulation_start_time(2) .eq. 9) .or. &
             (simulation_start_time(2) .eq. 11) ) then 
           simulation_start_time(3) = simulation_start_time(3) + 31
        else if((simulation_start_time(2) .eq. 5) .or. &
             (simulation_start_time(2) .eq. 7) .or. &
             (simulation_start_time(2) .eq. 10) .or. &
             (simulation_start_time(2) .eq. 12)) then
           simulation_start_time(3) = simulation_start_time(3) + 30
        else 
           simulation_start_time(3) = simulation_start_time(3) + 28
        endif
     endif
  endif



  ! check elapsed time
  if(program_current_time(3).eq.1) then
     if((program_current_time(2) .eq. 1 ) .or. &
          (program_current_time(2) .eq. 2) .or. &
          (program_current_time(2) .eq. 4) .or. &
          (program_current_time(2) .eq. 6) .or. &
          (program_current_time(2) .eq. 8) .or. &
          (program_current_time(2) .eq. 9) .or. &
          (program_current_time(2) .eq. 11) ) then 
        program_current_time(3) = program_current_time(3) + 31
     else if((program_current_time(2) .eq. 5) .or. &
          (program_current_time(2) .eq. 7) .or. &
          (program_current_time(2) .eq. 10) .or. &
          (program_current_time(2) .eq. 12)) then
        program_current_time(3) = program_current_time(3) + 30
     else 
        program_current_time(3) = program_current_time(3) + 28
     endif
  endif

  elapsed_time1 = (program_current_time(3) - program_start_time(3)) * 1440d0 &
       + (program_current_time(5) - program_start_time(5)) * 60d0 &
       + (program_current_time(6) - program_start_time(6)) * 1d0 &
       + (program_current_time(7) - program_start_time(7)) * 0.016666666666667d0 &
       + (program_current_time(8) - program_start_time(8)) * 0.0166666666666667d-3 
 
  elapsed_time2 = (program_current_time(3) - simulation_start_time(3)) * 1440d0 &
       + (program_current_time(5) - simulation_start_time(5)) * 60d0 &
       + (program_current_time(6) - simulation_start_time(6)) * 1d0 &
       + (program_current_time(7) - simulation_start_time(7)) * 0.016666666666667d0 &
       + (program_current_time(8) - simulation_start_time(8)) * 0.0166666666666667d-3 &
       - time_writing_field*0.01666666666666667d0


  remained_time = program_run_hours*60d0 - elapsed_time1

  ! check how many steps remained
  timestep_code_remain = target_timestep - timestep

  avg_time_per_step = dble(timestep_code_run)/elapsed_time2


  ! decide program stop or not
  if (timestep_code_run.eq.0) then
     target_code_timestep = 1000000
  else
     target_code_timestep = int(avg_time_per_step * remained_time) 
  endif
  

  !default flag
  stop_program_or_not = .false. ! program stop or not
  data_store_or_not = .false. ! write field or not
  stat_store_or_not = .false. ! write field or not
  


  ! Check if data need to be stored or not
  if (time.ge.target_data_time) data_store_or_not = .true. 
  !if (mod(timestep,10).eq.0) data_store_or_not = .true. 
  if (time.ge.target_stat_time) stat_store_or_not = .true. 



  ! decide if keep running program or not
  if ((target_code_timestep.le.40) .or. (timestep.eq.target_timestep)) then
     stop_program_or_not = .true.
     !data_store_or_not = .true.
  endif

  inquire(file='RUNNING', exist=file_exist)
  if(.not. file_exist) then
     stop_program_or_not = .true.
     print*, 'RUNNING deleted !'
  end if


  return
end subroutine check_remaining_time
