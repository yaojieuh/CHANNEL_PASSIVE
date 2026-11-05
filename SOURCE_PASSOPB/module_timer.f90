!--------------------------------------------------------------------
! This module contains variables to be used 
! when measuring average and standard deviation of elapsed time
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------


module timer
  use, intrinsic :: iso_c_binding 

  implicit none

  real(c_double) :: time_transpose_start
  real(c_double) :: time_transpose_end
  real(c_double) :: time_transform_start
  real(c_double) :: time_transform_end
  real(c_double) :: time_one_full_step_start
  real(c_double) :: time_one_full_step_end
  real(c_double) :: time_initial_step_start
  real(c_double) :: time_initial_step_end
  real(c_double) :: time_realspace_start
  real(c_double) :: time_realspace_end
  real(c_double) :: time_wavespace_start
  real(c_double) :: time_wavespace_end
  real(c_double) :: time_y_to_x_start
  real(c_double) :: time_y_to_x_end
  real(c_double) :: time_x_to_y_start
  real(c_double) :: time_x_to_y_end
  real(c_double) :: time_substep_start
  real(c_double) :: time_substep_end
  real(c_double) :: time_data_store_start
  real(c_double) :: time_data_store_end
  real(c_double) :: time_statistics_start
  real(c_double) :: time_statistics_end

  real(c_double) :: time_plan_start
  real(c_double) :: time_plan_end
 
  real(c_double) :: time_writing_field


  integer :: program_start_time(8)
  integer :: program_end_time(8)
  integer :: program_current_time(8)
  integer :: simulation_start_time(8)
  



end module timer
