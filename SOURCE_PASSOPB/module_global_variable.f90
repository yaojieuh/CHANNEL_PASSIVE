!--------------------------------------------------------------------
! This module contains variables to be used in general  
! throughout entire code
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------

module global_variable
  use omp_lib
  use, intrinsic :: iso_c_binding
  implicit none


  integer :: timestep, timestep_code_start, target_timestep, target_code_timestep
  integer :: iostep, iostep_safety, statstep
  integer :: write_u, write_v, write_w, write_g, write_phi, write_t
  real(C_DOUBLE) :: program_run_hours
  real(C_DOUBLE) :: target_stat_time, target_data_time, stat_time_interval, data_time_interval
  real(C_DOUBLE) :: nu, pr, bulk_vel_u, bulk_vel_w, bulk_t, u_gal, u_center, t_center
  real(C_DOUBLE) :: Re, wall_shear_stress, wall_theta_stress, forcing_factor
  integer :: abort_signal, delete_status
  integer :: omp_reorder
  integer :: omp_compute
  integer :: restart_file_exist, restart_option
  real(C_DOUBLE),parameter :: pi = 3.14159265358979323d0
  real(C_DOUBLE) :: Lx,Ly,Lz,dx,dz,dxr,dzr,LMFx,LMFy,LMFz,LMFxr,LMFzr,dy_factor
  character(len=255) :: filename, filename_old, command_line, restart_file_name, output_file_path
  character(len=255) :: filename_stat, time_of_saving, date_of_saving
  real(C_DOUBLE) :: velmax,CFL,CFLx,rescale,dt_factor
  real(C_DOUBLE), allocatable :: dyr(:) 
  logical :: stop_program_or_not, data_store_or_not, stat_store_or_not
end module global_variable
