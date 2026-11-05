!--------------------------------------------------------------------
! This subroutine write field data in HDF5 file
! Uses native parallel HDF5 through channel_hdf5_io
!
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------
subroutine data_store(Udata,RK_DATA)

  use channel_hdf5_io
  use global_variable
  use labs
  use ppdfftw
  use runge_kutta

  implicit none
  real(C_DOUBLE) :: Udata(2,ny,Yxsz,Yzsz,no_of_field)
  real(C_DOUBLE) :: RK_DATA(2,ny,Yxsz,Yzsz,2)
  type(io_handle) :: h
  character(len=255) :: str, str2, time_string!, time_of_saving, date_of_saving
  character(len=9) :: timestep_string
  character(len=4) :: year
  character(len=2) :: month, date, hr, min, sec
  integer, dimension(8) :: date_value


!!$  if(myrank .eq. 0) then
!!$     write(6,*) ""
!!$     write(6,'(a48,i8,f20.10)') " WRITING FIELDS : Creating File at timestep/time", timestep, time
!!$  endif


  ! native HDF5 handler setting
  write(time_string,'(es10.3)') time
  write(timestep_string,'(i9.9)') timestep
  filename=trim(adjustl(output_file_path))//'vel_step_'//timestep_string//'_time_'//trim(adjustl(time_string))//'.h5'
  call io_handle_initialize (h, MPI_COMM_WORLD)
  call io_file_create (h,filename, .true.)


  ! writing file information
  str="This file contains the information necessary to restart inhomogeneous channel flow. " &
       //"It features the instantaneous 3D velocity field and the problem definition."
  call date_and_time(VALUES=date_value)
  write(year,'(i4)')  date_value(1)
  write(month,'(i2)') date_value(2)
  write(date,'(i2)')  date_value(3)
  write(hr,'(i2)')    date_value(5)
  write(min,'(i2)')   date_value(6)
  write(sec,'(i2)')   date_value(7)
  date_of_saving=year//"/"//month//"/"//date
  time_of_saving=hr//":"//min//":"//sec
  call io_string_set(h,"/","Readme", str)
  call io_string_set(h,"/","Date of Saving (Year/Month/Day) : ", date_of_saving)
  call io_string_set(h,"/","Time of Saving (Hour:Minute:Second) : ", time_of_saving)



  ! write scalar values
  if (myrank.eq.0) then
     call io_line_establish (h, 1, 1, 1,ierr)
  else
     call io_line_establish (h, 1, 1, 0,ierr)
  endif
  call io_line_write_integer(h,"timestep", [INT(timestep)],0,"Current step in time")
  call io_line_write_double(h,"time",[DBLE(time)],0,"Simulation Time")
  call io_line_write_double(h,"nu",[DBLE(nu)],0,"Viscosity or 1/Re")
  call io_line_write_double(h,"Pr",[DBLE(pr)],0,"Prantal number")
  call io_line_write_double(h,"U_gal",[DBLE(U_gal)],0,"Galilean transformation velocity")
  call io_line_write_double(h,"U_avg",[DBLE(bulk_vel_u)],0,"Bulk(average) velocity U - used for forcing")
  call io_line_write_double(h,"W_avg",[DBLE(bulk_vel_w)],0,"Bulk(average) velocity W - used for forcing")
  call io_line_write_double(h,"T_avg",[DBLE(bulk_t)],0,"Bulk(average) passive scalar - used for forcing")
  call io_line_write_double(h,"Re",[DBLE(Re)],0,"Reynold Number (not Re_tau)")
  call io_line_write_double(h,"LMFx",[DBLE(LMFx)],0,"Length Multiplication Factor(streamwise direction) : Lx=LMFx*2pi")
  call io_line_write_double(h,"LMFy",[DBLE(LMFy)],0,"Length Multiplication Factor(wall normal direction) : Ly=LMFy*2")
  call io_line_write_double(h,"LMFz",[DBLE(LMFz)],0,"Length Multiplication Factor(spanwise direction) : Lz=LMFz*2pi")
  call io_line_write_double(h,"Lx",[DBLE(Lx)],0,"Length in streamwise direction : Lx")
  call io_line_write_double(h,"Ly",[DBLE(Ly)],0,"Length in wall normal direction : Ly")
  call io_line_write_double(h,"Lz",[DBLE(Lz)],0,"Length in spanwise direction : Lz")
  call io_line_write_double(h,"Stretch",[DBLE(stretch)],0,"Stretch Factor for non-uniform grid")
  call io_line_write_double(h,"dt",[DBLE(dt)],0,"Time advance interval")
  call io_line_write_integer(h,"nx",[INT(nx)],0,"Number of x modes")
  call io_line_write_integer(h,"ny",[INT(ny)],0,"Number of y points")
  call io_line_write_integer(h,"nz",[INT(nz)],0,"Number of z modes")
  call io_line_write_integer(h,"LABS_K",[INT(LABS_K)],0,"B-Spline Order")
  call io_line_write_integer(h,"LABS_BW",[INT(LABS_BW)],0,"Bandwidth of Matrix")

  ! write vector values
  if (myrank.eq.0) then
     call io_line_establish (h, int(ny), 1, int(ny), ierr)
     Udata(1,:,1,1,1) = Udata(1,:,1,1,1) + U_gal
  else
     call io_line_establish (h, int(ny), 1, 0, ierr)
  endif
  call io_line_write_double(h,"U_mean",[DBLE(Udata(1,:,1,1,1))],1,"U average")
  call io_line_write_double(h,"W_mean",[DBLE(Udata(1,:,1,1,3))],1,"W average")
  call io_line_write_double(h,"T_mean",[DBLE(Udata(1,:,1,1,7))],1,"T average")
  call io_line_write_double(h,"LABS_COL",[dble(LABS_Y(:))],1,"Collocation points in y")
  call io_line_write_double(h,"LABS_CENTER_LINE",[dble(LABS_CENTER_LINE(:))],1,&
       "Center Line Coeff: LABS_CENTER_LINE / e.g. U_c = LABS_CENTER_LINE * TRANSPOSE(D-^-1 * u)")
  call io_line_write_double(h,"LABS_INTEGRATE",[dble(LABS_INTEGRATE(:))],1,&
       "Integration Coefficient: LABS_INT / e.g. Int u dy = LABS_INT * TRANSPOSE(D0^-1 * u) = LABS_INT * TRANSPOSE(cu)")




  ! Write out system Field
  call io_field_establish (h, INT(ny), 1, INT(ny), INT(nx)/2, INT(Yxst)+1, INT(Yxsz), INT(nz), INT(Yzst)+1, INT(Yzsz), ierr)
  if (write_u.ne.0)    call io_field_writev_double (h, "u",  Udata(1,1,1,1,1),2)
  if (write_v.ne.0)    call io_field_writev_double (h, "v",  Udata(1,1,1,1,2),2)
  if (write_w.ne.0)    call io_field_writev_double (h, "w",  Udata(1,1,1,1,3),2)
  if (write_g.ne.0)    call io_field_writev_double (h, "g",  RK_data(1,1,1,1,1),2)
  if (write_phi.ne.0)  call io_field_writev_double (h, "phi",  RK_data(1,1,1,1,2),2)
  if (write_t.ne.0)    call io_field_writev_double (h, "theta",  Udata(1,1,1,1,7),2)
  if (myrank.eq.0) Udata(1,:,1,1,1) = Udata(1,:,1,1,1) - U_gal
  call mpi_barrier(MPI_COMM_WORLD,ierr)







  if(myrank.eq.0) then
     call io_line_establish(h,LABS_NY-LABS_K+1,1,LABS_NY-LABS_K+1,ierr)
  else
     call io_line_establish(h,LABS_NY-LABS_K+1,1,0,ierr)
  endif
  call io_line_write_double(h,"LABS_KNOT",[dble(LABS_Y_KNOT(:))],1,"KNOT points in y")


  if (myrank.eq.0) then
     call io_plane_establish(h,LABS_NY,1,LABS_NY,LABS_NY,1,LABS_NY,ierr)
  else
     call io_plane_establish(h,LABS_NY,1,0,LABS_NY,1,0,ierr)
  endif

  call io_plane_write_double(h, "LABS_D0_FULL", LABS_D0_FULL(1,1), 1,0, &
       "Zeroth order operator(FULL): D0 / e.g. u = D0*cu ")
  call io_plane_write_double(h, "LABS_D1_FULL", LABS_D1_FULL(1,1), 1,0, &
       "First order operator(FULL): D1 / e.g. dudy = D1*cu = D1*(D0^-1)*u")
  call io_plane_write_double(h, "LABS_D2_FULL", LABS_D2_FULL(1,1), 1,0, &
       "Second order operator(FULL): D2 / e.g. d2udy2 = D2*cu = D2*(D0^-1)*u")

  if (myrank.eq.0) then
     call io_plane_establish(h,LABS_BW*2+1,1,LABS_BW*2+1,LABS_NY,1,LABS_NY,ierr)
  else
     call io_plane_establish(h,LABS_BW*2+1,1,0,LABS_NY,1,0,ierr)
  endif

  call io_plane_write_double(h, "LABS_D0_BAND", LABS_D0_B(1,1,1), 1,0, &
       "Zeroth order operator(BAND): D0 / e.g. u = D0*cu ")
  call io_plane_write_double(h, "LABS_D1_BAND", LABS_D1_B(1,1,1), 1,0, &
       "First order operator(BAND): D1 / e.g. dudy = D1*cu = D1*(D0^-1)*u")
  call io_plane_write_double(h, "LABS_D2_BAND", LABS_D2_B(1,1,1), 1,0, &
       "Second order operator(BAND): D2 / e.g. d2udy2 = D2*cu = D2*(D0^-1)*u")






  call io_file_close (h)
  call io_handle_finalize (h)

!!$  if(myrank .eq. 0) write(6,*) "WRITING FIELDS : Closing File"

  return

end subroutine data_store

! subroutine to save final data
subroutine data_store_f(Udata,RK_DATA)

  use channel_hdf5_io
  use global_variable
  use labs
  use ppdfftw
  use runge_kutta

  implicit none
  real(C_DOUBLE) :: Udata(2,ny,Yxsz,Yzsz,no_of_field)
  real(C_DOUBLE) :: RK_DATA(2,ny,Yxsz,Yzsz,2)
  type(io_handle) :: h
  character(len=255) :: str, str2, time_string!, time_of_saving, date_of_saving
  character(len=9) :: timestep_string
  character(len=4) :: year
  character(len=2) :: month, date, hr, min, sec
  integer, dimension(8) :: date_value


  ! native HDF5 handler setting
  !write(time_string,'(es10.3)') time
  !write(timestep_string,'(i9.9)') timestep
  filename=trim(adjustl(output_file_path))//'vel_savef.h5'
  call io_handle_initialize (h, MPI_COMM_WORLD)
  call io_file_create (h,filename, .true.)


  ! writing file information
  str="This file contains the information necessary to restart inhomogeneous channel flow. " &
       //"It features the instantaneous 3D velocity field and the problem definition."
  call date_and_time(VALUES=date_value)
  write(year,'(i4)')  date_value(1)
  write(month,'(i2)') date_value(2)
  write(date,'(i2)')  date_value(3)
  write(hr,'(i2)')    date_value(5)
  write(min,'(i2)')   date_value(6)
  write(sec,'(i2)')   date_value(7)
  date_of_saving=year//"/"//month//"/"//date
  time_of_saving=hr//":"//min//":"//sec
  call io_string_set(h,"/","Readme", str)
  call io_string_set(h,"/","Date of Saving (Year/Month/Day) : ", date_of_saving)
  call io_string_set(h,"/","Time of Saving (Hour:Minute:Second) : ", time_of_saving)



  ! write scalar values
  if (myrank.eq.0) then
     call io_line_establish (h, 1, 1, 1,ierr)
  else
     call io_line_establish (h, 1, 1, 0,ierr)
  endif
  call io_line_write_integer(h,"timestep", [INT(timestep)],0,"Current step in time")
  call io_line_write_double(h,"time",[DBLE(time)],0,"Simulation Time")
  call io_line_write_double(h,"nu",[DBLE(nu)],0,"Viscosity or 1/Re")
  call io_line_write_double(h,"Pr",[DBLE(pr)],0,"Prantal number")
  call io_line_write_double(h,"U_gal",[DBLE(U_gal)],0,"Galilean transformation velocity")
  call io_line_write_double(h,"U_avg",[DBLE(bulk_vel_u)],0,"Bulk(average) velocity U - used for forcing")
  call io_line_write_double(h,"W_avg",[DBLE(bulk_vel_w)],0,"Bulk(average) velocity W - used for forcing")
  call io_line_write_double(h,"T_avg",[DBLE(bulk_t)],0,"Bulk(average) passive scalar - used for forcing")
  call io_line_write_double(h,"Re",[DBLE(Re)],0,"Reynold Number (not Re_tau)")
  call io_line_write_double(h,"LMFx",[DBLE(LMFx)],0,"Length Multiplication Factor(streamwise direction) : Lx=LMFx*2pi")
  call io_line_write_double(h,"LMFy",[DBLE(LMFy)],0,"Length Multiplication Factor(wall normal direction) : Ly=LMFy*2")
  call io_line_write_double(h,"LMFz",[DBLE(LMFz)],0,"Length Multiplication Factor(spanwise direction) : Lz=LMFz*2pi")
  call io_line_write_double(h,"Lx",[DBLE(Lx)],0,"Length in streamwise direction : Lx")
  call io_line_write_double(h,"Ly",[DBLE(Ly)],0,"Length in wall normal direction : Ly")
  call io_line_write_double(h,"Lz",[DBLE(Lz)],0,"Length in spanwise direction : Lz")
  call io_line_write_double(h,"Stretch",[DBLE(stretch)],0,"Stretch Factor for non-uniform grid")
  call io_line_write_double(h,"dt",[DBLE(dt)],0,"Time advance interval")
  call io_line_write_integer(h,"nx",[INT(nx)],0,"Number of x modes")
  call io_line_write_integer(h,"ny",[INT(ny)],0,"Number of y points")
  call io_line_write_integer(h,"nz",[INT(nz)],0,"Number of z modes")
  call io_line_write_integer(h,"LABS_K",[INT(LABS_K)],0,"B-Spline Order")
  call io_line_write_integer(h,"LABS_BW",[INT(LABS_BW)],0,"Bandwidth of Matrix")

  ! write vector values
  if (myrank.eq.0) then
     call io_line_establish (h, int(ny), 1, int(ny), ierr)
     Udata(1,:,1,1,1) = Udata(1,:,1,1,1) + U_gal
  else
     call io_line_establish (h, int(ny), 1, 0, ierr)
  endif
  call io_line_write_double(h,"U_mean",[DBLE(Udata(1,:,1,1,1))],1,"U average")
  call io_line_write_double(h,"W_mean",[DBLE(Udata(1,:,1,1,3))],1,"W average")
  call io_line_write_double(h,"T_mean",[DBLE(Udata(1,:,1,1,7))],1,"T average")
  call io_line_write_double(h,"LABS_COL",[dble(LABS_Y(:))],1,"Collocation points in y")
  call io_line_write_double(h,"LABS_CENTER_LINE",[dble(LABS_CENTER_LINE(:))],1,&
       "Center Line Coeff: LABS_CENTER_LINE / e.g. U_c = LABS_CENTER_LINE * TRANSPOSE(D-^-1 * u)")
  call io_line_write_double(h,"LABS_INTEGRATE",[dble(LABS_INTEGRATE(:))],1,&
       "Integration Coefficient: LABS_INT / e.g. Int u dy = LABS_INT * TRANSPOSE(D0^-1 * u) = LABS_INT * TRANSPOSE(cu)")




  ! Write out system Field
  call io_field_establish (h, INT(ny), 1, INT(ny), INT(nx)/2, INT(Yxst)+1, INT(Yxsz), INT(nz), INT(Yzst)+1, INT(Yzsz), ierr)
  if (write_u.ne.0)    call io_field_writev_double (h, "u",  Udata(1,1,1,1,1),2)
  if (write_v.ne.0)    call io_field_writev_double (h, "v",  Udata(1,1,1,1,2),2)
  if (write_w.ne.0)    call io_field_writev_double (h, "w",  Udata(1,1,1,1,3),2)
  if (write_g.ne.0)    call io_field_writev_double (h, "g",  RK_data(1,1,1,1,1),2)
  if (write_phi.ne.0)  call io_field_writev_double (h, "phi",  RK_data(1,1,1,1,2),2)
  if (write_t.ne.0)  call io_field_writev_double (h, "theta",  Udata(1,1,1,1,7),2)
  if (myrank.eq.0) Udata(1,:,1,1,1) = Udata(1,:,1,1,1) - U_gal
  call mpi_barrier(MPI_COMM_WORLD,ierr)







  if(myrank.eq.0) then
     call io_line_establish(h,LABS_NY-LABS_K+1,1,LABS_NY-LABS_K+1,ierr)
  else
     call io_line_establish(h,LABS_NY-LABS_K+1,1,0,ierr)
  endif
  call io_line_write_double(h,"LABS_KNOT",[dble(LABS_Y_KNOT(:))],1,"KNOT points in y")


  if (myrank.eq.0) then
     call io_plane_establish(h,LABS_NY,1,LABS_NY,LABS_NY,1,LABS_NY,ierr)
  else
     call io_plane_establish(h,LABS_NY,1,0,LABS_NY,1,0,ierr)
  endif

  call io_plane_write_double(h, "LABS_D0_FULL", LABS_D0_FULL(1,1), 1,0, &
       "Zeroth order operator(FULL): D0 / e.g. u = D0*cu ")
  call io_plane_write_double(h, "LABS_D1_FULL", LABS_D1_FULL(1,1), 1,0, &
       "First order operator(FULL): D1 / e.g. dudy = D1*cu = D1*(D0^-1)*u")
  call io_plane_write_double(h, "LABS_D2_FULL", LABS_D2_FULL(1,1), 1,0, &
       "Second order operator(FULL): D2 / e.g. d2udy2 = D2*cu = D2*(D0^-1)*u")

  if (myrank.eq.0) then
     call io_plane_establish(h,LABS_BW*2+1,1,LABS_BW*2+1,LABS_NY,1,LABS_NY,ierr)
  else
     call io_plane_establish(h,LABS_BW*2+1,1,0,LABS_NY,1,0,ierr)
  endif

  call io_plane_write_double(h, "LABS_D0_BAND", LABS_D0_B(1,1,1), 1,0, &
       "Zeroth order operator(BAND): D0 / e.g. u = D0*cu ")
  call io_plane_write_double(h, "LABS_D1_BAND", LABS_D1_B(1,1,1), 1,0, &
       "First order operator(BAND): D1 / e.g. dudy = D1*cu = D1*(D0^-1)*u")
  call io_plane_write_double(h, "LABS_D2_BAND", LABS_D2_B(1,1,1), 1,0, &
       "Second order operator(BAND): D2 / e.g. d2udy2 = D2*cu = D2*(D0^-1)*u")






  call io_file_close (h)
  call io_handle_finalize (h)

!!$  if(myrank .eq. 0) write(6,*) "WRITING FIELDS : Closing File"

  return

end subroutine data_store_f


subroutine print_labs_data
  use channel_hdf5_io
  use labs
  use ppdfftw
  implicit none

  type(io_handle) :: h

  call io_handle_initialize(h, MPI_COMM_WORLD)
  call io_file_create(h,trim(adjustl(output_file_path))//'LABS_info.h5',.true.)

  if(myrank.eq.0) then
     call io_line_establish(h,1,1,1,ierr)
  else
     call io_line_establish(h,1,1,0,ierr)
  endif
  call io_line_write_double(h,"stretch",[dble(stretch)],0,"Grid Stretch Factor")
  call io_line_write_double(h,"Ly",[dble(Ly)],0,"Length in Y direction")
  call io_line_write_integer(h,"Ny",[INT(LABS_NY)],0,"Number of grid in Y direction")
  call io_line_write_integer(h,"LABS_K",[INT(LABS_K)],0,"B-Spline Order")
  call io_line_write_integer(h,"LABS_BW",[INT(LABS_BW)],0,"Bandwidth of Matrix")


  if(myrank.eq.0) then
     call io_line_establish(h,LABS_NY,1,LABS_NY,ierr)
  else
     call io_line_establish(h,LABS_NY,1,0,ierr)
  endif
  call io_line_write_double(h,"LABS_COL",[dble(LABS_Y(:))],1,"Collocation points in y")
  call io_line_write_double(h,"LABS_CENTER_LINE",[dble(LABS_CENTER_LINE(:))],1,&
       "Center Line Coeff: LABS_CENTER_LINE / e.g. U_c = LABS_CENTER_LINE * TRANSPOSE(D-^-1 * u)")
  call io_line_write_double(h,"LABS_INTEGRATE",[dble(LABS_INTEGRATE(:))],1,&
       "Integration Coefficient: LABS_INT / e.g. Int u dy = LABS_INT * TRANSPOSE(D0^-1 * u) = LABS_INT * TRANSPOSE(cu)")

  if(myrank.eq.0) then
     call io_line_establish(h,LABS_NY-LABS_K+1,1,LABS_NY-LABS_K+1,ierr)
  else
     call io_line_establish(h,LABS_NY-LABS_K+1,1,0,ierr)
  endif
  call io_line_write_double(h,"LABS_KNOT",[dble(LABS_Y_KNOT(:))],1,"KNOT points in y")


  if (myrank.eq.0) then
     call io_plane_establish(h,LABS_NY,1,LABS_NY,LABS_NY,1,LABS_NY,ierr)
  else
     call io_plane_establish(h,LABS_NY,1,0,LABS_NY,1,0,ierr)
  endif

  call io_plane_write_double(h, "LABS_D0_FULL", LABS_D0_FULL(1,1), 1,0, &
       "Zeroth order operator(FULL): D0 / e.g. u = D0*cu ")
  call io_plane_write_double(h, "LABS_D1_FULL", LABS_D1_FULL(1,1), 1,0, &
       "First order operator(FULL): D1 / e.g. dudy = D1*cu = D1*(D0^-1)*u")
  call io_plane_write_double(h, "LABS_D2_FULL", LABS_D2_FULL(1,1), 1,0, &
       "Second order operator(FULL): D2 / e.g. d2udy2 = D2*cu = D2*(D0^-1)*u")

  if (myrank.eq.0) then
     call io_plane_establish(h,LABS_BW*2+1,1,LABS_BW*2+1,LABS_NY,1,LABS_NY,ierr)
  else
     call io_plane_establish(h,LABS_BW*2+1,1,0,LABS_NY,1,0,ierr)
  endif

  call io_plane_write_double(h, "LABS_D0_BAND", LABS_D0_B(1,1,1), 1,0, &
       "Zeroth order operator(BAND): D0 / e.g. u = D0*cu ")
  call io_plane_write_double(h, "LABS_D1_BAND", LABS_D1_B(1,1,1), 1,0, &
       "First order operator(BAND): D1 / e.g. dudy = D1*cu = D1*(D0^-1)*u")
  call io_plane_write_double(h, "LABS_D2_BAND", LABS_D2_B(1,1,1), 1,0, &
       "Second order operator(BAND): D2 / e.g. d2udy2 = D2*cu = D2*(D0^-1)*u")




  call io_file_close(h)
  call io_handle_finalize(h)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

!!$  deallocate(LABS_D0_FULL)
!!$  deallocate(LABS_D1_FULL)
!!$  deallocate(LABS_D2_FULL)
!!$  deallocate(LABS_D3_FULL)
!!$  deallocate(LABS_D4_FULL)
!!$  deallocate(LABS_Y_KNOT)


  return
end subroutine print_labs_data
