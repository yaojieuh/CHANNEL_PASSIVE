!-----------------------------------------------------------------------------------------------
! This routine calculate the 1st order statistics and store in HDF5 files
!
! Developer : M. K. Lee
!
! Last update : 2013.07.01
!-----------------------------------------------------------------------------------------------

subroutine get_statistics_1(Udata)

  use channel_hdf5_io
  use global_variable
  use labs
  use ppdfftw
  use runge_kutta

  implicit none
  real(C_DOUBLE) :: Udata(2,ny,Yxsz,Yzsz,no_of_field)
  type(io_handle) :: h
  character(len=255) :: str, str2, time_string !, filename_stat, time_of_saving, date_of_saving
  character(len=9) :: timestep_string
  character(len=4) :: year
  character(len=2) :: month, date, hr, min, sec
  integer, dimension(8) :: date_value
  integer :: i

  if (myrank.eq.0) then
     write(time_string,'(es10.3)') time
     write(timestep_string,'(i9.9)') timestep
     call date_and_time(VALUES=date_value)
     write(year,'(i4)') date_value(1)
     write(month,'(i2)') date_value(2)
     write(date,'(i2)') date_value(3)
     write(hr,'(i2)') date_value(5)
     write(min,'(i2)') date_value(6)

     filename_stat=trim(adjustl(output_file_path))//'stat_step_'//timestep_string//'_time_'//trim(adjustl(time_string))//'.h5'
     date_of_saving=year//"/"//month//"/"//date
     time_of_saving=hr//":"//min

     call io_handle_initialize (h, MPI_COMM_SELF)
     call io_file_create (h,filename_stat, .true.)
     str="This file contains the statistics of simulations. "
     call io_string_set(h,"/","Readme", str)
     call io_string_set(h,"/","Date of Saving (Year/Month/Day) : ", date_of_saving)
     call io_string_set(h,"/","Time of Saving (Hour:Minute) : ", time_of_saving)




     call io_line_establish (h, int(ny), 1, int(ny), ierr)
     call io_line_write_double(h,"LABS_COL",[dble(LABS_Y(:))],1,"Collocation points in y")
     call io_line_write_double(h,"LABS_CENTER_LINE",[dble(LABS_CENTER_LINE(:))],1,&
          "Center Line Coeff: LABS_CENTER_LINE / e.g. U_c = LABS_CENTER_LINE * TRANSPOSE(D-^-1 * u)")
     call io_line_write_double(h,"LABS_INTEGRATE",[dble(LABS_INTEGRATE(:))],1,&
          "Integration Coefficient: LABS_INT / e.g. Int u dy = LABS_INT * TRANSPOSE(D0^-1 * u) = LABS_INT * TRANSPOSE(cu)")
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  !---------------------------------------------------------------------------
  ! Calculate the Mean profile and its derivative
  ! Calculate U
  call get_U(Udata(1,1,1,1,1),stat_U,stat_buf_line_dummy,u_gal)
  call mpi_bcast(stat_U,labs_ny,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if (myrank.eq.0)   call io_line_write_double(h,"U", [DBLE(stat_U(:))],1,"Mean Velocity in X dir")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate W
  call get_U(Udata(1,1,1,1,3),stat_W,stat_buf_line_dummy,0d0)
  call mpi_bcast(stat_W,labs_ny,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if (myrank.eq.0)   call io_line_write_double(h,"W", [DBLE(stat_W(:))],1,"Mean Velocity in Z dir")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate T
  call get_U(Udata(1,1,1,1,7),stat_T,stat_buf_line_dummy,0d0)
  call mpi_bcast(stat_T,labs_ny,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if (myrank.eq.0)   call io_line_write_double(h,"Theta", [DBLE(stat_T(:))],1,"Mean passive scalar")
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  ! Calculate dUdy
  call get_dudy_u_center(Udata(1,1,1,1,1),stat_dUdy,stat_buf_line_dummy1, stat_buf_line_dummy2)
  if (myrank.eq.0)   call io_line_write_double(h,"dUdy", [DBLE(stat_dUdy(:))],1,"Mean Vel Gradient of U (=-<Omega_z>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate dWdy
  call get_dwdy(Udata(1,1,1,1,3),stat_dWdy,stat_buf_line_dummy1, stat_buf_line_dummy2)
  if (myrank.eq.0)   call io_line_write_double(h,"dWdy", [DBLE(stat_dWdy(:))],1,"Mean Vel Gradient of W (=<Omega_x>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate dTdy
  call get_dtdy_t_center(Udata(1,1,1,1,7),stat_dTdy,stat_buf_line_dummy1, stat_buf_line_dummy2)
  if (myrank.eq.0)   call io_line_write_double(h,"dTdy", [DBLE(stat_dTdy(:))],1,"Mean  Gradient of T")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  !---------------------------------------------------------------------------



  !---------------------------------------------------------------------------
  ! Calculate Turbulent Intensity
  ! Calculate <uu>
  call get_uv(Udata(1,1,1,1,1),Udata(1,1,1,1,1),stat_uu,stat_buf_line_dummy)
  if(myrank.eq.0) call io_line_write_double(h,"uu", [DBLE(stat_uu(:))],1,"<uu>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <vv>
  call get_uv(Udata(1,1,1,1,2),Udata(1,1,1,1,2),stat_vv,stat_buf_line_dummy)
  if(myrank.eq.0) call io_line_write_double(h,"vv", [DBLE(stat_vv(:))],1,"<vv>")
  if(myrank.eq.0) call io_line_write_double(h,"P", [DBLE(-stat_vv(:))],1,"Mean Pressure = -<vv>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <ww>
  call get_uv(Udata(1,1,1,1,3),Udata(1,1,1,1,3),stat_ww,stat_buf_line_dummy)
  if(myrank.eq.0) call io_line_write_double(h,"ww", [DBLE(stat_ww(:))],1,"<ww>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <TT>
  call get_uv(Udata(1,1,1,1,7),Udata(1,1,1,1,7),stat_tt,stat_buf_line_dummy)
  if(myrank.eq.0) call io_line_write_double(h,"tt", [DBLE(stat_tt(:))],1,"<tt>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  !---------------------------------------------------------------------------


  !---------------------------------------------------------------------------
  ! Calculate Turbulent Stress
  ! Calculate <uv>
  call get_uv(Udata(1,1,1,1,1),Udata(1,1,1,1,2),stat_uv,stat_buf_line_dummy)
  if(myrank.eq.0) call io_line_write_double(h,"uv", [DBLE(stat_uv(:))],1,"<uv>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <uw>
  call get_uv(Udata(1,1,1,1,1),Udata(1,1,1,1,3),stat_uw,stat_buf_line_dummy)
  if(myrank.eq.0) call io_line_write_double(h,"uw", [DBLE(stat_uw(:))],1,"<uw>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <vw>
  call get_uv(Udata(1,1,1,1,2),Udata(1,1,1,1,3),stat_vw,stat_buf_line_dummy)
  if(myrank.eq.0) call io_line_write_double(h,"vw", [DBLE(stat_vw(:))],1,"<vw>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <ut>
  call get_uv(Udata(1,1,1,1,1),Udata(1,1,1,1,7),stat_ut,stat_buf_line_dummy)
  if(myrank.eq.0) call io_line_write_double(h,"ut", [DBLE(stat_ut(:))],1,"<ut>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <wt>
  call get_uv(Udata(1,1,1,1,7),Udata(1,1,1,1,3),stat_wt,stat_buf_line_dummy)
  if(myrank.eq.0) call io_line_write_double(h,"wt", [DBLE(stat_wt(:))],1,"<wt>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <vt>
  call get_uv(Udata(1,1,1,1,2),Udata(1,1,1,1,7),stat_vt,stat_buf_line_dummy)
  if(myrank.eq.0) call io_line_write_double(h,"vt", [DBLE(stat_vt(:))],1,"<vt>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  !---------------------------------------------------------------------------



!!$  !---------------------------------------------------------------------------
!!$  ! Calculate turbulent production terms for turbulent transport equation
!!$  ! Calculate <uv>*dUdy
!!$  if (myrank.eq.0) then
!!$     do i=1,ny
!!$        stat_buf_line(i) = -2d0*stat_uv(i)*stat_dUdy(i)
!!$     enddo
!!$     call io_line_write_double(h,"Turb_Transport_P_uu", [DBLE(stat_buf_line(:))],1, &
!!$          "-2<uv>dUdy : Production term of Turbulent transport eq for uu")
!!$     ! Calculate <ww>*dWdy
!!$     do i=1,ny
!!$        stat_buf_line(i) = -2d0*stat_vw(i)*stat_dWdy(i)
!!$     enddo
!!$     call io_line_write_double(h,"Turb_Transport_P_ww", [DBLE(stat_buf_line(:))],1, &
!!$          "-2<vw>dWdy : Production term of Turbulent transport eq for ww")
!!$     ! Calculate <uv>*dUdy
!!$     do i=1,ny
!!$        stat_buf_line(i) = -stat_vv(i)*stat_dUdy(i)
!!$     enddo
!!$     call io_line_write_double(h,"Turb_Transport_P_uv", [DBLE(stat_buf_line(:))],1, &
!!$          "<vv>dUdy : Production term of Turbulent transport eq for uv")
!!$  endif
!!$  call mpi_barrier(MPI_COMM_WORLD,ierr)
!!$  !---------------------------------------------------------------------------




!!$  !---------------------------------------------------------------------------
!!$  ! Calculate viscous transport terms for turbulent transport equation
!!$  ! Calculate nu*d2<uu>dy2
!!$  if (myrank.eq.0) then
!!$     call get_turb_transport_D_ij(stat_uu,stat_buf_line, stat_buf_line_dummy1, stat_buf_line_dummy2)
!!$     call io_line_write_double(h,"Turb_Transport_D_uu", [DBLE(stat_buf_line(:))],1, &
!!$          "nu*d^2<uu>/dy^2 : Viscous transport term of Turbulent transport eq for uu")
!!$     ! Calculate nu*d2<uu>dy2
!!$     call get_turb_transport_D_ij(stat_vv,stat_buf_line, stat_buf_line_dummy1, stat_buf_line_dummy2)
!!$     call io_line_write_double(h,"Turb_Transport_D_vv", [DBLE(stat_buf_line(:))],1, &
!!$          "nu*d^2<vv>/dy^2 : Viscous transport term of Turbulent transport eq for vv")
!!$     ! Calculate nu*d2<uu>dy2
!!$     call get_turb_transport_D_ij(stat_ww,stat_buf_line, stat_buf_line_dummy1, stat_buf_line_dummy2)
!!$     call io_line_write_double(h,"Turb_Transport_D_ww", [DBLE(stat_buf_line(:))],1, &
!!$          "nu*d^2<ww>/dy^2 : Viscous transport term of Turbulent transport eq for ww")
!!$     ! Calculate nu*d2<uu>dy2
!!$     if (myrank.eq.0)   call get_turb_transport_D_ij(stat_uv,stat_buf_line, stat_buf_line_dummy1, stat_buf_line_dummy2)
!!$     call io_line_write_double(h,"Turb_Transport_D_uv", [DBLE(stat_buf_line(:))],1, &
!!$          "nu*d^2<uv>/dy^2 : Viscous transport term of Turbulent transport eq for uv")
!!$  endif
!!$  call mpi_barrier(MPI_COMM_WORLD,ierr)
!!$  !---------------------------------------------------------------------------



  !---------------------------------------------------------------------------
  ! Calculate Y-dir derivative of fields and save in stat_buf_field
  ! Calculate dudy
  call get_y_dir_deriv_of_field(Udata(1,1,1,1,1), stat_buf_field(1,1))
  ! Calculate dvdy
  call get_y_dir_deriv_of_field_with_PPE_BC(Udata(1,1,1,1,2), stat_buf_field(1,2), stat_PPE_Boundary(1,1,1))
  ! Calculate dwdy
  call get_y_dir_deriv_of_field(Udata(1,1,1,1,3), stat_buf_field(1,3))
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  !---------------------------------------------------------------------------



  !---------------------------------------------------------------------------
  ! Calculate viscous dissipation terms for turbulent transport equation
  ! Calculate -2*nu*(<dudx dudx> + <dudy dudy> + <dudz dudz>)
  call get_turb_transport_e_ij(Udata(1,1,1,1,1),Udata(1,1,1,1,1),stat_buf_field(1,1),stat_buf_field(1,1), &
       Udata(1,1,1,1,1),Udata(1,1,1,1,1),stat_buf_line,stat_buf_line_dummy)
  if (myrank.eq.0) call io_line_write_double(h,"2nududxidudxi", [DBLE(stat_buf_line(:))],1, &
       "2*nu*(<dudx dudx> + <dudy dudy> + <dudz dudz>)")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate -2*nu*(<dvdx dvdx> + <dvdy dvdy> + <dvdz dvdz>)
  call get_turb_transport_e_ij(Udata(1,1,1,1,2),Udata(1,1,1,1,2),stat_buf_field(1,2),stat_buf_field(1,2), &
       Udata(1,1,1,1,2),Udata(1,1,1,1,2),stat_buf_line,stat_buf_line_dummy)
  if (myrank.eq.0) call io_line_write_double(h,"2nudvdxidvdxi", [DBLE(stat_buf_line(:))],1, &
       "2*nu*(<dvdx dvdx> + <dvdy dvdy> + <dvdz dvdz>)")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate -2*nu*(<dwdx dwdx> + <dwdy dwdy> + <dwdz dwdz>)
  call get_turb_transport_e_ij(Udata(1,1,1,1,3),Udata(1,1,1,1,3),stat_buf_field(1,3),stat_buf_field(1,3), &
       Udata(1,1,1,1,3),Udata(1,1,1,1,3),stat_buf_line,stat_buf_line_dummy)
  if (myrank.eq.0) call io_line_write_double(h,"2nudwdxidwdxi", [DBLE(stat_buf_line(:))],1, &
       "2*nu*(<dwdx dwdx> + <dwdy dwdy> + <dwdz dwdz>)")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate -2*nu*(<dudx dvdx> + <dudy dvdy> + <dudz dvdz>)
  call get_turb_transport_e_ij(Udata(1,1,1,1,1),Udata(1,1,1,1,2),stat_buf_field(1,1),stat_buf_field(1,2), &
       Udata(1,1,1,1,1),Udata(1,1,1,1,2),stat_buf_line,stat_buf_line_dummy)
  if (myrank.eq.0) call io_line_write_double(h,"2nududxidvdxi", [DBLE(stat_buf_line(:))],1, &
       "2*nu*(<dudx dvdx> + <dudy dvdy> + <dudz dvdz>)")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  !---------------------------------------------------------------------------




  !---------------------------------------------------------------------------
  ! Calculate vorticity field
  call get_vorticity_field(Udata(1,1,1,1,1),Udata(1,1,1,1,2),Udata(1,1,1,1,3), &
       stat_buf_field(1,1),stat_buf_field(1,2),stat_buf_field(1,3))
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  !---------------------------------------------------------------------------


  !---------------------------------------------------------------------------
  ! Calculate vorticity variance and covariance
  ! Calculate <omega x omega x>
  call get_uv(stat_buf_field(1,1),stat_buf_field(1,1),stat_buf_line,stat_buf_line_dummy)
  if (myrank.eq.0)   call io_line_write_double(h,"omega_x_omega_x", [DBLE(stat_buf_line(:))],1,"<omega_x omega_x>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <omega y omega y>
  call get_uv(stat_buf_field(1,2),stat_buf_field(1,2),stat_buf_line,stat_buf_line_dummy)
  if (myrank.eq.0)   call io_line_write_double(h,"omega_y_omega_y", [DBLE(stat_buf_line(:))],1,"<omega_y omega_y>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <omega z omega z>
  call get_uv(stat_buf_field(1,3),stat_buf_field(1,3),stat_buf_line,stat_buf_line_dummy)
  if (myrank.eq.0)   call io_line_write_double(h,"omega_z_omega_z", [DBLE(stat_buf_line(:))],1,"<omega_z omega_z>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <omega x omega y>
  call get_uv(stat_buf_field(1,1),stat_buf_field(1,2),stat_buf_line,stat_buf_line_dummy)
  if (myrank.eq.0)   call io_line_write_double(h,"omega_x_omega_y", [DBLE(stat_buf_line(:))],1,"<omega_x omega_y>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <omega y omega z>
  call get_uv(stat_buf_field(1,2),stat_buf_field(1,3),stat_buf_line,stat_buf_line_dummy)
  if (myrank.eq.0)   call io_line_write_double(h,"omega_y_omega_z", [DBLE(stat_buf_line(:))],1,"<omega_y omega_z>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <omega x omega z>
  call get_uv(stat_buf_field(1,1),stat_buf_field(1,3),stat_buf_line,stat_buf_line_dummy)
  if (myrank.eq.0)   call io_line_write_double(h,"omega_x_omega_z", [DBLE(stat_buf_line(:))],1,"<omega_x omega_z>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! calculate uv vortiicty decompostion  added by jieyao
  ! Calculate <v omega z>
  call get_uv(Udata(1,1,1,1,2),stat_buf_field(1,3),stat_buf_line,stat_buf_line_dummy)
  if (myrank.eq.0)   call io_line_write_double(h,"v_omega_z", [DBLE(stat_buf_line(:))],1,"<v omega_z>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <w omega y>
  call get_uv(Udata(1,1,1,1,3),stat_buf_field(1,2),stat_buf_line,stat_buf_line_dummy)
  if (myrank.eq.0)   call io_line_write_double(h,"w_omega_y", [DBLE(stat_buf_line(:))],1,"<w omega_y>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  !---------------------------------------------------------------------------


  !---------------------------------------------------------------------------
  ! Copy vel data to buf
  call move_data_to_buf(Udata(1,1,1,1,1),Udata(1,1,1,1,2),Udata(1,1,1,1,3), &
       stat_buf_field(1,1),stat_buf_field(1,2),stat_buf_field(1,3))
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  !---------------------------------------------------------------------------




  if(myrank.eq.0) then
     call io_line_establish (h, 1, 1, 1, ierr)

     ! Write out step information
     call io_line_write_integer(h,"timestep", [INT(timestep)],0,"Current step in time")
     call io_line_write_double(h,"time",[DBLE(time)],0,"Simulation Time")
     ! Write out system parameters
     call io_line_write_double(h,"nu",[DBLE(nu)],0,"Viscosity or 1/Re")
     call io_line_write_double(h,"Pr",[DBLE(pr)],0,"Prantal number")
     call io_line_write_double(h,"Re",[DBLE(Re)],0,"Reynold Number (not Re_tau)")
     call io_line_write_double(h,"LMFx",[DBLE(LMFx)],0,"Length Multiplication Factor(streamwise direction) : Lx=LMFx*2pi")
     call io_line_write_double(h,"LMFy",[DBLE(LMFy)],0,"Length Multiplication Factor(wall normal direction) : Ly=LMFy*2")
     call io_line_write_double(h,"LMFz",[DBLE(LMFz)],0,"Length Multiplication Factor(spanwise direction) : Lz=LMFz*2pi")
     call io_line_write_double(h,"Lx",[DBLE(Lx)],0,"Length in streamwise direction : Lx")
     call io_line_write_double(h,"Ly",[DBLE(Ly)],0,"Length in wall normal direction : Ly")
     call io_line_write_double(h,"Lz",[DBLE(Lz)],0,"Length in spanwise direction : Lz")
     call io_line_write_double(h,"U_c",[DBLE(u_center)],0,"Center line velocity : U_c")
     call io_line_write_double(h,"T_c",[DBLE(t_center)],0,"Center line passive : T_c")
     call io_line_write_integer(h,"nx",[INT(nx)],0,"Number of x modes")
     call io_line_write_integer(h,"ny",[INT(ny)],0,"Number of y points")
     call io_line_write_integer(h,"nz",[INT(nz)],0,"Number of z modes")
     call io_line_write_integer(h,"LABS_K",[INT(LABS_K)],0,"B-Spline Order")
     call io_line_write_integer(h,"LABS_BW",[INT(LABS_BW)],0,"Bandwidth of Matrix")
     call io_line_write_double(h,"Stretch",[DBLE(stretch)],0,"Stretch Factor for non-uniform grid")
     call io_line_write_double(h,"dPdx",[DBLE(phi_dpdx(4))],0,"Instantaneous pressure gradient")
     call io_line_write_double(h,"dTdx",[DBLE(phi_dtdx(4))],0,"Instantaneous passive gradient")




     ! Write B-Spline Information
     call io_line_establish(h,LABS_NY-LABS_K+1,1,LABS_NY-LABS_K+1,ierr)

     call io_line_write_double(h,"LABS_KNOT",[dble(LABS_Y_KNOT(:))],1,"KNOT points in y")


     call io_plane_establish(h,LABS_NY,1,LABS_NY,LABS_NY,1,LABS_NY,ierr)

     call io_plane_write_double(h, "LABS_D0_FULL", [dble(LABS_D0_FULL(:,:))], 1,0, &
          "Zeroth order operator(FULL): D0 / e.g. u = D0*cu ")
     call io_plane_write_double(h, "LABS_D1_FULL", [dble(LABS_D1_FULL(:,:))], 1,0, &
          "First order operator(FULL): D1 / e.g. dudy = D1*cu = D1*(D0^-1)*u")
     call io_plane_write_double(h, "LABS_D2_FULL", [dble(LABS_D2_FULL(:,:))], 1,0, &
          "Second order operator(FULL): D2 / e.g. d2udy2 = D2*cu = D2*(D0^-1)*u")

     call io_plane_establish(h,LABS_BW*2+1,1,LABS_BW*2+1,LABS_NY,1,LABS_NY,ierr)

     call io_plane_write_double(h, "LABS_D0_BAND", [dble(LABS_D0_B(:,:,1))], 1,0, &
          "Zeroth order operator(BAND): D0 / e.g. u = D0*cu ")
     call io_plane_write_double(h, "LABS_D1_BAND", [dble(LABS_D1_B(:,:,1))], 1,0, &
          "First order operator(BAND): D1 / e.g. dudy = D1*cu = D1*(D0^-1)*u")
     call io_plane_write_double(h, "LABS_D2_BAND", [dble(LABS_D2_B(:,:,1))], 1,0, &
          "Second order operator(BAND): D2 / e.g. d2udy2 = D2*cu = D2*(D0^-1)*u")

  endif




  ! Finalize
  if(myrank.eq.0) call io_file_close(h)
  if(myrank.eq.0) call io_handle_finalize(h)
  call mpi_barrier(MPI_COMM_WORLD,ierr)


  return



end subroutine get_statistics_1







subroutine get_statistics_2(Udata)

  use channel_hdf5_io
  use global_variable
  use labs
  use ppdfftw
  use runge_kutta

  implicit none
  real(C_DOUBLE) :: Udata(2,ny,Yxsz,Yzsz,no_of_field)
  type(io_handle) :: h
  character(len=255) :: str, str2, time_string !, filename_stat, time_of_saving, date_of_saving
  character(len=9) :: timestep_string
  character(len=4) :: year
  character(len=2) :: month, date, hr, min, sec
  integer, dimension(8) :: date_value
  integer :: i


  if (myrank.eq.0) then
!!$     write(time_string,'(es10.3)') time
!!$     write(timestep_string,'(i9.9)') timestep
!!$     filename_stat=trim(adjustl(output_file_path))//'stat_step2_'//timestep_string//'_time_'//trim(adjustl(time_string))//'.h5'
!!$     write(*,*) filename_stat
     call io_handle_initialize (h, MPI_COMM_SELF)
     call io_file_open (h,filename_stat, .true.)
!!$     str="This file contains the statistics of simulations. "
!!$     call io_string_set(h,"/","Readme", str)
!!$     call io_string_set(h,"/","Date of Saving (Year/Month/Day) : ", date_of_saving)
!!$     call io_string_set(h,"/","Time of Saving (Hour:Minute) : ", time_of_saving)



     call io_line_establish (h, int(ny), 1, int(ny), ierr)


     !---------------------------------------------------------------------------
     ! Write skewness data which is calculated from realspace
     call io_line_write_double(h,"uuu", [DBLE(stat_skewness(:,1))],1,"Skewness of u")
     call io_line_write_double(h,"vvv", [DBLE(stat_skewness(:,2))],1,"Skewness of v")
     call io_line_write_double(h,"www", [DBLE(stat_skewness(:,3))],1,"Skewness of w")
     !---------------------------------------------------------------------------



     !---------------------------------------------------------------------------
     ! Write flatness data which is calculated from realspace
     call io_line_write_double(h,"uuuu", [DBLE(stat_flatness(:,1))],1,"Flatness of u")
     call io_line_write_double(h,"vvvv", [DBLE(stat_flatness(:,2))],1,"Flatness of v")
     call io_line_write_double(h,"wwww", [DBLE(stat_flatness(:,3))],1,"Flatness of w")
     !---------------------------------------------------------------------------

     !---------------------------------------------------------------------------
     ! Write high order statistic data which is calculated from realspace
     call io_line_write_double(h,"u6", [DBLE(stat_hs6(:,1))],1,"S6th of u")
     call io_line_write_double(h,"v6", [DBLE(stat_hs6(:,2))],1,"S6th of v")
     call io_line_write_double(h,"w6", [DBLE(stat_hs6(:,3))],1,"S6th of w")
     call io_line_write_double(h,"u8", [DBLE(stat_hs8(:,1))],1,"S8th of u")
     call io_line_write_double(h,"v8", [DBLE(stat_hs8(:,2))],1,"S8th of v")
     call io_line_write_double(h,"w8", [DBLE(stat_hs8(:,3))],1,"S8th of w")
     call io_line_write_double(h,"u10", [DBLE(stat_hs10(:,1))],1,"S10th of u")
     call io_line_write_double(h,"v10", [DBLE(stat_hs10(:,2))],1,"S10th of v")
     call io_line_write_double(h,"w10", [DBLE(stat_hs10(:,3))],1,"S10th of w")
     !---------------------------------------------------------------------------

     !---------------------------------------------------------------------------
     ! Calculate turbulent transport terms for turbulent transport equation
     ! Calculate -d<uuv>dy
     call io_line_write_double(h,"-duuvdy", [DBLE(stat_uuv(:))],1, "-d<uuv>dy")
     ! Calculate -d<vww>dy
     call io_line_write_double(h,"-dvvvdy", [DBLE(stat_vvv(:))],1, "-d<vvv>dy")
     ! Calculate -d<vvv>dy
     call io_line_write_double(h,"-dvwwdy", [DBLE(stat_vww(:))],1, "-d<vww>dy")
     ! Calculate -d<uvv>dy
     call io_line_write_double(h,"-duvvdy", [DBLE(stat_uvv(:))],1, "-d<uvv>dy")
     !---------------------------------------------------------------------------

  endif


  !---------------------------------------------------------------------------
  ! Calculate Pressure field and sttore in Udata(:,:,:,:,6)
  call get_pressure
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  !---------------------------------------------------------------------------


  !---------------------------------------------------------------------------
  ! Calculate Pressure-velocity covariance terms
  ! Calculate Pressure transport terms for turbulent transport equation
  ! Calculate <up>
  call get_uv(stat_buf_field(1,1),Udata(1,1,1,1,6),stat_buf_line,stat_buf_line_dummy)
  if(myrank.eq.0) then
     call io_line_write_double(h,"up", [DBLE(stat_buf_line(:))],1,"<up>")
     ! Calculate -d<up>dy
     do i=1,ny
        stat_buf_line_dummy1(i) = dcmplx(stat_buf_line(i),0d0)
     enddo
     call labs_bs(labs_d0_lu(:,:,1),stat_buf_line_dummy1)
     call labs_mat_mul(labs_d1_b(:,:,1),stat_buf_line_dummy1,stat_buf_line_dummy2)
     do i=1,ny
        stat_buf_line(i) = -real(stat_buf_line_dummy2(i))
     enddo
     call io_line_write_double(h,"-dpudy", [DBLE(stat_buf_line(:))],1,"-<pu>dy")
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <vp>
  call get_uv(stat_buf_field(1,2),Udata(1,1,1,1,6),stat_buf_line,stat_buf_line_dummy)
  if(myrank.eq.0) then
     call io_line_write_double(h,"vp", [DBLE(stat_buf_line(:))],1,"<vp>")
     ! Calculate -2d<vp>dy
     do i=1,ny
        stat_buf_line_dummy1(i) = dcmplx(stat_buf_line(i),0d0)
     enddo
     call labs_bs(labs_d0_lu(:,:,1),stat_buf_line_dummy1)
     call labs_mat_mul(labs_d1_b(:,:,1),stat_buf_line_dummy1,stat_buf_line_dummy2)
     do i=1,ny
        stat_buf_line(i) = -2d0*real(stat_buf_line_dummy2(i))
     enddo
     call io_line_write_double(h,"-2dpvdy", [DBLE(stat_buf_line(:))],1,"-2d<pv>dy")
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <wp>
  call get_uv(stat_buf_field(1,3),Udata(1,1,1,1,6),stat_buf_line,stat_buf_line_dummy)
  if(myrank.eq.0) call io_line_write_double(h,"wp", [DBLE(stat_buf_line(:))],1,"<wp>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! Calculate <pp>
  call get_uv(Udata(1,1,1,1,6),Udata(1,1,1,1,6),stat_buf_line,stat_buf_line_dummy)
  if(myrank.eq.0) call io_line_write_double(h,"pp", [DBLE(stat_buf_line(:))],1,"<pp>")
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  !---------------------------------------------------------------------------





  !---------------------------------------------------------------------------
  ! Calculate Pressure strain terms for turbulent transport equation
  call get_turb_transport_PI1_ij(Udata(1,1,1,1,6),stat_buf_field(1,1),stat_buf_field(1,2),stat_buf_field(1,3),&
       stat_buf_field(1,1),stat_buf_field(1,2),stat_buf_field(1,3), &
       stat_PI1_uu,stat_PI1_vv,stat_PI1_ww,stat_PI1_uv, stat_buf_line_dummy)
  if (myrank.eq.0) then
     call io_line_write_double(h,"2pdudx", [DBLE(stat_PI1_uu(:))],1,"<2p dudx>")
     call io_line_write_double(h,"2pdvdy", [DBLE(stat_PI1_vv(:))],1,"<2p dvdy>")
     call io_line_write_double(h,"2pdwdz", [DBLE(stat_PI1_ww(:))],1,"<2p dwdz>")
     call io_line_write_double(h,"pdudy+pdvdx", [DBLE(stat_PI1_uv(:))],1,"<p dudy> + <p dvdx>")
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  !---------------------------------------------------------------------------




  ! Finalize
  if(myrank.eq.0) call io_file_close(h)
  if(myrank.eq.0) call io_handle_finalize(h)
  call mpi_barrier(MPI_COMM_WORLD,ierr)


  return



end subroutine get_statistics_2





!---------------------------------------------------------------------
! This subroutine calculate <dU/dy> which is -<Omega_z>
!---------------------------------------------------------------------
subroutine get_dudy_u_center(U_data,STAT,dummy1,dummy2)

  use ppdfftw
  use labs

  implicit none

  complex(C_DOUBLE_COMPLEX), dimension(ny) :: U_data, dummy1, dummy2
  real(C_DOUBLE) :: STAT(ny)
  integer :: i,j,k

  STAT(:) = 0d0


  if ((myrank.eq.0).and.((Yzst.eq.0).and.(Yxst.eq.0))) then
     dummy1(:) = U_data(:)
     call labs_bs(LABS_D0_LU(:,:,1),dummy1)
     call labs_mat_mul(LABS_D1_B(:,:,1),dummy1,dummy2)

     u_center = U_gal
     do k=1,labs_ny
        u_center = u_center + real(dummy1(k))*labs_center_line(k)
     enddo

     do k=1,labs_ny
        STAT(k) = real(dummy2(k))
     enddo

  endif



  return
end subroutine get_dudy_u_center

!---------------------------------------------------------------------
! This subroutine calculate <dT/dy>
!---------------------------------------------------------------------
subroutine get_dtdy_t_center(T_data,STAT,dummy1,dummy2)

  use ppdfftw
  use labs

  implicit none

  complex(C_DOUBLE_COMPLEX), dimension(ny) :: T_data, dummy1, dummy2
  real(C_DOUBLE) :: STAT(ny)
  integer :: i,j,k

  STAT(:) = 0d0


  if ((myrank.eq.0).and.((Yzst.eq.0).and.(Yxst.eq.0))) then
     dummy1(:) = T_data(:)
     call labs_bs(LABS_D0_LU(:,:,1),dummy1)
     call labs_mat_mul(LABS_D1_B(:,:,1),dummy1,dummy2)

     !u_center = U_gal
     t_center = 0d0
     do k=1,labs_ny
        t_center =t_center + real(dummy1(k))*labs_center_line(k)
     enddo

     do k=1,labs_ny
        STAT(k) = real(dummy2(k))
     enddo

  endif



  return
end subroutine get_dtdy_t_center


!---------------------------------------------------------------------
! This subroutine calculate <dW/dy> which is -<Omega_x>
!---------------------------------------------------------------------
subroutine get_dwdy(W_data,STAT,dummy1,dummy2)

  use ppdfftw
  use labs

  implicit none

  complex(C_DOUBLE_COMPLEX), dimension(ny) :: W_data, dummy1, dummy2
  real(C_DOUBLE) :: STAT(ny)
  integer :: i,j,k

  STAT(:) = 0d0


  if ((myrank.eq.0).and.((Yzst.eq.0).and.(Yxst.eq.0))) then
     dummy1(:) = W_data(:)
     call labs_bs(LABS_D0_LU(:,:,1),dummy1)
     call labs_mat_mul(LABS_D1_B(:,:,1),dummy1,dummy2)

     do k=1,labs_ny
        STAT(k) = real(dummy2(k))
     enddo

  endif

  return
end subroutine get_dwdy


!---------------------------------------------------------------------
! This subroutine calculate y-dir derivative of a field
!---------------------------------------------------------------------
subroutine get_y_dir_deriv_of_field(Org_data, Tgt_data)
  use ppdfftw
  use labs
  implicit none
  integer :: i,j,k,l
  complex(C_DOUBLE_COMPLEX) :: Org_data(ny,yxsz*yzsz)
  complex(C_DOUBLE_COMPLEX) :: Tgt_data(ny,yxsz*yzsz)
  integer :: n_total
  n_total = yxsz*yzsz

  call omp_set_num_threads(omp_compute)
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,l)
  l = omp_get_thread_num() + 1
  !$OMP DO SCHEDULE(STATIC)
  do i=1,n_total
     buf_line_y_c1(:,l) = Org_data(:,i)
     call labs_bs(LABS_D0_LU(:,:,l), buf_line_y_c1(:,l))
     call labs_mat_mul(LABS_D1_B(:,:,l), buf_line_y_c1(:,l), Tgt_data(:,i))
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL

  return
end subroutine get_y_dir_deriv_of_field



!---------------------------------------------------------------------
! This subroutine calculate y-dir derivative of a field
! It also calculate the boundary condition of pressure poisson eq
!---------------------------------------------------------------------
subroutine get_y_dir_deriv_of_field_with_PPE_BC(Org_data, Tgt_data,PPE_BC)
  use ppdfftw
  use labs
  implicit none
  integer :: i,j,k,l
  complex(C_DOUBLE_COMPLEX) :: Org_data(ny,yxsz*yzsz)
  complex(C_DOUBLE_COMPLEX) :: Tgt_data(ny,yxsz*yzsz)
  complex(C_DOUBLE_COMPLEX) :: PPE_BC(2,yxsz*yzsz)
  integer :: n_total
  n_total = yxsz*yzsz

  call omp_set_num_threads(omp_compute)
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,l)
  l = omp_get_thread_num() + 1
  !$OMP DO SCHEDULE(STATIC)
  do i=1,n_total
     buf_line_y_c1(:,l) = Org_data(:,i)
     call labs_bs(LABS_D0_LU(:,:,l), buf_line_y_c1(:,l))
     call labs_mat_mul(LABS_D1_B(:,:,l), buf_line_y_c1(:,l), Tgt_data(:,i))
     PPE_BC(1,i) = nu*(LABS_D2_B(LABS_BW+1,1,l)*buf_line_y_c1(1,l) &
                     + LABS_D2_B(LABS_BW+2,1,l)*buf_line_y_c1(2,l) &
                     + LABS_D2_B(LABS_BW+3,1,l)*buf_line_y_c1(3,l))
     PPE_BC(2,i) = nu*(LABS_D2_B(LABS_BW-1,ny,l)*buf_line_y_c1(ny-2,l) &
                     + LABS_D2_B(LABS_BW  ,ny,l)*buf_line_y_c1(ny-1,l) &
                     + LABS_D2_B(LABS_BW+1,ny,l)*buf_line_y_c1(ny  ,l))
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL

  return
end subroutine get_y_dir_deriv_of_field_with_PPE_BC




!!$!---------------------------------------------------------------------
!!$! This subroutine calculate Viscous Transport term for turbulent stress eq
!!$!---------------------------------------------------------------------
!!$subroutine get_turb_transport_D_ij(U_data,STAT,dummy1,dummy2)
!!$
!!$  use ppdfftw
!!$  use labs
!!$
!!$  implicit none
!!$
!!$  complex(C_DOUBLE_COMPLEX), dimension(ny) :: dummy1, dummy2
!!$  real(C_DOUBLE) :: U_data(ny), STAT(ny)
!!$  integer :: i,j,k
!!$
!!$  STAT(:) = 0d0
!!$
!!$  do i=1,labs_ny
!!$     dummy1(i) = dcmplx(U_data(i)*nu,0d0)
!!$  enddo
!!$
!!$  call labs_bs(LABS_D0_LU(:,:,1),dummy1)
!!$  call labs_mat_mul(LABS_D2_B(:,:,1),dummy1,dummy2)
!!$
!!$
!!$  do i=1,labs_ny
!!$     STAT(i) = real(dummy2(i))
!!$  enddo
!!$
!!$  return
!!$end subroutine get_turb_transport_D_ij




!---------------------------------------------------------------------
! This subroutine calculate <u_i u_j> type statistics
!---------------------------------------------------------------------
subroutine get_turb_transport_e_ij(U1,U2,U3,U4,U5,U6,STAT,dummy)

  use ppdfftw
  use labs
  implicit none

  real(C_double), dimension(2,ny,yxsz,yzsz) :: U1, U2, U3, U4, U5, U6
  real(C_DOUBLE) :: STAT(ny), dummy(ny)
  integer :: i,j,k,l,size,xsz
  real(c_double) :: kz_wave(Yzsz)
  real(c_double) :: kx_wave(Yxsz)

  dummy(:) = 0d0
  size = yxsz * yzsz
  xsz = yxsz

  do j=1,Yzsz
     kz_wave(j) = LMFzr * kz(j)
  enddo
  do i=1,Yxsz
     kx_wave(i) = LMFxr * kx(i)
  enddo

  do l=1,size
     i=mod(l-1,xsz) + 1
     j = (l-1)/xsz + 1

     if ((Yxst.eq.0) .and. (i.eq.1)) then
!        if ((Yzst.ne.0) .or. (j.ne.1)) then
           do k=1,labs_ny
              dummy(k) = dummy(k) + 2d0*nu*kx_wave(i)*kx_wave(i)*(U1(1,k,i,j)*U2(1,k,i,j) + U1(2,k,i,j)*U2(2,k,i,j))&
                   + 2d0*nu*(U3(1,k,i,j)*U4(1,k,i,j) + U3(2,k,i,j)*U4(2,k,i,j)) &
                   + 2d0*nu*kz_wave(j)*kz_wave(j)*(U5(1,k,i,j)*U6(1,k,i,j) + U5(2,k,i,j)*U6(2,k,i,j))
           enddo
!        endif
     else
        do k=1,labs_ny
           dummy(k) = dummy(k) + 4d0*nu*kx_wave(i)*kx_wave(i)*(U1(1,k,i,j)*U2(1,k,i,j) + U1(2,k,i,j)*U2(2,k,i,j))&
                + 4d0*nu*(U3(1,k,i,j)*U4(1,k,i,j) + U3(2,k,i,j)*U4(2,k,i,j)) &
                + 4d0*nu*kz_wave(j)*kz_wave(j)*(U5(1,k,i,j)*U6(1,k,i,j) + U5(2,k,i,j)*U6(2,k,i,j))
        enddo
     endif
  enddo



  call mpi_reduce(dummy,STAT,labs_ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  return
end subroutine get_turb_transport_e_ij






subroutine get_turb_transport_PI1_ij(P,U1,U2,U3,U1C,U2C,U3C,PI1_uu,PI1_vv,PI1_ww,PI1_uv,dummy)

  use ppdfftw
  use labs
  implicit none

  real(C_double), dimension(2,ny,yxsz,yzsz) :: U1, U2, U3, P
  complex(C_double_complex), dimension(ny,yxsz,yzsz) :: U1C, U2C, U3C
  real(C_DOUBLE) :: PI1_uu(ny), PI1_vv(ny), PI1_ww(ny), PI1_uv(ny), dummy(ny)
  integer :: i,j,k,l,size,xsz
  real(c_double) :: kz_wave(Yzsz)
  real(c_double) :: kx_wave(Yxsz)

  dummy(:) = 0d0
  size = yxsz * yzsz
  xsz = yxsz

  do j=1,Yzsz
     kz_wave(j) = LMFzr * kz(j)
  enddo
  do i=1,Yxsz
     kx_wave(i) = LMFxr * kx(i)
  enddo

  PI1_uu(:) = 0d0
  PI1_vv(:) = 0d0
  PI1_ww(:) = 0d0
  PI1_uv(:) = 0d0

  do l=1,size
     i=mod(l-1,xsz) + 1
     j = (l-1)/xsz + 1

     if ((Yxst.eq.0) .and. (i.eq.1)) then
!        if ((Yzst.ne.0) .or. (j.ne.1)) then
           do k=1,labs_ny
              PI1_uu(k) = PI1_uu(k) + 2d0*kx_wave(i)*(U1(1,k,i,j)*P(2,k,i,j) - U1(2,k,i,j)*P(1,k,i,j))
              PI1_uv(k) = PI1_uv(k) +     kx_wave(i)*(U2(1,k,i,j)*P(2,k,i,j) - U2(2,k,i,j)*P(1,k,i,j))
              PI1_ww(k) = PI1_ww(k) + 2d0*kz_wave(j)*(U3(1,k,i,j)*P(2,k,i,j) - U3(2,k,i,j)*P(1,k,i,j))
           enddo
!        endif
     else
        do k=1,labs_ny
           PI1_uu(k) = PI1_uu(k) + 4d0*kx_wave(i)*(U1(1,k,i,j)*P(2,k,i,j) - U1(2,k,i,j)*P(1,k,i,j))
           PI1_uv(k) = PI1_uv(k) + 2d0*kx_wave(i)*(U2(1,k,i,j)*P(2,k,i,j) - U2(2,k,i,j)*P(1,k,i,j))
           PI1_ww(k) = PI1_ww(k) + 4d0*kz_wave(j)*(U3(1,k,i,j)*P(2,k,i,j) - U3(2,k,i,j)*P(1,k,i,j))
        enddo
     endif
  enddo




  do l=1,size
     i=mod(l-1,xsz) + 1
     j = (l-1)/xsz + 1
     call labs_bs(labs_d0_lu(:,:,1),U1C(:,i,j))
     call labs_mat_mul(labs_d1_B(:,:,1),U1C(:,i,j),U3C(:,i,j))
     call labs_bs(labs_d0_lu(:,:,1),U2C(:,i,j))
     call labs_mat_mul(labs_d1_B(:,:,1),U2C(:,i,j),U1C(:,i,j))
     if ((Yxst.eq.0) .and. (i.eq.1)) then
!        if ((Yzst.ne.0) .or. (j.ne.1)) then
           do k=1,labs_ny
              PI1_vv(k) = PI1_vv(k) + 2d0*(U1(1,k,i,j)*P(1,k,i,j) + U1(2,k,i,j)*P(2,k,i,j))
              PI1_uv(k) = PI1_uv(k) + (U3(1,k,i,j)*P(1,k,i,j) + U3(2,k,i,j)*P(2,k,i,j))
           enddo
!        endif
     else
        do k=1,labs_ny
           PI1_vv(k) = PI1_vv(k) + 4d0*(U1(1,k,i,j)*P(1,k,i,j) + U1(2,k,i,j)*P(2,k,i,j))
           PI1_uv(k) = PI1_uv(k) + 2d0*(U3(1,k,i,j)*P(1,k,i,j) + U3(2,k,i,j)*P(2,k,i,j))
        enddo
     endif
  enddo

  dummy(:) = PI1_uu(:)
  call mpi_reduce(dummy,PI1_uu,labs_ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  dummy(:) = PI1_vv(:)
  call mpi_reduce(dummy,PI1_vv,labs_ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  dummy(:) = PI1_ww(:)
  call mpi_reduce(dummy,PI1_ww,labs_ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  dummy(:) = PI1_uv(:)
  call mpi_reduce(dummy,PI1_uv,labs_ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  return
end subroutine get_turb_transport_PI1_ij






!---------------------------------------------------------------------
! This subroutine calculate vorticity field
!---------------------------------------------------------------------
subroutine get_vorticity_field(U1,U2,U3,vor1,vor2,vor3)

  use ppdfftw
  use labs
  implicit none

  complex(C_DOUBLE_COMPLEX), dimension(ny,Yxsz, Yzsz) :: U1, U2, U3, vor1, vor2, vor3
  real(C_DOUBLE) :: STAT(ny), dummy(ny)
  integer :: i,j,k,l,size,xsz
  real(c_double) :: kz_wave(Yzsz)
  real(c_double) :: kx_wave(Yxsz)
  complex(C_DOUBLE_COMPLEX) :: temp1, temp2, temp3
  dummy(:) = 0d0
  size = yxsz * yzsz
  xsz = yxsz

  do j=1,Yzsz
     kz_wave(j) = LMFzr * kz(j)
  enddo
  do i=1,Yxsz
     kx_wave(i) = LMFxr * kx(i)
  enddo

  do l=1,size
     i=mod(l-1,xsz) + 1
     j = (l-1)/xsz + 1
     do k=1,ny
        temp1 = vor1(k,i,j)
        temp3 = vor3(k,i,j)
        vor1(k,i,j) = temp3 - dcmplx(0d0,kz_wave(j))*U2(k,i,j)
        vor2(k,i,j) = dcmplx(0d0,kz_wave(j))*U1(k,i,j) - dcmplx(0d0,kx_wave(i))*U3(k,i,j)
        vor3(k,i,j) = dcmplx(0d0,kx_wave(i))*U2(k,i,j) - temp1
     enddo
  enddo

  return
end subroutine get_vorticity_field






!---------------------------------------------------------------------
! This subroutine calculate vorticity field
!---------------------------------------------------------------------
subroutine move_data_to_buf(U1,U2,U3,buf_data_1,buf_data_2,buf_data_3)

  use ppdfftw
  use labs
  implicit none

  complex(C_DOUBLE_COMPLEX), dimension(ny*Yxsz*Yzsz) :: U1, U2, U3, buf_data_1, buf_data_2, buf_data_3
  integer :: i,size

  buf_data_1(:) = U1(:)
  buf_data_2(:) = U2(:)
  buf_data_3(:) = U3(:)

  return
end subroutine move_data_to_buf






!---------------------------------------------------------------------
! This subroutine calculate <u_i u_j> type statistics
!---------------------------------------------------------------------
subroutine get_uv(U1,U2,STAT,dummy)

  use ppdfftw
  use labs
  implicit none

  real(C_DOUBLE), dimension(2,ny,Yxsz, Yzsz) :: U1, U2
  real(C_DOUBLE) :: STAT(ny), dummy(ny)
  integer :: i,j,k,l,size,xsz

  dummy(:) = 0d0
  size = yxsz * yzsz
  xsz = yxsz

  do l=1,size
     i=mod(l-1,xsz) + 1
     j = (l-1)/xsz + 1

     if ((Yxst.eq.0) .and. (i.eq.1)) then
!        if ((Yzst.ne.0) .or. (j.ne.1)) then
           do k=1,labs_ny
              dummy(k) = dummy(k) + U1(1,k,i,j)*U2(1,k,i,j) + U1(2,k,i,j)*U2(2,k,i,j)
           enddo
!        endif
     else
        do k=1,labs_ny
           dummy(k) = dummy(k) + 2d0*(U1(1,k,i,j)*U2(1,k,i,j) + U1(2,k,i,j)*U2(2,k,i,j))
        enddo
     endif

  enddo

  call mpi_reduce(dummy,STAT,labs_ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  return
end subroutine get_uv



!---------------------------------------------------------------------
! This subroutine calculate spatially averaged velocity <U> or <W>
!---------------------------------------------------------------------
subroutine get_u(U_data,STAT,dummy,u_galilean)

  use ppdfftw
  use labs
  implicit none

  complex(C_DOUBLE_COMPLEX), dimension(ny,Yxsz, Yzsz) :: U_data
  real(C_DOUBLE) :: STAT(ny), dummy(ny)
  real(C_DOUBLE) :: u_galilean
  integer :: i,j,k

  STAT(:) = 0d0


  if ((myrank.eq.0).and.((Yzst.eq.0).and.(Yxst.eq.0))) then
     do k=1,labs_ny
        STAT(k) = real(U_data(k,1,1)) + u_galilean
     enddo
  endif

  return
end subroutine get_u



