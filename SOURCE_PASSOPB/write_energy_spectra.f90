!-----------------------------------------------------------------------------------------------
! This routine calculate and write 1D energy (co)spectra
! 
! Developer : M. K. Lee
!
! Last update : 2013.07.01
!-----------------------------------------------------------------------------------------------

subroutine write_spectra(Udata)

  use esio
  use global_variable
  use labs
  use ppdfftw
  use runge_kutta

  implicit none
  real(C_DOUBLE) :: Udata(2,ny,Yxsz,Yzsz,no_of_field)
  type(esio_handle) :: h
  character(len=255) :: str, str2, time_string, filename_spectra!, time_of_saving, date_of_saving
  character(len=9) :: timestep_string
  character(len=4) :: year
  character(len=2) :: month, date, hr, min, sec
  integer, dimension(8) :: date_value

  !---------------------------------------------------------------------------
  ! Calculate Y-dir derivative of fields and save in stat_buf_field
  ! Calculate dudy
  call get_y_dir_deriv_of_field(Udata(1,1,1,1,1), stat_buf_field(1,1))
  ! Calculate dvdy
  !call get_y_dir_deriv_of_field_with_PPE_BC(Udata(1,1,1,1,2), stat_buf_field(1,2), stat_PPE_Boundary(1,1,1))
  ! Calculate dwdy
  call get_y_dir_deriv_of_field(Udata(1,1,1,1,3), stat_buf_field(1,3))
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  ! Calculate vorticity field
  call get_vorticity_field(Udata(1,1,1,1,1),Udata(1,1,1,1,2),Udata(1,1,1,1,3), &
       stat_buf_field(1,1),stat_buf_field(1,2),stat_buf_field(1,3))
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  !---------------------------------------------------------------------------

  write(time_string,'(es10.3)') time
  write(timestep_string,'(i9.9)') timestep
  call date_and_time(VALUES=date_value)
  write(year,'(i4)') date_value(1)
  write(month,'(i2)') date_value(2)
  write(date,'(i2)') date_value(3)
  write(hr,'(i2)') date_value(5)
  write(min,'(i2)') date_value(6)

  filename_spectra=trim(adjustl(output_file_path))//'energy_step_'//timestep_string//'_time_'//trim(adjustl(time_string))//'.h5'
  str="This file contains the energy spectra of simulations. " 
  date_of_saving=year//"/"//month//"/"//date
  time_of_saving=hr//":"//min


  call esio_handle_initialize (h, MPI_COMM_WORLD)
  call esio_file_create (h,filename_spectra, .true.)
  call esio_string_set(h,"/","Readme", str)
  call esio_string_set(h,"/","Date of Saving (Year/Month/Day) : ", date_of_saving)
  call esio_string_set(h,"/","Time of Saving (Hour:Minute) : ", time_of_saving)


  if (myrank.eq.0) then
     call esio_line_establish (h, int(ny), 1, int(ny), ierr)
  else
     call esio_line_establish (h, int(ny), 1, 0, ierr)
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call esio_line_write_double(h,"LABS_COL",[dble(LABS_Y(:))],1,"Collocation points in y")
  call esio_line_write_double(h,"LABS_CENTER_LINE",[dble(LABS_CENTER_LINE(:))],1,&
       "Center Line Coeff: LABS_CENTER_LINE / e.g. U_c = LABS_CENTER_LINE * TRANSPOSE(D-^-1 * u)")
  call esio_line_write_double(h,"LABS_INTEGRATE",[dble(LABS_INTEGRATE(:))],1,&
       "Integration Coefficient: LABS_INT / e.g. Int u dy = LABS_INT * TRANSPOSE(D0^-1 * u) = LABS_INT * TRANSPOSE(cu)")
  call mpi_barrier(MPI_COMM_WORLD,ierr)



  if (myrankB.eq.0) then
     call esio_plane_establish(h,int(ny),1,int(ny),int(nx/2),int(Yxst)+1,int(Yxsz),ierr)
  else
     call esio_plane_establish(h,int(ny),1,0,int(nx/2),int(Yxst)+1,0,ierr)
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call get_E_uu_x(Udata(1,1,1,1,1),Udata(1,1,1,1,1),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Euu_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(Udata(1,1,1,1,2),Udata(1,1,1,1,2),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Evv_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(Udata(1,1,1,1,3),Udata(1,1,1,1,3),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Eww_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(Udata(1,1,1,1,7),Udata(1,1,1,1,7),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Ett_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(Udata(1,1,1,1,1),Udata(1,1,1,1,2),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Euv_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(Udata(1,1,1,1,2),Udata(1,1,1,1,3),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Evw_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(Udata(1,1,1,1,1),Udata(1,1,1,1,3),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Euw_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(Udata(1,1,1,1,1),Udata(1,1,1,1,7),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Eut_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(Udata(1,1,1,1,2),Udata(1,1,1,1,7),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Evt_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(Udata(1,1,1,1,7),Udata(1,1,1,1,3),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Ewt_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)


  ! vortiicty spectrum added by Jie Yao
  call get_E_uu_x(stat_buf_field(1,1),stat_buf_field(1,1),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Eoxox_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(stat_buf_field(1,2),stat_buf_field(1,2),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Eoyoy_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(stat_buf_field(1,3),stat_buf_field(1,3),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Eozoz_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(stat_buf_field(1,1),stat_buf_field(1,2),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Eoxoy_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(stat_buf_field(1,1),stat_buf_field(1,3),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Eoxoz_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(stat_buf_field(1,2),stat_buf_field(1,3),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Eoyoz_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(Udata(1,1,1,1,2),stat_buf_field(1,3),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Evoz_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_x(Udata(1,1,1,1,3),stat_buf_field(1,2),stat_buf_plane_x_dummy(1,1),stat_buf_plane_x(1,1))
  call esio_plane_write_double(h, "Ewoy_kx", stat_buf_plane_x(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)


  if (myrankA.eq.0) then
     call esio_plane_establish(h,int(ny),1,int(ny),int(nz),int(Yzst)+1,int(Yzsz),ierr)
  else
     call esio_plane_establish(h,int(ny),1,0,int(nz),int(Yzst)+1,0,ierr)
  endif
  call get_E_uu_z(Udata(1,1,1,1,1),Udata(1,1,1,1,1),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Euu_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_z(Udata(1,1,1,1,2),Udata(1,1,1,1,2),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Evv_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_z(Udata(1,1,1,1,3),Udata(1,1,1,1,3),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Eww_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_z(Udata(1,1,1,1,7),Udata(1,1,1,1,7),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Ett_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)


  call get_E_uu_z(Udata(1,1,1,1,1),Udata(1,1,1,1,2),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Euv_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_z(Udata(1,1,1,1,2),Udata(1,1,1,1,3),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Evw_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_z(Udata(1,1,1,1,1),Udata(1,1,1,1,3),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Euw_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_z(Udata(1,1,1,1,1),Udata(1,1,1,1,7),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Eut_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_z(Udata(1,1,1,1,2),Udata(1,1,1,1,7),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Evt_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_z(Udata(1,1,1,1,3),Udata(1,1,1,1,7),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Ewt_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)


  ! vortiicty spectrum added by Jie Yao
  call get_E_uu_z(stat_buf_field(1,1),stat_buf_field(1,1),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Eoxox_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_z(stat_buf_field(1,2),stat_buf_field(1,2),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Eoyoy_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_z(stat_buf_field(1,1),stat_buf_field(1,3),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Eozoz_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_z(stat_buf_field(1,1),stat_buf_field(1,2),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Eoxoy_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_z(stat_buf_field(1,1),stat_buf_field(1,3),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Eoxoz_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_z(stat_buf_field(1,2),stat_buf_field(1,3),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Eoyoz_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_z(Udata(1,1,1,1,2),stat_buf_field(1,3),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Evoz_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  call get_E_uu_z(Udata(1,1,1,1,3),stat_buf_field(1,2),stat_buf_plane_z_dummy(1,1),stat_buf_plane_z(1,1))
  call esio_plane_write_double(h, "Ewoy_kz", stat_buf_plane_z(1,1),1)
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  !---------------------------------------------------------------------------
  ! Copy vel data to buf
  call move_data_to_buf(Udata(1,1,1,1,1),Udata(1,1,1,1,2),Udata(1,1,1,1,3), &
       stat_buf_field(1,1),stat_buf_field(1,2),stat_buf_field(1,3))
  call mpi_barrier(MPI_COMM_WORLD,ierr)


  if(myrank.eq.0) then
     call esio_line_establish (h, 1, 1, 1, ierr)
  else
     call esio_line_establish (h, 1, 1, 0, ierr)
  endif

  ! Write out step information
  call esio_line_write_integer(h,"timestep", [INT(timestep)],0,"Current step in time")
  call esio_line_write_double(h,"time",[DBLE(time)],0,"Simulation Time")
  ! Write out system parameters
  call esio_line_write_double(h,"nu",[DBLE(nu)],0,"Viscosity or 1/Re")
  call esio_line_write_double(h,"Pr",[DBLE(nu)],0,"Prantal number")
  call esio_line_write_double(h,"Re",[DBLE(Re)],0,"Reynold Number (not Re_tau)")
  call esio_line_write_double(h,"LMFx",[DBLE(LMFx)],0,"Length Multiplication Factor(streamwise direction) : Lx=LMFx*2pi")
  call esio_line_write_double(h,"LMFy",[DBLE(LMFy)],0,"Length Multiplication Factor(wall normal direction) : Ly=LMFy*2")
  call esio_line_write_double(h,"LMFz",[DBLE(LMFz)],0,"Length Multiplication Factor(spanwise direction) : Lz=LMFz*2pi")
  call esio_line_write_double(h,"Lx",[DBLE(Lx)],0,"Length in streamwise direction : Lx")
  call esio_line_write_double(h,"Ly",[DBLE(Ly)],0,"Length in wall normal direction : Ly")
  call esio_line_write_double(h,"Lz",[DBLE(Lz)],0,"Length in spanwise direction : Lz")
  call esio_line_write_integer(h,"nx",[INT(nx)],0,"Number of x modes")
  call esio_line_write_integer(h,"ny",[INT(ny)],0,"Number of y points")
  call esio_line_write_integer(h,"nz",[INT(nz)],0,"Number of z modes")
  call esio_line_write_integer(h,"LABS_K",[INT(LABS_K)],0,"B-Spline Order")
  call esio_line_write_integer(h,"LABS_BW",[INT(LABS_BW)],0,"Bandwidth of Matrix")
  call esio_line_write_double(h,"Stretch",[DBLE(stretch)],0,"Stretch Factor for non-uniform grid")


  if(myrank.eq.0) then
     call esio_line_establish(h,LABS_NY-LABS_K+1,1,LABS_NY-LABS_K+1,ierr)
  else
     call esio_line_establish(h,LABS_NY-LABS_K+1,1,0,ierr)
  endif
  call esio_line_write_double(h,"LABS_KNOT",[dble(LABS_Y_KNOT(:))],1,"KNOT points in y")


  if (myrank.eq.0) then
     call esio_plane_establish(h,LABS_NY,1,LABS_NY,LABS_NY,1,LABS_NY,ierr)
  else
     call esio_plane_establish(h,LABS_NY,1,0,LABS_NY,1,0,ierr)
  endif

  call esio_plane_write_double(h, "LABS_D0_FULL", LABS_D0_FULL(1,1), 1,0, &
       "Zeroth order operator(FULL): D0 / e.g. u = D0*cu ")
  call esio_plane_write_double(h, "LABS_D1_FULL", LABS_D1_FULL(1,1), 1,0, &
       "First order operator(FULL): D1 / e.g. dudy = D1*cu = D1*(D0^-1)*u")
  call esio_plane_write_double(h, "LABS_D2_FULL", LABS_D2_FULL(1,1), 1,0, &
       "Second order operator(FULL): D2 / e.g. d2udy2 = D2*cu = D2*(D0^-1)*u")

  if (myrank.eq.0) then
     call esio_plane_establish(h,LABS_BW*2+1,1,LABS_BW*2+1,LABS_NY,1,LABS_NY,ierr)
  else
     call esio_plane_establish(h,LABS_BW*2+1,1,0,LABS_NY,1,0,ierr)
  endif

  call esio_plane_write_double(h, "LABS_D0_BAND", LABS_D0_B(1,1,1), 1,0, &
       "Zeroth order operator(BAND): D0 / e.g. u = D0*cu ")
  call esio_plane_write_double(h, "LABS_D1_BAND", LABS_D1_B(1,1,1), 1,0, &
       "First order operator(BAND): D1 / e.g. dudy = D1*cu = D1*(D0^-1)*u")
  call esio_plane_write_double(h, "LABS_D2_BAND", LABS_D2_B(1,1,1), 1,0, &
       "Second order operator(BAND): D2 / e.g. d2udy2 = D2*cu = D2*(D0^-1)*u")


  call esio_file_close(h)
  call esio_handle_finalize(h)
  call mpi_barrier(MPI_COMM_WORLD,ierr)


  return



end subroutine write_spectra


!---------------------------------------------------------------------
! This subroutine calculate 1d energy spectra in x direction
! E_uu(kx) / E_vv(kx) / E_ww(kx)
!---------------------------------------------------------------------
subroutine get_E_uu_x(U_data1,U_data2,STAT_src,STAT_tgt)
  use ppdfftw
  use labs
  complex(C_DOUBLE_COMPLEX), dimension(ny,Yxsz,yzsz) :: U_data1,U_data2
  real(C_DOUBLE), dimension(ny,yxsz) :: STAT_src, STAT_tgt
  integer :: i,j,k,l,size,xsz

  STAT_src(:,:) = 0d0
  STAT_tgt(:,:) = 0d0
  size = yxsz*yzsz
  xsz = yxsz

  do l=1,size
     j=mod(l-1,xsz) + 1
     k = (l-1)/xsz + 1
     if ((Yxst.eq.0) .and. (j.eq.1)) then
!        if ((Yzst.ne.0) .or. (k.ne.1)) then
           do i=1,labs_ny
              STAT_src(i,j) = STAT_src(i,j) &
                   + real(U_data1(i,j,k))*real(U_data2(i,j,k)) + aimag(U_data1(i,j,k))*aimag(U_data2(i,j,k))
           enddo
!        endif
     else
        do i=1,labs_ny
           STAT_src(i,j) = STAT_src(i,j) &
                + 2d0*(real(U_data1(i,j,k))*real(U_data2(i,j,k)) + aimag(U_data1(i,j,k))*aimag(U_data2(i,j,k)))
        enddo
     endif
  enddo

  call mpi_reduce(STAT_src,STAT_tgt,int(ny*yxsz),MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM1DB,ierr)

  return
end subroutine get_E_uu_x


!---------------------------------------------------------------------
! This subroutine calculate 1d energy spectra in z direction
! E_uu(kz) / E_vv(kz) / E_ww(kz)
!---------------------------------------------------------------------
subroutine get_E_uu_z(U_data1,U_data2,STAT_src,STAT_tgt)
  use ppdfftw
  use labs
  complex(C_DOUBLE_COMPLEX), dimension(ny,Yxsz,yzsz) :: U_data1,U_data2
  real(C_DOUBLE), dimension(ny,yzsz) :: STAT_src, STAT_tgt
  integer :: i,j,k,l,size,xsz

  STAT_src(:,:) = 0d0
  STAT_tgt(:,:) = 0d0

  size = yxsz*yzsz
  xsz = yxsz

  do l=1,size
     j=mod(l-1,xsz) + 1
     k = (l-1)/xsz + 1
     if ((Yxst.eq.0) .and. (j.eq.1)) then
!        if ((Yzst.ne.0) .or. (k.ne.1)) then
           do i=1,labs_ny
              STAT_src(i,k) = STAT_src(i,k) &
                   + real(U_data1(i,j,k))*real(U_data2(i,j,k)) + aimag(U_data1(i,j,k))*aimag(U_data2(i,j,k))
           enddo
!        endif
     else
        do i=1,labs_ny
           STAT_src(i,k) = STAT_src(i,k) &
                + 2d0*(real(U_data1(i,j,k))*real(U_data2(i,j,k)) + aimag(U_data1(i,j,k))*aimag(U_data2(i,j,k)))
        enddo
     endif
  enddo

  call mpi_reduce(STAT_src,STAT_tgt,int(ny*yzsz),MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM1DA,ierr)
  return
end subroutine get_E_uu_z
