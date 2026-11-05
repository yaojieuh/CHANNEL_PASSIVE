! Native read/write verification; legacy fixtures were produced by real ESIO.
program io_roundtrip
  use channel_hdf5_io
  use mpi
  use iso_c_binding, only: c_double
  implicit none
  type(io_handle) :: h
  integer,parameter :: ny=3,nx=5,nz=4
  integer :: rank,np,e,xs,xn,zs,zn,i,j,k,c,p,f,iv(1),failures,total
  real(c_double),allocatable :: a(:),plane(:)
  real(c_double) :: profile(ny),scalar(1),want
  character(256) :: mode,filename,legacy
  character(5),parameter :: names(6)=[character(5)::'u','v','w','g','phi','theta']
  call MPI_Init(e)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,e)
  call MPI_Comm_size(MPI_COMM_WORLD,np,e)
  call get_command_argument(1,mode)
  call get_command_argument(2,filename)
  call get_command_argument(3,legacy)
  failures=0
  if(trim(mode)=='write') then
     xs=rank*nx/np; xn=(rank+1)*nx/np-xs
     zs=0; zn=nz
  else
     ! Read with a different axis decomposition, including empty ranks.
     xs=0; xn=nx
     zs=rank*nz/np; zn=(rank+1)*nz/np-zs
  endif
  allocate(a(max(1,2*ny*xn*zn)),plane(max(1,ny*xn)))
  profile=[11d0,22d0,33d0]
  call io_handle_initialize(h,MPI_COMM_WORLD)
  if(trim(mode)=='write') then
     call io_file_create(h,trim(filename),.true.)
     call io_string_set(h,'/','Readme','Channel I/O compatibility test')
     call io_line_establish(h,1,1,merge(1,0,rank==0),e)
     call io_line_write_integer(h,'timestep',[17],0,'Current step in time')
     call io_line_write_double(h,'time',[0.25d0],0,'Simulation time')
     call io_line_write_double(h,'dt',[0.125d0],0)
     call io_line_establish(h,ny,1,merge(ny,0,rank==0),e)
     call io_line_write_double(h,'U_mean',profile,1)
     call io_line_write_double(h,'W_mean',profile,1)
     call io_line_write_double(h,'T_mean',profile,1)
     call io_plane_establish(h,ny,1,ny,nx,xs+1,xn,e)
     do j=1,xn
        do i=1,ny
           plane(i+ny*(j-1))=i+100*(xs+j)
        enddo
     enddo
     call io_plane_write_double(h,'Euu_kx',plane,1,0,'Spectrum plane')
  else
     call io_file_open(h,trim(filename),.false.,e)
     ! All ranks read overlapping scalar metadata, as in the solver.
     call io_line_establish(h,1,1,1,e)
     call io_line_read_integer(h,'timestep',iv,0)
     if(iv(1)/=17) failures=failures+1
     call io_line_read_double(h,'time',scalar,0)
     if(scalar(1)/=0.25d0) failures=failures+1
     call io_line_read_double(h,'dt',scalar,0)
     if(scalar(1)/=0.125d0) failures=failures+1
     call io_line_establish(h,ny,1,merge(ny,0,rank==0),e)
     call io_line_read_double(h,'U_mean',profile,1)
     if(rank==0.and.any(profile/=[11d0,22d0,33d0])) failures=failures+1
     profile=0
     call io_line_read_double(h,'W_mean',profile,1)
     if(rank==0.and.any(profile/=[11d0,22d0,33d0])) failures=failures+1
  endif
  if(trim(mode)=='read'.and.trim(legacy)/='legacy') then
     profile=0
     call io_line_read_double(h,'T_mean',profile,1)
     if(rank==0.and.any(profile/=[11d0,22d0,33d0])) failures=failures+1
  endif
  call io_field_establish(h,ny,1,ny,nx,xs+1,xn,nz,zs+1,zn,e)
  do f=1,merge(5,6,trim(legacy)=="legacy")
     a=-999d0
     if(trim(mode)=='write') then
        p=0
        do k=1,zn
           do j=1,xn
              do i=1,ny
                 do c=1,2
                    p=p+1
                    a(p)=c+10*i+100*(xs+j)+1000*(zs+k)+10000*f
                 enddo
              enddo
           enddo
        enddo
        call io_field_writev_double(h,trim(names(f)),a,2)
     else
        call io_field_readv_double(h,trim(names(f)),a,2)
        p=0
        do k=1,zn
           do j=1,xn
              do i=1,ny
                 do c=1,2
                    p=p+1
                    want=c+10*i+100*(xs+j)+1000*(zs+k)+10000*f
                    if(a(p)/=want) failures=failures+1
                 enddo
              enddo
           enddo
        enddo
     endif
  enddo
  call io_file_close(h)
  call io_handle_finalize(h)
  ! Rank-zero-only statistics: reopen for update without truncation.
  if(trim(mode)=='write'.and.rank==0) then
     call io_handle_initialize(h,MPI_COMM_SELF)
     call io_file_open(h,trim(filename),.true.,e)
     call io_line_establish(h,1,1,1,e)
     call io_line_write_double(h,'appended',[42d0],0)
     call io_line_write_double(h,'appended',[43d0],0)
     call io_line_read_double(h,'time',scalar,0)
     if(scalar(1)/=0.25d0) failures=failures+1
     call io_file_close(h)
     call io_handle_finalize(h)
  endif
  call MPI_Allreduce(failures,total,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,e)
  if(total/=0) then
     if(rank==0) write(*,*) 'FAIL: mismatched values',total
     call MPI_Abort(MPI_COMM_WORLD,1,e)
  endif
  if(rank==0) write(*,*) 'PASS ',trim(mode),' ',trim(filename),' ranks=',np
  call MPI_Finalize(e)
end program
