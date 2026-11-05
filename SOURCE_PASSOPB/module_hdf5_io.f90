! Shared native I/O implementation adapted from the neighboring CHANNEL project.
! Parallel HDF5 I/O for the original channel solver. No external wrapper library is required.
! Fields retain the historical ESIO layout 0: a 3-D dataspace whose
! elements are HDF5 arrays of two doubles (real, imaginary), NOT a 4-D
! scalar dataset. Fortran dimensions (y,kx,kz) appear as (kz,kx,y) in C.
module channel_hdf5_io
  use hdf5
  use iso_c_binding, only: c_double, c_int, c_int64_t, c_ptr, c_loc, c_null_char
  use mpi, only: MPI_COMM_WORLD, MPI_COMM_NULL, MPI_INFO_NULL, MPI_INTEGER, &
       MPI_CHARACTER, MPI_Comm_dup, MPI_Comm_free, MPI_Comm_rank, MPI_Abort, MPI_Bcast
  implicit none
  private
  public :: io_handle, io_handle_initialize, io_handle_finalize, io_file_create, io_file_open, io_file_close
  public :: io_line_establish, io_plane_establish, io_field_establish, io_string_set
  public :: io_line_write_double, io_line_read_double, io_line_write_integer, io_line_read_integer
  public :: io_plane_write_double, io_field_writev_double, io_field_readv_double




  type :: decomposition
     integer(hsize_t) :: global(3)=0, start(3)=0, local(3)=0
     integer :: nd=0
  end type
  type :: io_handle
     private
     integer :: comm=MPI_COMM_NULL
     integer(hid_t) :: file=-1, dxpl=-1
     type(decomposition) :: line, plane, field
  end type
contains
  subroutine check(h, status, context)
    type(io_handle), intent(in) :: h
    integer, intent(in) :: status
    character(*), intent(in) :: context
    integer :: e,rank
    if (status >= 0) return
    rank=-1
    if(h%comm/=MPI_COMM_NULL) call MPI_Comm_rank(h%comm,rank,e)
    write(*,'(a,i0,a)') 'Channel HDF5 error on rank ',rank,': '//trim(context)
    ! Abort WORLD even for a rank-zero statistics handle; peers may be
    ! waiting in a solver collective outside this handle's communicator.
    call MPI_Abort(MPI_COMM_WORLD, 1, e)
    error stop 1
  end subroutine

  subroutine io_handle_initialize(h, comm)
    type(io_handle), intent(out) :: h
    integer, intent(in) :: comm
    integer :: e
    call MPI_Comm_dup(comm,h%comm,e)
    if(e /= 0) call check(h,-1,'MPI_Comm_dup')
    call h5open_f(e)
    call check(h,e,'initialize HDF5')
    call h5pcreate_f(H5P_DATASET_XFER_F,h%dxpl,e)
    call check(h,e,'create transfer properties')
    call h5pset_dxpl_mpio_f(h%dxpl,H5FD_MPIO_COLLECTIVE_F,e)
    call check(h,e,'enable collective I/O')
  end subroutine

  subroutine io_file_create(h, filename, overwrite)
    type(io_handle), intent(inout) :: h
    character(*), intent(in) :: filename
    logical, intent(in) :: overwrite
    integer(hid_t) :: fapl
    integer :: e, flags
    if(h%file >= 0) call check(h,-1,'file already open')
    call file_access(h,fapl)
    flags=H5F_ACC_EXCL_F
    if(overwrite) flags=H5F_ACC_TRUNC_F
    call h5fcreate_f(trim(filename),flags,h%file,e,access_prp=fapl)
    call check(h,e,'create '//trim(filename))
    call h5pclose_f(fapl,e)
    call check(h,e,'close file properties')
    call io_string_set(h,'/','Native_IO','POONGBACK channel native parallel HDF5 v1')
  end subroutine

  subroutine io_file_open(h, filename, readwrite, ierr)
    type(io_handle), intent(inout) :: h
    character(*), intent(in) :: filename
    logical, intent(in) :: readwrite
    integer, intent(out), optional :: ierr
    integer(hid_t) :: fapl
    integer :: e, flags
    if(h%file >= 0) call check(h,-1,'file already open')
    call file_access(h,fapl)
    flags=H5F_ACC_RDONLY_F
    if(readwrite) flags=H5F_ACC_RDWR_F
    call h5fopen_f(trim(filename),flags,h%file,e,access_prp=fapl)
    call check(h,e,'open '//trim(filename))
    call h5pclose_f(fapl,e)
    call check(h,e,'close file properties')
    if(present(ierr)) ierr=0
  end subroutine

  subroutine file_access(h,fapl)
    type(io_handle), intent(in) :: h
    integer(hid_t), intent(out) :: fapl
    integer :: e
    call h5pcreate_f(H5P_FILE_ACCESS_F,fapl,e)
    call check(h,e,'create file properties')
    call h5pset_fapl_mpio_f(fapl,h%comm,MPI_INFO_NULL,e)
    call check(h,e,'set MPI file access')
  end subroutine

  subroutine io_file_close(h)
    type(io_handle), intent(inout) :: h
    integer :: e
    if(h%file < 0) return
    call h5fclose_f(h%file,e)
    call check(h,e,'close file')
    h%file=-1
  end subroutine

  subroutine io_handle_finalize(h)
    type(io_handle), intent(inout) :: h
    integer :: e
    call io_file_close(h)
    if(h%dxpl >= 0) then
       call h5pclose_f(h%dxpl,e)
       call check(h,e,'close transfer properties')
       h%dxpl=-1
    endif
    if(h%comm /= MPI_COMM_NULL) call MPI_Comm_free(h%comm,e)
    ! Do not call h5close_f here: another handle may still be active.
  end subroutine

  subroutine establish(h,d,g,s,n)
    type(io_handle), intent(in) :: h
    type(decomposition), intent(out) :: d
    integer, intent(in) :: g(:),s(:),n(:)
    integer :: k
    k=size(g)
    if(any(g <= 0).or.any(s < 1).or.any(n < 0)) call check(h,-1,'invalid decomposition')
    d%nd=k
    d%global(:k)=int(g,hsize_t)
    d%start(:k)=int(s,hsize_t)-1_hsize_t
    d%local(:k)=int(n,hsize_t)
    if(any(d%start(:k)+d%local(:k)>d%global(:k))) call check(h,-1,'decomposition exceeds dataset')
  end subroutine

  subroutine io_line_establish(h,ag,as,an,ierr)
    type(io_handle), intent(inout) :: h
    integer,intent(in) :: ag,as,an
    integer,intent(out),optional :: ierr
    type(decomposition) :: d
    call establish(h,d,[ag],[as],[an])
    h%line=d
    if(present(ierr)) ierr=0
  end subroutine

  subroutine io_plane_establish(h,ag,as,an,bg,bs,bn,ierr)
    type(io_handle), intent(inout) :: h
    integer,intent(in) :: ag,as,an,bg,bs,bn
    integer,intent(out),optional :: ierr
    type(decomposition) :: d
    call establish(h,d,[ag,bg],[as,bs],[an,bn])
    h%plane=d
    if(present(ierr)) ierr=0
  end subroutine

  subroutine io_field_establish(h,ag,as,an,bg,bs,bn,cg,cs,cn,ierr)
    type(io_handle), intent(inout) :: h
    integer,intent(in) :: ag,as,an,bg,bs,bn,cg,cs,cn
    integer,intent(out),optional :: ierr
    type(decomposition) :: d
    call establish(h,d,[ag,bg,cg],[as,bs,cs],[an,bn,cn])
    h%field=d
    if(present(ierr)) ierr=0
  end subroutine

  subroutine io_string_set(h,path,name,value)
    type(io_handle),intent(in) :: h
    character(*),intent(in) :: path,name,value
    character(:),allocatable :: text
    integer :: e,rank,n
    integer(hid_t) :: obj,space,atype,attr
    integer(hsize_t) :: one(1)=[1_hsize_t]
    logical :: exists
    call MPI_Comm_rank(h%comm,rank,e)
    n=len_trim(value)
    call MPI_Bcast(n,1,MPI_INTEGER,0,h%comm,e)
    allocate(character(n+1)::text)
    if(rank==0) text=trim(value)//c_null_char
    call MPI_Bcast(text,n+1,MPI_CHARACTER,0,h%comm,e)
    call h5oopen_f(h%file,trim(path),obj,e)
    call check(h,e,'open attribute object '//path)
    call h5aexists_f(obj,trim(name),exists,e)
    call check(h,e,'query attribute '//name)
    if(exists) then
       call h5adelete_f(obj,trim(name),e)
       call check(h,e,'replace attribute '//name)
    endif
    call h5screate_f(H5S_SCALAR_F,space,e)
    call check(h,e,'attribute dataspace')
    call h5tcopy_f(H5T_C_S1,atype,e)
    call check(h,e,'attribute datatype')
    call h5tset_size_f(atype,int(n+1,size_t),e)
    call check(h,e,'attribute string size')
    call h5acreate_f(obj,trim(name),atype,space,attr,e)
    call check(h,e,'create attribute '//name)
    call h5awrite_f(attr,atype,text,one,e)
    call check(h,e,'write attribute '//name)
    call h5aclose_f(attr,e)
    call check(h,e,'close attribute')
    call h5tclose_f(atype,e)
    call check(h,e,'close attribute type')
    call h5sclose_f(space,e)
    call check(h,e,'close attribute space')
    call h5oclose_f(obj,e)
    call check(h,e,'close attribute object')
  end subroutine

  subroutine transfer(h,d,name,buffer,base,ncomponents,writing,comment)
    type(io_handle),intent(in) :: h
    type(decomposition),intent(in) :: d
    character(*),intent(in) :: name
    type(c_ptr),intent(in) :: buffer
    integer(hid_t),intent(in) :: base
    integer,intent(in) :: ncomponents
    logical,intent(in) :: writing
    character(*),intent(in),optional :: comment
    integer(hid_t) :: dtype,dset,fspace,mspace,stored,super
    type(c_ptr) :: data_pointer
    integer(hsize_t) :: dims(3),maxdims(3),adims(1),mdims(1)
    integer(size_t) :: bytes,expected_bytes
    integer :: e,nd,cls,expected_class,arank,layout
    logical :: exists
    if(d%nd==0.or.h%file<0) call check(h,-1,'uninitialized I/O for '//name)
    if(ncomponents<1) call check(h,-1,'invalid component count')
    if(ncomponents==1) then
       call h5tcopy_f(base,dtype,e)
    else
       call h5tarray_create_f(base,1,[int(ncomponents,hsize_t)],dtype,e)
    endif
    call check(h,e,'memory datatype '//name)
    call h5lexists_f(h%file,trim(name),exists,e)
    call check(h,e,'query dataset '//name)
    if(.not.exists) then
       if(.not.writing) call check(h,-1,'missing restart dataset '//name)
       call h5screate_simple_f(d%nd,d%global(:d%nd),fspace,e)
       call check(h,e,'create dataspace '//name)
       call h5dcreate_f(h%file,trim(name),dtype,fspace,dset,e)
       call check(h,e,'create dataset '//name)
       call h5sclose_f(fspace,e)
       call check(h,e,'close creation dataspace')
    else
       call h5dopen_f(h%file,trim(name),dset,e)
       call check(h,e,'open dataset '//name)
    endif
    ! Reject unexpected shapes/types before passing an unbounded C pointer.
    call h5dget_space_f(dset,fspace,e)
    call check(h,e,'get dataspace '//name)
    call h5sget_simple_extent_ndims_f(fspace,nd,e)
    call check(h,e,'get dataset rank '//name)
    layout=0
    if(d%nd==3) call legacy_layout(h,dset,d,ncomponents,layout)
    if(layout==2) then
       if(nd/=2) call check(h,-1,'unexpected flattened field rank '//name)
    else
       if(nd/=d%nd) call check(h,-1,'unexpected dataset rank '//name)
    endif
    call h5sget_simple_extent_dims_f(fspace,dims(:nd),maxdims(:nd),e)
    call check(h,e,'get dimensions '//name)
    if(layout==2) then
       if(any(dims(:2)/=[d%global(1),d%global(2)*d%global(3)])) &
            call check(h,-1,'flattened restart dimensions do not match input: '//name)
    else
       if(any(dims(:nd)/=d%global(:nd))) call check(h,-1,'restart dimensions do not match input: '//name)
    endif
    call h5dget_type_f(dset,stored,e)
    call check(h,e,'get datatype '//name)
    call h5tget_class_f(stored,cls,e)
    call check(h,e,'get datatype class '//name)
    if(ncomponents>1) then
       if(cls/=H5T_ARRAY_F) call check(h,-1,'expected array-valued field '//name)
       call h5tget_array_ndims_f(stored,arank,e)
       call check(h,e,'array rank '//name)
       if(arank/=1) call check(h,-1,'unexpected component rank '//name)
       call h5tget_array_dims_f(stored,adims,e)
       call check(h,e,'array dimensions '//name)
       if(adims(1)/=ncomponents) call check(h,-1,'unexpected component count '//name)
       call h5tget_super_f(stored,super,e)
       call check(h,e,'array base type '//name)
    else
       call h5tcopy_f(stored,super,e)
       call check(h,e,'scalar base type '//name)
    endif
    call h5tget_class_f(super,cls,e)
    call check(h,e,'base type class')
    call h5tget_class_f(base,expected_class,e)
    call check(h,e,'expected type class')
    call h5tget_size_f(super,bytes,e)
    call check(h,e,'base type size')
    call h5tget_size_f(base,expected_bytes,e)
    call check(h,e,'expected type size')
    if(cls/=expected_class.or.bytes/=expected_bytes) call check(h,-1,'unexpected numeric type '//name)
    call h5tclose_f(super,e)
    call check(h,e,'close base type')
    call h5tclose_f(stored,e)
    call check(h,e,'close stored type')
    mdims(1)=max(1_hsize_t,product(d%local(:d%nd)))
    call h5screate_simple_f(1,mdims,mspace,e)
    call check(h,e,'create memory space '//name)
    if(any(d%local(:d%nd)==0)) then
       call h5sselect_none_f(mspace,e)
       call check(h,e,'empty memory selection')
       call h5sselect_none_f(fspace,e)
       call check(h,e,'empty file selection')
    else if(layout==2) then
       ! Flattened (kz,kx) dimension: select the local kx block for
       ! each local kz without gathering or transposing the field.
       call h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F, &
            [d%start(1),d%start(2)+d%start(3)*d%global(2)], &
            [d%local(1),d%local(3)],e, &
            stride=[1_hsize_t,d%global(2)],block=[1_hsize_t,d%local(2)])
       call check(h,e,'select flattened field '//name)
    else
       call h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,d%start(:nd),d%local(:nd),e)
       call check(h,e,'select file region '//name)
    endif
    data_pointer=buffer
    if(writing) then
       call h5dwrite_f(dset,dtype,buffer,e,mem_space_id=mspace,file_space_id=fspace,xfer_prp=h%dxpl)
    else
       call h5dread_f(dset,dtype,data_pointer,e,mem_space_id=mspace,file_space_id=fspace,xfer_prp=h%dxpl)
    endif
    call check(h,e,'transfer '//name)
    if(writing.and.present(comment)) then
       if(len_trim(comment)>0) then
          call h5oset_comment_f(dset,trim(comment),e)
          call check(h,e,'dataset comment '//name)
       endif
    endif
    call h5sclose_f(mspace,e)
    call check(h,e,'close memory space')
    call h5sclose_f(fspace,e)
    call check(h,e,'close file space')
    call h5dclose_f(dset,e)
    call check(h,e,'close dataset')
    call h5tclose_f(dtype,e)
    call check(h,e,'close memory type')
  end subroutine

  subroutine legacy_layout(h,dset,d,ncomponents,layout)
    type(io_handle),intent(in) :: h
    integer(hid_t),intent(in) :: dset
    type(decomposition),intent(in) :: d
    integer,intent(in) :: ncomponents
    integer,intent(out) :: layout
    integer(hid_t) :: attr,space
    integer(hsize_t) :: dims(1),maxdims(1)
    integer :: e,nd,meta(8)
    logical :: exists
    layout=0
    call h5aexists_f(dset,'esio_field_metadata',exists,e)
    call check(h,e,'query legacy field layout')
    if(.not.exists) return
    call h5aopen_f(dset,'esio_field_metadata',attr,e)
    call check(h,e,'open legacy field metadata')
    call h5aget_space_f(attr,space,e)
    call check(h,e,'legacy metadata space')
    call h5sget_simple_extent_ndims_f(space,nd,e)
    call check(h,e,'legacy metadata rank')
    if(nd/=1) call check(h,-1,'malformed legacy metadata rank')
    call h5sget_simple_extent_dims_f(space,dims,maxdims,e)
    call check(h,e,'legacy metadata dimensions')
    if(dims(1)/=8) call check(h,-1,'malformed legacy metadata length')
    call h5aread_f(attr,H5T_NATIVE_INTEGER,meta,dims,e)
    call check(h,e,'read legacy metadata')
    layout=meta(4)
    if(layout<0.or.layout>2) call check(h,-1,'unsupported legacy field layout')
    if(any(int(meta(7:5:-1),hsize_t)/=d%global).or.meta(8)/=ncomponents) &
         call check(h,-1,'legacy field metadata does not match input')
    call h5sclose_f(space,e)
    call check(h,e,'close legacy metadata space')
    call h5aclose_f(attr,e)
    call check(h,e,'close legacy metadata')
  end subroutine

  subroutine contiguous(h,stride)
    type(io_handle),intent(in) :: h
    integer,intent(in) :: stride
    if(stride/=0.and.stride/=1) call check(h,-1,'Channel I/O expects contiguous buffers')
  end subroutine
  subroutine io_line_write_double(h,name,buffer,stride,comment)
    type(io_handle),intent(in) :: h
    character(*),intent(in) :: name
    real(c_double),target,intent(in) :: buffer(*)
    integer,intent(in) :: stride
    character(*),intent(in),optional :: comment
    call contiguous(h,stride)
    call transfer(h,h%line,name,c_loc(buffer(1)),H5T_NATIVE_DOUBLE,1,.true.,comment)
  end subroutine

  subroutine io_line_read_double(h,name,buffer,stride)
    type(io_handle),intent(in) :: h
    character(*),intent(in) :: name
    real(c_double),target,intent(inout) :: buffer(*)
    integer,intent(in) :: stride

    call contiguous(h,stride)
    call transfer(h,h%line,name,c_loc(buffer(1)),H5T_NATIVE_DOUBLE,1,.false.)
  end subroutine

  subroutine io_line_write_integer(h,name,buffer,stride,comment)
    type(io_handle),intent(in) :: h
    character(*),intent(in) :: name
    integer(c_int),target,intent(in) :: buffer(*)
    integer,intent(in) :: stride
    character(*),intent(in),optional :: comment
    call contiguous(h,stride)
    call transfer(h,h%line,name,c_loc(buffer(1)),H5T_NATIVE_INTEGER,1,.true.,comment)
  end subroutine

  subroutine io_line_read_integer(h,name,buffer,stride)
    type(io_handle),intent(in) :: h
    character(*),intent(in) :: name
    integer(c_int),target,intent(inout) :: buffer(*)
    integer,intent(in) :: stride

    call contiguous(h,stride)
    call transfer(h,h%line,name,c_loc(buffer(1)),H5T_NATIVE_INTEGER,1,.false.)
  end subroutine

  subroutine io_plane_write_double(h,name,buffer,astride,bstride,comment)
    type(io_handle),intent(in) :: h
    character(*),intent(in) :: name
    real(c_double),target,intent(in) :: buffer(*)
    integer,intent(in) :: astride
    integer,intent(in),optional :: bstride
    character(*),intent(in),optional :: comment
    call contiguous(h,astride)
    if(present(bstride)) then
       if(bstride/=0.and.int(bstride,hsize_t)/=h%plane%local(1)) &
            call check(h,-1,'Channel plane buffer must be contiguous')
    endif
    call transfer(h,h%plane,name,c_loc(buffer(1)),H5T_NATIVE_DOUBLE,1,.true.,comment)
  end subroutine

  subroutine io_field_writev_double(h,name,buffer,ncomponents)
    type(io_handle),intent(in) :: h
    character(*),intent(in) :: name
    real(c_double),target,intent(in) :: buffer(*)
    integer,intent(in) :: ncomponents
    call transfer(h,h%field,name,c_loc(buffer(1)),H5T_NATIVE_DOUBLE,ncomponents,.true.)
  end subroutine

  subroutine io_field_readv_double(h,name,buffer,ncomponents)
    type(io_handle),intent(in) :: h
    character(*),intent(in) :: name
    real(c_double),target,intent(inout) :: buffer(*)
    integer,intent(in) :: ncomponents
    call transfer(h,h%field,name,c_loc(buffer(1)),H5T_NATIVE_DOUBLE,ncomponents,.false.)
  end subroutine
end module channel_hdf5_io
