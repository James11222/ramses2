module output_part_module
contains
!#######################################################
!#######################################################
!#######################################################
!#######################################################
recursive subroutine r_output_part(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use amr_parameters, only: flen
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  character(LEN=flen)::filename,filename2
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_PART,pst%iUpper+1,input_size,output_size,input_array)
     call r_output_part(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     filename=transfer(input_array,filename)
     if(index(filename,'output')==0)then !backup dump
        if(pst%s%r%part)then
           filename2=TRIM(filename)//'part.'
           call backup_part(pst%s%r,pst%s%g,pst%s%p,filename2)
        endif
        if(pst%s%r%star)then
           filename2=TRIM(filename)//'star.'
           call backup_part(pst%s%r,pst%s%g,pst%s%star,filename2)
        endif
        if(pst%s%r%sink)then
           filename2=TRIM(filename)//'sink.'
           call backup_part(pst%s%r,pst%s%g,pst%s%sink,filename2)
        endif
#ifdef CORE
        if(pst%s%r%core)then
           filename2=TRIM(filename)//'core.'
           call backup_part(pst%s%r,pst%s%g,pst%s%core,filename2)
        endif
#endif
        if(pst%s%r%tree)then
           filename2=TRIM(filename)//'tree.'
           call backup_part(pst%s%r,pst%s%g,pst%s%tree,filename2)
        endif
        if(pst%s%r%trac)then
           filename2=TRIM(filename)//'trac.'
           call backup_part(pst%s%r,pst%s%g,pst%s%trac,filename2)
        endif
        if(pst%s%r%dust)then
           filename2=TRIM(filename)//'dust.'
           call backup_part(pst%s%r,pst%s%g,pst%s%dust,filename2)
        endif
     else !regular output dump
        !only output dark matter particle data if requested
        if(pst%s%r%part .and. pst%s%r%output_part)then
           filename2=TRIM(filename)//'part.'
           call output_part(pst%s,pst%s%p,filename2)
        endif
        if(pst%s%r%star)then
           filename2=TRIM(filename)//'star.'
           call output_part(pst%s,pst%s%star,filename2)
        endif
        if(pst%s%r%sink)then
           filename2=TRIM(filename)//'sink.'
           call output_part(pst%s,pst%s%sink,filename2)
        endif
#ifdef CORE
        if(pst%s%r%core)then
           filename2=TRIM(filename)//'core.'
           call output_part(pst%s,pst%s%core,filename2)
        endif
#endif
        if(pst%s%r%tree)then
           filename2=TRIM(filename)//'tree.'
           call output_part(pst%s,pst%s%tree,filename2)
        endif
        if(pst%s%r%trac)then
           filename2=TRIM(filename)//'trac.'
           call output_part(pst%s,pst%s%trac,filename2)
        endif
        if(pst%s%r%dust)then
           filename2=TRIM(filename)//'dust.'
           call output_part(pst%s,pst%s%dust,filename2)
        endif
     endif
  endif

end subroutine r_output_part
!#######################################################
!#######################################################
!#######################################################
!#######################################################
subroutine output_part(s,p,filename)
  use amr_parameters, only: ndim,i8b,flen
  use ramses_commons, only: ramses_t,open_part_file,close_part_file
  use pm_commons, only: part_t
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  character(LEN=flen)::filename
  !-----------------------------------
  ! Output part data to file
  !-----------------------------------
  integer::i,idim,ilun,ivar
  integer(kind=8),dimension(1:p%nvaralloc+1)::nskip
  real(kind=4),allocatable,dimension(:)::xsp
  integer(i8b),allocatable,dimension(:)::ii8
  integer,allocatable,dimension(:)::ll

  associate(r=>s%r,g=>s%g)

  call open_part_file(s,p,filename,nskip,ilun)

  allocate(xsp(1:p%npart))

  ! Write position
  do idim=1,ndim
     do i=1,p%npart
        xsp(i)=p%xp(i,idim)
     end do
     write(ilun,POS=nskip(idim))
     write(ilun)xsp
  end do

  ! Write velocity
  do idim=1,ndim
     do i=1,p%npart
        xsp(i)=p%vp(i,idim)
     end do
     write(ilun,POS=nskip(idim+ndim))
     write(ilun)xsp
  end do
  ivar=2*ndim

  ! Write mass
  do i=1,p%npart
     xsp(i)=p%mp(i)
  end do
  ivar=ivar+1
  write(ilun,POS=nskip(ivar))
  write(ilun)xsp

  ! Write potential (optional)
  if(allocated(p%phip))then
     do i=1,p%npart
        xsp(i)=p%phip(i)
     end do
     ivar=ivar+1
     write(ilun,POS=nskip(ivar))
     write(ilun)xsp
  endif

  ! Write metallicity
  if(allocated(p%zp))then
     do i=1,p%npart
        xsp(i)=p%zp(i)
     end do
     ivar=ivar+1
     write(ilun,POS=nskip(ivar))
     write(ilun)xsp
  endif

  ! Write acceleration
  if(allocated(p%fp))then
     do idim=1,ndim
        do i=1,p%npart
           xsp(i)=p%fp(i,idim)
        end do
        ivar=ivar+1
        write(ilun,POS=nskip(ivar))
        write(ilun)xsp
     end do
  endif

  ! Write angular momentum
  if(allocated(p%jp))then
     do idim=1,ndim
        do i=1,p%npart
           xsp(i)=p%jp(i,idim)
        end do
        ivar=ivar+1
        write(ilun,POS=nskip(ivar))
        write(ilun)xsp
     end do
  endif

  ! Write birth time
  if(allocated(p%tp))then
     do i=1,p%npart
        xsp(i)=p%tp(i)
     end do
     ivar=ivar+1
     write(ilun,POS=nskip(ivar))
     write(ilun)xsp
  endif

  ! Write merging time
  if(allocated(p%tm))then
     do i=1,p%npart
        xsp(i)=p%tm(i)
     end do
     ivar=ivar+1
     write(ilun,POS=nskip(ivar))
     write(ilun)xsp
  endif

#ifdef CORE
  ! Write specific internal energy
  if(allocated(p%up))then
     do i=1,p%npart
        xsp(i)=p%up(i)
     end do
     ivar=ivar+1
     write(ilun,POS=nskip(ivar))
     write(ilun)xsp
  endif
#endif

  ! Write size
  if(allocated(p%size))then
     do i=1,p%npart
        xsp(i)=p%size(i)
     end do
     ivar=ivar+1
     write(ilun,POS=nskip(ivar))
     write(ilun)xsp
  endif

  ! Write charge
  if(allocated(p%charge))then
     do i=1,p%npart
        xsp(i)=p%charge(i)
     end do
     ivar=ivar+1
     write(ilun,POS=nskip(ivar))
     write(ilun)xsp
  endif

  deallocate(xsp)

  allocate(ll(1:p%npart))

  ! Write level
  do i=1,p%npart
     ll(i)=p%levelp(i)
  end do
  ivar=ivar+1
  write(ilun,POS=nskip(ivar))
  write(ilun)ll

  deallocate(ll)

  allocate(ii8(1:p%npart))

  ! Write identity
  do i=1,p%npart
     ii8(i)=p%idp(i)
  end do
  ivar=ivar+1
  write(ilun,POS=nskip(ivar))
  write(ilun)ii8

  ! Write merging identity
  if(allocated(p%idm))then
     do i=1,p%npart
        ii8(i)=p%idm(i)
     end do
     ivar=ivar+1
     write(ilun,POS=nskip(ivar))
     write(ilun)ii8
  endif

  ! Write tracking identity
  if(allocated(p%idt))then
     do i=1,p%npart
        ii8(i)=p%idt(i)
     end do
     ivar=ivar+1
     write(ilun,POS=nskip(ivar))
     write(ilun)ii8
  endif

  deallocate(ii8)

  call close_part_file(s,p,filename,nskip,ilun)

  end associate

end subroutine output_part
!#######################################################
!#######################################################
!#######################################################
!#######################################################
subroutine backup_part(r,g,p,filename)
  use amr_parameters, only: ndim, dp, i8b, flen
  use amr_commons, only: run_t, global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  character(LEN=flen)::filename

  integer::i,idim,ilun,ierr
  character(LEN=flen)::fileloc
  character(LEN=5)::nchar
  real(dp),allocatable,dimension(:)::xdp
  integer(i8b),allocatable,dimension(:)::ii8
  integer,allocatable,dimension(:)::ll
  logical::file_exist

  ilun=10

  call title(g%myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  inquire(file=fileloc, exist=file_exist)
  if (file_exist) then
     open(unit=ilun,file=fileloc,iostat=ierr)
     close(ilun,status="delete")
  end if
  open(unit=ilun,file=TRIM(fileloc),access="stream",action="write",form='unformatted')
  rewind(ilun)

  ! Write header
  write(ilun)ndim
  write(ilun)p%npart

  allocate(xdp(1:p%npart))

  ! Write position
  do idim=1,ndim
     do i=1,p%npart
        xdp(i)=p%xp(i,idim)
     end do
     write(ilun)xdp
  end do

  ! Write velocity
  do idim=1,ndim
     do i=1,p%npart
        xdp(i)=p%vp(i,idim)
     end do
     write(ilun)xdp
  end do

  ! Write mass
  do i=1,p%npart
     xdp(i)=p%mp(i)
  end do
  write(ilun)xdp

  ! Write metallicity
  if(allocated(p%zp))then
     do i=1,p%npart
        xdp(i)=p%zp(i)
     end do
     write(ilun)xdp
  endif

  ! Write acceleration
  if(allocated(p%fp))then
     do idim=1,ndim
        do i=1,p%npart
           xdp(i)=p%fp(i,idim)
        end do
        write(ilun)xdp
     end do
  endif

  ! Write angular momentum
  if(allocated(p%jp))then
     do idim=1,ndim
        do i=1,p%npart
           xdp(i)=p%jp(i,idim)
        end do
        write(ilun)xdp
     end do
  endif

  ! Write birth time
  if(allocated(p%tp))then
     do i=1,p%npart
        xdp(i)=p%tp(i)
     end do
     write(ilun)xdp
  endif

  ! Write merging time
  if(allocated(p%tm))then
     do i=1,p%npart
        xdp(i)=p%tm(i)
     end do
     write(ilun)xdp
  endif

#ifdef CORE
  ! Write specific internal energy
  if(allocated(p%up))then
     do i=1,p%npart
        xdp(i)=p%up(i)
     end do
     write(ilun)xdp
  endif
#endif

  ! Write size
  if(allocated(p%size))then
     do i=1,p%npart
        xdp(i)=p%size(i)
     end do
     write(ilun)xdp
  endif

  ! Write charge
  if(allocated(p%charge))then
     do i=1,p%npart
        xdp(i)=p%charge(i)
     end do
     write(ilun)xdp
  endif

  deallocate(xdp)

  allocate(ll(1:p%npart))

  ! Write level
  do i=1,p%npart
     ll(i)=p%levelp(i)
  end do
  write(ilun)ll

  deallocate(ll)

  allocate(ii8(1:p%npart))

  ! Write identity
  do i=1,p%npart
     ii8(i)=p%idp(i)
  end do
  write(ilun)ii8

  ! Write merging identity
  if(allocated(p%idm))then
     do i=1,p%npart
        ii8(i)=p%idm(i)
     end do
     write(ilun)ii8
  endif

  ! Write tracking identity
  if(allocated(p%idt))then
     do i=1,p%npart
        ii8(i)=p%idt(i)
     end do
     write(ilun)ii8
  endif

  deallocate(ii8)

  close(ilun)

end subroutine backup_part
!#######################################################
!#######################################################
!#######################################################
!#######################################################
end module output_part_module
