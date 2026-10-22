module output_cr_module
contains
!###################################################
!###################################################
!###################################################
!###################################################
recursive subroutine r_output_cr(pst,input_array,input_size,output_array,output_size)
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

  character(LEN=flen)::filename
  integer::rID
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_CR,pst%iUpper+1,input_size,output_size,input_array)
     call r_output_cr(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     filename=transfer(input_array,filename)
     if(index(filename,'output')==0)then
        call backup_cr(pst%s%r,pst%s%g,pst%s%m,pst%s%mdl,filename)
     else
        call output_cr(pst%s,filename)
     endif
  endif

end subroutine r_output_cr
!###################################################
!###################################################
!###################################################
!###################################################
subroutine output_cr(s,filename)
  use amr_parameters, only: ndim,twotondim, flen
  use cr_parameters, only: ncruvar, ncrgrp
  use ramses_commons, only: ramses_t, open_file, close_file
  implicit none
  type(ramses_t)::s
  character(LEN=flen)::filename
  !-----------------------------------
  ! Output CR data in file
  !-----------------------------------
  integer::ilevel,igrid,ilun
  integer(kind=8),dimension(s%r%levelmin:s%r%nlevelmax)::nskip
  real(kind=4),dimension(1:twotondim,1:ncruvar)::qout
#ifdef DO_CR
  real(kind=8),dimension(1:twotondim,1:ncruvar)::cruold
#endif
  logical::overflow_reported=.false.

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)
  call open_file(s,filename,nskip,ilun)
  overflow_reported=.false.
  do ilevel=r%levelmin,r%nlevelmax
     write(ilun,POS=nskip(ilevel))
     do igrid=m%head(ilevel),m%tail(ilevel)
#ifdef DO_CR
        cruold=m%cruold(:,:,igrid)
        qout=real(cruold,kind=4)
#endif
        if(maxval(dble(abs(qout))).gt.1d31 .and. s%g%myid==1 .and. .not. overflow_reported) then
            print*,'The CR variables have very high values and are overflowing in the outputs'
            overflow_reported=.true.
        endif
        write(ilun)qout
     end do
  end do
  call close_file(s,filename,nskip,ilun)
  end associate

end subroutine output_cr
!###################################################
!###################################################
!###################################################
!###################################################
subroutine backup_cr(r,g,m,mdl,filename)
  use amr_parameters, only: ndim,flen
  use cr_parameters, only: ncruvar
  use amr_commons, only: run_t,global_t,mesh_t
  use mdl_module
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(mdl_t)::mdl
  character(LEN=flen)::filename

  integer::ilevel,igrid,ilun,ierr
  character(LEN=5)::nchar
  character(LEN=flen)::fileloc
  logical::file_exist
  ilun=10+mdl_core(mdl)
  call title(g%myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  inquire(file=fileloc,exist=file_exist)
  if (file_exist) then
     open(unit=ilun,file=fileloc,iostat=ierr)
     close(ilun,status="delete")
  end if
  open(unit=ilun,file=fileloc,access="stream",action="write",form='unformatted')
  write(ilun)ndim
  write(ilun)ncruvar
  write(ilun)r%levelmin
  write(ilun)r%nlevelmax
  do ilevel=r%levelmin,r%nlevelmax
     write(ilun)m%noct(ilevel)
  enddo
  do ilevel=r%levelmin,r%nlevelmax
     do igrid=m%head(ilevel),m%tail(ilevel)
#ifdef DO_CR
        write(ilun)m%cruold(:,:,igrid)
#endif
     end do
  enddo
  close(ilun)
end subroutine backup_cr
!########################################################################
!########################################################################
!########################################################################
!########################################################################
subroutine file_descriptor_cr(r,filename)
  use amr_parameters, only: ndim,flen
  use cr_parameters, only: ncrgrp
  use amr_commons, only: run_t
  implicit none
  type(run_t)::r
  character(LEN=flen)::filename
  character(len=1), dimension(1:3), parameter :: dim_keys = ["x", "y", "z"]

  character(LEN=flen)::fileloc
  integer::igrp,idim,ilun

  if(r%verbose)write(*,*)'Entering file_descriptor_cr'

  ilun=11

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write variable names in backup file
  write(ilun,'("nvar        =",I11)')ncrgrp*ndim
  do igrp = 1, ncrgrp
     do idim = 1, ndim
        write(ilun,'("variable #",I2,": cr_flux_", i0.2, "_", a)')idim+(igrp-1)*ndim, igrp, dim_keys(idim)
     end do
  end do
  close(ilun)

end subroutine file_descriptor_cr

!########################################################################
!########################################################################
!########################################################################
!########################################################################
subroutine output_crinfo(r, g, filename)

! Output cr information into info_cr_XXXXX.txt
!------------------------------------------------------------------------
  use amr_parameters, only: flen
  use amr_commons, only: run_t, global_t
  use cr_parameters, only: ncruvar, ncrgrp
  use constants, only: c_cgs
  implicit none
  type(run_t)::r
  type(global_t)::g
  character(LEN=flen)::filename, fileloc
  integer :: ilun
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
!------------------------------------------------------------------------
  if (r%verbose) write(*,*)'Entering output_crinfo'

  ilun=11

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write run parameters
  write(ilun,'("ncruvar     = ", I11)') ncruvar
  write(ilun,'("ncrgrp      = ", I11)') ncrgrp


  ! Write physical parameters
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  write(ilun,'("cr_c_frac   = ", 100(E15.7))') g%cr_c(r%levelmin:r%nlevelmax)*scale_t/c_cgs
  write(ilun,*)

  ! Write photon group properties
  !call write_group_props(r, .false., ilun)

  close(ilun)

end subroutine output_crInfo

!########################################################################
!########################################################################
!########################################################################
!########################################################################
subroutine write_cr_group_props(r, update, lun)

! Write CR group properties to file or std output.
! lun => File identifier (use 6 for std. output)
!------------------------------------------------------------------------
  use amr_commons, only: run_t
  use cr_parameters, only: ncrgrp
  implicit none
  type(run_t)::r
  logical :: update
  integer :: ip, lun
!------------------------------------------------------------------------
  if (.not. update) then
     write(lun,*) 'CR group properties=------------------------------ '
  else
     write(lun,*) 'CR properties have been changed to--------------- '
  end if
  !write(lun, 901) r%group_L0(:)
  !write(lun, 902) r%group_L1(:)
  do ip = 1, ncrgrp
     !write(lun, 907) ip
     !write(lun, 904) r%group_egy(ip)
     !write(lun, 905) r%group_csn(ip,:)
     !write(lun, 906) r%group_cse(ip,:)
  end do
  write (lun,*) '=-----------------------------------------------------'

!901 format ('  groupL0  [eV]  = ', 20f12.3)
!902 format ('  groupL1  [eV]  = ', 20f12.3)
!903 format ('  spec2group     = ', 20I12)
!904 format ('  egy      [eV]  = ', 1pe12.3)
!905 format ('  csn    [cm^2]  = ', 20(1pe12.3))
!906 format ('  cse    [cm^2]  = ', 20(1pe12.3))
!907 format ('  --=Group', I2)

end subroutine write_cr_group_props

end module output_cr_module
