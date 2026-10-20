module input_part_gadget_module

  type :: out_input_gadget_t
     real(kind=8)::mass_star
     real(kind=8)::mass_gas
     real(kind=8)::mass_halo
     integer(kind=8)::nstar
     integer(kind=8)::ngas
     integer(kind=8)::nhalo
  end type out_input_gadget_t

contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_part_gadget(pst)
  use mdl_module
  use amr_parameters, only: ndim
  use ramses_commons, only: pst_t
  use output_amr_module, only: input_header
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch particles
  ! from a Gadget file generated for example by DICE.
  !--------------------------------------------------------------------
  integer::dummy
  type(out_input_gadget_t)::output

  if(pst%s%r%verbose)write(*,*)'Entering input_part_gadget'

  ! Call recursive slave routine
  call r_input_part_gadget(pst,dummy,1,output,12)

  write(*,*)'Total mass in gas=',output%mass_gas
  write(*,*)'Total mass in dark matter=',output%mass_halo
  write(*,*)'Total mass in stars=',output%mass_star

  write(*,*)'Total number of particles in gas=',output%ngas
  write(*,*)'Total number of particles in dark matter=',output%nhalo
  write(*,*)'Total number of particles in stars=',output%nstar

  if(pst%s%r%part)pst%s%p%npart_tot=output%nhalo
  if(pst%s%r%hydro)pst%s%gas%npart_tot=output%ngas
  if(pst%s%r%star)pst%s%star%npart_tot=output%nstar

  pst%s%g%mass_star_tot=output%mass_star

end subroutine m_input_part_gadget
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_input_part_gadget(pst,dummy,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::dummy
  type(out_input_gadget_t)::output,next_output
  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to read and dispatch
  ! particles from a Gadget file.
  !--------------------------------------------------------------------
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INPUT_PART_GADGET,pst%iUpper+1,input_size,output_size,dummy)
     call r_input_part_gadget(pst%pLower,dummy,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass_star=output%mass_star+next_output%mass_star
     output%mass_gas=output%mass_gas+next_output%mass_gas
     output%mass_halo=output%mass_halo+next_output%mass_halo
     output%nstar=output%nstar+next_output%nstar
     output%ngas=output%ngas+next_output%ngas
     output%nhalo=output%nhalo+next_output%nhalo
  else
     call input_part_gadget(pst%s,output%mass_star,output%mass_gas,output%mass_halo,output%nstar,output%ngas,output%nhalo)
  endif

end subroutine r_input_part_gadget
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_part_gadget(s,mstar,mgas,mhalo,npart_star,npart_gas,npart_halo)
  use mdl_module, only: mdl_abort
  use amr_parameters, only: ndim,flen
  use ramses_commons, only: ramses_t
  use gadgetreadfilemod, only: gadgetheadertype
  use init_part_module, only: allocate_gas
  implicit none
  type(ramses_t)::s
  real(kind=8)::mstar,mgas,mhalo
  integer(kind=8)::npart_star,npart_gas,npart_halo
  !------------------------------------------------------------
  ! Read particles positions and velocities from a Ramses
  ! restart file and allocate particle-based arrays.
  !------------------------------------------------------------
  integer::ipart,i,j,k,idim,nrest,id
  integer::npart,nhalo,ngas,nstar,ihalo,igas,istar
  integer::type_index,kpart
  integer::dummy_int,blck_size,jump_blck,blck_cnt,stat
  integer::head_blck,pos_blck,vel_blck,id_blck,mass_blck,u_blck,metal_blck,age_blck
  integer::head_size,pos_size,vel_size,id_size,mass_size,u_size,metal_size,age_size
  type(gadgetheadertype)::header
  character(LEN=flen)::filename
  character(LEN=4)::blck_name
  character(LEN=12)::type_str
  real(kind=4)::dummy_real,x(1:3),v(1:3),z,ipbar
  logical::file_exist,skip(6)

  associate(r=>s%r,g=>s%g,m=>s%m,p=>s%p,gas=>s%gas,mdl=>s%mdl,star=>s%star)

  ! Reading header of the Gadget file
  ! Check the length of the joined initfile/ic_file path
  if(len_trim(r%initfile(r%levelmin))+1+len_trim(r%ic_file)>flen)then
     write(*,'(A,I0,A)')'initfile/ic_file path is longer than ',flen, &
          & ' characters; increase flen in amr_parameters.f90'
     call mdl_abort(mdl)
  endif
  filename = TRIM(r%initfile(r%levelmin))//'/'//TRIM(r%ic_file)
  INQUIRE(FILE=filename,EXIST=file_exist)
  if(.not.file_exist) then
     write(*,*) TRIM(filename)," not found"
     call mdl_abort(mdl)
  endif

  if((r%ic_format.ne.'Gadget1').and.(r%ic_format.ne.'Gadget2')) then
     write(*,*) 'Specify a valid IC file format [ic_format=Gadget1/Gadget2]'
     call mdl_abort(mdl)
  endif

  if(g%myid==1)write(*,'(A50)')"__________________________________________________"
  if(g%myid==1)write(*,'(A12,A)') " Opening -> ",filename
  OPEN(unit=1,file=filename,status='old',action='read',form='unformatted',access="stream")

  ! Init block address
  head_blck  = -1
  pos_blck   = -1
  vel_blck   = -1
  id_blck    = -1
  u_blck     = -1
  mass_blck  = -1
  metal_blck = -1
  age_blck   = -1

  ! Find all the data blocks in a Gadget1 file
  if(r%ic_format .eq. 'Gadget1') then
     ! Init block counter
     jump_blck = 1
     blck_cnt = 1
     do while(.true.)
        ! Reading data block header
        read(1,POS=jump_blck,iostat=stat) blck_size
        if(stat /= 0) exit
        ! Saving data block positions
        if(blck_cnt .eq. 1) then
           head_blck  = jump_blck+sizeof(blck_size)
           head_size  = blck_size
           if(g%myid==1)write(*,*)blck_cnt,blck_size
        endif
        if(blck_cnt .eq. 2) then
           pos_blck   = jump_blck+sizeof(blck_size)
           pos_size   = blck_size/(3*sizeof(dummy_real))
           if(g%myid==1)write(*,*)blck_cnt,blck_size
        endif
        if(blck_cnt .eq. 3) then
           vel_blck   = jump_blck+sizeof(blck_size)
           vel_size   = blck_size/(3*sizeof(dummy_real))
           if(g%myid==1)write(*,*)blck_cnt,blck_size
        endif
        if(blck_cnt .eq. 4) then
           id_blck    = jump_blck+sizeof(blck_size)
           id_size    = blck_size/sizeof(dummy_int)
           if(g%myid==1)write(*,*)blck_cnt,blck_size
        endif
        if(blck_cnt .eq. 5) then
           u_blck     = jump_blck+sizeof(blck_size)
           u_size     = blck_size/sizeof(dummy_real)
           if(g%myid==1)write(*,*)blck_cnt,blck_size
        endif
        if(blck_cnt .eq. 6) then
           mass_blck  = jump_blck+sizeof(blck_size)
           mass_size  = blck_size/sizeof(dummy_real)
           if(g%myid==1)write(*,*)blck_cnt,blck_size
        endif
        if(blck_cnt .eq. 7) then
           metal_blck = jump_blck+sizeof(blck_size)
           metal_size = blck_size/sizeof(dummy_real)
           if(g%myid==1)write(*,*)blck_cnt,blck_size
        endif
        if(blck_cnt .eq. 8) then
           age_blck   = jump_blck+sizeof(blck_size)
           age_size   = blck_size/sizeof(dummy_real)
           if(g%myid==1)write(*,*)blck_cnt,blck_size
        endif
        jump_blck = jump_blck+blck_size+2*sizeof(dummy_int)
        blck_cnt = blck_cnt+1
     enddo
  endif

  ! Find all the data blocks in a Gadget2 file
  if(r%ic_format .eq. 'Gadget2') then
     ! Init block counter
     jump_blck = 1
     do while(.true.)
        ! Reading data block header
        read(1,POS=jump_blck,iostat=stat) dummy_int
        if(stat /= 0) exit
        read(1,POS=jump_blck+sizeof(dummy_int),iostat=stat) blck_name
        if(stat /= 0) exit
        read(1,POS=jump_blck+sizeof(dummy_int)+sizeof(blck_name),iostat=stat) dummy_int
        if(stat /= 0) exit
        read(1,POS=jump_blck+2*sizeof(dummy_int)+sizeof(blck_name),iostat=stat) dummy_int
        if(stat /= 0) exit
        read(1,POS=jump_blck+3*sizeof(dummy_int)+sizeof(blck_name),iostat=stat) blck_size
        if(stat /= 0) exit
        ! Saving data block positions
        if(blck_name .eq. r%ic_head_name) then
           head_blck  = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
           head_size  = blck_size
           if(g%myid==1)write(*,*) '-> Found ',blck_name,' block ',head_blck,head_size
        endif
        if(blck_name .eq. r%ic_pos_name) then
           pos_blck   = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
           pos_size   = blck_size/(3*sizeof(dummy_real))
           if(g%myid==1)write(*,*) '-> Found ',blck_name,' block ',pos_blck,pos_size
        endif
        if(blck_name .eq. r%ic_vel_name) then
           vel_blck  = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
           vel_size  = blck_size/(3*sizeof(dummy_real))
           if(g%myid==1)write(*,*) '-> Found ',blck_name,' block ',vel_blck,vel_size
        endif
        if(blck_name .eq. r%ic_id_name) then
           id_blck    = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
           id_size    = blck_size/sizeof(dummy_int)
           if(g%myid==1)write(*,*) '-> Found ',blck_name,' block ',id_blck,id_size
        endif
        if(blck_name .eq. r%ic_mass_name) then
           mass_blck  = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
           mass_size  = blck_size/sizeof(dummy_real)
           if(g%myid==1)write(*,*) '-> Found ',blck_name,' block ',mass_blck,mass_size
        endif
        if(blck_name .eq. r%ic_u_name) then
           u_blck     = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
           u_size     = blck_size/sizeof(dummy_real)
           if(g%myid==1)write(*,*) '-> Found ',blck_name,' block ',u_blck,u_size
        endif
        if(blck_name .eq. r%ic_metal_name) then
           metal_blck = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
           metal_size = blck_size/sizeof(dummy_real)
           if(g%myid==1)write(*,*) '-> Found ',blck_name,' block ',metal_blck,metal_size
        endif
        if(blck_name .eq. r%ic_age_name) then
           age_blck   = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
           age_size   = blck_size/sizeof(dummy_real)
           if(g%myid==1)write(*,*) '-> Found ',blck_name,' block ',age_blck,age_size
        endif
        jump_blck = jump_blck+blck_size+sizeof(blck_name)+5*sizeof(dummy_int)
     enddo
  endif

  if((head_blck.eq.-1).or.(pos_blck.eq.-1).or.(vel_blck.eq.-1)) then
     write(*,*) 'Gadget file does not contain handful data'
     call mdl_abort(mdl)
  endif

  if(head_size.ne.256) then
     write(*,*) 'Gadget header is not 256 bytes'
     call mdl_abort(mdl)
  endif

  ! Read Gadget file header
  ! Byte swapping doesn't appear to work if you just do READ(1)header
  READ(1,POS=head_blck)header%npart,header%mass,header%time,header%redshift, &
       header%flag_sfr,header%flag_feedback,header%nparttotal, &
       header%flag_cooling,header%numfiles,header%boxsize, &
       header%omega0,header%omegalambda,header%hubbleparam, &
       header%flag_stellarage,header%flag_metals,header%totalhighword, &
       header%flag_entropy_instead_u, header%flag_doubleprecision, &
       header%flag_ic_info, header%lpt_scalingfactor

  npart = sum(header%npart) ! All particles in the file
  ngas = header%npart(1) ! Type 0 particles are gas particles
  nhalo = header%npart(2) ! Type 1 particles are dark matter particles
  nstar = sum(header%npart(3:5)) ! Type 2, 3, 4 particles are star particles

  if(g%myid==1)write(*,*)"Found ",npart," particles in total"

  ! Check if user asked to skip some particle types
  ! For example if gas is initialized from condinit
  do k=0,5
     write(type_str,*)k
     skip(k+1)=.false.
     do j=1,6
        if(r%ic_skip_type(j).eq.k) skip(k+1)=.true.
     enddo
     if(.not.skip(k+1).and.header%npart(k+1)>0)then
        if(g%myid==1)then
           if(header%mass(k+1)>0)then
              write(*,*)"----> ",header%npart(k+1)," type "//TRIM(ADJUSTL(type_str))// &
                   & " particles with header mass ",header%mass(k+1)
           else
              write(*,*)"----> ",header%npart(k+1)," type "//TRIM(ADJUSTL(type_str))// &
                   & " particles"
           endif
        endif
     endif
  end do

  if((pos_size.ne.npart).or.(vel_size.ne.npart)) then
     write(*,*) 'POS =',pos_size
     write(*,*) 'VEL =',vel_size
     write(*,*) 'Number of particles does not correspond to block sizes'
     call mdl_abort(mdl)
  endif

  ! Read gas particles
  if(r%hydro.and.ngas>0.and..not.skip(1))then
     if(g%myid==1)write(*,'(A)') " Reading gas particles..."
     gas%npart=ngas/g%ncpu ! Particles are split evenly across processors
     igas=(g%myid-1)*gas%npart
     nrest=ngas-gas%npart*g%ncpu
     if(g%myid.LE.nrest)then
        gas%npart=gas%npart+1
        igas=igas+(g%myid-1)
     else
        igas=igas+nrest
     endif
     call allocate_gas(r,g,gas) ! Allocate sph-like gas particles
     ipbar = 0.
     do i = 1,gas%npart
        kpart = igas+i
        read(1,POS=pos_blck+3*sizeof(dummy_real)*(kpart-1))x
        gas%xp(i,1:3)=x
        read(1,POS=vel_blck+3*sizeof(dummy_real)*(kpart-1))v
        gas%vp(i,1:3)=v
        if(mass_blck.ne.-1) then
           kpart = igas+i
           read(1,POS=mass_blck+sizeof(dummy_real)*(kpart-1))z
           gas%mp(i)=z
        endif
        if(metal_blck.ne.-1) then
           kpart = igas+i
           read(1,POS=metal_blck+sizeof(dummy_real)*(kpart-1))z
           gas%zp(i)=z
        endif
        if(u_blck.ne.-1) then
           kpart = igas+i
           read(1,POS=u_blck+sizeof(dummy_real)*(kpart-1))z
           gas%up(i)=z
        endif
        if(i.ge.ipbar*(gas%npart/49.0))then
           if(g%myid==1)write(*,'(A1)',advance='no') "_"
           ipbar=ipbar+1.0
        endif
     end do
     if(g%myid==1)write(*,'(A1)') " "
  endif

  ! Read halo particles
  if(r%part.and.nhalo>0.and..not.skip(2))then
     if(g%myid==1)write(*,'(A)') " Reading halo particles..."
     p%npart=nhalo/g%ncpu
     p%npart_tot=nhalo
     ihalo=(g%myid-1)*p%npart
     nrest=nhalo-p%npart*g%ncpu
     if(g%myid.LE.nrest)then
        p%npart=p%npart+1
        ihalo=ihalo+(g%myid-1)
     else
        ihalo=ihalo+nrest
     endif
     ipbar=0.
     do i=1,p%npart
        kpart = ngas+ihalo+i
        read(1,POS=pos_blck+3*sizeof(dummy_real)*(kpart-1))x
        p%xp(i,1:3)=x
        read(1,POS=vel_blck+3*sizeof(dummy_real)*(kpart-1))v
        p%vp(i,1:3)=v
        if(mass_blck.ne.-1) then
           kpart = ngas+ihalo+i
           read(1,POS=mass_blck+sizeof(dummy_real)*(kpart-1))z
           p%mp(i)=z
        endif
        if(id_blck.ne.-1) then
           kpart = ngas+ihalo+i
           read(1,POS=id_blck+sizeof(dummy_int)*(kpart-1))id
           p%idp(i)=id
        endif
        if(i.ge.ipbar*(p%npart/49.0))then
           if(g%myid==1)write(*,'(A1)',advance='no') "_"
           ipbar=ipbar+1.0
        endif
     end do
     if(g%myid==1)write(*,'(A1)') " "
  endif

  ! Read star particles
  if(r%star.and.nstar>0.and..not.skip(3))then
     if(g%myid==1)write(*,'(A)') " Reading star particles..."
     star%npart=nstar/g%ncpu
     star%npart_tot=nstar
     istar=(g%myid-1)*star%npart
     nrest=nstar-star%npart*g%ncpu
     if(g%myid.LE.nrest)then
        star%npart=star%npart+1
        istar=istar+(g%myid-1)
     else
        istar=istar+nrest
     endif
     do i=1,star%npart
        kpart = ngas+nhalo+istar+i
        read(1,POS=pos_blck+3*sizeof(dummy_real)*(kpart-1))x
        star%xp(i,1:3)=x
        read(1,POS=vel_blck+3*sizeof(dummy_real)*(kpart-1))v
        star%vp(i,1:3)=v
        if(mass_blck.ne.-1) then
           kpart = ngas+nhalo+istar+i
           read(1,POS=mass_blck+sizeof(dummy_real)*(kpart-1))z
           star%mp(i)=z
        endif
        if(id_blck.ne.-1) then
           kpart = ngas+nhalo+istar+i
           read(1,POS=id_blck+sizeof(dummy_int)*(kpart-1))id
           star%idp(i)=id
        endif
        if(metal_blck.ne.-1) then
           kpart = ngas+nhalo+istar+i
           read(1,POS=metal_blck+sizeof(dummy_real)*(kpart-1))z
           star%zp(i)=z
        endif
        if(age_blck.ne.-1) then
           kpart = istar+i
           read(1,POS=age_blck+sizeof(dummy_real)*(kpart-1))z
           star%tp(i)=z
        endif
        if(i.ge.ipbar*(star%npart/49.0))then
           if(g%myid==1)write(*,'(A1)',advance='no') "_"
           ipbar=ipbar+1.0
        endif
     end do
     if(g%myid==1)write(*,'(A1)') " "
  endif

  if(g%myid==1)write(*,'(A)')' '//TRIM(r%ic_format)//' file successfully loaded'
  if(g%myid==1)write(*,'(A50)')"__________________________________________________"

  ! Use header mass if missing mass block
  if(mass_blck.eq.-1)then
     if(r%part)p%mp(1:p%npart)=header%mass(2)
     if(r%hydro)gas%mp(1:gas%npart)=header%mass(1)
     if(r%star)star%mp(1:star%npart)=header%mass(3)
  endif

  ! Set particle level to levelmin
  if(r%part)p%levelp(1:p%npart)=r%levelmin
  if(r%hydro)gas%levelp(1:gas%npart)=r%levelmin
  if(r%star)star%levelp(1:star%npart)=r%levelmin

  ! Compute ids if missing id block
  if(id_blck.eq.-1)then
     if(r%part)then
        do i=1,p%npart
           p%idp(i)=ihalo+i
        end do
     endif
     if(r%star)then
        do i=1,star%npart
           star%idp(i)=istar+i
        end do
     end if
  endif

  ! Put default metallicity if missing metal block
  if(metal_blck.eq.-1)then
     if(r%hydro)gas%zp(1:gas%npart)=0.02*r%z_ave
     if(r%star)star%zp(1:star%npart)=0.02*r%z_ave
  endif

  ! Rescale from Gadget code units to Ramses code units
  if(r%part)call rescale_gadget(r,g,p)
  if(r%hydro)call rescale_gadget(r,g,gas)
  if(r%star)call rescale_gadget(r,g,star)

  ! Trim particles that are outside of the Ramses box
  if(r%part)call trim_box(r,g,p)
  if(r%hydro)call trim_box(r,g,gas)
  if(r%star)call trim_box(r,g,star)

  ! Compute total mass in star, gas and dark matter
  mstar=0.; mhalo=0.; mgas=0.
  if(r%part)mhalo=sum(p%mp(1:p%npart))
  if(r%hydro)mgas=sum(gas%mp(1:gas%npart))
  if(r%star)mstar=sum(star%mp(1:star%npart))

  ! Compute number of particles in star, gas and dark matter
  npart_star=0; npart_halo=0; npart_gas=0
  if(r%part)npart_halo=p%npart
  if(r%hydro)npart_gas=gas%npart
  if(r%star)npart_star=star%npart

  ! Put all particles inside levelmin
  if(r%part)call init_levelmin(r,p)
  if(r%hydro)call init_levelmin(r,gas)
  if(r%star)call init_levelmin(r,star)

  end associate

end subroutine input_part_gadget
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine rescale_gadget(r,g,p)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t, global_t
  use pm_commons, only: part_t
  type(part_t)::p
  type(run_t)::r
  type(global_t)::g
  ! This routine convert particle variables from Gadget code units
  ! to Ramses code units
  integer::i
  real(kind=8)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,scale_m
  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m = scale_d*scale_l**3
  do i=1,p%npart
     p%xp(i,1:ndim)=p%xp(i,1:ndim)*(r%gadget_scale_l/scale_l)+r%box_size(1:ndim)/2.0d0
     p%vp(i,1:ndim)=p%vp(i,1:ndim)*(r%gadget_scale_v/scale_v)
     p%mp(i)=p%mp(i)*(r%gadget_scale_m/scale_m)
  end do
  if(allocated(p%up))then
     do i=1,p%npart
        p%up(i)=p%up(i)*(r%gadget_scale_v/scale_v)**2
     end do
  endif
  if(allocated(p%tp))then
     do i=1,p%npart
        p%tp(i)=p%tp(i)*(r%gadget_scale_t/scale_t)
     end do
  endif
end subroutine rescale_gadget
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine trim_box(r,g,p)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t, global_t
  use pm_commons, only: part_t
  type(part_t)::p
  type(run_t)::r
  type(global_t)::g
  ! This routine remove particles that are outside of the Ramses box.
  ! Surviving particles are compacted.
  integer::ipart,i,idim
  logical::ok_part
  ipart=0
  do i=1,p%npart
     ok_part=.true.
     do idim=1,ndim
        ok_part=ok_part.and.p%xp(i,idim)>0.and.p%xp(i,idim)<r%box_size(idim)
     end do
     if(ok_part)then
        ipart=ipart+1
        p%sortp(ipart)=i
     endif
  end do
  do i=1,ipart
     p%xp(i,1:ndim)=p%xp(p%sortp(i),1:ndim)
     p%vp(i,1:ndim)=p%vp(p%sortp(i),1:ndim)
     p%mp(i)=p%mp(p%sortp(i))
  end do
  if(allocated(p%idp))then
     do i=1,ipart
        p%idp(i)=p%idp(p%sortp(i))
     end do
  end if
  if(allocated(p%zp))then
     do i=1,ipart
        p%zp(i)=p%zp(p%sortp(i))
     end do
  end if
  if(allocated(p%tp))then
     do i=1,ipart
        p%tp(i)=p%tp(p%sortp(i))
     end do
  end if
  if(allocated(p%up))then
     do i=1,ipart
        p%up(i)=p%up(p%sortp(i))
     end do
  end if
  ! Set variables to zero for trimmed particles
  do i=ipart+1,p%npart
     p%xp(i,1:ndim)=0.0
     p%vp(i,1:ndim)=0.0
     p%mp(i)=0.0
  end do
  if(allocated(p%idp))then
     do i=ipart+1,p%npart
        p%idp(i)=0
     end do
  endif
  if(allocated(p%zp))then
     do i=ipart+1,p%npart
        p%zp(i)=0.0
     end do
  endif
  if(allocated(p%tp))then
     do i=ipart+1,p%npart
        p%tp(i)=0.0
     end do
  endif
  if(allocated(p%up))then
     do i=ipart+1,p%npart
        p%up(i)=0.0
     end do
  endif
  p%npart=ipart
end subroutine trim_box
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_levelmin(r,p)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t
  use pm_commons, only: part_t
  type(run_t)::r
  type(part_t)::p
  ! This routine puts all particles at the coarsest level.
  p%headp=p%npart+1
  p%tailp=p%npart
  p%headp(r%levelmin)=1
  p%tailp(r%levelmin)=p%npart
  if(ANY(.not.r%periodic(1:ndim)))then
     p%headp(r%levelmin-1)=1
     p%tailp(r%levelmin-1)=0
  endif
end subroutine init_levelmin
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module input_part_gadget_module
