module output_clump_module
  use clump_merger_module
contains
!###################################################
!###################################################
!###################################################
!###################################################
recursive subroutine r_output_clump(pst,input_array,input_size,output_array,output_size)
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

  character(LEN=flen)::filename,fileloc
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_CLUMP,pst%iUpper+1,input_size,output_size,input_array)
     call r_output_clump(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     filename=transfer(input_array,filename)
     ! Write clump properties
     if(pst%s%r%output_clump)then
        call output_clump_properties(pst%s,filename)
     endif
     ! Write density, peak id and halo id for AMR cells
     if(pst%s%r%output_peak_grid)then
        fileloc=TRIM(filename)//'peak_grid.'
        call output_clump_field(pst%s,fileloc)
     endif
     ! Write peak id and halo id for dark matter particles
     if(pst%s%r%output_peak_part.and.pst%s%r%part)then
        fileloc=TRIM(filename)//'peak_part.'
        call output_peak_part(pst%s,pst%s%p   ,fileloc)
     endif
     ! Write peak id and halo for star particles
     if(pst%s%r%output_peak_star.and.pst%s%r%star)then
        fileloc=TRIM(filename)//'peak_star.'
        call output_peak_part(pst%s,pst%s%star,fileloc)
     endif
     ! Write peak id and halo for sink particles
     if(pst%s%r%output_peak_sink.and.pst%s%r%sink)then
        fileloc=TRIM(filename)//'peak_sink.'
        call output_peak_part(pst%s,pst%s%sink,fileloc)
     endif
     ! Write peak id and halo for tree particles
     if(pst%s%r%output_peak_tree.and.pst%s%r%tree)then
        fileloc=TRIM(filename)//'peak_tree.'
        call output_peak_part(pst%s,pst%s%tree,fileloc)
     endif
  endif

end subroutine r_output_clump
!###################################################
!###################################################
!###################################################
!###################################################
subroutine output_clump_properties(s,filename)
  use amr_parameters, only: flen, nbin, ndim
  use ramses_commons, only: ramses_t,open_file,close_file
  use clfind_commons
  implicit none
  type(ramses_t)::s
  character(LEN=flen)::filename
  !----------------------------------------------------------------
  ! This routine output the clump properties for each processor
  !----------------------------------------------------------------
  integer::ilun,j,ind,igrid,hid
  character(LEN=flen)::fileloc
  integer(kind=8),dimension(s%r%levelmin:s%r%nlevelmax)::nskip
  real(kind=8)::pi,grav,rad,mass,r200c,rmax,concentration,purity
  real(kind=8),dimension(1:nbin)::mbin
  logical::part_clump,ok_profile

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c)

  pi=ACOS(-1d0)
  grav=1d0
  if(s%r%cosmo)grav=3d0/8d0/pi*s%g%omega_m*s%g%aexp

  ! The particle count criterion only applies to particle-based clumps
  part_clump=(r%rho_type_clump>=1.and.r%rho_type_clump<=3)
  ! Mass profiles are only computed for dark matter clumps with halo merging
  ok_profile=(c%saddle_threshold>0.and.r%part.and.r%rho_type_clump==1)

  ! Write clump file
  fileloc=TRIM(filename)//'clump.'

  call open_file(s,fileloc,nskip,ilun)

  do j=1,c%npeak
     if (    c%relevance(j) > c%relevance_threshold.AND. &
          & c%clump_mass(j) > c%mass_threshold.AND. &
          & (c%npart(j) > 0 .or. .not. part_clump))then

        if(ok_profile)then
           ! Get cumulative mass profile
           mbin=c%mass_bin(j,1:nbin)
           ! Get clump tidal radius
           rad=c%clump_rad(j)
           ! Get clump particle mass
           mass=c%mass_bin(j,nbin)
           ! Comnpute r200, rmax and concentration
           call halo_mass_def(s,mbin,rad,r200c,rmax,concentration)
           ! Compute clump purity
           purity=c%npart(j)*g%mp_min/mass
        else
           ! No mass profile: use the patch mass and skip the purity cut
           mass=c%clump_mass(j)
           r200c=0d0; rmax=0d0; concentration=0d0
           purity=huge(purity)
        endif
        ! Without halo merging, each clump is its own halo
        if(c%saddle_threshold>0)then
           hid=c%ind_halo(j)
        else
           hid=c%ind_final(j)
        endif

        ! Output clump properties to file
        if(purity>c%purity_threshold)then
           write(ilun,'(I10,4(1X,I10),15(1X,1PE18.9E2))')&
                 c%ind_final(j)&
                ,c%new_peak(j)&
                ,hid&
                ,c%n_cells(j)&
                ,c%npart(j)&
                ,c%min_dens(j)&
                ,c%max_dens(j)&
                ,c%tidal_dens(j)&
                ,c%clump_mass(j)/c%clump_vol(j)&
                ,c%clump_mass(j)&
                ,mass&
                ,r200c&
                ,rmax&
                ,concentration&
                ,c%peak_pos(j,1),c%peak_pos(j,2),c%peak_pos(j,3)&
                ,c%peak_vel(j,1),c%peak_vel(j,2),c%peak_vel(j,3)
        endif

     end if
  end do

  call close_file(s,fileloc,nskip,ilun)

  end associate

end subroutine output_clump_properties
!###################################################
!###################################################
!###################################################
!###################################################
subroutine halo_mass_def(s,mbin,rad,r200c,rmax,c)
  use amr_parameters, only: nbin
  use ramses_commons, only: ramses_t
  type(ramses_t)::s
  real(kind=8),dimension(1:nbin)::mbin
  real(kind=8)::rad,rmax,r200c,c,deltamax,cmin,cmax
  real(kind=8)::pi,G,delta,deltaold,rbin,vcirc,volbin,alpha,omega_mz
  real(kind=8)::d200c,vmax,v200c,c0,cl,cr,err,const,dr
  integer::i,imax
  ! Constants
  pi=ACOS(-1.0D0)
  G=1d0
  if(s%r%cosmo)G=3d0/8d0/pi*s%g%omega_m*s%g%aexp
  if(s%r%cosmo)omega_mz = (s%g%omega_m/s%g%aexp**3) / (s%g%omega_m/s%g%aexp**3 + s%g%omega_l) !include redshift dependence

  d200c = 200d0/omega_mz !200c is 200 times the critical density, not the mean density
  ! Find densest bin
  deltamax=0
  imax=1
  ! Loop over radial bins
  dr=2d0*rad/dble(nbin)
  do i=1,nbin
     if(mbin(i)==0)cycle
     rbin=dble(i)*dr
     volbin=4d0/3d0*pi*rbin**3
     delta=mbin(i)/volbin
     if(delta>deltamax)then
        deltamax=delta
        imax=i
     endif
  end do
  ! Initializations
  vmax=0
  r200c=0
  rmax=0
  delta=0
  ! Loop over radial bins
  do i=imax,nbin
     if(mbin(i)==0)cycle
     rbin=dble(i)*dr
     vcirc=sqrt(G*mbin(i)/rbin)
     if(vcirc>vmax)then
        vmax=vcirc
        rmax=rbin
     endif
     volbin=4d0/3d0*pi*rbin**3
     deltaold=delta
     delta=mbin(i)/volbin
     if(delta<=d200c.and.r200c==0)then
        if(deltaold>0)then
           alpha=log(d200c/delta)/log(deltaold/delta)
           r200c=rbin*(dble(i-1)/dble(i))**alpha
        else
           r200c=rbin
        endif
     endif
  end do
  ! Compute quantities at r200c
  if(delta>d200c.and.r200c==0)then
     alpha=log(d200c/delta)/log(deltaold/delta)
     r200c=2d0*rad*(dble(nbin-1)/dble(nbin))**alpha
  endif
  c0=2.1626
  v200c=sqrt(G*4d0*pi/3d0*d200c*r200c**2)
  ! Find concentration parameter
  vmax=max(vmax,v200c)
  const=MIN(dble(vmax**2/v200c**2),4d0)*c0/(log(1d0+c0)-c0/(1d0+c0))
  ! Find large root
  cl=c0
  cr=100d0
  err=1d0
  i=0
  do while(abs(err)>1d-3.and.i<100)
     c=0.5*(cl+cr)
     err=c/(log(1d0+c)-c/(1d0+c))-const
     if(err>0)then
        cr=c
     else
        cl=c
     endif
     i=i+1
  end do
  cmax=c
  ! Find small root
  cl=0.1
  cr=c0
  err=1d0
  i=0
  do while(abs(err)>1d-3.and.i<100)
     c=0.5*(cl+cr)
     err=c/(log(1d0+c)-c/(1d0+c))-const
     if(err>0)then
        cl=c
     else
        cr=c
     endif
     i=i+1
  end do
  cmin=c
  ! Choose the right root
  if(rmax>r200c)then
     c=cmin
  else
     c=cmax
  endif
!!$  write(*,*)c,r200c,v200c,rmax,vmax
!!$  do i=1,nbin
!!$     if(mbin(i)==0)cycle
!!$     rbin=dble(i)*dr
!!$     vcirc=sqrt(G*mbin(i)/rbin)
!!$     volbin=4d0/3d0*pi*rbin**3
!!$     delta=mbin(i)/volbin
!!$     write(*,*)rbin,mbin(i),vcirc,delta
!!$  end do
end subroutine halo_mass_def
!###################################################
!###################################################
!###################################################
!###################################################
subroutine output_clump_field(s,filename)
  use amr_parameters, only: ndim,twotondim,flen
  use ramses_commons, only: ramses_t,open_file,close_file
  implicit none
  type(ramses_t)::s
  character(LEN=flen)::filename
  !----------------------------------------------------------------
  ! This routine output the peak patch fields for each processor
  !----------------------------------------------------------------
  integer::ilevel,igrid,ilun
  integer(kind=8),dimension(s%r%levelmin:s%r%nlevelmax)::nskip
  real(kind=4),dimension(1:twotondim)::flg1,flg2
  real(kind=4),dimension(1:twotondim)::rho

  associate(g=>s%g,r=>s%r,m=>s%m,mdl=>s%mdl)

#ifdef GRAV

  call open_file(s,filename,nskip,ilun)

  do ilevel=r%levelmin,r%nlevelmax
     write(ilun,POS=nskip(ilevel))
     do igrid=m%head(ilevel),m%tail(ilevel)
        rho(1:twotondim)=real(m%rho(1:twotondim,igrid),kind=4)
        flg1=real(m%flag1(1:twotondim,igrid),kind=4)
        flg2=real(m%flag2(1:twotondim,igrid),kind=4)
        write(ilun)flg1
        write(ilun)flg2
        write(ilun)rho
     end do
  enddo

  call close_file(s,filename,nskip,ilun)

#endif

  end associate

end subroutine output_clump_field
!###################################################
!###################################################
!###################################################
!###################################################
subroutine file_descriptor_clump(r,filename)
  use amr_parameters, only: ndim,flen
  use amr_commons, only: run_t
  implicit none
  type(run_t)::r
  character(LEN=flen)::filename

  character(LEN=flen)::fileloc
  integer::ivar,ilun

  if(r%verbose)write(*,*)'Entering file_descriptor_clump'

  ilun=11
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  write(ilun,'("nvar        =",I11)')3

  ivar=1
  write(ilun,'("variable #",I2,": halo patch ID")')ivar

  ivar=2
  write(ilun,'("variable #",I2,": peak patch ID")')ivar

  ivar=3
  write(ilun,'("variable #",I2,": density")')ivar

  close(ilun)

end subroutine file_descriptor_clump
!#######################################################
!#######################################################
!#######################################################
!#######################################################
subroutine output_peak_part(s,p,filename)
  use amr_parameters, only: ndim,i8b,flen
  use ramses_commons, only: ramses_t,open_part_file,close_part_file
  use pm_commons, only: part_t
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  character(LEN=flen)::filename
  !------------------------------------
  ! Output peak id of particles to file
  !------------------------------------
  integer::i,idim,ilun,ivar
  integer(kind=8),dimension(1:p%nvaralloc+1)::nskip
  integer,allocatable,dimension(:)::id

  associate(r=>s%r,g=>s%g)

  call open_part_file(s,p,filename,nskip,ilun)

  allocate(id(1:p%npart))

  ! Write halo id
  do i=1,p%npart
     id(i)=p%hid(i)
  end do
  write(ilun,POS=nskip(1))
  write(ilun)id

  ! Write peak id
  do i=1,p%npart
     id(i)=p%pid(i)
  end do
  write(ilun,POS=nskip(2))
  write(ilun)id

  deallocate(id)

  call close_part_file(s,p,filename,nskip,ilun)

  end associate

end subroutine output_peak_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_peak_header(r,g,p,filename)
  use amr_parameters, only: flen
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  character(LEN=flen)::filename

  ! Local variables
  integer::ilun
  character(LEN=flen)::fileloc

  if(r%verbose)write(*,*)'Entering output_peak_header'

  ilun=10!+g%myid

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write header information
  write(ilun,*)'Total number of particles'
  write(ilun,*)p%npart_tot

  write(ilun,*)'Total number of files'
  write(ilun,*)r%nfile

  ! Keep track of what particle fields are present
  write(ilun,*)'Particle fields'
  write(ilun,'(a)',advance='no')'halo_id '
  write(ilun,'(a)',advance='no')'peak_id '
  close(ilun)

end subroutine output_peak_header
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module output_clump_module
