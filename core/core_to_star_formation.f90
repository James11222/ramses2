module core_to_star_formation_module
  implicit none

  ! Input: level. Output: star mass formed.
  type :: in_core_to_star_t
     integer::ilevel
  end type in_core_to_star_t
  type :: out_core_to_star_t
     real(kind=8)::mass
  end type out_core_to_star_t

contains

recursive subroutine r_core_to_star_formation(pst,input,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  use star_formation_common_module, only: finish_star_formation
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  type(in_core_to_star_t)::input
  type(out_core_to_star_t)::output,next_output
  integer::ilevel,rID

  ilevel=input%ilevel
  if(pst%nLower>0)then
     rID=mdl_send_request(pst%s%mdl,MDL_CORE_TO_STAR_FORMATION,pst%iUpper+1,input_size,output_size,input)
     call r_core_to_star_formation(pst%pLower,input,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
  else
     call star_formation_from_core(pst%s,pst%s%r,pst%s%g,pst%s%star, &
          pst%s%core,ilevel,output%mass)
     call finish_star_formation(pst%s)
  endif
end subroutine r_core_to_star_formation

!###########################################################
subroutine star_formation_from_core(s,r,g,star,core,ilevel,mstar_loc)
  use constants
  use amr_parameters, only: ndim
  use amr_commons, only: run_t, global_t
  use core_sf_support_module, only: core_top_level
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use mdl_module, only: mdl_abort
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  type(run_t)::r
  type(global_t)::g
  type(part_t)::star
  type(part_t)::core
  integer::ilevel
  real(kind=8)::mstar_loc
  !-------------------------------------------------------------------
  ! Convert core mass into stars in units of core_to_star_mass.
  ! The gas was already removed by core accretion.
  !-------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:g%ncpu)::nstar_cpu_tot
#endif
  integer(kind=8),dimension(0:g%ncpu)::nstar_cum
  integer,dimension(1:g%ncpu)::nstar_cpu
  integer::i,idim,nstar_loc,ibirth
  integer,allocatable::birth_core_index(:)
  real(kind=8),allocatable::birth_mass(:)
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_msun
  real(kind=8)::mstar_unit,mnew

  mstar_loc=0.0d0

  if(ilevel/=core_top_level(r,g))return

  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun=scale_l**3*scale_d/M_sun
  mstar_unit=r%core_to_star_mass/scale_msun

  ! Pick the cores that form a star, so the new ids can be numbered first.
  nstar_loc=0
  allocate(birth_core_index(max(1,core%npart)),birth_mass(max(1,core%npart)))
  birth_core_index=0
  do i=1,core%npart
     mnew=mstar_unit*dble(int(core%mp(i)/mstar_unit))
     if(mnew>0d0)then
        nstar_loc=nstar_loc+1
        birth_core_index(nstar_loc)=i
        birth_mass(nstar_loc)=mnew
     endif
  end do

  nstar_cpu=0
  nstar_cpu(g%myid)=nstar_loc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nstar_cpu,nstar_cpu_tot,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nstar_cpu=nstar_cpu_tot
#endif
  nstar_cum=0
  do i=1,g%ncpu
     nstar_cum(i)=nstar_cum(i-1)+int(nstar_cpu(i),kind=8)
  end do

  do ibirth=1,nstar_loc
     i=birth_core_index(ibirth)
     mnew=birth_mass(ibirth)
     mstar_loc=mstar_loc+mnew
     star%npart=star%npart+1
     if(star%npart>r%nstarmax)then
        write(*,*)'Not enough memory for star particle'
        write(*,*)'Increase nstarmax in the namelist'
        call mdl_abort(s%mdl)
     endif
     do idim=1,ndim
        star%xp(star%npart,idim)=core%xp(i,idim)
        star%vp(star%npart,idim)=core%vp(i,idim)
     end do
     star%mp(star%npart)=mnew
     star%tp(star%npart)=g%texp
     star%zp(star%npart)=core%zp(i)
     star%levelp(star%npart)=ilevel
     core%mp(i)=max(0d0,core%mp(i)-mnew)
  end do

  star%tailp(r%nlevelmax)=star%tailp(r%nlevelmax)+nstar_loc

  do i=1,nstar_loc
     star%idp(star%npart-nstar_loc+i)=star%npart_tot+nstar_cum(g%myid-1)+i
  end do
  star%npart_tot=star%npart_tot+nstar_cum(g%ncpu)
  deallocate(birth_core_index,birth_mass)

end subroutine star_formation_from_core

end module core_to_star_formation_module
