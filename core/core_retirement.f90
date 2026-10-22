module core_retirement_module
  implicit none

  ! Input: coarse step. Output: core mass turned into stars.
  type :: in_core_retire_t
     integer::nstep_coarse
  end type in_core_retire_t
  type :: out_core_retire_t
     real(kind=8)::mass
  end type out_core_retire_t

contains

recursive subroutine r_core_retirement(pst,input,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  type(in_core_retire_t)::input
  type(out_core_retire_t)::output,next_output
  integer::rID

  if(pst%nLower>0)then
     rID=mdl_send_request(pst%s%mdl,MDL_CORE_RETIREMENT,pst%iUpper+1,input_size,output_size,input)
     call r_core_retirement(pst%pLower,input,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
  else
     call core_retirement(pst%s,pst%s%r,pst%s%g,pst%s%star,pst%s%core, &
          input%nstep_coarse,output%mass)
  endif
end subroutine r_core_retirement

!###########################################################
subroutine core_retirement(s,r,g,star,core,nstep_coarse,mstar_loc)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t, global_t
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use mdl_module, only: mdl_abort
  use star_formation_common_module, only: finish_star_formation
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  type(run_t)::r
  type(global_t)::g
  type(part_t)::star
  type(part_t)::core
  integer::nstep_coarse
  real(kind=8)::mstar_loc
  !-------------------------------------------------------------------
  ! Turn the cores that have not accreted for core_zero_acc_intervals
  ! coarse steps into stars. Their rows stay with zero mass, as merged
  ! sinks do.
  !-------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:g%ncpu)::nstar_cpu_tot
#endif
  integer(kind=8),dimension(0:g%ncpu)::nstar_cum
  integer,dimension(1:g%ncpu)::nstar_cpu
  integer::i,nstar_loc
  logical,allocatable::retire(:)

  mstar_loc=0.0d0

  ! nacc is the step of the last accretion, or of the formation
  nstar_loc=0
  allocate(retire(max(1,core%npart)))
  do i=1,core%npart
     retire(i)=core%mp(i)>0d0.and.nstep_coarse-core%nacc(i)>r%core_zero_acc_intervals
     if(retire(i))nstar_loc=nstar_loc+1
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

  ! Append the retiring cores to the stars
  if(nstar_loc>0)then
     do i=1,core%npart
        if(.not.retire(i))cycle
        star%npart=star%npart+1
        if(star%npart>r%nstarmax)then
           write(*,*)'Not enough memory for star particle'
           write(*,*)'Increase nstarmax in the namelist'
           call mdl_abort(s%mdl)
        endif
        star%xp(star%npart,1:ndim)=core%xp(i,1:ndim)
        star%vp(star%npart,1:ndim)=core%vp(i,1:ndim)
        star%mp(star%npart)=core%mp(i)
        star%tp(star%npart)=g%texp
        star%zp(star%npart)=core%zp(i)
        star%levelp(star%npart)=core%levelp(i)
        mstar_loc=mstar_loc+core%mp(i)
        core%mp(i)=0d0
     end do
     star%tailp(r%nlevelmax)=star%tailp(r%nlevelmax)+nstar_loc
     do i=1,nstar_loc
        star%idp(star%npart-nstar_loc+i)=star%npart_tot+nstar_cum(g%myid-1)+i
     end do
  endif
  deallocate(retire)

  if(nstar_cum(g%ncpu)>0)then
     star%npart_tot=star%npart_tot+nstar_cum(g%ncpu)
     call finish_star_formation(s)
  endif

end subroutine core_retirement

end module core_retirement_module
