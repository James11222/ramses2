module core_step_module
  use ramses_commons, only: pst_t
  use rho_fine_module, only: m_rho_fine
  use core_formation_module, only: m_core_formation,m_core_formation_site
  use core_evolution_module, only: r_core_evolution,in_core_acc_t,out_core_acc_t, &
       CORE_ACC_PREPARE,CORE_ACC_FINALIZE
  use core_retirement_module, only: r_core_retirement,in_core_retire_t,out_core_retire_t
  use core_merger_module, only: r_core_merger
  implicit none
  private

  integer,parameter,public :: CORE_PHASE_RETIRE_BOUNDARY=1
  integer,parameter,public :: CORE_PHASE_FORMATION=2
  integer,parameter,public :: CORE_PHASE_ACC_PREPARE=3
  integer,parameter,public :: CORE_PHASE_ACC_FINALIZE=4
  integer,parameter,public :: CORE_PHASE_MERGER=5

  ! Accretion input carried from ACC_PREPARE to ACC_FINALIZE
  type,public :: core_step_state_t
     private
     type(in_core_acc_t) :: input_acc
  end type core_step_state_t

  public :: m_core_step

contains

subroutine m_core_step(pst,state,phase,ilevel,icount)
  type(pst_t),intent(inout)::pst
  type(core_step_state_t),intent(inout)::state
  integer,intent(in)::phase
  integer,intent(in),optional::ilevel,icount
  integer::level,count

  level=0
  count=0
  if(present(ilevel))level=ilevel
  if(present(icount))count=icount

  select case(phase)
  case(CORE_PHASE_RETIRE_BOUNDARY)
     call core_retirement_phase(pst,level,count)
  case(CORE_PHASE_FORMATION)
     call core_formation_phase(pst,level)
  case(CORE_PHASE_ACC_PREPARE)
     call core_accretion_prepare(pst,state,level)
  case(CORE_PHASE_ACC_FINALIZE)
     call core_accretion_finalize(pst,state)
  case(CORE_PHASE_MERGER)
     call core_merger_phase(pst,level)
  end select
end subroutine m_core_step

subroutine core_retirement_phase(pst,ilevel,icount)
  type(pst_t),intent(inout)::pst
  integer,intent(in)::ilevel,icount
  type(in_core_retire_t)::input
  type(out_core_retire_t)::output
  ! Turn the cores that stopped accreting into stars, once per coarse step

  associate(r=>pst%s%r,g=>pst%s%g)
  if(.not.r%core.or..not.r%core_retire_if_no_acc)return
  if(ilevel/=r%levelmin.or.icount/=1)return
  call m_timer('core - retirement','start')
  input%nstep_coarse=g%nstep_coarse
  call r_core_retirement(pst,input,storage_size(input)/32, &
       output,storage_size(output)/32)
  if(output%mass>0d0)then
     g%mass_star_tot=g%mass_star_tot+output%mass
     g%mass_sink_tot=g%mass_sink_tot-output%mass
  endif
  end associate
end subroutine core_retirement_phase

subroutine core_formation_phase(pst,ilevel)
  type(pst_t),intent(inout)::pst
  integer,intent(in)::ilevel

  associate(r=>pst%s%r)
  if(.not.r%core.or.ilevel/=r%levelmin.or..not.r%core_form)return
  call m_timer('core - formation','start')
  call m_core_formation_site(pst)
  call m_core_formation(pst)
#ifdef GRAV
  ! Restore m%rho and m%nref, used by the clump finder, from the debited gas
  ! and the new cores, and sort the new cores into their level buckets
  call m_rho_fine(pst,r%levelmin,0)
#endif
  end associate
end subroutine core_formation_phase

subroutine core_merger_phase(pst,ilevel)
  type(pst_t),intent(inout)::pst
  integer,intent(in)::ilevel

  associate(r=>pst%s%r)
  if(.not.r%core.or..not.r%core_merge)return
  call m_timer('core - merger','start')
  call r_core_merger(pst,ilevel,1)
  end associate
end subroutine core_merger_phase

subroutine core_accretion_prepare(pst,state,ilevel)
  type(pst_t),intent(inout)::pst
  type(core_step_state_t),intent(inout)::state
  integer,intent(in)::ilevel
  type(out_core_acc_t)::output
  ! Build the accretion request of the cores of this level before the hydro step

  associate(r=>pst%s%r,g=>pst%s%g)
  if(.not.r%core.or.r%core_accretion_type<=0)return
  call m_timer('core - prepare','start')
  state%input_acc%action=CORE_ACC_PREPARE
  state%input_acc%ilevel=ilevel
  state%input_acc%nstep_coarse=g%nstep_coarse
  call r_core_evolution(pst,state%input_acc,storage_size(state%input_acc)/32, &
       output,storage_size(output)/32)
  end associate
end subroutine core_accretion_prepare

subroutine core_accretion_finalize(pst,state)
  type(pst_t),intent(inout)::pst
  type(core_step_state_t),intent(inout)::state
  type(out_core_acc_t)::output
  ! Apply the accretion after cooling

  associate(r=>pst%s%r,g=>pst%s%g)
  if(.not.r%core.or.r%core_accretion_type<=0)return
  call m_timer('core - accretion','start')
  state%input_acc%action=CORE_ACC_FINALIZE
  call r_core_evolution(pst,state%input_acc,storage_size(state%input_acc)/32, &
       output,storage_size(output)/32)
  g%mass_sink_tot=g%mass_sink_tot+output%mass
  end associate
end subroutine core_accretion_finalize

end module core_step_module
