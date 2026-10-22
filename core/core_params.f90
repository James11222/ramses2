! Namelist parameters of the core module and the routine that reads them.
module core_params_module
  implicit none

  integer(kind=8)::ncoretot=0
  logical::core=.false. ! Cores and their formation activated
  integer::ncoremax=0
  real(kind=8)::core_to_star_mass=100d0
  integer::rho_type_core=4 ! 1: DM, 2: stars, 3: sinks, 4: gas
  real(kind=8)::core_relevance_threshold=2
  real(kind=8)::core_density_threshold=-1 ! Formation density gate in H/cc (proper), like n_star
  ! Clump-finder threshold precedence: _LP_frac, then _nH, then core_clump_density_threshold
  real(kind=8)::core_clump_density_threshold=-1 ! Clump-finder threshold in H/cc; <=0: core_density_threshold/10
  real(kind=8)::core_clump_density_threshold_nH=-1 ! Clump-finder threshold in H/cc at _zref, scaled as (aexp_ref/aexp)**2
  real(kind=8)::core_clump_density_threshold_zref=3 ! Reference redshift of _nH
  real(kind=8)::core_clump_density_threshold_LP_frac=-1 ! Clump-finder threshold as a fraction of the LP density at dx(nlevelmax)
  real(kind=8)::core_clump_density_threshold_LP_T=100 ! Temperature (K) of that LP density
  real(kind=8)::core_clump_density_threshold_LP_mu=1.22 ! Mean molecular weight of that LP density (not mu_mol)
  logical::core_use_LP_criteria=.false. ! Form cores above the LP density instead of core_density_threshold
  real(kind=8)::core_saddle_threshold=-1
  real(kind=8)::core_mass_threshold=0
  real(kind=8)::core_fraction_threshold=2d0
  real(kind=8)::core_mass4_threshold=0d0
  logical::core_use_clump_prop_criteria=.false.
  logical::core_use_mass4_criteria=.false.
  logical::core_use_energy_criteria=.false. ! Require the clump to be gravitationally bound
  logical::core_use_contracting_criteria=.false.
  logical::core_use_contracting_all_axes=.false. ! Require contraction along all principal axes
  logical::core_include_radiation_support=.true.
  logical::core_include_magnetic_support=.true.
  real(kind=8)::core_seed_mass=0d0 ! Seed core mass in Msun; 0: no seed (cores from the initial conditions)
  logical::core_form=.false.
  logical::core_to_star=.false. ! Turn core mass into stars at the finest level; cores reach it with core_refine
  logical::core_refine=.false.
  logical::output_core_clump=.false. ! Write core_clump.* catalogue
  logical::core_merge=.false. ! Merge cores with overlapping accretion zones
  integer::core_merger_search_order=4 ! Core merger search half-width in dx_min
  real(kind=8)::core_merging_timescale=1000d0 ! Core merger age gate in yr
  integer::core_accretion_type = 2 ! 0: None, 2: flux accretion
  real(kind=8)::core_accretion_boost = 1.0d0 ! Boost for core accretion
  integer::core_b_spline_order = 2 ! Core accretion sampling: 2=CIC, 3=TSC, 4=PCS (order**ndim samples)
  logical::core_use_rho_inf = .true. ! Extrapolate density at infinity from Bondi solution
  logical::core_use_bondi_lambda = .true. ! Use Bondi sonic-point lambda factor
  logical::core_use_bondi_for_subsonic_accretion = .false. ! Use Bondi instead of flux rate for subsonic accretion
  logical::core_use_fa_fact_B14 = .false. ! Apply Bleuler+14 flux-accretion density-threshold correction
  logical::core_momentum_conserving = .true. ! Update core momentum with accreted gas
  logical::core_fix_mass = .false. ! Keep core mass fixed during accretion
  logical::core_retire_if_no_acc = .false. ! Retire a core that has stopped accreting
  integer::core_zero_acc_intervals = 2 ! Coarse steps without accretion before a core retires

contains

  subroutine read_core_params(lun,s,star,nml_ok)
    use ramses_commons, only: ramses_t
    implicit none
    integer,intent(in)::lun
    type(ramses_t)::s
    logical,intent(in)::star
    logical,intent(inout)::nml_ok
  namelist/core_params/core,ncoremax,ncoretot,core_form,core_to_star &
       & ,core_refine,output_core_clump,core_merge &
       & ,core_merger_search_order,core_merging_timescale &
       & ,rho_type_core,core_seed_mass,core_to_star_mass &
       & ,core_relevance_threshold,core_density_threshold &
       & ,core_clump_density_threshold,core_use_LP_criteria &
       & ,core_clump_density_threshold_nH,core_clump_density_threshold_zref &
       & ,core_clump_density_threshold_LP_frac &
       & ,core_clump_density_threshold_LP_T,core_clump_density_threshold_LP_mu &
       & ,core_saddle_threshold,core_mass_threshold &
       & ,core_fraction_threshold &
       & ,core_use_clump_prop_criteria,core_use_mass4_criteria &
       & ,core_use_energy_criteria &
       & ,core_use_contracting_criteria,core_use_contracting_all_axes &
       & ,core_include_radiation_support,core_include_magnetic_support &
       & ,core_mass4_threshold
  namelist/core_accretion_params/core_accretion_type &
       & ,core_accretion_boost,core_b_spline_order &
       & ,core_use_rho_inf,core_use_bondi_lambda &
       & ,core_use_bondi_for_subsonic_accretion,core_use_fa_fact_B14 &
       & ,core_momentum_conserving,core_fix_mass &
       & ,core_retire_if_no_acc,core_zero_acc_intervals

  rewind(lun)
  read(lun,NML=core_params,END=1121)
1121 continue
  rewind(lun)
  read(lun,NML=core_accretion_params,END=1122)
1122 continue
  rewind(lun)

  ! Compute maximum number of core particles
  if(ncoremax==0)then
     ncoremax=int(ncoretot/int(s%g%ncpu,kind=8),kind=4)
     if(ncoremax==0)core=.false.
  endif

  if(rho_type_core<1 .or. rho_type_core>4)then
     write(*,*)'Error in the namelist:'
     write(*,*)'rho_type_core should be 1, 2, 3, or 4'
     nml_ok=.false.
  end if
  if(core_form .and. core_seed_mass<=0d0)then
     write(*,*)'Error in the namelist:'
     write(*,*)'Choose core_seed_mass to be >0 before enabling core formation'
     nml_ok=.false.
  end if
  if(core_to_star .and. core_to_star_mass<=0d0)then
     write(*,*)'Error in the namelist:'
     write(*,*)'core_to_star_mass should be >0 when core_to_star is enabled'
     nml_ok=.false.
  end if
  ! LP criteria need a clump-finder threshold.
  if(core_form .and. core_use_LP_criteria .and. &
       core_clump_density_threshold<=0d0 .and. &
       core_clump_density_threshold_nH<=0d0 .and. &
       core_clump_density_threshold_LP_frac<=0d0 .and. &
       core_density_threshold<=0d0)then
     write(*,*)'Error in the namelist:'
     write(*,*)'core_use_LP_criteria needs a collection threshold in H/cc:'
     write(*,*)'set core_clump_density_threshold (or leave'
     write(*,*)'core_density_threshold positive to derive it)'
     nml_ok=.false.
  end if
  if(core_clump_density_threshold_nH>0d0 .and. &
       core_clump_density_threshold_zref<=-1d0)then
     write(*,*)'Error in the namelist:'
     write(*,*)'core_clump_density_threshold_zref should be > -1'
     nml_ok=.false.
  end if
  if(core_clump_density_threshold_LP_frac>0d0)then
     if(core_clump_density_threshold_LP_T<=0d0)then
        write(*,*)'Error in the namelist:'
        write(*,*)'core_clump_density_threshold_LP_T should be > 0'
        nml_ok=.false.
     end if
     if(core_clump_density_threshold_LP_mu<=0d0)then
        write(*,*)'Error in the namelist:'
        write(*,*)'core_clump_density_threshold_LP_mu should be > 0'
        nml_ok=.false.
     end if
  end if

  ! These criteria act only with core_use_clump_prop_criteria or core_use_contracting_all_axes.
  if(core_form.and. &
       .not.(core_use_clump_prop_criteria.or. &
             core_use_contracting_all_axes).and. &
       (core_use_mass4_criteria.or.core_use_energy_criteria.or. &
        core_use_contracting_criteria))then
     if(s%g%myid==1)then
        write(*,*)'Warning: the mass4, energy and contracting criteria have no effect without'
        write(*,*)'core_use_clump_prop_criteria or core_use_contracting_all_axes'
     endif
  endif

  if(core.and.core_accretion_type/=0.and.core_accretion_type/=2)then
     write(*,*)'Error in the namelist:'
     write(*,*)'core_accretion_type should be 0 or 2'
     nml_ok=.false.
  endif
  if(core.and.(core_b_spline_order<2.or.core_b_spline_order>4))then
     write(*,*)'Error in the namelist:'
     write(*,*)'core_b_spline_order should be 2, 3, or 4'
     nml_ok=.false.
  endif
  if(core.and.core_retire_if_no_acc.and..not.star)then
     write(*,*)'Error in the namelist:'
     write(*,*)'core_retire_if_no_acc requires star particles'
     nml_ok=.false.
  endif

  s%r%core=core
  s%r%ncoremax=ncoremax
  s%r%core_to_star_mass=core_to_star_mass

  s%r%rho_type_core=rho_type_core
  s%r%core_relevance_threshold=core_relevance_threshold
  s%r%core_density_threshold=core_density_threshold
  s%r%core_clump_density_threshold=core_clump_density_threshold
  s%r%core_clump_density_threshold_nH=core_clump_density_threshold_nH
  s%r%core_clump_density_threshold_zref=core_clump_density_threshold_zref
  s%r%core_clump_density_threshold_LP_frac=core_clump_density_threshold_LP_frac
  s%r%core_clump_density_threshold_LP_T=core_clump_density_threshold_LP_T
  s%r%core_clump_density_threshold_LP_mu=core_clump_density_threshold_LP_mu
  s%r%core_use_LP_criteria=core_use_LP_criteria
  s%r%core_saddle_threshold=core_saddle_threshold
  s%r%core_mass_threshold=core_mass_threshold
  s%r%core_fraction_threshold=core_fraction_threshold
  s%r%core_mass4_threshold=core_mass4_threshold
  s%r%core_use_clump_prop_criteria=core_use_clump_prop_criteria
  s%r%core_use_mass4_criteria=core_use_mass4_criteria
  s%r%core_use_energy_criteria=core_use_energy_criteria
  s%r%core_use_contracting_criteria=core_use_contracting_criteria
  s%r%core_use_contracting_all_axes=core_use_contracting_all_axes
  s%r%core_include_radiation_support=core_include_radiation_support
  s%r%core_include_magnetic_support=core_include_magnetic_support
  s%r%core_seed_mass = core_seed_mass
  s%r%core_accretion_type = core_accretion_type
  s%r%core_accretion_boost = core_accretion_boost
  s%r%core_b_spline_order = core_b_spline_order
  s%r%core_use_rho_inf = core_use_rho_inf
  s%r%core_use_bondi_lambda = core_use_bondi_lambda
  s%r%core_use_bondi_for_subsonic_accretion = core_use_bondi_for_subsonic_accretion
  s%r%core_use_fa_fact_B14 = core_use_fa_fact_B14
  s%r%core_momentum_conserving = core_momentum_conserving
  s%r%core_fix_mass = core_fix_mass
  s%r%core_retire_if_no_acc = core_retire_if_no_acc
  s%r%core_zero_acc_intervals = core_zero_acc_intervals
  s%r%core_form = core_form
  s%r%core_to_star = core_to_star
  s%r%core_refine = core_refine
  s%r%output_core_clump = output_core_clump
  s%r%core_merge = core_merge
  s%r%core_merger_search_order = core_merger_search_order
  s%r%core_merging_timescale = core_merging_timescale

  end subroutine read_core_params

end module core_params_module
