module core_sf_support_module
  implicit none
contains

  pure integer function core_top_level(r,g) result(ilev)
    ! Finest level the mesh may reach now: nlevelmax, lowered in
    ! cosmology by aexp_lock_refine as in the refinement flags.
    use amr_commons, only: run_t,global_t
    type(run_t),intent(in)::r
    type(global_t),intent(in)::g
    ilev = r%nlevelmax
    if(r%cosmo.and.r%aexp_lock_refine>0d0) &
         ilev = min(r%nlevelmax, max(g%nlevelmax_part+4, &
                r%nlevelmax+int(log(g%aexp/r%aexp_lock_refine/2.0d0)/log(2d0))))
  end function core_top_level

end module core_sf_support_module
