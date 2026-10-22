#ifndef HYDRO
#error "CORE=1 requires HYDRO=1"
#endif
#if NDIM != 3
#error "CORE=1 requires NDIM=3"
#endif
! Core accretion: ACC_PREPARE builds the request of each core before the
! hydro step, ACC_FINALIZE caps it with the gas after cooling and applies it.
module core_evolution_module
  use amr_parameters, only: ndim, twotondim, dp
  use hydro_parameters, only: nvar, nener
  use sink_evolution_module, only: sink_B_spline_weights_CIC, sink_B_spline_weights_TSC, sink_B_spline_weights_PCS

  integer,parameter::CORE_ACC_PREPARE=1
  integer,parameter::CORE_ACC_FINALIZE=2

  type in_core_acc_t
     sequence
     integer::action=0
     integer::ilevel=0
     integer::nstep_coarse=0
  end type in_core_acc_t

  type out_core_acc_t
     sequence
     real(kind=8)::mass=0d0
  end type out_core_acc_t

  ! One gas cell a core accretes from
  type core_acc_contrib_t
     integer(kind=8),dimension(0:ndim)::key
     integer::ind=0
     integer::ipart=0
     real(kind=8)::d_req=0d0
     real(kind=8),dimension(1:nvar)::q=0d0
     real(kind=8),dimension(1:ndim)::x_rel=0d0
     real(kind=8),dimension(1:ndim)::v_rel=0d0
  end type core_acc_contrib_t

  ! Summed requests and owner scales per cell, on the mesh and cache grids
  real(dp),allocatable,dimension(:,:)::core_acc_req,core_acc_scale

  ! Owner cells touched by this call
  integer,allocatable,dimension(:)::core_owner_touch_igrid,core_owner_touch_ind
  integer::core_owner_touch_n=0

  ! Gas cells the cores accrete from, built by PREPARE
  type(core_acc_contrib_t),allocatable,dimension(:)::core_acc_contrib
  integer::core_acc_ncontrib=0

contains
  !##############################################################################
  subroutine append_owner_touch(igrid,ind)
    implicit none
    integer,intent(in)::igrid,ind
    integer,allocatable,dimension(:)::tmp_i,tmp_j
    integer::n
    if(.not.allocated(core_owner_touch_igrid))then
       allocate(core_owner_touch_igrid(1024),core_owner_touch_ind(1024))
    endif
    if(core_owner_touch_n==size(core_owner_touch_igrid))then
       n=core_owner_touch_n
       allocate(tmp_i(2*n),tmp_j(2*n))
       tmp_i(1:n)=core_owner_touch_igrid(1:n)
       tmp_j(1:n)=core_owner_touch_ind(1:n)
       call move_alloc(tmp_i,core_owner_touch_igrid)
       call move_alloc(tmp_j,core_owner_touch_ind)
    endif
    core_owner_touch_n=core_owner_touch_n+1
    core_owner_touch_igrid(core_owner_touch_n)=igrid
    core_owner_touch_ind(core_owner_touch_n)=ind
  end subroutine append_owner_touch
  !##############################################################################
  subroutine append_contrib(item)
    implicit none
    type(core_acc_contrib_t),intent(in)::item
    type(core_acc_contrib_t),allocatable::tmp(:)
    if(.not.allocated(core_acc_contrib))allocate(core_acc_contrib(1024))
    if(core_acc_ncontrib==size(core_acc_contrib))then
       allocate(tmp(2*core_acc_ncontrib))
       tmp(1:core_acc_ncontrib)=core_acc_contrib(1:core_acc_ncontrib)
       call move_alloc(tmp,core_acc_contrib)
    endif
    core_acc_ncontrib=core_acc_ncontrib+1
    core_acc_contrib(core_acc_ncontrib)=item
  end subroutine append_contrib
  !##############################################################################
  subroutine init_flush_core_req(mesh,igrid,hash_key)
    use amr_commons, only: mesh_t
    implicit none
    type(mesh_t)::mesh
    integer::igrid
    integer(kind=8),dimension(0:ndim)::hash_key
    mesh%grid(igrid)%lev=int(hash_key(0),kind=4)
    mesh%grid(igrid)%ckey(1:ndim)=int(hash_key(1:ndim),kind=4)
    core_acc_req(1:twotondim,igrid)=0d0
  end subroutine init_flush_core_req
  !##############################################################################
  subroutine pack_flush_core_req(mesh,igrid,msg_size,msg_array)
    use amr_commons, only: mesh_t
    use cache_commons, only: msg_small_realdp
    implicit none
    type(mesh_t)::mesh
    integer::igrid,msg_size
    integer,dimension(1:msg_size),optional::msg_array
    type(msg_small_realdp)::msg
    msg%realdp(1:twotondim)=core_acc_req(1:twotondim,igrid)
    msg_array=transfer(msg,msg_array)
  end subroutine pack_flush_core_req
  !##############################################################################
  subroutine unpack_flush_core_req(mesh,igrid,msg_size,msg_array,hash_key)
    use amr_commons, only: mesh_t
    use cache_commons, only: msg_small_realdp
    implicit none
    type(mesh_t)::mesh
    integer::igrid,msg_size
    integer,dimension(1:msg_size),optional::msg_array
    integer(kind=8),dimension(0:ndim)::hash_key
    type(msg_small_realdp)::msg
    integer::ind
    mesh%grid(igrid)%lev=int(hash_key(0),kind=4)
    mesh%grid(igrid)%ckey(1:ndim)=int(hash_key(1:ndim),kind=4)
    msg=transfer(msg_array,msg)
    do ind=1,twotondim
       if(msg%realdp(ind)<=0d0)cycle
       if(core_acc_req(ind,igrid)==0d0)call append_owner_touch(igrid,ind)
       core_acc_req(ind,igrid)=core_acc_req(ind,igrid)+msg%realdp(ind)
    end do
  end subroutine unpack_flush_core_req
  !##############################################################################
  subroutine pack_fetch_core_scale(mesh,igrid,msg_size,msg_array)
    use amr_commons, only: mesh_t
    use cache_commons, only: msg_small_realdp
    implicit none
    type(mesh_t)::mesh
    integer::igrid,msg_size
    integer,dimension(1:msg_size),optional::msg_array
    type(msg_small_realdp)::msg
    msg%realdp(1:twotondim)=core_acc_scale(1:twotondim,igrid)
    msg_array=transfer(msg,msg_array)
  end subroutine pack_fetch_core_scale
  !##############################################################################
  subroutine unpack_fetch_core_scale(mesh,igrid,msg_size,msg_array,hash_key)
    use amr_commons, only: mesh_t
    use cache_commons, only: msg_small_realdp
    implicit none
    type(mesh_t)::mesh
    integer::igrid,msg_size
    integer,dimension(1:msg_size),optional::msg_array
    integer(kind=8),dimension(0:ndim)::hash_key
    type(msg_small_realdp)::msg
    mesh%grid(igrid)%lev=int(hash_key(0),kind=4)
    mesh%grid(igrid)%ckey(1:ndim)=int(hash_key(1:ndim),kind=4)
    msg=transfer(msg_array,msg)
    core_acc_scale(1:twotondim,igrid)=msg%realdp(1:twotondim)
  end subroutine unpack_fetch_core_scale
  !##############################################################################
  !##############################################################################
  !##############################################################################
  !##############################################################################
  recursive subroutine r_core_evolution(pst,input,input_size,output,output_size)
    use mdl_module
    use ramses_commons, only: pst_t
    use mdl_parameters
    implicit none
    type(pst_t)::pst
    integer,VALUE::input_size
    integer::output_size
    type(in_core_acc_t)::input
    type(out_core_acc_t)::output,next_output
    integer::rID

    if(pst%nLower>0)then
       rID = mdl_send_request(pst%s%mdl,MDL_CORE_EVOLUTION,pst%iUpper+1,input_size,output_size,input)
       call r_core_evolution(pst%pLower,input,input_size,output,output_size)
       call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
       output%mass=output%mass+next_output%mass
    else
       output%mass=0d0
       select case(input%action)
       case(CORE_ACC_PREPARE)
          call core_prepare(pst%s,pst%s%core,input)
       case(CORE_ACC_FINALIZE)
          call core_finalize(pst%s,pst%s%core,input,output)
       end select
    endif

  end subroutine r_core_evolution
  !##############################################################################
  !##############################################################################
  !##############################################################################
  !##############################################################################
  subroutine core_prepare(s,p,input)
    use constants, only: twopi
    use ramses_commons, only: ramses_t
    use pm_commons, only: part_t
    use cache_commons
    use cache
    use marshal, only: pack_fetch_refine,unpack_fetch_refine
    use boundaries, only: init_bound_refine
    implicit none
    type(ramses_t)::s
    type(part_t)::p
    type(in_core_acc_t),intent(in)::input
    !------------------------------------------------------------------
    ! Build the accretion request of every core of this level from the
    ! gas before the hydro step.
    !------------------------------------------------------------------
    integer::ilevel,ipart,n_acc_sample
    real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,factG
    real(kind=8)::dx_loc,vol_loc,lambda_sonic
    type(msg_large_realdp)::dummy_large_realdp

    associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

    ilevel=input%ilevel
    core_acc_ncontrib=0

    call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

    factG=1d0
    if(r%cosmo)factG=3d0/4d0/twopi*g%omega_m*g%aexp

    ! Bondi sonic constant
    if(r%core_use_bondi_lambda)then
       if(abs(r%gamma - 1d0).le.0.01d0)then
          lambda_sonic = 0.25d0*exp(1.5d0)
       else if(abs(r%gamma - 5d0/3d0).le.0.01d0)then
          lambda_sonic = 0.25d0
       else
          lambda_sonic = 0.5d0**((r%gamma + 1d0)/(2d0*(r%gamma - 1d0))) &
               & * (0.25d0*(5d0-3d0*r%gamma))**(-(5d0-3d0*r%gamma)/(2d0*(r%gamma - 1d0)))
       end if
    else
       lambda_sonic = 1d0
    end if

    dx_loc=r%boxlen/2**ilevel
    vol_loc=dx_loc**ndim
    n_acc_sample = int(r%core_b_spline_order**ndim)

    call open_cache(mdl, m, pack_size=storage_size(dummy_large_realdp)/32, &
         pack=pack_fetch_refine, unpack=unpack_fetch_refine, bound=init_bound_refine)
    do ipart = 1, p%npart
       if(p%mp(ipart)<=0d0)cycle
       if(p%levelp(ipart)/=ilevel)cycle
       call core_flux_accretion_request(s,p,ilevel,ipart,dx_loc,vol_loc,n_acc_sample, &
            scale_nH,factG,lambda_sonic)
    end do
    call close_cache(mdl)

    end associate
  end subroutine core_prepare
  !##############################################################################
  !##############################################################################
  !##############################################################################
  !##############################################################################
  subroutine core_finalize(s,p,input,output)
    use ramses_commons, only: ramses_t
    use pm_commons, only: part_t, cross
    use nbors_utils
    use cache_commons
    use cache
    use marshal, only: pack_fetch_refine,unpack_fetch_refine
    implicit none
    type(ramses_t),intent(inout)::s
    type(part_t),intent(inout)::p
    type(in_core_acc_t),intent(in)::input
    type(out_core_acc_t),intent(inout)::output
    !------------------------------------------------------------------
    ! Cap the requests with the gas after cooling, scale them so that the
    ! cores sharing a cell take at most 0.75(rho-smallr) of it, and move
    ! the accreted gas from uold to the cores.
    !------------------------------------------------------------------
    type(msg_large_realdp)::dummy_large_realdp
    type(msg_small_realdp)::dummy_small_realdp
    integer::icontrib,igridn,ind,idim,itouch,ipart,n
#if NENER>0
    integer::irad
    real(kind=8),dimension(1:nener)::erad
#endif
    real(kind=8)::dx_loc,vol_loc,rho,cap,req,scale,d_acc,u_acc,mold,mnew
#ifdef MHD
    real(kind=8)::bx,by,bz,emag
#endif
    real(kind=8),allocatable,dimension(:)::macc,uacc,zacc
    real(kind=8),allocatable,dimension(:,:)::xacc,pacc,lacc

    associate(r=>s%r,m=>s%m,mdl=>s%mdl)
    dx_loc=r%boxlen/2**input%ilevel; vol_loc=dx_loc**ndim
    if(.not.allocated(core_acc_req))then
       allocate(core_acc_req(1:twotondim,1:m%ngridmax+m%ncachemax))
       allocate(core_acc_scale(1:twotondim,1:m%ngridmax+m%ncachemax))
       core_acc_req=0d0; core_acc_scale=1d0
    endif

    ! Cap each request with the gas of its cell
    call open_cache(mdl,m,pack_size=storage_size(dummy_large_realdp)/32, &
         pack=pack_fetch_refine,unpack=unpack_fetch_refine)
    do icontrib=1,core_acc_ncontrib
       call get_grid(s,core_acc_contrib(icontrib)%key,igridn,flush_cache=.false.,fetch_cache=.true.)
       ind=core_acc_contrib(icontrib)%ind
       rho=m%uold(ind,1,igridn)
       core_acc_contrib(icontrib)%q(1:nvar)=m%uold(ind,1:nvar,igridn)/rho
       ! Thermal energy only: magnetic and non-thermal energies stay on the grid
#ifdef MHD
       bx=0.5d0*(m%bold(ind,1,igridn)+m%bold(ind,4,igridn))
       by=0.5d0*(m%bold(ind,2,igridn)+m%bold(ind,5,igridn))
       bz=0.5d0*(m%bold(ind,3,igridn)+m%bold(ind,6,igridn))
       emag=0.5d0*(bx**2+by**2+bz**2)
       core_acc_contrib(icontrib)%q(5)=core_acc_contrib(icontrib)%q(5)-emag/rho
#endif
#if NENER>0
       do irad=1,nener
          core_acc_contrib(icontrib)%q(5)=core_acc_contrib(icontrib)%q(5)-core_acc_contrib(icontrib)%q(5+irad)
       end do
#endif
       core_acc_contrib(icontrib)%v_rel(1:ndim)= &
            core_acc_contrib(icontrib)%q(2:ndim+1)-p%vp(core_acc_contrib(icontrib)%ipart,1:ndim)
       cap=0.75d0*max(rho-r%smallr,0d0)
       core_acc_contrib(icontrib)%d_req=min(core_acc_contrib(icontrib)%d_req,cap)
    enddo
    call close_cache(mdl)

    ! Sum the capped requests on the owner cells
    call open_cache(mdl,m,pack_size=storage_size(dummy_small_realdp)/32, &
         init=init_flush_core_req,flush=pack_flush_core_req,combine=unpack_flush_core_req)
    do icontrib=1,core_acc_ncontrib
       if(core_acc_contrib(icontrib)%d_req<=0d0)cycle
       call get_grid(s,core_acc_contrib(icontrib)%key,igridn,flush_cache=.true.,fetch_cache=.false.)
       ind=core_acc_contrib(icontrib)%ind
       if(core_acc_req(ind,igridn)==0d0.and.igridn<=m%ngridmax)call append_owner_touch(igridn,ind)
       core_acc_req(ind,igridn)=core_acc_req(ind,igridn)+core_acc_contrib(icontrib)%d_req
    enddo
    call close_cache(mdl)

    ! Scale the requests of each owner cell to the shared cap and remove the gas
    do itouch=1,core_owner_touch_n
       igridn=core_owner_touch_igrid(itouch); ind=core_owner_touch_ind(itouch)
       rho=m%uold(ind,1,igridn); req=core_acc_req(ind,igridn)
       cap=0.75d0*max(rho-r%smallr,0d0)
       scale=min(1d0,cap/req)
       core_acc_scale(ind,igridn)=scale
       d_acc=scale*req
#if NENER>0
       erad(1:nener)=m%uold(ind,6:5+nener,igridn)
#endif
       m%uold(ind,1:nvar,igridn)=m%uold(ind,1:nvar,igridn)*(1d0-d_acc/rho)
       ! Put the magnetic and non-thermal energies back into the cell
#ifdef MHD
       bx=0.5d0*(m%bold(ind,1,igridn)+m%bold(ind,4,igridn))
       by=0.5d0*(m%bold(ind,2,igridn)+m%bold(ind,5,igridn))
       bz=0.5d0*(m%bold(ind,3,igridn)+m%bold(ind,6,igridn))
       emag=0.5d0*(bx**2+by**2+bz**2)
       m%uold(ind,5,igridn)=m%uold(ind,5,igridn)+d_acc/rho*emag
#endif
#if NENER>0
       do irad=1,nener
          m%uold(ind,5,igridn)=m%uold(ind,5,igridn)+d_acc/rho*erad(irad)
          m%uold(ind,5+irad,igridn)=erad(irad)
       end do
#endif
    enddo

    ! Fetch the owner scales and sum what each core accretes
    call open_cache(mdl,m,pack_size=storage_size(dummy_small_realdp)/32, &
         pack=pack_fetch_core_scale,unpack=unpack_fetch_core_scale)
    n=max(1,p%npart)
    allocate(macc(n),uacc(n),zacc(n),xacc(ndim,n),pacc(ndim,n),lacc(ndim,n))
    macc=0d0; uacc=0d0; zacc=0d0; xacc=0d0; pacc=0d0; lacc=0d0
    do icontrib=1,core_acc_ncontrib
       call get_grid(s,core_acc_contrib(icontrib)%key,igridn,flush_cache=.false.,fetch_cache=.true.)
       ind=core_acc_contrib(icontrib)%ind
       d_acc=core_acc_scale(ind,igridn)*core_acc_contrib(icontrib)%d_req
       if(d_acc<=0d0)cycle
       ipart=core_acc_contrib(icontrib)%ipart
       macc(ipart)=macc(ipart)+d_acc*vol_loc
       u_acc=core_acc_contrib(icontrib)%q(5)-0.5d0*sum(core_acc_contrib(icontrib)%q(2:ndim+1)**2)
       u_acc=max(u_acc,r%smallc**2/(r%gamma-1d0))
       uacc(ipart)=uacc(ipart)+d_acc*vol_loc*u_acc
       if(r%metal)zacc(ipart)=zacc(ipart)+d_acc*vol_loc*core_acc_contrib(icontrib)%q(r%imetal)
       if(r%core_momentum_conserving)then
          xacc(:,ipart)=xacc(:,ipart)+d_acc*core_acc_contrib(icontrib)%x_rel*vol_loc*dx_loc
          pacc(:,ipart)=pacc(:,ipart)+d_acc*core_acc_contrib(icontrib)%v_rel*vol_loc
          lacc(:,ipart)=lacc(:,ipart)+ &
               d_acc*cross(core_acc_contrib(icontrib)%x_rel,core_acc_contrib(icontrib)%v_rel)*vol_loc*dx_loc
       endif
    enddo
    call close_cache(mdl)

    ! Add the accreted gas to the cores
    do ipart=1,p%npart
       if(macc(ipart)==0d0)cycle
       mold=p%mp(ipart); mnew=mold+macc(ipart)
       if(.not.p%static)then
          p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+xacc(:,ipart)/mnew
          do idim=1,ndim
             if(r%periodic(idim))then
                if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
                if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
             endif
          end do
       endif
       p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+pacc(:,ipart)/mnew
       p%jp(ipart,1:ndim)=(mold*p%jp(ipart,1:ndim)+lacc(:,ipart))/mnew
       p%up(ipart)=(mold*p%up(ipart)+uacc(ipart))/mnew
       if(r%metal)p%zp(ipart)=(mold*p%zp(ipart)+zacc(ipart))/mnew
       ! core_fix_mass keeps the core mass; the gas is still removed
       if(.not.r%core_fix_mass)p%mp(ipart)=mnew
       ! Coarse step of this accretion
       p%nacc(ipart)=input%nstep_coarse
    enddo
    output%mass=sum(macc)

    ! Leave the owner cells ready for the next call
    do itouch=1,core_owner_touch_n
       igridn=core_owner_touch_igrid(itouch); ind=core_owner_touch_ind(itouch)
       core_acc_req(ind,igridn)=0d0; core_acc_scale(ind,igridn)=1d0
    enddo
    core_owner_touch_n=0

    end associate
  end subroutine core_finalize

  !##############################################################################
  !##############################################################################
  !##############################################################################
  !##############################################################################
  subroutine core_flux_accretion_request(s,p,ilevel,ipart,dx_loc,vol_loc,n_acc_sample, &
       scale_nH,factG,lambda_sonic)
    use constants, only: pi
    use ramses_commons, only: ramses_t
    use pm_commons, only: part_t
    use nbors_utils
    implicit none
    type(ramses_t)::s
    type(part_t)::p
    integer::ilevel,ipart,n_acc_sample
    real(kind=8)::dx_loc,vol_loc,scale_nH,factG,lambda_sonic
    !------------------------------------------------------------------
    ! Flux accretion rate of one core (Bleuler et al. 2014), with the
    ! Bondi rate below the sonic radius, spread over the B-spline cells.
    !------------------------------------------------------------------
    real(kind=8),dimension(1:ndim,1:n_acc_sample)::x_acc_sample
    integer,dimension(1:ndim,1:n_acc_sample)::ckeynei
    integer,dimension(1:ndim)::ckey_div
    real(kind=8),dimension(1:n_acc_sample)::vol
    real(kind=8),dimension(1:ndim)::xcen,xnei,x_rel,vv,vel_gas
    real(kind=8),dimension(1:nvar)::fluid_var
    integer(kind=8),dimension(0:ndim)::hash_nbor
    integer::j,idim,icelln,igridn
#if NENER>0
    integer::irad
#endif
    real(kind=8)::weight,weight_sum,eff_cells,acc_volume,r_acc
    real(kind=8)::rho_sum,eth_sum,ethermal,rho_gas,c2,v2,v_bondi,boost,d_star
    real(kind=8)::m_gas,r2_core,rho_inf,dMdt_flux,dMBondi_overdt,fa_fact
    real(kind=8)::dMdt_cap,sonic_radius
    real(kind=8)::d,e,d_acc,div_cell,div_sum,div_right,div_left
    type(core_acc_contrib_t)::item
#ifdef MHD
    real(kind=8)::bx,by,bz,emag
#endif

    associate(r=>s%r,g=>s%g,m=>s%m)

    hash_nbor(0) = ilevel+1
    x_acc_sample=0d0; ckeynei=0; vol=0d0

    xcen(1:ndim) = (p%xp(ipart,1:ndim)+m%skip(1:ndim)) / dx_loc

    if      (r%core_b_spline_order==2)then
       call sink_B_spline_weights_CIC(s,xcen(1:ndim),x_acc_sample,ckeynei,vol,ilevel)
    else if (r%core_b_spline_order==3)then
       call sink_B_spline_weights_TSC(s,xcen(1:ndim),x_acc_sample,ckeynei,vol,ilevel)
    else if (r%core_b_spline_order==4)then
       call sink_B_spline_weights_PCS(s,xcen(1:ndim),x_acc_sample,ckeynei,vol,ilevel)
    end if

    ! No accretion from cells outside the box
    do j = 1, n_acc_sample
       do idim = 1, ndim
          if(.not.r%periodic(idim))then
             if(ckeynei(idim,j)<m%box_ckey_min(idim,ilevel+1))vol(j)=0d0
             if(ckeynei(idim,j)>=m%box_ckey_max(idim,ilevel+1))vol(j)=0d0
          endif
       end do
    end do

    rho_sum=0d0; eth_sum=0d0; vel_gas=0d0; div_sum=0d0; weight_sum=0d0

    do j = 1, n_acc_sample
       weight = vol(j)
       if(weight<=0d0)cycle
       hash_nbor(1:ndim) = int(ckeynei(1:ndim,j),kind=8)
       call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.false.,fetch_cache=.true.)
       if(igridn==0)cycle
       fluid_var(1:nvar) = m%uold(icelln,1:nvar,igridn)
       weight_sum = weight_sum + weight

       d = max(dble(fluid_var(1)),r%smallr)
       vv(1:ndim) = fluid_var(2:ndim+1)/d
       e = fluid_var(5)/d
       ! Thermal energy only: magnetic and non-thermal energies stay on the grid
#ifdef MHD
       bx=0.5d0*(m%bold(icelln,1,igridn) + m%bold(icelln,4,igridn))
       by=0.5d0*(m%bold(icelln,2,igridn) + m%bold(icelln,5,igridn))
       bz=0.5d0*(m%bold(icelln,3,igridn) + m%bold(icelln,6,igridn))
       emag=0.5d0*(bx**2+by**2+bz**2)
       e = e - emag/d
#endif
#if NENER>0
       do irad=1,nener
          e = e - fluid_var(5+irad)/d
       end do
#endif
       ethermal = e - 0.5d0*sum(vv(:)**2)

       rho_sum = rho_sum + weight*d
       vel_gas(1:ndim) = vel_gas(1:ndim) + weight*d*vv(1:ndim)
       eth_sum = eth_sum + weight*d*ethermal

       div_cell = 0d0
       do idim = 1, ndim
          ! Wrap the shifted key across periodic faces
          ckey_div(1:ndim) = ckeynei(1:ndim,j)
          ckey_div(idim) = ckey_div(idim) + 1
          if(r%periodic(idim))then
             if(ckey_div(idim)>=m%box_ckey_max(idim,ilevel+1)) &
                  ckey_div(idim)=ckey_div(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
          endif
          hash_nbor(1:ndim) = int(ckey_div(1:ndim),kind=8)
          call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.false.,fetch_cache=.true.)
          if(igridn/=0)then
             div_right = (m%uold(icelln,1+idim,igridn) - max(dble(m%uold(icelln,1,igridn)),r%smallr)*p%vp(ipart,idim))/dx_loc
          else
             div_right = 0d0
          end if

          ckey_div(1:ndim) = ckeynei(1:ndim,j)
          ckey_div(idim) = ckey_div(idim) - 1
          if(r%periodic(idim))then
             if(ckey_div(idim)< m%box_ckey_min(idim,ilevel+1)) &
                  ckey_div(idim)=ckey_div(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
          endif
          hash_nbor(1:ndim) = int(ckey_div(1:ndim),kind=8)
          call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.false.,fetch_cache=.true.)
          if(igridn/=0)then
             div_left = (m%uold(icelln,1+idim,igridn) - max(dble(m%uold(icelln,1,igridn)),r%smallr)*p%vp(ipart,idim))/dx_loc
          else
             div_left = 0d0
          end if
          ! Centred difference; 1/dx is already in div_right and div_left
          div_cell = div_cell + (div_right - div_left)/2.0d0
       end do
       div_sum = div_sum + weight*div_cell
    end do

    if(weight_sum<=0d0)return

    rho_gas = rho_sum/weight_sum
    vel_gas(1:ndim) = vel_gas(1:ndim)/rho_sum
    ethermal = eth_sum/rho_sum
    ! Ideal-gas sound speed squared
    c2 = max(r%gamma*(r%gamma-1d0)*ethermal,r%smallc**2)
    v2 = sum((vel_gas(1:ndim)-p%vp(ipart,1:ndim))**2)

    if(r%core_accretion_boost<0d0 .and. r%n_star>0d0)then
       d_star = r%n_star/scale_nH
       boost = max((rho_gas/d_star)**2,1d0)
    else
       boost = max(abs(r%core_accretion_boost),1d0)
    end if
    v_bondi = sqrt(c2 + v2) * boost**(-1d0/3d0)

    eff_cells = weight_sum*dble(n_acc_sample)
    acc_volume = eff_cells*vol_loc
    m_gas = rho_gas*acc_volume
    r_acc = 0.5d0*dble(r%core_b_spline_order)*dx_loc
    r2_core = (factG*(p%mp(ipart)+m_gas)/v_bondi**2)**2

    if(r%core_use_rho_inf)then
       rho_inf = rho_gas/(bondi_alpha(r_acc/sqrt(r2_core + tiny(0d0))))
    else
       rho_inf = rho_gas
    end if
    dMBondi_overdt = 4d0*pi*rho_inf*r2_core*v_bondi*lambda_sonic

    dMdt_flux = -(div_sum/weight_sum)*acc_volume
    if(r%core_use_fa_fact_B14.and.r%core_density_threshold>0d0)then
       ! The threshold is in H/cc and rho_gas in code units
       fa_fact = 0.1d0*(log10(rho_gas) &
            - log10(r%core_density_threshold/scale_nH)) + 1d0
       dMdt_flux = dMdt_flux*fa_fact
    end if
    sonic_radius = 0.5d0*factG*p%mp(ipart)/c2
    if(r%core_use_bondi_for_subsonic_accretion)then
       if(sonic_radius < r_acc)then
          dMdt_flux = dMBondi_overdt
       end if
    end if

    dMdt_flux = max(0d0,dMdt_flux)
    dMdt_cap = 0.75d0*m_gas/g%dtnew(ilevel)
    dMdt_flux = min(dMdt_flux,dMdt_cap)
    if(dMdt_flux<=0d0)return

    do j = 1, n_acc_sample
       xnei(1:ndim) = x_acc_sample(1:ndim,j)
       x_rel(1:ndim) = xnei(1:ndim) - xcen(1:ndim)
       hash_nbor(1:ndim) = int(ckeynei(1:ndim,j),kind=8)
       call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.false.,fetch_cache=.true.)
       if(igridn==0)cycle
       if(m%grid(igridn)%refined(icelln))cycle
       weight = vol(j)
       if(weight<=0d0)cycle

       d = max(dble(m%uold(icelln,1,igridn)),r%smallr)
       d_acc = dMdt_flux*g%dtnew(ilevel)*(weight/weight_sum)/vol_loc
       d_acc = d_acc*d/rho_gas
       if(d_acc<=0d0)cycle

       ! Record the cell by its parent oct
       item%key(0)=hash_nbor(0)-1
       item%key(1:ndim)=hash_nbor(1:ndim)/2
       item%ind=icelln
       item%ipart=ipart
       item%d_req=d_acc
       item%x_rel(1:ndim)=x_rel(1:ndim)
       call append_contrib(item)
    end do

    end associate
  contains
    ! Bondi density factor alpha(x)
    real(kind=8) function bondi_alpha(x)
      implicit none
      real(kind=8)::x
      real(kind=8), parameter :: XMIN=0.01d0, XMAX=2.0d0
      integer, parameter :: NTABLE=51
      real(kind=8)::lambda_c, xtable, xtablep1, alpha_exp
      integer::idx
      ! Table of alpha values for x from 0.01 to 2.0 with uniform
      ! logarithmic spacing
      real(kind=8), parameter, dimension(NTABLE) :: alphatable = (/ &
           820.254d0, 701.882d0, 600.752d0, 514.341d0, 440.497d0, 377.381d0, 323.427d0, &
           277.295d0, 237.845d0, 204.1d0, 175.23d0, 150.524d0, 129.377d0, 111.27d0, 95.7613d0, &
           82.4745d0, 71.0869d0, 61.3237d0, 52.9498d0, 45.7644d0, 39.5963d0, 34.2989d0, &
           29.7471d0, 25.8338d0, 22.4676d0, 19.5705d0, 17.0755d0, 14.9254d0, 13.0714d0, &
           11.4717d0, 10.0903d0, 8.89675d0, 7.86467d0, 6.97159d0, 6.19825d0, 5.52812d0, &
           4.94699d0, 4.44279d0, 4.00497d0, 3.6246d0, 3.29395d0, 3.00637d0, 2.75612d0, &
           2.53827d0, 2.34854d0, 2.18322d0, 2.03912d0, 1.91344d0, 1.80378d0, 1.70804d0, &
           1.62439d0 /)

      lambda_c = exp(1.5d0) / 4d0
      if (x .le. XMIN) then
         bondi_alpha = lambda_c / sqrt(2d0 * x**ndim)
      else if (x .ge. XMAX) then
         bondi_alpha = exp(1d0/x)
      else
         idx = floor((NTABLE-1) * log(x/XMIN) / log(XMAX/XMIN))
         xtable = exp(log(XMIN) + idx*log(XMAX/XMIN)/(NTABLE-1))
         xtablep1 = exp(log(XMIN) + (idx+1)*log(XMAX/XMIN)/(NTABLE-1d0))
         alpha_exp = log(x/xtable) / log(xtablep1/xtable)
         bondi_alpha = alphatable(idx+1) * (alphatable(idx+2)/alphatable(idx+1))**alpha_exp
      end if
    end function bondi_alpha

  end subroutine core_flux_accretion_request

end module core_evolution_module
