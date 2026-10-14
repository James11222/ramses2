module sink_merger_module
#ifndef WITHOUTMPI
  use mpi
#endif
  use constants
  use marshal, only: pack_fetch_flag, unpack_fetch_flag
  use sink_evolution_module, only: sink_B_spline_weights_TSC
  use nbors_utils, only: get_parent_cell
#ifdef _CUDA
  use sink_device, only: gpu_sink_merger
#endif

  integer, parameter :: ncollision_max = 1000

  real(kind=8), dimension(1:3, 1:ncollision_max) :: pos1, pos2
  real(kind=8), dimension(1:3, 1:ncollision_max) :: vel1, vel2
  real(kind=8), dimension(1:ncollision_max) :: mass1, mass2
  integer, dimension(1:ncollision_max) :: id1, id2
  integer, dimension(1:ncollision_max) :: all_id1, all_id2
  logical, dimension(1:ncollision_max) :: exists1, exists2
  integer, dimension(1:ncollision_max) :: owner1, owner2
  integer, dimension(1:ncollision_max) :: local_id1, local_id2
  integer, dimension(1:ncollision_max) :: flag_duplicate

contains

  !==============================================================================
  ! Recursive function
  !==============================================================================
  recursive subroutine r_sink_merger(pst,ilevel,input_size)
    use mdl_module
    use ramses_commons, only: pst_t
    use mdl_parameters
    implicit none
    type(pst_t)::pst
    integer,VALUE::input_size
    integer::ilevel

    integer::rID

    if(pst%nLower>0)then
       rID = mdl_send_request(pst%s%mdl,MDL_SINK_MERGER,pst%iUpper+1,input_size,0,ilevel)
       call r_sink_merger(pst%pLower,ilevel,input_size)
       call mdl_get_reply(pst%s%mdl,rID,0)
    else
#ifdef _CUDA
       call gpu_sink_merger(pst%s, ilevel)
#else
       call sink_merger(pst,ilevel)
#endif
    endif

  end subroutine r_sink_merger

  !==============================================================================
  ! Sink merging main routine
  !==============================================================================
  subroutine sink_merger(pst, ilevel)
    use ramses_commons, only: pst_t, ramses_t
    use pm_commons, only: part_t
    use mdl_parameters
    use amr_parameters, only: ndim
    implicit none
    type(pst_t)::pst
    integer::ilevel

    integer::n_count_local, n_count_total, n_valid_mergers

    if(pst%s%g%myid==1.and.pst%s%r%verbose) write(*,*) 'Entering sink merger...'

    ! Step 1: Deposit sink IDs to grid
    call sink_id_deposition(pst%s, pst%s%sink, ilevel)

    ! Step 2: Collect local collisions
    call collect_collisions(pst%s, pst%s%sink, ilevel, n_count_local)
    if(pst%s%g%myid==1.and.pst%s%r%verbose) write(*,*) 'Proc:', pst%s%g%myid
    if(pst%s%g%myid==1.and.pst%s%r%verbose) write(*,*) 'Local collision count:', n_count_local

    ! Step 3: Gather all collisions across MPI ranks
    call gather_all_collisions(pst%s, n_count_local, n_count_total)
    if(n_count_total == 0) then
       if(pst%s%g%myid==1.and.pst%s%r%verbose) write(*,*) 'No collisions detected'
       return
    endif
    if(pst%s%g%myid==1.and.pst%s%r%verbose) write(*,*) 'Total collision count:', n_count_total

    ! Step 4: Deduplicate collision list
    call deduplicate_collisions(n_count_total)

    ! Step 5: Gather sink data with MPI
    call gather_sink_data(pst%s, pst%s%sink, n_count_total)

    ! Step 6: Apply filtering using merging criteria
    call filter_collision(pst%s, ilevel, n_count_total, n_valid_mergers)
    if(pst%s%g%myid==1.and.pst%s%r%verbose) write(*,*) 'Valid mergers after filtering:', n_valid_mergers

    ! Step 7: Break merger chains (defer conflicting pairs to next call)
    call break_merger_chains(n_count_total, n_valid_mergers)
    if(pst%s%r%verbose .and. pst%s%g%myid == 1) write(*,*) 'Valid mergers after chain breaking:', n_valid_mergers

    ! Step 8: Execute mergers
    if(n_valid_mergers > 0) then
       call execute_mergers(pst%s, pst%s%sink, n_count_total)
    endif

    if(pst%s%g%myid==1.and.pst%s%r%verbose) write(*,*) 'Sink merger process complete, executed:', n_valid_mergers

  end subroutine sink_merger

  !==============================================================================
  ! Sink ID deposition
  !==============================================================================
  subroutine sink_id_deposition(s,p,ilevel)
    use amr_parameters, only: ndim, twotondim, threetondim
    use ramses_commons, only: ramses_t
    use pm_commons, only: part_t
    use cache_commons, only: msg_int4
    use cache
    implicit none
    type(ramses_t)::s
    type(part_t)::p
    integer::ilevel

    real(kind=8),dimension(1:ndim)::xcen
    real(kind=8)::dx_loc
    integer(kind=8),dimension(0:ndim)::hash_nbor
    integer::ipart,icelln,igridn,j

    real(kind=8),dimension(1:ndim,1:threetondim)::xBHnei
    integer,dimension(1:ndim,1:threetondim)::ckeynei
    real(kind=8),dimension(1:threetondim)::vol
    real(kind=8)::weight
    integer::ind,igrid,new_id
    type(msg_int4)::dummy_int4

    associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

    dx_loc = s%r%boxlen/2**ilevel

    ! Set flag1 to max possible index
    do igrid=m%head(ilevel),m%tail(ilevel)
       do ind=1,twotondim
          m%flag1(ind,igrid)=huge(1)
       end do
    end do

    call open_cache(mdl, m, pack_size=storage_size(dummy_int4)/32, &
         init=init_flush_idsinkmin, flush=pack_flush_idsinkmin, combine=unpack_flush_idsinkmin)

    hash_nbor(0) = ilevel+1
    do ipart = p%headp(ilevel), p%tailp(ilevel)

       ! Skip zero-mass sinks
       if(p%mp(ipart) <= 0.0d0) cycle

       ! Compute TSC weights
       do j=1,ndim
          xcen(j) = (p%xp(ipart,j)+m%skip(j))/dx_loc
       end do
       call sink_B_spline_weights_TSC(s,xcen(1:ndim),xBHnei,ckeynei,vol,ilevel)

       do j = 1, threetondim
          if(vol(j) <= 0.0d0) cycle

          hash_nbor(1:ndim) = ckeynei(1:ndim,j)
          call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.true.,fetch_cache=.false.)

          if(igridn==0)cycle

          weight = vol(j)
          ! flag1 is 32-bit: make sure the sink id fits
          if(p%idp(ipart) > int(huge(1),kind=8))then
             write(*,*)'sink_id_deposition: sink id exceeds 32-bit range'
             call mdl_abort(mdl)
          endif
          new_id = int(p%idp(ipart), kind=4)
          if(new_id < m%flag1(icelln,igridn)) then
             m%flag1(icelln,igridn) = new_id
          endif
       end do
    end do

    call close_cache(mdl)

    end associate

  end subroutine sink_id_deposition

  subroutine init_flush_idsinkmin(mesh,igrid,hash_key)
    use amr_parameters, only: ndim, twotondim
    use amr_commons, only: mesh_t
    use cache_commons, only: msg_int4
    implicit none
    type(mesh_t)::mesh
    integer::igrid
    integer(kind=8),dimension(0:ndim)::hash_key
    integer::ind

    mesh%grid(igrid)%lev=hash_key(0)
    mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)

    do ind=1,twotondim
       mesh%flag1(ind,igrid)=huge(1)
    end do
  end subroutine init_flush_idsinkmin

  subroutine pack_flush_idsinkmin(mesh,igrid,msg_size,msg_array)
    use amr_parameters, only: twotondim
    use amr_commons, only: mesh_t
    use cache_commons, only: msg_int4
    implicit none
    type(mesh_t)::mesh
    integer::igrid
    integer::msg_size
    integer,dimension(1:msg_size),optional::msg_array
    integer::ind
    type(msg_int4)::msg

    do ind=1,twotondim
       msg%int4(ind)=mesh%flag1(ind,igrid)
    end do

    msg_array=transfer(msg,msg_array)
  end subroutine pack_flush_idsinkmin

  subroutine unpack_flush_idsinkmin(mesh,igrid,msg_size,msg_array,hash_key)
    use amr_parameters, only: ndim,twotondim
    use amr_commons, only: mesh_t
    use cache_commons, only: msg_int4
    implicit none
    type(mesh_t)::mesh
    integer::igrid
    integer::msg_size
    integer,dimension(1:msg_size),optional::msg_array
    integer(kind=8),dimension(0:ndim)::hash_key
    integer::ind
    type(msg_int4)::msg

    msg=transfer(msg_array,msg)

    do ind=1,twotondim
       if(msg%int4(ind) < mesh%flag1(ind,igrid)) then
          mesh%flag1(ind,igrid)=msg%int4(ind)
       endif
    end do
  end subroutine unpack_flush_idsinkmin

  !==============================================================================
  ! Collect collision pairs
  !==============================================================================
  subroutine collect_collisions(s,p,ilevel,n_local)
    use amr_parameters, only: ndim, twotondim, threetondim
    use ramses_commons, only: ramses_t
    use pm_commons, only: part_t
    use cache_commons, only: msg_int4
    use cache
    implicit none
    type(ramses_t)::s
    type(part_t)::p
    integer::ilevel,n_local

    real(kind=8),dimension(1:ndim)::xcen
    real(kind=8)::dx_loc
    integer(kind=8),dimension(0:ndim)::hash_nbor
    integer::ipart,icelln,igridn,j,my_id,neighbor_id,neighbor_id_min

    real(kind=8),dimension(1:ndim,1:threetondim)::xBHnei
    integer,dimension(1:ndim,1:threetondim)::ckeynei
    real(kind=8),dimension(1:threetondim)::vol
    type(msg_int4)::dummy_int4

    associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

    if(r%verbose .and. g%myid == 1)write(*,*)'Collecting collision pairs (using TSC)...'

    n_local = 0

    dx_loc = s%r%boxlen/2**ilevel

    hash_nbor(0) = ilevel+1
    call open_cache(mdl, m, pack_size=storage_size(dummy_int4)/32, &
         pack=pack_fetch_flag, unpack=unpack_fetch_flag)

    do ipart = p%headp(ilevel), p%tailp(ilevel)
       my_id = p%idp(ipart)

       ! Skip zero-mass sinks
       if(p%mp(ipart) <= 0.0d0) cycle

       do j=1,ndim
          xcen(j) = (p%xp(ipart,j)+m%skip(j))/dx_loc
       end do

       call sink_B_spline_weights_TSC(s,xcen(1:ndim),xBHnei,ckeynei,vol,ilevel)

       neighbor_id_min = huge(1)
       do j = 1, threetondim
          if(vol(j) <= 0.0d0) cycle

          hash_nbor(1:ndim) = ckeynei(1:ndim,j)
          call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.false.,fetch_cache=.true.)

          if(igridn>0) then
             neighbor_id = m%flag1(icelln,igridn)

             if(neighbor_id < my_id) then
                if(neighbor_id < neighbor_id_min) then
                   neighbor_id_min = neighbor_id
                endif
             endif
          endif
       end do

       if(neighbor_id_min < my_id) then
          n_local = n_local + 1
          if(n_local <= ncollision_max) then
             if(r%verbose)write(*,*)'Sink collision pair:', neighbor_id_min, my_id
             id1(n_local) = neighbor_id_min
             id2(n_local) = my_id
          endif
       endif
    end do

    call close_cache(mdl)

    if(r%verbose .and. g%myid == 1)write(*,*)'Collected',n_local,'collision pairs'

    end associate

  end subroutine collect_collisions

  !==============================================================================
  ! Gather all collisions across MPI ranks
  !==============================================================================
  subroutine gather_all_collisions(s, n_local, n_total)
    use ramses_commons, only: ramses_t
    use mdl_module, only: mdl_abort
    type(ramses_t) :: s
    integer::n_local
    integer::n_total
#ifndef WITHOUTMPI
    integer::i, ierr
    integer,dimension(1:s%g%ncpu)::recvcounts, displs
#endif
    n_total = n_local
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(n_local, n_total, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    if(n_total == 0) return
    if(n_total > ncollision_max) then
       write(*,*)'Too many collisions'
       write(*,*)'Increase ncollision_max'
       call mdl_abort(s%mdl)
    endif
#ifndef WITHOUTMPI
    call MPI_ALLGATHER(n_local, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    displs(1) = 0
    do i = 2, s%g%ncpu
       displs(i) = displs(i-1) + recvcounts(i-1)
    end do
#endif
    all_id1(1:n_local) = id1(1:n_local)
    all_id2(1:n_local) = id2(1:n_local)
#ifndef WITHOUTMPI
    call MPI_ALLGATHERV(id1, n_local, MPI_INTEGER, all_id1, recvcounts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHERV(id2, n_local, MPI_INTEGER, all_id2, recvcounts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#endif
  end subroutine gather_all_collisions

  !==============================================================================
  ! Deduplicate collisions to get unique collision pairs
  !==============================================================================
  subroutine deduplicate_collisions(n_pairs)
    integer::n_pairs
    integer::i, idx

    ! 1. Sort the pairs by all_id1 (so duplicate all_id1 values are adjacent)
    call sort_collision_pairs(n_pairs)

    ! 2. Mark subsequent duplicate all_id1 entries as 0
    flag_duplicate(1) = 1
    do i = 2, n_pairs
       if (all_id1(i) == all_id1(i-1)) then
          flag_duplicate(i) = 0
       else
          flag_duplicate(i) = 1
       endif
    end do

    ! 3. Compress all_id1 and all_id2 to remove the marked duplicate entries
    idx = 0
    do i = 1, n_pairs
       if (flag_duplicate(i) == 1) then
          idx = idx + 1
          all_id1(idx) = all_id1(i)
          all_id2(idx) = all_id2(i)
       endif
    end do
    n_pairs = idx
  end subroutine deduplicate_collisions

  !==============================================================================
  ! Gather sink physical data
  !==============================================================================
  subroutine gather_sink_data(s, p, n_pairs)
    use ramses_commons, only: ramses_t
    use pm_commons, only: part_t
    implicit none
    type(ramses_t)::s
    type(part_t)::p
    integer::n_pairs

    integer::i, ipart, ierr

    ! 1. Initialize flat module arrays
    mass1(1:n_pairs) = 0.0d0
    mass2(1:n_pairs) = 0.0d0
    pos1(1:3, 1:n_pairs) = 0.0d0
    pos2(1:3, 1:n_pairs) = 0.0d0
    vel1(1:3, 1:n_pairs) = 0.0d0
    vel2(1:3, 1:n_pairs) = 0.0d0
    owner1(1:n_pairs) = 0
    owner2(1:n_pairs) = 0
    local_id1(1:n_pairs) = 0
    local_id2(1:n_pairs) = 0
    exists1(1:n_pairs) = .false.
    exists2(1:n_pairs) = .false.
#if NDIM==3
    ! 2. Fill properties
    do ipart = 1, p%npart
       do i = 1, n_pairs
          if(all_id1(i) == p%idp(ipart)) then
             mass1(i) = p%mp(ipart)
             pos1(1:3, i) = p%xp(ipart, 1:3)
             vel1(1:3, i) = p%vp(ipart, 1:3)
             owner1(i) = s%g%myid
             local_id1(i) = ipart
          endif
          if(all_id2(i) == p%idp(ipart)) then
             mass2(i) = p%mp(ipart)
             pos2(1:3, i) = p%xp(ipart, 1:3)
             vel2(1:3, i) = p%vp(ipart, 1:3)
             owner2(i) = s%g%myid
             local_id2(i) = ipart
          endif
       end do
    end do
#endif
#ifndef WITHOUTMPI
    ! 3. Combine properties globally
    call MPI_ALLREDUCE(MPI_IN_PLACE, mass1, n_pairs, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, mass2, n_pairs, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    call MPI_ALLREDUCE(MPI_IN_PLACE, pos1, 3*n_pairs, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, pos2, 3*n_pairs, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    call MPI_ALLREDUCE(MPI_IN_PLACE, vel1, 3*n_pairs, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, vel2, 3*n_pairs, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    call MPI_ALLREDUCE(MPI_IN_PLACE, owner1, n_pairs, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, owner2, n_pairs, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif

    ! 4. Set exists flag
    do i = 1, n_pairs
       if(mass1(i) > 0.0d0) exists1(i) = .true.
       if(mass2(i) > 0.0d0) exists2(i) = .true.
    end do

  end subroutine gather_sink_data

  !==============================================================================
  ! Filter collisions based on adopted merging criteria
  !==============================================================================
  subroutine filter_collision(s, ilevel, n_total, n_filtered)
    use ramses_commons, only: ramses_t
    implicit none
    type(ramses_t)::s
    integer::ilevel, n_total, n_filtered

    integer::i, idim
    real(kind=8)::dx_loc, factG
    real(kind=8)::separation, dvx, dvy, dvz, vel_squared
    real(kind=8)::total_mass, binding_criteria
    real(kind=8),dimension(1:3)::dpos

    dx_loc = s%r%boxlen/2**ilevel
    factG = 1.0d0
    if(s%r%cosmo) factG = 3.0d0/4.0d0/twopi * s%g%omega_m * s%g%aexp

    ! Apply merger criteria directly using gathered data
    do i = 1, n_total
       ! Check if both sinks exist and have mass
       if (.not. exists1(i) .or. .not. exists2(i) .or. &
            & mass1(i) <= 0.0d0 .or. mass2(i) <= 0.0d0) then
          all_id1(i) = 0
          all_id2(i) = 0
          cycle
       endif

       ! Physical check inlined
       ! Use the minimum-image convention across periodic boundaries
       do idim = 1, 3
          dpos(idim) = pos1(idim, i) - pos2(idim, i)
          if(s%r%periodic(idim))then
             if(dpos(idim) >  0.5d0*s%r%box_size(idim)) dpos(idim) = dpos(idim) - s%r%box_size(idim)
             if(dpos(idim) < -0.5d0*s%r%box_size(idim)) dpos(idim) = dpos(idim) + s%r%box_size(idim)
          endif
       end do
       separation = norm2(dpos)

       if (separation < dx_loc) then
          dvx = vel1(1, i) - vel2(1, i)
          dvy = vel1(2, i) - vel2(2, i)
          dvz = vel1(3, i) - vel2(3, i)
          vel_squared = dvx*dvx + dvy*dvy + dvz*dvz

          total_mass = mass1(i) + mass2(i)
          binding_criteria = (factG * total_mass) / dx_loc * (1.0d0 - (separation/dx_loc)**2)

          if (vel_squared >= binding_criteria) then
             all_id1(i) = 0
             all_id2(i) = 0
          endif
       else
          all_id1(i) = 0
          all_id2(i) = 0
       endif
    end do

    ! Count filtered pairs
    n_filtered = 0
    do i = 1, n_total
       if (all_id1(i) /= 0) n_filtered = n_filtered + 1
    end do

  end subroutine filter_collision

  !==============================================================================
  ! Break merger chains. Conflicting pairs are invalidated.
  !==============================================================================
  subroutine break_merger_chains(n_total, n_valid)
    implicit none
    integer::n_total, n_valid

    integer::i, k

    do i = 1, n_total
       if(all_id1(i) == 0) cycle
       do k = 1, i-1
          if(all_id1(k) /= 0 .and. all_id2(k) == all_id1(i)) then
             all_id1(i) = 0
             all_id2(i) = 0
             exit
          endif
       end do
    end do

    n_valid = 0
    do i = 1, n_total
       if (all_id1(i) /= 0) n_valid = n_valid + 1
    end do

  end subroutine break_merger_chains

  !==============================================================================
  ! Execute valid mergers (set masses to zero for merged sinks)
  !==============================================================================
  subroutine execute_mergers(s, p, n_total)
    use ramses_commons, only: ramses_t
    use pm_commons, only: part_t
    implicit none
    type(ramses_t)::s
    type(part_t)::p
    integer::n_total

    integer::i, j, myrank, idim
    integer::id_keep, id_delete
    real(kind=8)::m1, m2, total_mass
    real(kind=8),dimension(1:3)::com_position, com_velocity, rel

    myrank = s%g%myid

    do i = 1, n_total
       id_keep = all_id1(i)
       id_delete = all_id2(i)

       ! Skip invalid pairs
       if(id_keep == 0) cycle

       ! Calculate merged properties
       m1 = mass1(i)
       m2 = mass2(i)
       total_mass = m1 + m2

       ! Compute the center of mass using the minimum-image offset of
       ! sink 2 relative to sink 1, then wrap back into the box
       do idim = 1, 3
          rel(idim) = pos2(idim, i) - pos1(idim, i)
          if(s%r%periodic(idim))then
             if(rel(idim) >  0.5d0*s%r%box_size(idim)) rel(idim) = rel(idim) - s%r%box_size(idim)
             if(rel(idim) < -0.5d0*s%r%box_size(idim)) rel(idim) = rel(idim) + s%r%box_size(idim)
          endif
          com_position(idim) = pos1(idim, i) + m2 * rel(idim) / total_mass
          if(s%r%periodic(idim))then
             if(com_position(idim) <  0.0d0             ) com_position(idim) = com_position(idim) + s%r%box_size(idim)
             if(com_position(idim) >= s%r%box_size(idim)) com_position(idim) = com_position(idim) - s%r%box_size(idim)
          endif
       end do
       com_velocity(1:3) = (m1 * vel1(1:3, i) + m2 * vel2(1:3, i)) / total_mass

#if NDIM==3
       ! Each CPU only updates the sinks it owns
       if(owner1(i) == myrank) then
          j = local_id1(i)
          if(j > 0) then
             p%mp(j) = total_mass
             p%xp(j,1:3) = com_position(1:3)
             p%vp(j,1:3) = com_velocity(1:3)
          endif
       endif
       if(owner2(i) == myrank) then
          j = local_id2(i)
          if(j > 0) then
             p%mp(j) = 0.0d0
             p%xp(j,1:3) = com_position(1:3)
             p%vp(j,1:3) = com_velocity(1:3)
             p%tm(j)     = s%g%texp
             p%idm(j)    = id_keep
             write(*, '(" Sink ",I0," merged into sink ",I0)') id_delete, id_keep
          endif
       endif
#endif

    end do
  end subroutine execute_mergers

  !==============================================================================
  ! Sort collision pairs by id1, then id2 (bubble sort)
  !==============================================================================
  subroutine sort_collision_pairs(n)
    implicit none
    integer::n

    integer::i, j
    integer::temp_id1, temp_id2
    logical::swapped

    do i = n, 1, -1
       swapped = .false.
       do j = 1, i-1
          if(all_id1(j) > all_id1(j+1) .or. &
               & (all_id1(j) == all_id1(j+1) .and. all_id2(j) > all_id2(j+1))) then
             temp_id1 = all_id1(j)
             temp_id2 = all_id2(j)
             all_id1(j) = all_id1(j+1)
             all_id2(j) = all_id2(j+1)
             all_id1(j+1) = temp_id1
             all_id2(j+1) = temp_id2
             swapped = .true.
          endif
       end do
       if(.not.swapped) exit
    end do

  end subroutine sort_collision_pairs

end module sink_merger_module
