#if NDIM != 3
#error "CORE=1 requires NDIM=3"
#endif
module core_merger_module
#ifndef WITHOUTMPI
  use mpi
#endif
  use amr_parameters, only: ndim, twotondim, i8b
  use amr_commons
  use nbors_utils, only: get_parent_cell

  ! Core ids a cell keeps for the candidate search
  integer,parameter::ncore_marker_slot=4

  type :: core_collision_pair_t
     integer :: id1
     integer :: id2
  end type core_collision_pair_t

  type :: core_data_t
     integer :: marker
     integer(i8b) :: id
     real(kind=8) :: mass
     real(kind=8), dimension(1:3) :: position
     real(kind=8), dimension(1:3) :: velocity
     real(kind=8), dimension(1:3) :: angular_momentum
     real(kind=8) :: birth_time
     real(kind=8) :: metallicity
     real(kind=8) :: thermal_energy
     logical :: exists
     logical :: changed
     integer :: cpu_owner
  end type core_data_t

  ! Cache message with the marker ids of the cells of one oct
  type :: msg_core_marker
     integer(kind=4),dimension(1:ncore_marker_slot,1:twotondim)::mark
  end type msg_core_marker

  ! Sparse per-oct store of up to ncore_marker_slot marker ids per cell.
  type core_marker_oct_t
     integer::igrid=0
     integer(kind=8)::key(0:ndim)=0_8
     integer::mark(ncore_marker_slot,twotondim)=huge(1)
  end type core_marker_oct_t
  type(core_marker_oct_t),allocatable::marker_oct(:)
  integer,allocatable::marker_hash(:)
  integer::marker_n=0,marker_cap=0

contains

  !==============================================================================
  ! Sparse marker store (open-addressing hash on igrid)
  !==============================================================================
  subroutine core_marker_store_clear()
    implicit none
    if(allocated(marker_oct))deallocate(marker_oct)
    if(allocated(marker_hash))deallocate(marker_hash)
    marker_n=0
    marker_cap=0
  end subroutine core_marker_store_clear

  integer function core_marker_bucket(igrid,nbucket) result(bucket)
    implicit none
    integer,intent(in)::igrid,nbucket
    bucket=int(modulo(int(igrid,kind=8)*104729_8,int(nbucket,kind=8)))+1
  end function core_marker_bucket

  integer function core_marker_find(igrid) result(slot)
    implicit none
    integer,intent(in)::igrid
    integer::bucket,entry
    slot=0
    if(.not.allocated(marker_hash))return
    bucket=core_marker_bucket(igrid,size(marker_hash))
    do
       entry=marker_hash(bucket)
       if(entry==0)return
       if(marker_oct(entry)%igrid==igrid)then
          slot=entry
          return
       endif
       bucket=modulo(bucket,size(marker_hash))+1
    enddo
  end function core_marker_find

  subroutine core_marker_grow()
    implicit none
    type(core_marker_oct_t),allocatable::octs(:)
    integer,allocatable::buckets(:)
    integer::newcap,i,bucket
    newcap=max(64,2*marker_cap)
    allocate(octs(newcap))
    allocate(buckets(2*newcap))
    buckets=0
    if(marker_n>0)octs(1:marker_n)=marker_oct(1:marker_n)
    do i=1,marker_n
       bucket=core_marker_bucket(octs(i)%igrid,size(buckets))
       do while(buckets(bucket)/=0)
          bucket=modulo(bucket,size(buckets))+1
       enddo
       buckets(bucket)=i
    enddo
    call move_alloc(octs,marker_oct)
    call move_alloc(buckets,marker_hash)
    marker_cap=newcap
  end subroutine core_marker_grow

  integer function core_marker_entry(igrid,key,create) result(slot)
    implicit none
    integer,intent(in)::igrid
    integer(kind=8),intent(in)::key(0:ndim)
    logical,intent(in)::create
    integer::bucket
    slot=core_marker_find(igrid)
    if(slot>0)then
       if(any(marker_oct(slot)%key/=key))then
          marker_oct(slot)%key=key
          marker_oct(slot)%mark=huge(1)
       endif
       return
    endif
    if(.not.create)return
    if(marker_n==marker_cap)call core_marker_grow()
    marker_n=marker_n+1
    slot=marker_n
    marker_oct(slot)%igrid=igrid
    marker_oct(slot)%key=key
    marker_oct(slot)%mark=huge(1)
    bucket=core_marker_bucket(igrid,size(marker_hash))
    do while(marker_hash(bucket)/=0)
       bucket=modulo(bucket,size(marker_hash))+1
    enddo
    marker_hash(bucket)=slot
  end function core_marker_entry

  subroutine core_marker_reset(igrid,key)
    implicit none
    integer,intent(in)::igrid
    integer(kind=8),intent(in)::key(0:ndim)
    integer::slot
    slot=core_marker_find(igrid)
    if(slot==0)return
    marker_oct(slot)%key=key
    marker_oct(slot)%mark=huge(1)
  end subroutine core_marker_reset

  subroutine core_marker_read(igrid,key,mark)
    implicit none
    integer,intent(in)::igrid
    integer(kind=8),intent(in)::key(0:ndim)
    integer,intent(out)::mark(ncore_marker_slot,twotondim)
    integer::slot
    mark=huge(1)
    slot=core_marker_entry(igrid,key,.false.)
    if(slot>0)mark=marker_oct(slot)%mark
  end subroutine core_marker_read

  subroutine core_marker_mesh_key(mesh,igrid,key)
    implicit none
    type(mesh_t),intent(in)::mesh
    integer,intent(in)::igrid
    integer(kind=8),intent(out)::key(0:ndim)
    key(0)=mesh%grid(igrid)%lev
    key(1:ndim)=mesh%grid(igrid)%ckey(1:ndim)
  end subroutine core_marker_mesh_key

  subroutine core_marker_cell_insert(mesh,igrid,ind,id)
    implicit none
    type(mesh_t),intent(in)::mesh
    integer,intent(in)::igrid,ind,id
    integer(kind=8)::key(0:ndim)
    integer::slot
    call core_marker_mesh_key(mesh,igrid,key)
    slot=core_marker_entry(igrid,key,.true.)
    call core_marker_insert(marker_oct(slot)%mark(:,ind),id)
  end subroutine core_marker_cell_insert

  subroutine core_marker_mesh_read(mesh,igrid,mark)
    implicit none
    type(mesh_t),intent(in)::mesh
    integer,intent(in)::igrid
    integer,intent(out)::mark(ncore_marker_slot,twotondim)
    integer(kind=8)::key(0:ndim)
    call core_marker_mesh_key(mesh,igrid,key)
    call core_marker_read(igrid,key,mark)
  end subroutine core_marker_mesh_read

  !==============================================================================
  ! Recursive function
  !==============================================================================
  recursive subroutine r_core_merger(pst,ilevel,input_size)
    use mdl_module
    use ramses_commons, only: pst_t
    use mdl_parameters
    implicit none
    type(pst_t)::pst
    integer,VALUE::input_size
    integer::ilevel
    integer::rID

    if(pst%nLower>0)then
       rID = mdl_send_request(pst%s%mdl,MDL_CORE_MERGER,pst%iUpper+1,input_size,0,ilevel)
       call r_core_merger(pst%pLower,ilevel,input_size)
       call mdl_get_reply(pst%s%mdl,rID,0)
    else
       call core_merger(pst,ilevel)
    endif

  end subroutine r_core_merger

  !==============================================================================
  ! Core merging main routine
  !==============================================================================
  subroutine core_merger(pst, ilevel)
    use constants, only: yr2sec
    use ramses_commons, only: pst_t
    implicit none
    type(pst_t)::pst
    integer::ilevel

    type(core_collision_pair_t),dimension(:),allocatable::local_collision_pairs
    integer,dimension(:),allocatable::all_id1, all_id2
    integer,dimension(:),allocatable::unique_ids
    type(core_data_t),dimension(:),allocatable::global_core_data
    integer::nnei, search_order
    integer::n_count_local, n_count_total, n_unique
    integer::n_accepted
    integer(kind=8)::pair_capacity
    real(kind=8)::dx_loc, merge_radius, t_merge
    real(kind=8)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

    if(pst%s%r%core_merging_timescale <= 0.0d0)return

    search_order = pst%s%r%core_merger_search_order

    dx_loc = pst%s%r%boxlen/2**ilevel
    merge_radius = dble(search_order)*dx_loc
    call units(pst%s%r,pst%s%g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    ! tp and texp are proper times in 1/H0: convert yr with scale_t/aexp**2.
    t_merge = pst%s%r%core_merging_timescale*yr2sec &
         /(scale_t/pst%s%g%aexp**2)
    nnei = (2*search_order+1)**ndim

    call core_id_deposition(pst%s, pst%s%core, ilevel, dx_loc, nnei, search_order)
    ! Each core at this level pairs once with each smaller id in its stencil:
    ! at most nnei*ncore_marker_slot ids, and fewer than npart_tot.
    pair_capacity = max(1_8, int(pst%s%core%tailp(ilevel)-pst%s%core%headp(ilevel)+1,kind=8)* &
         min(int(nnei*ncore_marker_slot,kind=8),pst%s%core%npart_tot-1))
    allocate(local_collision_pairs(pair_capacity))
    call collect_core_collisions(pst%s, pst%s%core, ilevel, dx_loc, nnei, search_order, &
         local_collision_pairs, n_count_local)
    call core_marker_store_clear()
    call gather_all_core_collisions(pst%s%g, local_collision_pairs, n_count_local, all_id1, all_id2, n_count_total)
    if(n_count_total == 0) then
       deallocate(local_collision_pairs)
       return
    endif

    call sort_core_collision_pairs(all_id1, all_id2, n_count_total)

    call build_unique_core_list(all_id1, all_id2, n_count_total, unique_ids, n_unique)

    allocate(global_core_data(n_unique))
    call gather_all_core_data(pst%s%g, pst%s%core, ilevel, unique_ids, n_unique, global_core_data)

    call merge_core_pairs(all_id1,all_id2,n_count_total,merge_radius,t_merge, &
         pst%s%g%texp,pst%s%r%box_size,pst%s%r%periodic,global_core_data,n_unique, &
         n_accepted)

    if(n_accepted>0)then
       call update_merged_cores(pst%s%g,pst%s%core,global_core_data,n_unique)
    endif

    deallocate(local_collision_pairs, all_id1, all_id2, unique_ids, global_core_data)

  end subroutine core_merger

  !==============================================================================
  ! Keys of the (2*nsearch+1)**3 cells around x; periodic keys wrap and cells
  ! outside the box get vol=0
  !==============================================================================
  subroutine core_enlarged_stencil(s,x,nsearch,nnei,ckey,vol,ilevel)
    use ramses_commons, only: ramses_t
    implicit none
    type(ramses_t)::s
    real(kind=8),dimension(1:ndim)::x
    integer,intent(in)::nsearch,nnei,ilevel
    integer,dimension(1:ndim,1:nnei)::ckey
    real(kind=8),dimension(1:nnei)::vol
    integer::idim,idx,ii,jj,kk
    integer,dimension(1:ndim)::center,key
    logical::valid

    associate(r=>s%r,m=>s%m)

    idx=0
    ckey=0
    vol=0d0

    do idim=1,ndim
       center(idim)=int(x(idim))
    end do

    do kk=-nsearch,nsearch
       do jj=-nsearch,nsearch
          do ii=-nsearch,nsearch
             idx=idx+1
             key(1)=center(1)+ii
             key(2)=center(2)+jj
             key(3)=center(3)+kk
             valid=.true.
             do idim=1,ndim
                if(r%periodic(idim))then
                   do while(key(idim)<m%box_ckey_min(idim,ilevel+1))
                      key(idim)=key(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
                   end do
                   do while(key(idim)>=m%box_ckey_max(idim,ilevel+1))
                      key(idim)=key(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
                   end do
                else
                   if(key(idim)<m%box_ckey_min(idim,ilevel+1))valid=.false.
                   if(key(idim)>=m%box_ckey_max(idim,ilevel+1))valid=.false.
                endif
             end do
             ckey(1:ndim,idx)=key(1:ndim)
             if(valid)vol(idx)=1d0
          end do
       end do
    end do

    end associate
  end subroutine core_enlarged_stencil

  !==============================================================================
  ! First marker id of this rank: ids are dense over ranks per level
  !==============================================================================
  subroutine core_marker_offset(g,p,ilevel,offset)
    use pm_commons, only: part_t
    implicit none
    type(global_t)::g
    type(part_t)::p
    integer,intent(in)::ilevel
    integer,intent(out)::offset
    integer::ipart,n_local_active
#ifndef WITHOUTMPI
    integer::ierr
#endif

    n_local_active = 0
    do ipart = p%headp(ilevel), p%tailp(ilevel)
       if(p%mp(ipart) > 0.0d0)n_local_active = n_local_active + 1
    end do

#ifndef WITHOUTMPI
    call MPI_EXSCAN(n_local_active, offset, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    if(g%myid == 1)offset = 0
#else
    offset = 0
#endif
  end subroutine core_marker_offset

  !==============================================================================
  ! Core marker deposition
  !==============================================================================
  subroutine core_id_deposition(s,p,ilevel,dx_loc,nnei,search_order)
    use ramses_commons, only: ramses_t
    use pm_commons, only: part_t
    use cache
    implicit none
    type(ramses_t)::s
    type(part_t)::p
    integer::ilevel,nnei,search_order
    real(kind=8)::dx_loc

    real(kind=8),dimension(1:ndim)::xcen
    integer(kind=8),dimension(0:ndim)::hash_nbor
    integer::ipart,icelln,igridn,j
    integer,dimension(1:ndim,1:nnei)::ckeynei
    real(kind=8),dimension(1:nnei)::vol
    integer::marker_id
    type(msg_core_marker)::dummy_marker

    associate(r=>s%r,m=>s%m,mdl=>s%mdl)

    call open_cache(mdl, m, pack_size=storage_size(dummy_marker)/32, &
         init=init_flush_core_marker, flush=pack_flush_core_marker, combine=unpack_flush_core_marker)

    hash_nbor(0) = ilevel+1
    call core_marker_offset(s%g,p,ilevel,marker_id)
    do ipart = p%headp(ilevel), p%tailp(ilevel)
       if(p%mp(ipart) <= 0.0d0) cycle
       marker_id = marker_id + 1

       do j=1,ndim
          xcen(j) = (p%xp(ipart,j)+m%skip(j))/dx_loc
       end do

       call core_enlarged_stencil(s,xcen(1:ndim),search_order,nnei,ckeynei,vol,ilevel)

       do j = 1, nnei
          if(vol(j) <= 0.0d0) cycle
          hash_nbor(1:ndim) = ckeynei(1:ndim,j)
          call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.true.,fetch_cache=.false.)
          if(igridn==0)cycle
          call core_marker_cell_insert(m,igridn,icelln,marker_id)
       end do
    end do

    call close_cache(mdl)

    end associate
  end subroutine core_id_deposition

  ! Keep the ncore_marker_slot smallest ids per cell, ascending, padded with huge(1).
  subroutine core_marker_insert(slot,id)
    implicit none
    integer,dimension(1:ncore_marker_slot)::slot
    integer::id
    integer::k,pos

    if(id>=slot(ncore_marker_slot))return
    do k=1,ncore_marker_slot
       if(slot(k)==id)return
    end do
    pos=ncore_marker_slot
    do k=1,ncore_marker_slot
       if(id<slot(k))then
          pos=k
          exit
       endif
    end do
    do k=ncore_marker_slot,pos+1,-1
       slot(k)=slot(k-1)
    end do
    slot(pos)=id
  end subroutine core_marker_insert

  ! Combine a ghost copy's ids with the owner's.
  subroutine core_marker_merge(dst,src)
    implicit none
    integer,dimension(1:ncore_marker_slot)::dst,src
    integer::k

    do k=1,ncore_marker_slot
       if(src(k)==huge(1))exit
       call core_marker_insert(dst,src(k))
    end do
  end subroutine core_marker_merge

  !==============================================================================
  ! Cache routines for the marker store
  !==============================================================================
  subroutine init_flush_core_marker(mesh,igrid,hash_key)
    use amr_commons, only: mesh_t
    implicit none
    type(mesh_t)::mesh
    integer::igrid
    integer(kind=8),dimension(0:ndim)::hash_key
    mesh%grid(igrid)%lev=hash_key(0)
    mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
    call core_marker_reset(igrid,hash_key)
  end subroutine init_flush_core_marker

  subroutine pack_flush_core_marker(mesh,igrid,msg_size,msg_array)
    use amr_commons, only: mesh_t
    implicit none
    type(mesh_t)::mesh
    integer::igrid,msg_size
    integer,dimension(1:msg_size),optional::msg_array
    type(msg_core_marker)::msg

    call core_marker_mesh_read(mesh,igrid,msg%mark)
    msg_array=transfer(msg,msg_array)
  end subroutine pack_flush_core_marker

  subroutine unpack_flush_core_marker(mesh,igrid,msg_size,msg_array,hash_key)
    use amr_commons, only: mesh_t
    implicit none
    type(mesh_t)::mesh
    integer::igrid,msg_size
    integer,dimension(1:msg_size),optional::msg_array
    integer(kind=8),dimension(0:ndim)::hash_key
    integer::ind,slot
    integer(kind=8)::key(0:ndim)
    type(msg_core_marker)::msg

    msg=transfer(msg_array,msg)
    if(all(msg%mark==huge(1)))return
    call core_marker_mesh_key(mesh,igrid,key)
    slot=core_marker_entry(igrid,key,.true.)
    do ind=1,twotondim
       call core_marker_merge(marker_oct(slot)%mark(1:ncore_marker_slot,ind), &
            msg%mark(1:ncore_marker_slot,ind))
    end do
  end subroutine unpack_flush_core_marker

  subroutine unpack_fetch_core_marker(mesh,igrid,msg_size,msg_array,hash_key)
    use amr_commons, only: mesh_t
    implicit none
    type(mesh_t)::mesh
    integer::igrid,msg_size
    integer,dimension(1:msg_size),optional::msg_array
    integer(kind=8),dimension(0:ndim)::hash_key
    integer::slot
    type(msg_core_marker)::msg

    mesh%grid(igrid)%lev=hash_key(0)
    mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
    msg=transfer(msg_array,msg)
    if(all(msg%mark==huge(1)))then
       call core_marker_reset(igrid,hash_key)
       return
    endif
    slot=core_marker_entry(igrid,hash_key,.true.)
    marker_oct(slot)%mark=msg%mark
  end subroutine unpack_fetch_core_marker

  !==============================================================================
  ! Gather all collision pairs across MPI ranks
  !==============================================================================
  subroutine gather_all_core_collisions(g, local_pairs, n_local, all_id1, all_id2, n_total)
    implicit none
    type(global_t)::g
    type(core_collision_pair_t),dimension(:)::local_pairs
    integer::n_local
    integer,dimension(:),allocatable::all_id1, all_id2
    integer::n_total

    integer::i
#ifndef WITHOUTMPI
    integer::ierr
    integer,dimension(1:g%ncpu)::recvcounts, displs, counts2
    integer,dimension(:),allocatable::temp_pairs, all_pairs

    call MPI_ALLGATHER(n_local, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    n_total = sum(recvcounts)
    if(n_total == 0)return

    ! Gather the pairs interleaved in one call.
    counts2 = 2*recvcounts
    displs(1) = 0
    do i = 2, g%ncpu
       displs(i) = displs(i-1) + counts2(i-1)
    end do

    allocate(temp_pairs(max(1,2*n_local)))
    do i = 1, n_local
       temp_pairs(2*i-1) = local_pairs(i)%id1
       temp_pairs(2*i)   = local_pairs(i)%id2
    end do

    allocate(all_pairs(2*n_total), all_id1(n_total), all_id2(n_total))
    call MPI_ALLGATHERV(temp_pairs, 2*n_local, MPI_INTEGER, &
         all_pairs, counts2, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    do i = 1, n_total
       all_id1(i) = all_pairs(2*i-1)
       all_id2(i) = all_pairs(2*i)
    end do

    deallocate(temp_pairs, all_pairs)
#else
    n_total = n_local
    if(n_total == 0)return
    allocate(all_id1(n_total), all_id2(n_total))
    do i = 1, n_total
       all_id1(i) = local_pairs(i)%id1
       all_id2(i) = local_pairs(i)%id2
    end do
#endif
  end subroutine gather_all_core_collisions

  !==============================================================================
  ! List of the markers that appear in a pair
  !==============================================================================
  subroutine build_unique_core_list(all_id1, all_id2, n_pairs, unique_ids, n_unique)
    implicit none
    integer,dimension(:)::all_id1, all_id2
    integer::n_pairs
    integer,dimension(:),allocatable::unique_ids
    integer::n_unique
    integer::i, j
    logical::found
    integer,dimension(:),allocatable::temp_ids

    allocate(temp_ids(2*n_pairs))
    n_unique = 0
    do i = 1, n_pairs
       found = .false.
       do j = 1, n_unique
          if(temp_ids(j) == all_id1(i)) then
             found = .true.
             exit
          endif
       end do
       if(.not. found) then
          n_unique = n_unique + 1
          temp_ids(n_unique) = all_id1(i)
       endif

       found = .false.
       do j = 1, n_unique
          if(temp_ids(j) == all_id2(i)) then
             found = .true.
             exit
          endif
       end do
       if(.not. found) then
          n_unique = n_unique + 1
          temp_ids(n_unique) = all_id2(i)
       endif
    end do

    allocate(unique_ids(n_unique))
    unique_ids(1:n_unique) = temp_ids(1:n_unique)
    deallocate(temp_ids)
  end subroutine build_unique_core_list

  !==============================================================================
  ! Gather the data of every listed core on every rank
  !==============================================================================
  subroutine gather_all_core_data(g, p, ilevel, unique_ids, n_unique, cores)
    use pm_commons, only: part_t
    implicit none
    type(global_t)::g
    type(part_t)::p
    integer,intent(in)::ilevel
    integer,dimension(:)::unique_ids
    integer::n_unique
    type(core_data_t),dimension(:)::cores

    integer::i, ipart, base
    integer::marker_id
#ifndef WITHOUTMPI
    integer::ierr, mpi_id_type
#endif
    ! Per core: m, x(3), v(3), j(3), tp, zp, up; id, owner myid.
    integer,parameter::NR=13, NI=2
    real(kind=8),dimension(:),allocatable::rbuf
    integer(i8b),dimension(:),allocatable::ibuf

    allocate(rbuf(NR*max(1,n_unique)), ibuf(NI*max(1,n_unique)))
    rbuf = 0.0d0
    ibuf = 0_i8b

    ! Marker ids are per level: walk the level the deposition walked.
    call core_marker_offset(g,p,ilevel,marker_id)
    do ipart = p%headp(ilevel), p%tailp(ilevel)
       if(p%mp(ipart) <= 0.0d0) cycle
       marker_id = marker_id + 1
       do i = 1, n_unique
          if(unique_ids(i) == marker_id) then
             base = NR*(i-1)
             rbuf(base+1) = p%mp(ipart)
             rbuf(base+2:base+4) = p%xp(ipart, 1:3)
             rbuf(base+5:base+7) = p%vp(ipart, 1:3)
             rbuf(base+8:base+10) = p%jp(ipart, 1:3)
             rbuf(base+11) = p%tp(ipart)
             rbuf(base+12) = p%zp(ipart)
             rbuf(base+13) = p%up(ipart)
             ibuf(NI*(i-1)+1) = p%idp(ipart)
             ibuf(NI*(i-1)+2) = int(g%myid,kind=i8b)
             exit
          endif
       end do
    end do

    ! Only the owner fills a core's slots, so a sum over ranks gathers it.
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(MPI_IN_PLACE, rbuf, NR*n_unique, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_COMM_WORLD, ierr)
    mpi_id_type=MPI_INTEGER4
    if(i8b==8)mpi_id_type=MPI_INTEGER8
    call MPI_ALLREDUCE(MPI_IN_PLACE, ibuf, NI*n_unique, mpi_id_type, &
         MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    do i = 1, n_unique
       base = NR*(i-1)
       cores(i)%marker = unique_ids(i)
       cores(i)%id = ibuf(NI*(i-1)+1)
       cores(i)%cpu_owner = int(ibuf(NI*(i-1)+2))
       cores(i)%exists = .true.
       cores(i)%changed = .false.
       cores(i)%mass = rbuf(base+1)
       cores(i)%position(1:3) = rbuf(base+2:base+4)
       cores(i)%velocity(1:3) = rbuf(base+5:base+7)
       cores(i)%angular_momentum(1:3) = rbuf(base+8:base+10)
       cores(i)%birth_time = rbuf(base+11)
       cores(i)%metallicity = rbuf(base+12)
       cores(i)%thermal_energy = rbuf(base+13)
    end do

    deallocate(rbuf, ibuf)
  end subroutine gather_all_core_data

  !==============================================================================
  ! Merge the accepted pairs in order on the replicated data
  !==============================================================================
  subroutine merge_core_pairs(all_id1,all_id2,n_pairs,merge_radius,t_merge,t_now, &
       box_size,periodic,cores,n_cores,n_accepted)
    use pm_commons, only: cross
    implicit none
    integer,dimension(:)::all_id1,all_id2
    integer::n_pairs,n_cores,n_accepted
    real(kind=8)::merge_radius,t_merge,t_now
    real(kind=8),dimension(1:3)::box_size
    logical,dimension(1:ndim)::periodic
    type(core_data_t),dimension(:)::cores
    integer::i,j,icore1,icore2
    real(kind=8)::mass1,mass2,total_mass,znew,unew,tpnew
    real(kind=8),dimension(1:3)::com_position,com_unwrapped,momentum,vnew,jnew,ltot
    real(kind=8),dimension(1:3)::pos_keep,pos_delete,dpos
    logical::should_merge

    n_accepted=0

    do i=1,n_pairs
       icore1=0
       icore2=0
       do j=1,n_cores
          if(cores(j)%marker==all_id1(i))icore1=j
          if(cores(j)%marker==all_id2(i))icore2=j
       end do

       if(.not.cores(icore1)%exists .or. .not.cores(icore2)%exists)cycle

       should_merge=check_core_merger_criteria(cores(icore1),cores(icore2), &
            merge_radius,t_merge,t_now,box_size,periodic)
       if(.not.should_merge)cycle

       mass1=cores(icore1)%mass
       mass2=cores(icore2)%mass
       total_mass=mass1+mass2

       pos_keep(1:3)=cores(icore1)%position(1:3)
       dpos(1:3)=cores(icore2)%position(1:3)-pos_keep(1:3)
       call core_minimum_image_delta(dpos,box_size,periodic)
       pos_delete(1:3)=pos_keep(1:3)+dpos(1:3)
       com_unwrapped(1:3)=(mass1*pos_keep(1:3)+mass2*pos_delete(1:3))/total_mass
       com_position(1:3)=com_unwrapped(1:3)
       call core_wrap_position(com_position,box_size,periodic)
       momentum(1:3)=mass1*cores(icore1)%velocity(1:3)+ &
            mass2*cores(icore2)%velocity(1:3)
       vnew(1:3)=momentum(1:3)/total_mass
       unew=(mass1*cores(icore1)%thermal_energy+ &
            mass2*cores(icore2)%thermal_energy)/total_mass
       ! Mass-weighted metallicity.
       znew=(mass1*cores(icore1)%metallicity+ &
            mass2*cores(icore2)%metallicity)/total_mass
       tpnew=min(cores(icore1)%birth_time,cores(icore2)%birth_time)
       ltot(1:3)=mass1*cores(icore1)%angular_momentum(1:3)+ &
            mass2*cores(icore2)%angular_momentum(1:3)+ &
            mass1*cross(pos_keep(1:3)-com_unwrapped(1:3), &
                        cores(icore1)%velocity(1:3)-vnew(1:3))+ &
            mass2*cross(pos_delete(1:3)-com_unwrapped(1:3), &
                        cores(icore2)%velocity(1:3)-vnew(1:3))
       jnew(1:3)=ltot(1:3)/total_mass

       cores(icore1)%mass=total_mass
       cores(icore1)%position(1:3)=com_position(1:3)
       cores(icore1)%velocity(1:3)=vnew(1:3)
       cores(icore1)%angular_momentum(1:3)=jnew(1:3)
       cores(icore1)%birth_time=tpnew
       cores(icore1)%metallicity=znew
       cores(icore1)%thermal_energy=unew
       cores(icore1)%changed=.true.

       cores(icore2)%mass=0.0d0
       cores(icore2)%position(1:3)=com_position(1:3)
       cores(icore2)%velocity(1:3)=vnew(1:3)
       cores(icore2)%angular_momentum(1:3)=jnew(1:3)
       cores(icore2)%birth_time=tpnew
       cores(icore2)%metallicity=znew
       cores(icore2)%thermal_energy=unew
       cores(icore2)%exists=.false.
       cores(icore2)%changed=.true.
       n_accepted=n_accepted+1
    end do
  end subroutine merge_core_pairs

  !==============================================================================
  ! Each owner writes back its changed cores
  !==============================================================================
  subroutine update_merged_cores(g,p,cores,n_cores)
    use pm_commons, only: part_t
    implicit none
    type(global_t)::g
    type(part_t)::p
    type(core_data_t),dimension(:)::cores
    integer::n_cores
    integer::i,j

    do i=1,n_cores
       if(.not.cores(i)%changed)cycle
       if(cores(i)%cpu_owner/=g%myid)cycle
       do j=1,p%npart
          if(p%idp(j)==cores(i)%id)then
             p%mp(j)=cores(i)%mass
             p%xp(j,1:3)=cores(i)%position(1:3)
             p%vp(j,1:3)=cores(i)%velocity(1:3)
             p%jp(j,1:3)=cores(i)%angular_momentum(1:3)
             p%zp(j)=cores(i)%metallicity
             p%up(j)=cores(i)%thermal_energy
             p%tp(j)=cores(i)%birth_time
             ! The merged core counts as having accreted at this step
             if(cores(i)%mass>0d0)p%nacc(j)=g%nstep_coarse
             exit
          endif
       end do
    end do
  end subroutine update_merged_cores

  !==============================================================================
  ! Collect local collision pairs
  !==============================================================================
  subroutine collect_core_collisions(s,p,ilevel,dx_loc,nnei,search_order,all_pairs,n_local)
    use ramses_commons, only: ramses_t
    use pm_commons, only: part_t
    use cache
    implicit none
    type(ramses_t)::s
    type(part_t)::p
    integer::ilevel,nnei,search_order,n_local
    real(kind=8)::dx_loc
    type(core_collision_pair_t),dimension(:)::all_pairs

    real(kind=8),dimension(1:ndim)::xcen
    integer(kind=8),dimension(0:ndim)::hash_nbor
    integer::ipart,icelln,igridn,j,k,my_id,neighbor_id
    integer::islot,nnear
    integer::cell_mark(ncore_marker_slot,twotondim)
    integer,dimension(1:nnei*ncore_marker_slot)::near_id
    integer,dimension(1:ndim,1:nnei)::ckeynei
    real(kind=8),dimension(1:nnei)::vol
    type(msg_core_marker)::dummy_marker

    associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

    n_local = 0
    hash_nbor(0) = ilevel+1
    call open_cache(mdl, m, pack_size=storage_size(dummy_marker)/32, &
         pack=pack_flush_core_marker, unpack=unpack_fetch_core_marker)

    call core_marker_offset(g,p,ilevel,my_id)
    do ipart = p%headp(ilevel), p%tailp(ilevel)
       if(p%mp(ipart) <= 0.0d0) cycle
       my_id = my_id + 1

       do j=1,ndim
          xcen(j) = (p%xp(ipart,j)+m%skip(j))/dx_loc
       end do

       call core_enlarged_stencil(s,xcen(1:ndim),search_order,nnei,ckeynei,vol,ilevel)

       ! Ids below my_id in the stencil cells, ascending and each once
       nnear = 0
       do j = 1, nnei
          if(vol(j) <= 0.0d0) cycle
          hash_nbor(1:ndim) = ckeynei(1:ndim,j)
          call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.false.,fetch_cache=.true.)

          if(igridn>0) then
             call core_marker_mesh_read(m,igridn,cell_mark)
             ! Slots are ascending: stop at the first id not below my_id.
             do islot = 1, ncore_marker_slot
                neighbor_id = cell_mark(islot,icelln)
                if(neighbor_id >= my_id) exit
                call core_id_add(near_id,nnear,neighbor_id)
             end do
          endif
       end do

       do k = 1, nnear
          n_local = n_local + 1
          all_pairs(n_local)%id1 = near_id(k)
          all_pairs(n_local)%id2 = my_id
       end do
    end do

    call close_cache(mdl)
    end associate
  end subroutine collect_core_collisions

  ! Add id to the ascending list ids(1:n) unless it is there already.
  subroutine core_id_add(ids,n,id)
    implicit none
    integer,dimension(:)::ids
    integer::n,id
    integer::k

    do k=1,n
       if(ids(k)==id)return
       if(ids(k)>id)exit
    end do
    ids(k+1:n+1)=ids(k:n)
    ids(k)=id
    n=n+1
  end subroutine core_id_add

  !==============================================================================
  ! Periodic helpers
  !==============================================================================
  ! Minimum-image offset along the periodic directions
  subroutine core_minimum_image_delta(dpos,box_size,periodic)
    implicit none
    real(kind=8),dimension(1:3)::dpos
    real(kind=8),dimension(1:3)::box_size
    logical,dimension(1:ndim)::periodic
    integer::idim

    do idim=1,ndim
       if(periodic(idim))dpos(idim)=dpos(idim)-dnint(dpos(idim)/box_size(idim))*box_size(idim)
    end do
  end subroutine core_minimum_image_delta

  subroutine core_wrap_position(xpos,box_size,periodic)
    implicit none
    real(kind=8),dimension(1:3),intent(inout)::xpos
    real(kind=8),dimension(1:3),intent(in)::box_size
    logical,dimension(1:ndim),intent(in)::periodic
    integer::idim

    do idim=1,ndim
       if(periodic(idim))xpos(idim)=modulo(xpos(idim),box_size(idim))
    end do
  end subroutine core_wrap_position

  !==============================================================================
  ! Merging criteria: closer than merge_radius, both younger than t_merge
  !==============================================================================
  logical function check_core_merger_criteria(core1,core2,merge_radius,t_merge,t_now,box_size,periodic)
    implicit none
    type(core_data_t)::core1,core2
    real(kind=8)::merge_radius,t_merge,t_now
    real(kind=8),dimension(1:3)::box_size
    logical,dimension(1:ndim)::periodic
    real(kind=8)::separation,age1,age2
    real(kind=8),dimension(1:3)::dpos

    check_core_merger_criteria = .false.

    dpos(1:3) = core1%position(1:3) - core2%position(1:3)
    call core_minimum_image_delta(dpos,box_size,periodic)
    separation = sqrt(sum(dpos(1:ndim)**2))

    age1 = t_now - core1%birth_time
    age2 = t_now - core2%birth_time

    if(age1 < t_merge .and. age2 < t_merge .and. separation < merge_radius) &
         check_core_merger_criteria = .true.
  end function check_core_merger_criteria

  !==============================================================================
  ! Sort collision pairs by id1, then id2; no pair appears twice
  !==============================================================================
  subroutine sort_core_collision_pairs(all_id1, all_id2, n)
    implicit none
    integer::n
    integer,dimension(1:n)::all_id1, all_id2
    integer::i, j

    call quick_sort_int_int(all_id1(1), all_id2(1), n)
    ! Sort id2 within each run of equal id1, carrying the equal id1 along
    i = 1
    do while(i <= n)
       j = i
       do while(j < n)
          if(all_id1(j+1) /= all_id1(i)) exit
          j = j + 1
       end do
       call quick_sort_int_int(all_id2(i), all_id1(i), j-i+1)
       i = j + 1
    end do
  end subroutine sort_core_collision_pairs

end module core_merger_module
