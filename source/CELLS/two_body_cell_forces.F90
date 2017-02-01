Module two_body_cell_forces_module

Implicit none

Private
Public :: two_body_cell_forces

Contains


Subroutine two_body_cell_forces                        &
           (rcut,rlnk,rvdw,rmet,pdplnc,keyens,    &
           alpha,epsq,keyfce,nstfce,lbook,megfrz, &
           lrdf,nstrdf,leql,nsteql,nstep,         &
           elrc,virlrc,elrcm,vlrcm,               &
           engcpe,vircpe,engsrp,virsrp,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly 4 subroutine for calculating interatomic forces and rdf
! using the cell list method
!
!NYI ntpvdw > 0 ------- switch for vdw potential calculation
!NYI ntpmet > 0 ------- switch for metal local density and potentials
!                    calculation
!
!
! ELECTROSTATICS KEYS
!
!NYI keyfce = 0 ------ no electrostatics
!NYI keyfce = 2 ------ Ewald sum (ewald_spme,ewald_real,ewald_excl)
!NYI keyfce = 4 ------ distance dependent dielectric potential (coul_dddp)
!NYI keyfce = 6 ------ coulombic 1/r potential (coul_cp)
!NYI keyfce = 8 ------ force-shifted and damped coulombic potential (coul_fscp)
!NYI keyfce =10 ------ reaction field and damped coulombic potential (coul_rfp)
!NYI keyfce =12 ------ direct space Poisson solver (possion_module)
!
!NYI nstfce - the rate at which the k-space contributions of SPME are
!          refreshed.  Once every 1 <= nstfce <= 7 steps.
!
! copyright - daresbury laboratory
! author    - a.b.g.chalk november 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
#ifdef __OPENMP
  Use comms_module,  Only : idnode,mxnode,gsum,mxthreads
  Use config_module, Only : cell,volm,sumchg,natms,list,xxx,yyy,zzz,fxx,fyy,fzz
  Use rdf_module,    Only : ncfrdf,rdf, block_size, l_block, l_jack
#else
  Use comms_module,  Only : idnode,mxnode,gsum
  Use config_module, Only : cell,volm,sumchg,natms,list,xxx,yyy,zzz
  Use rdf_module,    Only : ncfrdf, block_size, l_block, l_jack
#endif
  Use setup_module
  Use vnl_module,    Only : l_vnl
  Use site_module,    Only : ntpatm,unqatm
  Use ewald_module
  Use mpoles_module,  Only : induce
  Use poisson_module, Only : poisson_forces,poisson_excl_forces,poisson_frzn_forces
  Use vdw_module,     Only : ntpvdw, ld_vdw, ls_vdw
  Use metal_module,   Only : ntpmet
  Use kim_module
  Use block_averages_module
  Use cell_list_module
  Use sort_cell_module
  Use cell_ewald_real_forces_module2
  Use cell_rdf_collect_module
  Use cell_vdw_forces_module2

!  Use mpi ! remove me soon
#ifdef CHRONO                                                                                           
#ifdef SERIAL
  Use mpi_module
#else
  Use mpi
#endif
#endif

#ifdef __OPENMP
  Use omp_lib
#endif

Implicit None

  Logical,                                  Intent( In    ) :: lrdf,leql
  Logical, Intent(InOut) :: lbook
  Integer,                                  Intent( In    ) :: keyens,        &
                                                               keyfce,nstfce, &
                                                               megfrz,nstrdf, &
                                                               nsteql,nstep
  Real( Kind = wp ),                        Intent( In    ) :: rcut,rlnk,rvdw,rmet, &
                                                               pdplnc,alpha,epsq
  Real( Kind = wp ),                        Intent( In    ) :: elrc,virlrc
  Real( Kind = wp ), Dimension( 0:mxatyp ), Intent( InOut ) :: elrcm,vlrcm
  Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe, &
                                                               engsrp,virsrp
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
  Logical,           Save :: new_nz    = .true.
  Real( Kind = wp ), Save :: factor_nz = 0.0_wp

  Logical           :: safe = .true., l_do_rdf
  Integer           :: fail,i,j,k,limit, id, p, id2
  Real( Kind = wp ) :: engcpe_rc,vircpe_rc, &
                       engcpe_fr,vircpe_fr, &
                       engcpe_nz,vircpe_nz,          vircpe_dt, &
                       engden,virden,engmet,virmet,             &
                       engkim,virkim,             &
                       engacc,viracc,tmp,buffer(0:17), engacc2, viracc2

#ifdef __OPENMP
  Real( Kind = wp ), Allocatable, Dimension(:), Save   :: engvdw, virvdw,&
  engcpe_rl, vircpe_rl, engcpe_ex, vircpe_ex
  Real( Kind = wp ), Allocatable, Dimension(:,:), Save :: l_stress
  Integer, Save :: nr_threads
  Logical, Save :: first = .false.
  Integer       :: tid
#else
  Real( Kind = wp ) :: engvdw, virvdw, engcpe_rl, vircpe_rl, engcpe_ex,&
  vircpe_ex
#endif


!!TODO Add more variables

#ifdef __OPENMP
if(.not. first) then
!$omp parallel default(none) shared(nr_threads, first)
!$omp single
nr_threads = omp_get_num_threads()
first = .true.
!$omp end single
!$omp end parallel
Allocate(engvdw(0:nr_threads-1), virvdw(0:nr_threads-1), l_stress(1:9,0:(nr_threads-1))) !TODO status
Allocate(engcpe_rl(0:nr_threads-1), vircpe_rl(0:nr_threads-1),engcpe_ex(0:nr_threads-1),vircpe_ex(0:nr_threads-1))
!l_stress(1:9,0) = stress(1:9)
Do i=1,9
  l_stress(i,0) = stress(i)
End Do
end if
#ifdef TASK_TIMERS
timestep = timestep + 1
step_start = omp_get_wtime()
#endif
#endif

!  lbook = .false.
!Do we need to do rdf in this step?
  l_do_rdf = (lrdf .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,nstrdf) == 0)


  If (keyfce == 2 .or. keyfce == 12) Call ewald_check(keyens,megfrz,nsteql,nstfce,nstep)


! initialise energy and virial accumulators

  engkim = 0.0_wp
  virkim = 0.0_wp


  engden    = 0.0_wp
  virden    = 0.0_wp

  engmet    = 0.0_wp
  virmet    = 0.0_wp

  engvdw    = 0.0_wp
  virvdw    = 0.0_wp

  engsrp    = 0.0_wp
  virsrp    = 0.0_wp

  engcpe_rc = 0.0_wp
  vircpe_rc = 0.0_wp

  engcpe_rl = 0.0_wp
  vircpe_rl = 0.0_wp

  engcpe_ex = 0.0_wp
  vircpe_ex = 0.0_wp

  engcpe_fr = 0.0_wp
  vircpe_fr = 0.0_wp

  engcpe_nz = 0.0_wp
  vircpe_nz = 0.0_wp

  vircpe_dt = 0.0_wp

  engcpe    = 0.0_wp
  vircpe    = 0.0_wp

  engkim = 0.0_wp
  virkim = 0.0_wp


!Set up non-bonded interaction cell lists (aka pseudo verlet lists) if needed, also computed neighbour cells.

!Believe this works.
  If ((.not. induce) .and. l_vnl) Call cell_lists(rcut, rlnk, rvdw, rmet, pdplnc, lbook, megfrz)

! Update halo data positions in the xxt/yyt/zzt arrays. Should work
  Call perform_halo_mapping(xxx, yyy, zzz, natms+1, nlast)

!Reset temp forces array to 0. Should work.
  Call reset_force_store()

!TODO If we redid cell lists then create the sorted indices arrays (this should be tasks, so maybe later). Currently done in serial.
 If ((.not. induce) .and. l_vnl) then
!$omp parallel do
  do i=1, number_cells
    Call sort_cell(cell_list(i))
  end do
!$omp end parallel do
End If

!TODO NYI - Calculate KIM interactions.

!TODO NYI - If ntpmet > 0 calculate local density in metals
!TODO NYI - If keyfce == 2 .and. l_fce then calculate coulombic forces, Ewald sum- fourier contribution
!TODO NYI - Might be we don't need to change this.
  If(keyfce == 2 .and. l_fce) Then
    If(mximpl > 0) Then
    If(mxompl <= 2) Then
        !print *, NYI
        Call ewald_spme_mforces_d(alpha, epsq, engcpe_rc, vircpe_rc, stress)
    Else
        !print *, NYI
        Call ewald_spme_mforces(alpha, epsq, engcpe_rc, vircpe_rc, stress)
    End If
    Else
      Call ewald_spme_forces(alpha, epsq, engcpe_rc, vircpe_rc, stress)
    End If
  End If
t_start = omp_get_wtime()


!TODO NYI - Calculate metal_forces
!TODO NYI - Calculate vdw_forces
!Make sure we update engvdw in some way after each task.
!$omp parallel default(none) shared(ld_vdw, ls_vdw,x_cells, y_cells, z_cells,mximpl,keyfce) &
!$omp shared(cell_list,engvdw,virvdw, rvdw, rcut, l_stress, neighbours, number_cells, l_do_rdf) &
!$omp shared(engcpe_rl, vircpe_rl, engcpe_ex, vircpe_ex, alpha, epsq) &
!$omp private (i,j,k,p,id,engacc,viracc, id2, tid)
if(ld_vdw) then
if(ls_vdw) then
!ld & ls
!$omp do collapse(3) 
do i=0, x_cells+1
  do j=0, y_cells+1
    do k=0, z_cells+1
      id = 1+i+(x_cells+2)*(j+(y_cells+2)*k)
      If(.not. cell_list(id)%is_halo) then

#ifdef __OPENMP
      !$omp task default(none) shared(rvdw, cell_list, l_stress,engvdw,virvdw) private(engacc, viracc) firstprivate(id, id2) private(tid) untied
       call cell_self_vdw_ld_ls_interaction(rvdw,cell_list(id), engacc, viracc,l_stress)
tid = omp_get_thread_num()
       engvdw(tid) = engvdw(tid) + engacc
       virvdw(tid) = virvdw(tid) + viracc
       !$omp end task
#else
       call cell_self_vdw_ld_ls_interaction(rvdw,cell_list(id), engacc, viracc,stress)
       engvdw = engvdw + engacc
       virvdw = virvdw + viracc
#endif

       end if
      do p = 1,13
        id2 = neighbours(p, id)
        if(id2 <= 0 .or. id2 > number_cells) CYCLE
        if(cell_list(id)%is_halo) then
          if(.not. cell_list(id2)%is_halo) then
            !If this is a halo cell, we interact it with all of its non-halo neighbours and only store updates on the domain cells.
      
#ifdef __OPENMP
!TODO envdw needs to be a array of size num_thrads
      !$omp task default(none) shared(rvdw, cell_list,l_stress,engvdw,virvdw) private(engacc,viracc) firstprivate(p) firstprivate(id, id2) private(tid) untied
            call cell_pair_vdw_ld_ls_interaction(rvdw,cell_list(id2), cell_list(id), p+13, engacc, viracc, l_stress)
            tid = omp_get_thread_num()
            engvdw(tid) = engvdw(tid) + engacc
            virvdw(tid) = virvdw(tid) + viracc
            !$omp end task
#else
            call cell_pair_vdw_ld_ls_interaction(rvdw,cell_list(id2), cell_list(id), p+13, engacc, viracc, stress)
            engvdw = engvdw + engacc
            virvdw = virvdw + viracc
#endif
      

          if(l_do_rdf) then
            !$omp task default(none) shared(rcut, cell_list) firstprivate(p) firstprivate(id, id2) private(tid) untied
            Call pair_rdf_collect(rcut, cell_list(id2), cell_list(id), p+13)
            !$omp end task
           end if
         end if
        else
          if(cell_list(id2)%is_halo) then
            !If cell_list id2 is a halo cell then don't store forces on that cell.

#ifdef __OPENMP
!TODO envdw needs to be a array of size num_thrads
            !$omp task default(none) shared(rvdw, cell_list,l_stress,engvdw,virvdw) firstprivate(p) private(engacc, viracc) firstprivate(id, id2) private(tid) untied
            call cell_pair_vdw_ld_ls_interaction(rvdw,cell_list(id), cell_list(id2),p, engacc, viracc, l_stress)
            tid = omp_get_thread_num()
            engvdw(tid) = engvdw(tid) + engacc
            virvdw(tid) = virvdw(tid) + viracc
            !$omp end task
#else
            call cell_pair_vdw_ld_ls_interaction(rvdw,cell_list(id), cell_list(id2),p, engacc, viracc, stress)
            engvdw = engvdw + engacc
            virvdw = virvdw + viracc
#endif

          !TODO Do we need to collect rdf here?
          else
            !cell_list id2 is not a halo cell, store forces on both cells (remove if statement).

#ifdef __OPENMP
          !TODO envdw needs to be a array of size num_thrads
            !$omp task default(none) shared(rvdw, cell_list,engvdw,virvdw,l_stress) firstprivate(p) private(engacc, viracc) firstprivate(id, id2) private(tid) untied 
            call cell_local_pair_vdw_ld_ls_interaction(rvdw,cell_list(id), cell_list(id2), p, engacc, viracc, l_stress)
            tid = omp_get_thread_num()
            engvdw(tid) = engvdw(tid) + engacc
            virvdw(tid) = virvdw(tid) + viracc
            !$omp end task
#else
            call cell_local_pair_vdw_ld_ls_interaction(rvdw,cell_list(id), cell_list(id2), p, engacc, viracc, stress)
            engvdw = engvdw + engacc
            virvdw = virvdw + viracc
#endif

            if(l_do_rdf) then
              !$omp task default(none) shared(rcut, cell_list) firstprivate(p) firstprivate(id, id2) private(tid) untied
              Call pair_rdf_collect(rcut, cell_list(id), cell_list(id2), p)
              !$omp end task
            end if
          end if
        end if
      end do !p
      If(l_do_rdf .and. (.not. cell_list(id)%is_halo)) then
        !$omp task default(none) shared(rcut, cell_list) firstprivate(id, id2) private(tid) untied
        Call self_rdf_collect(rcut, cell_list(id)) 
        !$omp end task
      end if
    end do !k
  end do !j
end do !i
!$omp end do nowait
else !not ls_vdw
!ld & not ls
!$omp do collapse(3)
do i=0, x_cells+1
  do j=0, y_cells+1
    do k=0, z_cells+1
      id = 1+i+(x_cells+2)*(j+(y_cells+2)*k)
      If(.not. cell_list(id)%is_halo) then

#ifdef __OPENMP
          !TODO envdw needs to be a array of size num_thrads
      !$omp task default(none) shared(rvdw, cell_list, l_stress,engvdw,virvdw) private(engacc, viracc) firstprivate(id, id2) private(tid) untied
       call cell_self_vdw_ld_non_ls_interaction(rvdw,cell_list(id), engacc, viracc,l_stress)
            tid = omp_get_thread_num()
            engvdw(tid) = engvdw(tid) + engacc
            virvdw(tid) = virvdw(tid) + viracc
      !$omp end task
#else
       call cell_self_vdw_ld_non_ls_interaction(rvdw,cell_list(id), engacc, viracc,stress)
       engvdw = engvdw + engacc
       virvdw = virvdw + viracc
#endif

       end if
      do p = 1,13
        id2 = neighbours(p, id)
        if(id2 <= 0 .or. id2 > number_cells) CYCLE
        if(cell_list(id)%is_halo) then
          if(.not. cell_list(id2)%is_halo) then
            !If this is a halo cell, we interact it with all of its non-halo neighbours and only store updates on the domain cells.

#ifdef __OPENMP
      !$omp task default(none) shared(rvdw, cell_list, l_stress,engvdw,virvdw) private(engacc, viracc) firstprivate(p) firstprivate(id, id2) private(tid) untied
            call cell_pair_vdw_ld_non_ls_interaction(rvdw,cell_list(id2), cell_list(id), p+13, engacc, viracc, l_stress)
            tid = omp_get_thread_num()
            engvdw(tid) = engvdw(tid) + engacc
            virvdw(tid) = virvdw(tid) + viracc
      !$omp end task
#else
            call cell_pair_vdw_ld_non_ls_interaction(rvdw,cell_list(id2), cell_list(id), p+13, engacc, viracc, stress)
            engvdw = engvdw + engacc
            virvdw = virvdw + viracc
#endif

          if(l_do_rdf) then
           !$omp task default(none) shared(rcut, cell_list) firstprivate(p) firstprivate(id, id2) private(tid) untied
            Call pair_rdf_collect(rcut, cell_list(id2), cell_list(id), p+13)
           !$omp end task
          end if
         end if
        else
          if(cell_list(id2)%is_halo) then
            !If cell_list id2 is a halo cell then don't store forces on that cell.

#ifdef __OPENMP
          !TODO envdw needs to be a array of size num_thrads
      !$omp task default(none) shared(rvdw, cell_list, l_stress,engvdw,virvdw) private(engacc, viracc) firstprivate(p) firstprivate(id, id2) private(tid) untied
            call cell_pair_vdw_ld_non_ls_interaction(rvdw,cell_list(id), cell_list(id2),p, engacc, viracc, l_stress)
            tid = omp_get_thread_num()
            engvdw(tid) = engvdw(tid) + engacc
            virvdw(tid) = virvdw(tid) + viracc
          !$omp end task      
#else
            call cell_pair_vdw_ld_non_ls_interaction(rvdw,cell_list(id), cell_list(id2),p, engacc, viracc, stress)
            engvdw = engvdw + engacc
            virvdw = virvdw + viracc
#endif      

          !TODO Do we need to collect rdf here?
          else
            !cell_list id2 is not a halo cell, store forces on both cells (remove if statement).

#ifdef __OPENMP
          !TODO envdw needs to be a array of size num_thrads
      !$omp task default(none) shared(rvdw, cell_list, l_stress,engvdw,virvdw) private(engacc, viracc) firstprivate(p) firstprivate(id, id2) private(tid) untied
            call cell_local_pair_vdw_ld_non_ls_interaction(rvdw,cell_list(id), cell_list(id2), p, engacc, viracc, l_stress)
            tid = omp_get_thread_num()
            engvdw(tid) = engvdw(tid) + engacc
            virvdw(tid) = virvdw(tid) + viracc
            !$omp end task
#else
            call cell_local_pair_vdw_ld_non_ls_interaction(rvdw,cell_list(id), cell_list(id2), p, engacc, viracc, stress)
            engvdw = engvdw + engacc
            virvdw = virvdw + viracc
#endif

            if(l_do_rdf) then
              !$omp task default(none) shared(rcut, cell_list) firstprivate(p) firstprivate(id, id2) private(tid) untied
              Call pair_rdf_collect(rcut, cell_list(id), cell_list(id2), p)
              !$omp end task
            end if
          end if
        end if
      end do !p
      If(l_do_rdf .and. (.not. cell_list(id)%is_halo)) then 
        !$omp task default(none) shared(rcut, cell_list) firstprivate(id, id2) private(tid) untied
        Call self_rdf_collect(rcut, cell_list(id)) 
        !$omp end task
      end if
    end do !k
  end do !j
end do !i
!$omp end do nowait
end if

else !no ld_vdw
if(ls_vdw) then
!no ld, ls
!$omp do collapse(3) 
do i=0, x_cells+1
  do j=0, y_cells+1
    do k=0, z_cells+1
      id = 1+i+(x_cells+2)*(j+(y_cells+2)*k)
      If(.not. cell_list(id)%is_halo) then
#ifdef __OPENMP
      !$omp task default(none) shared(rvdw, cell_list, l_stress,engvdw,virvdw) private(engacc, viracc) firstprivate(id, id2) private(tid) untied
       call cell_self_vdw_non_ld_ls_interaction(rvdw,cell_list(id), engacc, viracc,l_stress)
          !TODO envdw needs to be a array of size num_thrads
            tid = omp_get_thread_num()
            engvdw(tid) = engvdw(tid) + engacc
            virvdw(tid) = virvdw(tid) + viracc
#else
       call cell_self_vdw_non_ld_ls_interaction(rvdw,cell_list(id), engacc, viracc, stress)
       engvdw = engvdw + engacc
       virvdw = virvdw + viracc
#endif      
      !$omp end task
       end if
      do p = 1,13
        id2 = neighbours(p, id)
        if(id2 <= 0 .or. id2 > number_cells) CYCLE
        if(cell_list(id)%is_halo) then
          if(.not. cell_list(id2)%is_halo) then
            !If this is a halo cell, we interact it with all of its non-halo neighbours and only store updates on the domain cells.

#ifdef __OPENMP
          !TODO envdw needs to be a array of size num_thrads
      !$omp task default(none) shared(rvdw, cell_list, l_stress,engvdw,virvdw) private(engacc, viracc) firstprivate(p) firstprivate(id, id2) private(tid) untied
            call cell_pair_vdw_non_ld_ls_interaction(rvdw,cell_list(id2), cell_list(id), p+13, engacc, viracc, l_stress)
            tid = omp_get_thread_num()
            engvdw(tid) = engvdw(tid) + engacc
            virvdw(tid) = virvdw(tid) + viracc
      !$omp end task
#else
            call cell_pair_vdw_non_ld_ls_interaction(rvdw,cell_list(id2), cell_list(id), p+13, engacc, viracc, stress)
            engvdw = engvdw + engacc
            virvdw = virvdw + viracc
#endif

          if(l_do_rdf) then
            !$omp task default(none) shared(rcut, cell_list) firstprivate(p) firstprivate(id, id2) private(tid) untied
            Call pair_rdf_collect(rcut, cell_list(id2), cell_list(id), p+13)
            !$omp end task
          end if
        end if
        else
          if(cell_list(id2)%is_halo) then
            !If cell_list id2 is a halo cell then don't store forces on that cell.

#ifdef __OPENMP
          !TODO envdw needs to be a array of size num_thrads
      !$omp task default(none) shared(rvdw, cell_list, l_stress,engvdw,virvdw) private(engacc, viracc) firstprivate(p) firstprivate(id, id2) private(tid) untied
            call cell_pair_vdw_non_ld_ls_interaction(rvdw,cell_list(id), cell_list(id2),p, engacc, viracc, l_stress)
            tid = omp_get_thread_num()
            engvdw(tid) = engvdw(tid) + engacc
            virvdw(tid) = virvdw(tid) + viracc
          !$omp end task
#else
            call cell_pair_vdw_non_ld_ls_interaction(rvdw,cell_list(id), cell_list(id2),p, engacc, viracc, stress)
            engvdw = engvdw + engacc
            virvdw = virvdw + viracc
#endif
          !TODO Do we need to collect rdf here?

          else
            !cell_list id2 is not a halo cell, store forces on both cells (remove if statement).

#ifdef __OPENMP
          !TODO envdw needs to be a array of size num_thrads
      !$omp task default(none) shared(rvdw, cell_list, l_stress,engvdw,virvdw) private(engacc, viracc) firstprivate(p) firstprivate(id, id2) private(tid) untied
            call cell_local_pair_vdw_non_ld_ls_interaction(rvdw,cell_list(id), cell_list(id2), p, engacc, viracc, l_stress)
            tid = omp_get_thread_num()
            engvdw(tid) = engvdw(tid) + engacc
            virvdw(tid) = virvdw(tid) + viracc
        !$omp end task
#else
            call cell_local_pair_vdw_non_ld_ls_interaction(rvdw,cell_list(id), cell_list(id2), p, engacc, viracc, stress)
            engvdw = engvdw + engacc
            virvdw = virvdw + viracc
#endif

            if(l_do_rdf) then
              !$omp task default(none) shared(rcut, cell_list) firstprivate(p) firstprivate(id, id2) private(tid) untied
              Call pair_rdf_collect(rcut, cell_list(id), cell_list(id2), p)
              !$omp end task
            end if
          end if
        end if
      end do !p
      If(l_do_rdf .and. (.not. cell_list(id)%is_halo)) then
        !$omp task default(none) shared(rcut, cell_list) firstprivate(id, id2) private(tid) untied
        Call self_rdf_collect(rcut, cell_list(id)) 
        !$omp end task
      end if
    end do !k
  end do !j
end do !i
!omp end do nowait

else !no ls_vdw
!no ls, no ld
!$omp do collapse(3) 
do i=0, x_cells+1
  do j=0, y_cells+1
    do k=0, z_cells+1
      id = 1+i+(x_cells+2)*(j+(y_cells+2)*k)
      If(.not. cell_list(id)%is_halo) then

#ifdef __OPENMP
          !TODO envdw needs to be a array of size num_thrads
      !$omp task default(none) shared(rvdw, cell_list, l_stress,engvdw,virvdw) private(engacc, viracc)  firstprivate(id, id2) private(tid) untied
       call cell_self_vdw_non_ld_non_ls_interaction(rvdw,cell_list(id), engacc, viracc,l_stress)
            tid = omp_get_thread_num()
            engvdw(tid) = engvdw(tid) + engacc
            virvdw(tid) = virvdw(tid) + viracc
      !$omp end task
#else
       call cell_self_vdw_non_ld_non_ls_interaction(rvdw,cell_list(id), engacc, viracc,stress)
       engvdw = engvdw + engacc
       virvdw = virvdw + viracc
#endif

       end if
      do p = 1,13
        id2 = neighbours(p, id)
        if(id2 <= 0 .or. id2 > number_cells) CYCLE
        if(cell_list(id)%is_halo) then
          if(.not. cell_list(id2)%is_halo) then
            !If this is a halo cell, we interact it with all of its non-halo neighbours and only store updates on the domain cells.

#ifdef __OPENMP
          !TODO envdw needs to be a array of size num_thrads
      !$omp task default(none) shared(rvdw, cell_list, l_stress,engvdw,virvdw) private(engacc, viracc) firstprivate(p) firstprivate(id, id2) private(tid) untied
            call cell_pair_vdw_non_ld_non_ls_interaction(rvdw,cell_list(id2), cell_list(id), p+13, engacc, viracc, l_stress)
            tid = omp_get_thread_num()
            engvdw(tid) = engvdw(tid) + engacc
            virvdw(tid) = virvdw(tid) + viracc
      !$omp end task
#else
            call cell_pair_vdw_non_ld_non_ls_interaction(rvdw,cell_list(id2), cell_list(id), p+13, engacc, viracc, stress)
            engvdw = engvdw + engacc
            virvdw = virvdw + viracc
#endif

          if(l_do_rdf) then
            !$omp task default(none) shared(rcut, cell_list) firstprivate(p) firstprivate(id, id2) private(tid) untied
            Call pair_rdf_collect(rcut, cell_list(id2), cell_list(id), p+13)
            !$omp end task
          end if
          end if
        else
          if(cell_list(id2)%is_halo) then
            !If cell_list id2 is a halo cell then don't store forces on that cell.

#ifdef __OPENMP
          !TODO envdw needs to be a array of size num_thrads
      !$omp task default(none) shared(rvdw, cell_list, l_stress,engvdw,virvdw) private(engacc, viracc) firstprivate(p) firstprivate(id, id2) private(tid) untied
            call cell_pair_vdw_non_ld_non_ls_interaction(rvdw,cell_list(id), cell_list(id2),p, engacc, viracc, l_stress)
            tid = omp_get_thread_num()
            engvdw(tid) = engvdw(tid) + engacc
            virvdw(tid) = virvdw(tid) + viracc
      !$omp end task
#else
            call cell_pair_vdw_non_ld_non_ls_interaction(rvdw,cell_list(id), cell_list(id2),p, engacc, viracc, stress)
            engvdw = engvdw + engacc
            virvdw = virvdw + viracc
#endif

          !TODO Do we need to collect rdf here?
          else
            !cell_list id2 is not a halo cell, store forces on both cells (remove if statement).

#ifdef __OPENMP
          !TODO envdw needs to be a array of size num_thrads
      !$omp task default(none) shared(rvdw, cell_list, l_stress,engvdw,virvdw) private(engacc, viracc) firstprivate(p) firstprivate(id, id2) private(tid) untied
            call cell_local_pair_vdw_non_ld_non_ls_interaction(rvdw,cell_list(id), cell_list(id2), p, engacc, viracc, l_stress)
            tid = omp_get_thread_num()
            engvdw(tid) = engvdw(tid) + engacc
            virvdw(tid) = virvdw(tid) + viracc
      !$omp end task
#else
            call cell_local_pair_vdw_non_ld_non_ls_interaction(rvdw,cell_list(id), cell_list(id2), p, engacc, viracc, stress)
            engvdw = engvdw + engacc
            virvdw = virvdw + viracc
#endif

            if(l_do_rdf) then 
              !$omp task shared(rcut, cell_list) firstprivate(p) firstprivate(id, id2) private(tid) untied
              Call pair_rdf_collect(rcut, cell_list(id), cell_list(id2), p)
              !$omp end task
            end if
          end if
        end if
      end do !p
      If(l_do_rdf .and. (.not. cell_list(id)%is_halo)) then
        !$omp task shared(rcut, cell_list) firstprivate(id, id2) private(tid) untied
        Call self_rdf_collect(rcut, cell_list(id)) 
        !$omp end task
      end if
    end do !k
  end do !j
end do !i
!$omp end do nowait
end if !ls_vdw
end if !if ld_vdw



!TODO NYI - Coulombic contributions (mximpl > 0 - not done in previous OpenMP version unless keyfce==2.
!!!!!!!!!!!!!!!!!!!!!!!!!
!COULOMBIC CONTRIBUTIONS!
!!!!!!!!!!!!!!!!!!!!!!!!!

If (mximpl > 0) Then
  !!TODO NYI MUILTPOLAR ATOMIC SITES
  print *, "Multipolar atomic sites are not yet implemented"
Else
!calculate coulombic forces, Ewald sum - real space contribution
  If (keyfce ==2) Then
!$omp do collapse(2)
    do i=0, x_cells+1
      do j=0, y_cells+1
        do k=0, z_cells+1
          id = 1+i+(x_cells+2)*(j+(y_cells+2)*k)
          if(.not. cell_list(id)%is_halo) then
!            call self_ewald_real_forces(rcut, alpha, epsq, engacc, viracc, stress,cell_list(id))
#ifdef __OPENMP
!TODO envdw needs to be a array of size num_thrads
      !$omp task default(none) shared(rcut, cell_list,l_stress,engcpe_rl,vircpe_rl, engcpe_ex, vircpe_ex) private(engacc,viracc,engacc2,viracc2) firstprivate(id) private(tid) shared(epsq,alpha)
            call self_ewald_real_forces(rcut, alpha, epsq, engacc, viracc, l_stress, cell_list(id), engacc2, viracc2)
            tid = omp_get_thread_num()
            engcpe_rl(tid) = engcpe_rl(tid) + engacc
            vircpe_rl(tid) = vircpe_rl(tid) + viracc
            engcpe_ex(tid) = engcpe_ex(tid) + engacc2
            vircpe_ex(tid) = vircpe_ex(tid) + viracc2
      !$omp end task
#else
            call self_ewald_real_forces(rcut, alpha, epsq, engacc, viracc, stress,cell_list(id), engacc2, viracc2)

            engcpe_rl=engcpe_rl+engacc
            vircpe_rl=vircpe_rl+viracc
            engcpe_ex=engcpe_ex+engacc2
            vircpe_ex=vircpe_ex+viracc2
#endif
          end if
          do p = 1,13
            id2 = neighbours(p,id)
            if(id2 < 0 .or. id2 > number_cells ) CYCLE
            if(cell_list(id)%is_halo) then
              if(.not. cell_list(id2)%is_halo) then
#ifdef __OPENMP
      !$omp task default(none) shared(rcut, cell_list,l_stress,engcpe_rl,vircpe_rl, engcpe_ex, vircpe_ex) private(engacc,viracc,engacc2,viracc2) firstprivate(p) firstprivate(id,id2) private(tid) shared(epsq,alpha)
                  call pair_ewald_real_Forces(rcut, alpha, epsq, engacc, viracc, l_stress, cell_list(id2), cell_list(id), p+13, engacc2, viracc2)
            tid = omp_get_thread_num()
            engcpe_rl(tid) = engcpe_rl(tid) + engacc
            vircpe_rl(tid) = vircpe_rl(tid) + viracc
            engcpe_ex(tid) = engcpe_ex(tid) + engacc2
            vircpe_ex(tid) = vircpe_ex(tid) + viracc2
      !$omp end task
#else
                  call pair_ewald_real_Forces(rcut, alpha, epsq, engacc, viracc, stress, cell_list(id2), cell_list(id), p+13, engacc2, viracc2)
                  engcpe_rl=engcpe_rl+engacc
                  vircpe_rl=vircpe_rl+viracc
                  engcpe_ex=engcpe_ex+engacc2
                  vircpe_ex=vircpe_ex+viracc2
#endif
              end if
            else
             if(cell_list(id2)%is_halo) then
               !If cell_list id2 is a halo cell then don't store forces on that cell.
#ifdef __OPENMP
      !$omp task default(none) shared(rcut, cell_list,l_stress,engcpe_rl,vircpe_rl, engcpe_ex, vircpe_ex) private(engacc,viracc,engacc2,viracc2) firstprivate(p) firstprivate(id,id2) private(tid)  shared(epsq,alpha)
                  call pair_ewald_real_Forces(rcut, alpha, epsq, engacc, viracc, l_stress, cell_list(id), cell_list(id2), p, engacc2, viracc2)
            tid = omp_get_thread_num()
            engcpe_rl(tid) = engcpe_rl(tid) + engacc
            vircpe_rl(tid) = vircpe_rl(tid) + viracc
            engcpe_ex(tid) = engcpe_ex(tid) + engacc2
            vircpe_ex(tid) = vircpe_ex(tid) + viracc2
      !$omp end task
#else
               call pair_ewald_real_forces(rcut, alpha, epsq, engacc, viracc, stress, cell_list(id), cell_list(id2), p, engacc2, viracc2)
               engcpe_rl=engcpe_rl+engacc
               vircpe_rl=vircpe_rl+viracc
               engcpe_ex=engcpe_ex+engacc2
               vircpe_ex=vircpe_ex+viracc2
#endif
             else
#ifdef __OPENMP
      !$omp task default(none) shared(rcut, cell_list,l_stress,engcpe_rl,vircpe_rl, engcpe_ex, vircpe_ex) private(engacc,viracc,engacc2,viracc2) firstprivate(p) firstprivate(id,id2) private(tid)  shared(epsq,alpha)
               !cell_list id2 is not a halo cell, store forces on both cells (remove if statement).
               call local_pair_ewald_real_forces(rcut, alpha, epsq, engacc, viracc, l_stress, cell_list(id), cell_list(id2), p, engacc2, viracc2)
            tid = omp_get_thread_num()
            engcpe_rl(tid) = engcpe_rl(tid) + engacc
            vircpe_rl(tid) = vircpe_rl(tid) + viracc
            engcpe_ex(tid) = engcpe_ex(tid) + engacc2
            vircpe_ex(tid) = vircpe_ex(tid) + viracc2
      !$omp end task
#else
               call local_pair_ewald_real_forces(rcut, alpha, epsq, engacc, viracc, stress, cell_list(id), cell_list(id2), p, engacc2, viracc2)
               engcpe_rl=engcpe_rl+engacc
               vircpe_rl=vircpe_rl+viracc
               engcpe_ex=engcpe_ex+engacc2
               vircpe_ex=vircpe_ex+viracc2
#endif
             end if
            end if
          end do !p
        end do !k
      end do !j
    end do !i
!$omp end do nowait  
  Else If(keyfce >0) then
    !!TODO NYI Implement other types
    print *, "only ewald_real_forces implemented, keyfce=", keyfce, "requested"
  End If
End If

!$omp end parallel
#ifdef __OPENMP
do i = 1, nr_threads-1
  engvdw(0) = engvdw(0) + engvdw(i)
  virvdw(0) = virvdw(0) + virvdw(i)
  engcpe_rl(0) = engcpe_rl(0) + engcpe_rl(i)
  vircpe_rl(0) = vircpe_rl(0) + vircpe_rl(i)
  engcpe_ex(0) = engcpe_ex(0) + engcpe_ex(i)
  vircpe_ex(0) = vircpe_ex(0) + vircpe_ex(i)
  l_stress(1:9,0) = l_stress(1:9,0) + l_stress(1:9,i)
end do
do i = 1,9
  stress(i) = l_stress(i,0)
end do
#endif

t_total = t_total + (omp_get_wtime() - t_start)

!TODO NYI - Poisson solver alternative to Ewald (never done with OpenMP)

!TODO NYI - Reverse excluded interactions. (Note this also needs doing for ewald forces somehow, NYI)
!TODO NYI - Reset frozen particles forces to 0.


!TODO NYI - Frozen pairs corrections to coulombic forces.

  If ( keyfce == 2 .or. keyfce == 12) Then
     If (l_fce) Then

! TODO frozen pairs corrections to coulombic forces

        If (megfrz /= 0) Then
           If (keyfce == 2) Then ! Ewald
              If (mximpl > 0) Then
                 Call ewald_frzn_mforces(rcut,alpha,epsq,engcpe_fr,vircpe_fr,stress)
              Else
                 Call ewald_frzn_forces(rcut,alpha,epsq,engcpe_fr,vircpe_fr,stress)
              End If
           Else !If (keyfce == 12) Then ! Poisson Solver
              Call poisson_frzn_forces(rcut,epsq,engcpe_fr,vircpe_fr,stress)
           End If
        End If

     Else

! Refresh all Ewald k-space contributions

        Call ewald_refresh(engcpe_rc,vircpe_rc,engcpe_fr,vircpe_fr,stress)

     End If

! non-zero total system charge correction (for the whole system)
! ( Fuchs, Proc. R. Soc., A, 151, (585),1935 )

     If (Abs(sumchg) > 1.0e-6_wp) Then
        If (new_nz) Then
           new_nz = .false.
           factor_nz = -0.5_wp * (pi*r4pie0/epsq) * (sumchg/alpha)**2
        End If

        engcpe_nz=factor_nz/volm
        vircpe_nz=-3.0_wp*engcpe_nz
     End If
  End If


!TODO NYI - Update particle forces to global stores
  do i=1,x_cells
    do j=1, y_cells
      do k = 1,z_cells
        id = 1+i+(x_cells+2)*(j+(y_cells+2)*k)
        Call update_forces(cell_list(id), fxx, fyy, fzz)
      end do
    end do
  end do

!TODO NYI - Line 698 to 771
  If (mxnode > 1) Then
     buffer( 0) = tmp 
     buffer( 1) = engkim
     buffer( 2) = virkim
     buffer( 3) = engden
     buffer( 4) = virden
     buffer( 5) = engmet
     buffer( 6) = virmet
#ifdef __OPENMP
     buffer( 7) = engvdw(0)
     buffer( 8) = virvdw(0)
#else     
     buffer( 7) = engvdw
     buffer( 8) = virvdw
#endif
     buffer( 9) = engcpe_rc
     buffer(10) = vircpe_rc
#ifdef __OPENMP
     buffer(11) = engcpe_rl(0)
     buffer(12) = vircpe_rl(0)
     buffer(13) = engcpe_ex(0)
     buffer(14) = vircpe_ex(0)
#else
     buffer(11) = engcpe_rl
     buffer(12) = vircpe_rl
     buffer(13) = engcpe_ex
     buffer(14) = vircpe_ex
#endif
     buffer(15) = engcpe_fr
     buffer(16) = vircpe_fr
     buffer(17) = vircpe_dt

     Call gsum(buffer(0:17))

     tmp       = buffer( 0)
     engkim    = buffer( 1)
     virkim    = buffer( 2)
     engden    = buffer( 3)
     virden    = buffer( 4)
     engmet    = buffer( 5)
     virmet    = buffer( 6)
#ifdef __OPENMP
     engvdw(0) = buffer( 7)
     virvdw(0) = buffer( 8)
#else     
     engvdw    = buffer( 7)
     virvdw    = buffer( 8)
#endif
     engcpe_rc = buffer( 9)
     vircpe_rc = buffer(10)
#ifdef __OPENMP
     engcpe_rl(0) = buffer(11)
     vircpe_rl(0) = buffer(12)
     engcpe_ex(0) = buffer(13)
     vircpe_ex(0) = buffer(14)
#else
     engcpe_rl = buffer(11)
     vircpe_rl = buffer(12)
     engcpe_ex = buffer(13)
     vircpe_ex = buffer(14)

#endif
     engcpe_fr = buffer(15)
     vircpe_fr = buffer(16)
     vircpe_dt = buffer(17)
  End If

  safe=(tmp < 0.5_wp)
  If (.not.safe) Call error(505)

!TODO NYI - Line 773 to End
! Self-interaction is constant for the default charges only SPME

  If (mxnode > 1 .and. keyfce == 2) Then ! Sum it up for multipolar SPME
     If (mximpl > 0 .and. mxompl <= 2) Call gsum(engsic)
  End If

! Globalise coulombic contributions: cpe
#ifdef __OPENMP
  engcpe = engcpe_rc + engcpe_rl(0) + engcpe_ex(0) + engcpe_fr + engcpe_nz
  vircpe = vircpe_rc + vircpe_rl(0) + vircpe_ex(0) + vircpe_fr + vircpe_nz + vircpe_dt
#else
  engcpe = engcpe_rc + engcpe_rl + engcpe_ex + engcpe_fr + engcpe_nz
  vircpe = vircpe_rc + vircpe_rl + vircpe_ex + vircpe_fr + vircpe_nz + vircpe_dt
#endif


! Add non-zero total system charge correction to
! diagonal terms of stress tensor (per node)

  tmp = - vircpe_nz/(3.0_wp*Real(mxnode,wp))
  stress(1) = stress(1) + tmp
  stress(5) = stress(5) + tmp
  stress(9) = stress(9) + tmp

! Globalise short-range, KIM and metal interactions with
! their long-range corrections contributions: srp
#ifdef __OPENMP
  engsrp = engkim + (engden + engmet + elrcm(0)) + (engvdw(0) + elrc)
  virsrp = virkim + (virden + virmet + vlrcm(0)) + (virvdw(0) + virlrc)
#else
  engsrp = engkim + (engden + engmet + elrcm(0)) + (engvdw + elrc)
  virsrp = virkim + (virden + virmet + vlrcm(0)) + (virvdw + virlrc)
#endif

! Add long-range corrections to diagonal terms of stress tensor (per node)
#ifdef __OPENMP
#ifdef TASK_TIMERS
if(timestep == 25) then
  Open(Unit=nrdfdt, File='TASKTIMERS', Status='replace')
do i=1,number_cells
  do j=1,27
    if(cell_list(i)%tid(j) >= 0) then
      write(nrdfdt,'(2i10, 2f25.15, 3i10)'), 1, cell_list(i)%tid(j),cell_list(i)%start(j)-step_start, cell_list(i)%finish(j)-step_start, cell_list(i)%cell_id, cell_list(neighbours(j,i))%cell_id, cell_list(i)%nr_yields(j)
    end if
  end do
  iF(cell_list(i)%tid(28) >= 0) then
    write(nrdfdt,'(2i10, 2f25.15, 3i10)'), 2, cell_list(i)%tid(28),cell_list(i)%start(28)-step_start, cell_list(i)%finish(28)-step_start, cell_list(i)%cell_id, cell_list(i)%cell_id, cell_list(i)%nr_yields(28)
  end if
end do


Close(Unit=nrdfdt)
!TODO OUTPUT
end if
#endif
#endif


  tmp = - (virlrc+vlrcm(0))/(3.0_wp*Real(mxnode,wp))
  stress(1) = stress(1) + tmp
  stress(5) = stress(5) + tmp
  stress(9) = stress(9) + tmp

End Subroutine two_body_cell_forces

End Module two_body_cell_forces_module
