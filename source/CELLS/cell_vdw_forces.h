#ifdef __GFORTRAN__
#define PASTE2(x,y) x/**/y
#define PASTE(x,y) PASTE2(PASTE2(x,_),y)
#else
#define PASTE(x,y) x##_##y
#endif


#define _local_pair(N) PASTE(cell_local_pair_vdw,N)
#define _pair(N) PASTE(cell_pair_vdw,N)
#define _self(N) PASTE(cell_self_vdw,N)


Subroutine _local_pair(FUNC)(rvdw,cell_i, cell_j, neighbour_number, engvdw, virvdw, stress)

!Cell_i is the left cell (loop from count to 1), cell_j is right cell (loop from 1 to count). However
!since we look over axis+13 for cell_j we nede to loop also from count to 1 and multiply the sorted
!values by -1 when performing j-i<r_c

#ifdef __OPENMP
  Real(Kind = wp), Intent(Out) :: engvdw, virvdw
  Real(Kind = wp), Dimension(:, :), Allocatable, Intent(InOut) :: stress
  Integer :: tid
  Logical :: locked
#else
    Real(Kind = wp), Intent(Out) :: engvdw, virvdw
    Real(Kind = wp), Dimension(1:9), Intent(InOut) :: stress
#endif
  Real(Kind = wp), Intent(In) :: rvdw
  type(Cells), Intent(InOut) :: cell_i, cell_j
  Integer, Intent(In) :: neighbour_number
  Logical, Save :: newjob = .true.
  Real(Kind = wp), Save :: dlrpot, rdr

    Integer           :: mm,idi,ai,aj,jatm,key,k,l,ityp,n,m,i,j, pi, pj
    Integer :: itype ,jtype
    Real( Kind = wp ) :: rrr,rsq,ppp,gamma,eng,            &
                         r0,r0rn,r0rm,r_6,sor6,            &
                         rho,a,b,c,d,e0,kk,                &
                         nr,mr,rc,sig,eps,alpha,beta,      &
                         fix,fiy,fiz,fx,fy,fz,             &
                         gk,gk1,gk2,vk,vk1,vk2,t1,t2,t3,t, &
                         strs1,strs2,strs3,strs5,strs6,strs9
    Real( Kind = wp ) :: irrr, spos_i, spos_j, sqvdw, dx, dy, dz
    Integer :: j_loop_end, temp_end, dir2
  Integer, Dimension(1:table_size) :: excl_table
  Integer :: diff, maxdiff
  Logical :: excluded
  !$omp threadprivate(newjob, dlrpot, rdr)
#ifdef __OPENMP
  locked = .false.
  if(cell_i%cell_id < cell_j%cell_id) then
    do while(.not. locked)
      if(omp_test_lock(cell_i%cell_lock)) then
          if(omp_test_lock(cell_j%cell_lock)) then
          locked = .true.
          else
          Call omp_unset_lock(cell_i%cell_lock)
          end if
      end if
      if(.not. locked) then
#ifdef TASK_TIMERS
        !$omp atomic
        cell_i%nr_yields(neighbour_number) = cell_i%nr_yields(neighbour_number)+1
#endif
        !$omp taskyield
      end if
    end do
  else
    do while(.not. locked)
      if(omp_test_lock(cell_j%cell_lock)) then
        if(omp_test_lock(cell_i%cell_lock)) then
          locked = .true.
        else
          Call omp_unset_lock(cell_j%cell_lock)
        end if
      end if
      if(.not.locked ) then
#ifdef TASK_TIMERS
        !$omp atomic
        cell_i%nr_yields(neighbour_number) = cell_i%nr_yields(neighbour_number)+1
#endif
        !$omp taskyield
      end if
    end do
  endif

  tid = omp_get_thread_num()
#ifdef TASK_TIMERS
cell_i%tid(neighbour_number) = tid
cell_i%start(neighbour_number) = omp_get_wtime()
#endif
#endif


    If (newjob) Then
       newjob = .false.

       dlrpot = rvdw/Real(mxgvdw-4,wp)
       rdr    = 1.0_wp/dlrpot

  ! set cutoff condition
    End If
  maxdiff = 0
    sqvdw = rvdw*rvdw
  ! initialise potential energy and virial

    engvdw=0.0_wp
  virvdw=0.0_wp

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

dir2 = MOD(neighbour_number + 12 , 26)+1
j_loop_end = 1
temp_end = j_loop_end
!TODO NYI This will (probably) not vectorise - probably want a blocking algorithm of some sorts here.
Do i = cell_i%cell_count, 1, -1

  pi = cell_i%sorted_indices(i, neighbour_number)
  ai = ltype(new_to_old_map(pi))
  idi = ltg(new_to_old_map(pi))
  fix = temp_forces_x(pi)
  fiy = temp_forces_y(pi)
  fiz = temp_forces_z(pi)
  spos_i = cell_i%vector_positions(i, neighbour_number)

  !Load the exclusion table for this particle
  if(table_size > 0) then
  excl_table(1:table_size) = exclusion_table(table_position(pi):table_position(pi)+table_size-1)
  maxdiff = excl_maxdiff(pi)
  end if 

  Do j= cell_j%cell_count, j_loop_end, -1
    spos_j = cell_j%vector_positions(j, dir2) * -1.0_wp
  If( Abs(spos_j - spos_i) < shell_cutoff) then
    !Particle might be in range so we need to check actual distance
    pj = cell_j%sorted_indices(j, dir2)
    dx = xxt(pi) - xxt(pj)
    dy = yyt(pi) - yyt(pj)
    dz = zzt(pi) - zzt(pj)
    rsq = dx ** 2 + dy ** 2 + dz ** 2

    jatm = new_to_old_map(pj)
    aj = ltype(jatm)

    If (ai > aj) Then
      key=ai*(ai-1)/2 + aj
    Else
      key=aj*(aj-1)/2 + ai
    End If
     k=lstvdw(key)
    itype = ltpvdw(k)

    If(rsq < sqvdw) Then
      diff = ltg(jatm) - idi
      if(abs(diff) <= maxdiff) then
        diff  = diff + maxdiff
        excluded = iand(excl_table(rshift(diff, 5)+1), lshift(1, iand(diff, 31) )) /= 0
        If(excluded) CYCLE
      !Don't need to mark others as not excluded currently as we use CYCLE to avoid exclusions.
      end if
      !Particles are in range, need to interact
      rrr = sqrt(rsq)
      irrr = 1.0_wp/rrr

      Call FUNC( irrr, rrr, rdr, eng, gamma, itype,rvdw ,k)

      fx = gamma*dx
      fy = gamma*dy
      fz = gamma*dz

      fix = fix + fx
      fiy = fiy + fy
      fiz = fiz + fz

      temp_forces_x(pj) = temp_forces_x(pj) - fx
      temp_forces_y(pj) = temp_forces_y(pj) - fy
      temp_forces_z(pj) = temp_forces_z(pj) - fz

      !store the forces on pj.
      !add interaction energy
        engvdw = engvdw + eng
        !add virial
        virvdw = virvdw - gamma*rrr*rrr

        !add stress tensor
        strs1 = strs1 + dx*fx
        strs2 = strs2 + dx*fy
        strs3 = strs3 + dx*fz
        strs5 = strs5 + dy*fy
        strs6 = strs6 + dy*fz
        strs9 = strs9 + dz*fz

    End If !Particles are in range.

    Else
      temp_end = j
      Exit
    End If

  End Do !Loop over j
  j_loop_end = temp_end

!Update forces on pi
temp_forces_x(pi) = fix
temp_forces_y(pi) = fiy
temp_forces_z(pi) = fiz

End Do !Loop over cell i

!complete stress tensor
#ifdef __OPENMP
stress(1, tid) = stress(1, tid) + strs1
stress(2, tid) = stress(2, tid) + strs2
stress(3, tid) = stress(3, tid) + strs3
stress(4, tid) = stress(4, tid) + strs2
stress(5, tid) = stress(5, tid) + strs5
stress(6, tid) = stress(6, tid) + strs6
stress(7, tid) = stress(7, tid) + strs3
stress(8, tid) = stress(8, tid) + strs6
stress(9, tid) = stress(9, tid) + strs9
Call omp_unset_lock(cell_i%cell_lock)
Call omp_unset_lock(cell_j%cell_lock)
#ifdef TASK_TIMERS
cell_i%finish(neighbour_number) = omp_get_wtime()
#endif
#else
stress(1) = stress(1) + strs1
stress(2) = stress(2) + strs2
stress(3) = stress(3) + strs3
stress(4) = stress(4) + strs2
stress(5) = stress(5) + strs5
stress(6) = stress(6) + strs6
stress(7) = stress(7) + strs3
stress(8) = stress(8) + strs6
stress(9) = stress(9) + strs9
#endif

End Subroutine _local_pair(FUNC)


Subroutine _pair(FUNC)(rvdw,cell_i, cell_j, neighbour_number, engvdw, virvdw, stress)
!Cell_i is the left cell (loop from count to 1), cell_j is right cell (loop from 1 to count). However
!since we look over axis+13 for cell_j we nede to loop also from count to 1 and multiply the sorted
!values by -1 when performing j-i<r_c

#ifdef __OPENMP
  Real(Kind = wp), Intent(Out) :: engvdw, virvdw
  Real(Kind = wp), Dimension(:, :), Allocatable, Intent(InOut) :: stress
  Integer :: tid
  Logical :: locked
#else
  Real(Kind = wp), Intent(Out) :: engvdw, virvdw
  Real(Kind = wp), Dimension(1:9), Intent(InOut) :: stress
#endif
Real(Kind = wp), Intent(In) :: rvdw
type(Cells), Intent(InOut) :: cell_i, cell_j
Integer, Intent(In) :: neighbour_number
Logical, Save :: newjob = .true.
Real(Kind = wp), Save :: dlrpot, rdr

  Integer           :: mm,idi,ai,aj,jatm,key,k,l,ityp,n,m,i,j, pi, pj
  Integer :: itype ,jtype
  Real( Kind = wp ) :: rrr,rsq,ppp,gamma,eng,            &
                       r0,r0rn,r0rm,r_6,sor6,            &
                       rho,a,b,c,d,e0,kk,                &
                       nr,mr,rc,sig,eps,alpha,beta,      &
                       fix,fiy,fiz,fx,fy,fz,             &
                       gk,gk1,gk2,vk,vk1,vk2,t1,t2,t3,t, &
                       strs1,strs2,strs3,strs5,strs6,strs9
  Real( Kind = wp ) :: irrr, spos_i, spos_j, sqvdw, dx, dy, dz
  Integer :: j_loop_end, temp_end, dir2
Integer, Dimension(1:table_size) :: excl_table
Integer :: diff, maxdiff
Logical :: excluded


!$omp threadprivate(newjob, dlrpot, rdr)
#ifdef __OPENMP
locked = .false.
if(cell_i%is_halo) then
  do while(.not. locked)
    if(omp_test_lock(cell_j%cell_lock)) then
          locked = .true.
    end if
    if(.not. locked) then
#ifdef TASK_TIMERS
        !$omp atomic
        cell_i%nr_yields(neighbour_number) = cell_i%nr_yields(neighbour_number)+1
#endif
      !$omp taskyield
    end if
  end do
else
  do while(.not. locked)
     if(omp_test_lock(cell_i%cell_lock)) then
       locked = .true.
      end if
    if(.not.locked ) then
#ifdef TASK_TIMERS
        !$omp atomic
        cell_i%nr_yields(neighbour_number) = cell_i%nr_yields(neighbour_number)+1
#endif
      !$omp taskyield
    end if
  end do
end if

tid = omp_get_thread_num()
#ifdef TASK_TIMERS
cell_i%tid(neighbour_number) = tid
cell_i%start(neighbour_number) = omp_get_wtime()
#endif
#endif
  If (newjob) Then
     newjob = .false.

     dlrpot = rvdw/Real(mxgvdw-4,wp)
     rdr    = 1.0_wp/dlrpot

! set cutoff condition
  End If

  sqvdw = rvdw*rvdw

  maxdiff = 0

! initialise potential energy and virial

  engvdw=0.0_wp
  virvdw=0.0_wp

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

j_loop_end = 1
temp_end = j_loop_end
dir2 = MOD(neighbour_number + 12 , 26)+1
!TODO NYI This will (probably) not vectorise - probably want a blocking algorithm of some sorts here.
Do i = cell_i%cell_count, 1, -1

  pi = cell_i%sorted_indices(i, neighbour_number)
  ai = ltype(new_to_old_map(pi))
  idi = ltg(new_to_old_map(pi))
  fix = temp_forces_x(pi)
  fiy = temp_forces_y(pi)
  fiz = temp_forces_z(pi)
  spos_i = cell_i%vector_positions(i, neighbour_number)

  !Load the exclusion table for this particle
  if(table_size > 0) then
  excl_table(1:table_size) = exclusion_table(table_position(pi):table_position(pi)+table_size-1)
  maxdiff = excl_maxdiff(pi)
  end if

  Do j= cell_j%cell_count, j_loop_end, -1

    spos_j = cell_j%vector_positions(j, dir2) * -1.0_wp
  If( Abs(spos_j - spos_i) < shell_cutoff) then
    !Particle might be in range so we need to check actual distance

    pj = cell_j%sorted_indices(j, dir2)
    dx = xxt(pi) - xxt(pj)
    dy = yyt(pi) - yyt(pj)
    dz = zzt(pi) - zzt(pj)
    rsq = dx ** 2 + dy ** 2 + dz ** 2

    jatm = new_to_old_map(pj)
    aj = ltype(jatm)

    If (ai > aj) Then
      key=ai*(ai-1)/2 + aj
    Else
      key=aj*(aj-1)/2 + ai
    End If
     k=lstvdw(key)
    itype = ltpvdw(k)

    If(rsq < sqvdw) Then
      diff = ltg(jatm) - idi
      if(abs(diff) <= maxdiff) then
        diff  = diff + maxdiff
        excluded = iand(excl_table(rshift(diff, 5)+1), lshift(1, iand(diff, 31) )) /= 0
        If(excluded) CYCLE
      end if
      !Particles are in range, need to interact
      rrr = sqrt(rsq)
      irrr = 1.0_wp/rrr
      !TODO Intieract particles
      Call FUNC( irrr, rrr, rdr, eng, gamma, itype,rvdw,k )

      fx = gamma*dx
      fy = gamma*dy
      fz = gamma*dz

    !Using kahan summation algorithm to avoid numerical inaccuracy
      fix = fix + fx
      fiy = fiy + fy
      fiz = fiz + fz

      !TODO Does this ever evaluate to true?
      if( idi < ltg(jatm)) then

        !add interaction energy
        engvdw = engvdw + eng

        !add virial
        virvdw = virvdw - gamma*rrr*rrr

        !add stress tensor
        strs1 = strs1 + dx*fx
        strs2 = strs2 + dx*fy
        strs3 = strs3 + dx*fz
        strs5 = strs5 + dy*fy
        strs6 = strs6 + dy*fz
        strs9 = strs9 + dz*fz

      end if !Store more values if required
    End If !Particles are in range.

    Else
      temp_end = j
      Exit
    End If

  End Do !Loop over j

  j_loop_end = temp_end

  !Update forces on pi
temp_forces_x(pi) = fix
temp_forces_y(pi) = fiy
temp_forces_z(pi) = fiz


End Do !Loop over cell i

!complete stress tensor
#ifdef __OPENMP
stress(1, tid) = stress(1, tid) + strs1
stress(2, tid) = stress(2, tid) + strs2
stress(3, tid) = stress(3, tid) + strs3
stress(4, tid) = stress(4, tid) + strs2
stress(5, tid) = stress(5, tid) + strs5
stress(6, tid) = stress(6, tid) + strs6
stress(7, tid) = stress(7, tid) + strs3
stress(8, tid) = stress(8, tid) + strs6
stress(9, tid) = stress(9, tid) + strs9
if(.not. cell_i%is_halo) Call omp_unset_lock(cell_i%cell_lock)
if(.not. cell_j%is_halo) Call omp_unset_lock(cell_j%cell_lock)
#ifdef TASK_TIMERS
cell_i%finish(neighbour_number) = omp_get_wtime()
#endif
#else
  stress(1) = stress(1) + strs1
  stress(2) = stress(2) + strs2
  stress(3) = stress(3) + strs3
  stress(4) = stress(4) + strs2
  stress(5) = stress(5) + strs5
  stress(6) = stress(6) + strs6
  stress(7) = stress(7) + strs3
  stress(8) = stress(8) + strs6
  stress(9) = stress(9) + strs9
#endif


End Subroutine _pair(FUNC)


Subroutine _self(FUNC)(rvdw, cell_i, engvdw, virvdw, stress)
Implicit None
#ifdef __OPENMP
  Real(Kind = wp), Intent(Out) :: engvdw, virvdw
  Real(Kind = wp), Dimension(:, :), Allocatable, Intent(InOut) :: stress
  Integer :: tid
  Logical :: locked
#else
  Real(Kind = wp), Intent(Out) :: engvdw, virvdw
  Real(Kind = wp), Dimension(1:9), Intent(InOut) :: stress
#endif
Real(Kind = wp), Intent(In) :: rvdw
type(Cells), Intent(InOut) :: cell_i
Logical, Save :: newjob = .true.
Real(Kind = wp), Save :: dlrpot, rdr

  Integer           :: mm,idi,ai,aj,jatm,key,k,l,ityp,n,m,i,j, pi, pj
  Integer :: itype ,jtype, loop_start, loop_end
  Real( Kind = wp ) :: rrr,rsq,ppp,gamma,eng,            &
                       r0,r0rn,r0rm,r_6,sor6,            &
                       rho,a,b,c,d,e0,kk,                &
                       nr,mr,rc,sig,eps,alpha,beta,      &
                       fix,fiy,fiz,fx,fy,fz,             &
                       gk,gk1,gk2,vk,vk1,vk2,t1,t2,t3,t, &
                       strs1,strs2,strs3,strs5,strs6,strs9
  Real( Kind = wp ) :: irrr, spos_i, spos_j, sqvdw, dx, dy, dz
Integer, Dimension(1:table_size) :: excl_table
Integer :: diff, maxdiff
Logical :: excluded


!$omp threadprivate(newjob, dlrpot, rdr)
#ifdef __OPENMP
locked = .false.
  do while(.not. locked)
    locked = omp_test_lock(cell_i%cell_lock)
    if(.not. locked) then
#ifdef TASK_TIMERS
        !$omp atomic
        cell_i%nr_yields(28) = cell_i%nr_yields(28)+1
#endif
      !$omp taskyield
    end if
  end do

tid = omp_get_thread_num()
#ifdef TASK_TIMERS
cell_i%tid(28) = tid
cell_i%start(28) = omp_get_wtime()
#endif
#endif
  If (newjob) Then
     newjob = .false.

     dlrpot = rvdw/Real(mxgvdw-4,wp)
     rdr    = 1.0_wp/dlrpot

! set cutoff condition
  End If

  sqvdw = rvdw*rvdw
  maxdiff = 0
! initialise potential energy and virial

  engvdw=0.0_wp
  virvdw=0.0_wp

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

loop_start = cell_i%cell_start_index
loop_end = loop_start + cell_i%cell_count -1
!double loop over the particles in the cell. No need for sorted lists.
do i = loop_start, loop_end
  !Load details on particle i.
  pi = i
  ai = ltype(new_to_old_map(pi))
  idi = ltg(new_to_old_map(pi))
  fix = temp_forces_x(pi)
  fiy = temp_forces_y(pi)
  fiz = temp_forces_z(pi)

!Load the exclusion table for this particle
  if(table_size > 0) then
  excl_table(1:table_size) = exclusion_table(table_position(pi):table_position(pi)+table_size-1)
  maxdiff = excl_maxdiff(pi)
  end if

  !Loop over the remainder of the particles.
  do j = i+1, loop_end
    pj = j
    dx = xxt(pi) - xxt(pj)
    dy = yyt(pi) - yyt(pj)
    dz = zzt(pi) - zzt(pj)
    rsq = dx **2  + dy ** 2 + dz ** 2

    jatm = new_to_old_map(pj)
    aj = ltype(jatm)

    If (ai > aj) Then
      key=ai*(ai-1)/2 + aj
    Else
      key=aj*(aj-1)/2 + ai
    End If
     k=lstvdw(key)
    itype = ltpvdw(k)
   If(rsq < sqvdw) Then
      diff = ltg(jatm) - idi
      if(abs(diff) <= maxdiff) then
        diff  = diff + maxdiff
        excluded = iand(excl_table(rshift(diff, 5)+1), lshift(1, iand(diff, 31) )) /= 0
        If(excluded) CYCLE
      end if
      !Particles are in range, need to interact
      rrr = sqrt(rsq)
      irrr = 1.0_wp/rrr

      Call FUNC( irrr, rrr, rdr, eng, gamma, itype,rvdw,k )

      fx = gamma*dx
      fy = gamma*dy
      fz = gamma*dz

    !Using kahan summation algorithm to avoid numerical inaccuracy
      fix = fix + fx
      fiy = fiy + fy
      fiz = fiz + fz


      temp_forces_x(pj) = temp_forces_x(pj) - fx
      temp_forces_y(pj) = temp_forces_y(pj) - fy
      temp_forces_z(pj) = temp_forces_z(pj) - fz


      !j is local store the forces - we only do self interaction for non-halo cells.

        !add interaction energy
        engvdw = engvdw + eng

        !add virial
        virvdw = virvdw - gamma*rrr*rrr

        !add stress tensor
        strs1 = strs1 + dx*fx
        strs2 = strs2 + dx*fy
        strs3 = strs3 + dx*fz
        strs5 = strs5 + dy*fy
        strs6 = strs6 + dy*fz
        strs9 = strs9 + dz*fz

    End If !Particles are in range.



  end do
!Update forces on pi
temp_forces_x(pi) = fix
temp_forces_y(pi) = fiy
temp_forces_z(pi) = fiz

end do

!complete stress tensor
#ifdef __OPENMP
stress(1, tid) = stress(1, tid) + strs1
stress(2, tid) = stress(2, tid) + strs2
stress(3, tid) = stress(3, tid) + strs3
stress(4, tid) = stress(4, tid) + strs2
stress(5, tid) = stress(5, tid) + strs5
stress(6, tid) = stress(6, tid) + strs6
stress(7, tid) = stress(7, tid) + strs3
stress(8, tid) = stress(8, tid) + strs6
stress(9, tid) = stress(9, tid) + strs9
Call omp_unset_lock(cell_i%cell_lock)
#ifdef TASK_TIMERS
cell_i%finish(28) = omp_get_wtime()
#endif
#else
  stress(1) = stress(1) + strs1
  stress(2) = stress(2) + strs2
  stress(3) = stress(3) + strs3
  stress(4) = stress(4) + strs2
  stress(5) = stress(5) + strs5
  stress(6) = stress(6) + strs6
  stress(7) = stress(7) + strs3
  stress(8) = stress(8) + strs6
  stress(9) = stress(9) + strs9
#endif



End Subroutine _self(FUNC)
