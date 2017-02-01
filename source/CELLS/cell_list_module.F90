Module cell_list_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring variables required for the construction and
! storage of cell lists and sorted indices
!
! copyright - daresbury laboratory
! author    - a.b.g.chalk november 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use Kinds_f90
#ifdef __OPENMP
Use omp_lib
#endif
Implicit None

type Cells
  Real(Kind=wp), Dimension(3)  :: cell_centre !Store cell centre.
  Real(Kind=wp), Dimension(9)  :: cell_shape !Store cell shape (as vectors) - parallelpiped.
  Integer, Dimension(3) :: cell_position ! Store cell "position" - this is the integer value in the cell grid
  Integer             :: cell_start_index, cell_count !Cell start index and number of particles in it
  Integer             :: cell_id !Cell ID (used to compute position)
  Integer, Dimension(:,:), Allocatable :: sorted_indices !Arrays to store the sorted indices. These are not 1 to count, but rather cell_start_index to +cell_count
  Real(Kind=wp), Dimension(:,:), Allocatable :: vector_positions
  Real(Kind=wp), Dimension(:,:), Allocatable :: neighbour_vectors
  Integer, Dimension(:), Allocatable :: exclusion_list_key, exclusion_list_value
  Integer :: exclusion_count
  Logical :: is_halo

#ifdef __OPENMP
  Integer (Kind=omp_lock_kind) :: cell_lock
#ifdef TASK_TIMERS
  Integer, Dimension(1:28) :: tid
  Double Precision, Dimension(1:28) :: start, finish
  Integer, Dimension(1:28) :: nr_yields
#endif
#endif

!TODO NYI - Cell lock.
end type Cells

type (Cells), Allocatable, Dimension(:), Save:: cell_list
Integer , Save :: number_cells !Note that this is the number cells -1 since it doesn't count the "0 cell"
Real(Kind = wp) , Parameter :: half_plus = 0.5_wp
Real(Kind = wp), Dimension(:), Allocatable, Save :: xxt, yyt, zzt
Integer, Allocatable, Save, Dimension(:) :: domain_map !Domain map stores new index from original position.
Integer, Allocatable, Save, Dimension(:) :: new_to_old_map !Stores a mapping of new positions to original positions.
Real(Kind = wp), Dimension(:), Allocatable, Save :: temp_forces_x, temp_forces_y, temp_forces_z
Integer, Allocatable, Save, Dimension(:,:) :: neighbours
Integer, Save :: x_cells, y_cells, z_cells
Real(Kind = wp) :: shell_cutoff
Integer, Save :: excl_index
Real(Kind = wp), Save :: ewald_c, vir_c
! x_cells, y_cells, z_cells is number of non-halo cells
Integer, Save, Allocatable, Dimension(:) :: exclusion_table
Integer, Save, Allocatable,Dimension(:) :: table_position
Integer, Save, Allocatable, Dimension(:) :: excl_maxdiff
Integer, Save :: table_size
#ifdef __OPENMP
Double Precision :: t_start, t_total
Integer :: timestep = 0
Double Precision :: step_start
#endif

!We need to know the size of the halos for this algorithm, so we'll store them here
Real(Kind = wp), Save :: pos_x, pos_y, pos_z, neg_x, neg_y, neg_z

Contains

!#ifdef TODO
subroutine compute_exclusions2(lbook)
use config_module
Implicit None
Logical, Intent(in) :: lbook
Integer :: i,j,maxdiff, global_maxdiff, idi, idj
Integer :: ierr,a_size, diff, ind

if(.not. lbook) return

If(Allocated(exclusion_table)) then
  Deallocate(exclusion_table, table_position, excl_maxdiff, stat=ierr)
If(ierr /= 0) print *, "failed allocation"
End if

Allocate(table_position(1:nlast), excl_maxdiff(1:nlast), stat=ierr)
If(ierr /= 0) print *, "failed allocation"
table_position = 0
excl_maxdiff = 0

!First we need to compute the maxdiff between the global ids of the excluded
!pairs.
global_maxdiff = 0
do i=1, natms
  maxdiff = 0
  do j=1, lexatm(0,i)
    idi = ltg(i)
    idj = ltg(lexatm(j,i))
    if(idi - idj /= i - idj) print *, "i /= idi"
    maxdiff = max(maxdiff, abs(idi-idj))
  end do
  global_maxdiff = max(maxdiff, global_maxdiff)
  excl_maxdiff(domain_map(i)) = maxdiff
end do

!We need 2 int per 32 maxdiff (since store as a bit mask, and can be positive
!and negative direction)
!
a_size = 2* ((global_maxdiff)/32 +1)
table_size = a_size

Allocate(exclusion_table(1:(nlast*a_size)), stat = ierr)
if(ierr /=0 ) print *, "failed allocation"

exclusion_table = 0

!Create the exclusion table
do idi=1, natms
  i = domain_map(idi)
  !The start index in the table is new index position * a_size
  table_position(i) = 1+ (i-1)*a_size
  maxdiff = excl_maxdiff(i)
  do j=1, lexatm(0,idi)
    diff = ltg(lexatm(j, idi)) - idi
    !Store the eclusion in j-i+maxdiff
    diff = diff + maxdiff
    ind = diff/32
    diff = IAND(diff, 31)
    exclusion_table(table_position(i)+ind) = IOR(exclusion_table(table_position(i)+ind),LSHIFT( 1 , diff))
  end do
end do
end subroutine compute_exclusions2
!#endif

subroutine compute_exclusions(lbook)
use config_module
Implicit None
logical, intent(in) :: lbook
Integer :: i, j, counts, part, z, k, v

if(lbook) then
do i= 1, number_cells
counts = 0
if(cell_list(i)%is_halo) then
cell_list(i)%exclusion_count = 0
CYCLE
end if
!Loop over the cell and compute how many exclusions exist.
  do j = cell_list(i)%cell_start_index, cell_list(i)%cell_start_index+cell_list(i)%cell_count-1
    part = new_to_old_map(j)
    counts = counts + lexatm(0,part)
  end do

  call init_exclusion_list(cell_list(i), counts)
  z = 0
  do j = cell_list(i)%cell_start_index, cell_list(i)%cell_start_index+cell_list(i)%cell_count-1
    part = new_to_old_map(j) 
      do k = 1, lexatm(0,part)
        z = z + 1
!Calculate one-sided exclusions. Write to key, interact with value
        cell_list(i)%exclusion_list_key(z) = j
        cell_list(i)%exclusion_list_value(z) = domain_map(lexatm(k,part))
        excl_index = max(excl_index, abs(part - lexatm(k,part) ) )
      end do  
  end do
end do
else
  do i = 1, number_cells
    cell_list(i)%exclusion_count = 0
  end do
end if

end subroutine compute_exclusions

subroutine init_exclusion_list(this, list_size)
  Implicit none
  type(cells), Intent(InOut) :: this
  Integer, Intent(in) :: list_size
  Integer :: ierr

  If(allocated(this%exclusion_list_key)) then
    deallocate(this%exclusion_list_key, this%exclusion_list_value, stat = ierr)
  end if
  this%exclusion_count = list_size
  Allocate(this%exclusion_list_key(1:list_size),this%exclusion_list_value(1:list_size), stat = ierr)
end subroutine init_exclusion_list

subroutine compute_neighbours(nlx, nly, nlz, imcon, px, py, pz)

Implicit None
Integer, Intent(In) :: nlx, nly, nlz !nlx, nly, nlz say how many cells in each dimension
Integer, Intent(In) :: px, py, pz ! px, py, pz say how many ranks in each dimension
Integer, Intent(In) :: imcon !imcon tells us the periodicness of the box
Integer :: ierr
Integer :: x, y, z, id
Integer :: nx, nxx, ny, nz, nyy, nzz !Store the values of the neighbouring x/y/z position
!Each cell has at most 26 neighbours
If(Allocated(neighbours)) Then
Deallocate(neighbours, stat=ierr)
!If(ierr /= 0) error(0)
End If


Allocate(neighbours(1:26, 1:number_cells), stat = ierr )

neighbours = -1

!We may have some periodicity to work out.
!TODO If non-periodic then we should only have a halo neighbour if there is a halo there.
!TODO Negative value means no neighbour in this dimension
  Do z = 0, nlz+1
    Do y = 0, nly+1
      Do x = 0, nlx+1
          nx = Merge(x+1, number_cells+1, x < nlx+1)
          nxx = Merge(x-1, -1 * (nly+2) * (nlx+2) * number_cells -1, x>0)
          ny = Merge(y+1, number_cells+1, y < nly+1)
          nyy = Merge(y-1, -1* (nly+2) * number_cells -1, y>0)
          nz = Merge(z+1, number_cells+1, z < nlz+1)
          nzz = Merge(z-1, number_cells * -1, z > 0 )

        !Lets compute the neighbours now we have worked out any potential periodicity requirements.
        id = 1+x + (nlx+2)*(y+(nly+2)*z)
        neighbours(1, id) = 1 + x + (nlx+2)*((ny)+(nly+2)*z)! 0 1 0  !! Reverse is 0 -1 0
        neighbours(2, id) = 1 + x + (nlx+2)*(y+(nly+2)*(nz))! 0 0 1  !! Reverse is 0 0 -1
        neighbours(3, id) = 1 + x + (nlx+2)*((ny)+(nly+2)*(nz))! 0 1 1  !! Reverse is 0 -1 -1
        neighbours(4, id) = 1 + x + (nlx+2)*((ny)+(nly+2)*(nzz))! 0 1 -1 !! Reverse is 0 -1 1
        neighbours(5, id) = 1 + nx + (nlx+2)*(y + (nly+2)*(z))! 1 0 0  !! Reverse is -1 0 0
        neighbours(6, id) = 1 + nx + (nlx+2)*(y + (nly+2)*(nz))! 1 0 1  !! Reverse is -1 0 -1
        neighbours(7, id) = 1 + nx + (nlx+2)*((ny) + (nly+2)*z)! 1 1 0  !! Reverse is -1 -1 0
        neighbours(8, id) = 1 + nx + (nlx+2)*(y + (nly+2)*(nzz))! 1 0 -1 !! Reverse is -1 0 1
        neighbours(9, id) = 1 + nx + (nlx+2)*((nyy)+(nly+2)*(z))! 1 -1 0 !! Reverse is -1 1 0
        neighbours(10, id) = 1 + nx + (nlx+2)*((ny)+(nly+2)*(nz))! 1 1 1 !! Reverse is -1 -1 -1
        neighbours(11, id) = 1 + nx + (nlx+2)*((nyy)+(nly+2)*(nzz))! 1 -1 -1 !! Reverse is -1 1 1
        neighbours(12, id) = 1 + nx + (nlx+2)*((nyy)+(nly+2)*(nz))! 1 -1 1 !! Reverse is -1 1 -1
        neighbours(13, id) = 1 + nx + (nlx+2)*((ny)+(nly+2)*(nzz))! 1 1 -1 !! Reverse is -1 -1 1
      !Symmetric interactions mean cells only have 13 neighbours (plus their self interaction)
      !These are the 13 cells in "positive" directions.
      !TODO Due to sorting complications, it is easier to store all 26 neighbours.
      !Neighbours(i-13,id) is the reverse of neighbours(i,id) for 14-26.
        neighbours(14, id) = 1 + x + (nlx+2)*((nyy)+(nly+2)*z)! 0 -1 0
        neighbours(15, id) = 1 + x + (nlx+2)*(y+(nly+2)*(nzz))! 0 0 -1
        neighbours(16, id) = 1 + x + (nlx+2)*((nyy)+(nly+2)*(nzz)) ! 0 -1 -1
        neighbours(17, id) = 1 + x + (nlx+2)*((nyy)+(nly+2)*(nz)) ! 0 -1 1 
        neighbours(18, id) = 1 + nxx + (nlx+2)*(y+(nly+2)*z) ! -1 0 0
        neighbours(19, id) = 1 + nxx + (nlx+2)*(y+(nly+2)*(nzz)) ! -1 0 -1
        neighbours(20, id) = 1 + nxx + (nlx+2)*((nyy)+(nly+2)*z) ! -1 -1 0
        neighbours(21, id) = 1 + nxx + (nlx+2)*(y+(nly+2)*nz) ! -1 0 1
        neighbours(22, id) = 1 + nxx + (nlx+2)*(ny+(nly+2)*z) ! -1 1 0
        neighbours(23, id) = 1 + nxx + (nlx+2)*(nyy+(nly+2)*nzz) ! -1 -1 -1
        neighbours(24, id) = 1 + nxx + (nlx+2)*(ny+(nly+2)*nz) ! -1 1 1
        neighbours(25, id) = 1 + nxx + (nlx+2)*(ny+(nly+2)*nzz) ! -1 1 -1
        neighbours(26, id) = 1 + nxx + (nlx+2)*(nyy+(nly+2)*nz) ! -1 -1 1
      End Do
    End Do
  End Do
end subroutine compute_neighbours

subroutine create_force_store(count)

Implicit None
Integer, Intent(in) :: count
Integer :: ierr
ierr = 0

If(Allocated(temp_forces_x)) Then
Deallocate(temp_forces_x,temp_forces_y, temp_forces_z ,stat = ierr)
!If(ierr /=0 ) error(0)
End If


Allocate(temp_forces_x(1:count), temp_forces_y(1:count), temp_forces_z(1:count),stat = ierr)
!If(ierr /=0 ) error(0)

end subroutine create_force_store

subroutine reset_force_store()

Implicit None

temp_forces_x = 0.0_wp
temp_forces_y = 0.0_wp
temp_forces_z = 0.0_wp
ewald_c = 0.0_wp
vir_c = 0.0_wp
end subroutine

subroutine create_sorted_stores(count)

Implicit None
Integer, Intent(in) :: count
Integer :: ierr
ierr = 0

If(Allocated(xxt) ) Then
Deallocate(xxt, yyt, zzt,domain_map,new_to_old_map, stat = ierr)
!If(ierr /=0 ) error(0)
End If

Allocate(xxt(1:count), yyt(1:count), zzt(1:count),domain_map(1:count),new_to_old_map(1:count),stat = ierr)
!TODO REMVOE THIS
new_to_old_map(1:count) = -1
domain_map(1:count) = -1
!If(ierr /=0 ) error(0)
end subroutine create_sorted_stores

subroutine update_forces(cell_i, global_fx, global_fy, global_fz)
Use config_module
Implicit None
Type(cells), Intent(in) :: cell_i
Integer :: i, loop_end, pi
Real(Kind = wp), Dimension(:), Intent(inout) :: global_fx, global_fy, global_fz

loop_end = cell_i%cell_start_index + cell_i%cell_count -1

do i=cell_i%cell_start_index, loop_end
pi = new_to_old_map(i)
global_fx(pi) = global_fx(pi) +  temp_forces_x(i)
global_fy(pi) = global_fy(pi) +  temp_forces_y(i)
global_fz(pi) = global_fz(pi) +  temp_forces_z(i)
end do

end subroutine update_forces

!Places the updated halo particle data into the cells containing the particle.
!halo_start should be natms+1
subroutine perform_halo_mapping(xxx, yyy, zzz, halo_start, halo_count)
Implicit none
Real(Kind = wp), Dimension(:), Intent(in) :: xxx, yyy, zzz
Integer, Intent(in) :: halo_start, halo_count
Integer :: i, pos

do i=halo_start, halo_count
  pos = domain_map(i)
  xxt(pos) = xxx(i)
  yyt(pos) = yyy(i)
  zzt(pos) = zzz(i)
  
end do


end subroutine perform_halo_mapping

subroutine create_cells(cell_count) 

Implicit None

Integer, Intent(In) :: cell_count
Integer :: ierr

ierr = 0

If(Allocated(cell_list)) Then
Deallocate(cell_list, stat = ierr)
!if(ierr /= 0) call error(108)
End If

number_cells = cell_count
Allocate(cell_list(0:cell_count), stat = ierr)
!If(ierr /= 0) call error(108)


end subroutine create_cells

subroutine cell_init(this, centre, cshape, x,y,z, nlx, nly, halo)
Implicit None
type (Cells), Intent(InOut) :: this
Integer, Intent(In):: nlx, nly, x, y, z
Real(Kind=wp), Intent(In), Dimension(3) :: centre
Real(Kind=wp), Intent(In), dimension(9) :: cshape
Logical, Intent(In) :: halo

this%cell_centre = centre
this%cell_position(1) = x
this%cell_position(2) = y
this%cell_position(3) = z
this%cell_id = 1+x + (nlx+2)*(y+(nly+2)*z)
!if(this%cell_id == 17) print *, x, y, z
this%cell_shape = cshape
this%is_halo = halo

#ifdef __OPENMP
Call omp_init_lock(this%cell_lock)
#ifdef TASK_TIMERS
this%tid = -1
this%nr_yields = 0
#endif
#endif

end subroutine

subroutine cell_set_parts(this, c_count, start)
Implicit None
type (Cells), Intent(InOut) :: this
Integer, Intent(In) :: c_count, start
Integer :: ierr
ierr = 0


this%cell_count = c_count
this%cell_start_index = start
!Due to halo cells being potentially different sizes we store all 26 directions of sort.
Allocate(this%sorted_indices(1:c_count, 1:26), stat = ierr)
!If(ierr /= 0) call error(108)
Allocate(this%vector_positions(1:c_count, 1:26), stat = ierr)
Allocate(this%neighbour_vectors(1:3, 1:26), stat = ierr)
!If(ierr /= 0) call error(108)


end subroutine cell_set_parts 

subroutine domain_particle_to_cell(rcell, xxx, yyy, zzz, xdc, ydc, zdc, jx, jy, jz, icell, nlx, nly, nlx1s, nlx0e, nly1s, nly0e, nlz1s, nlz0e)
Implicit None
Real(Kind=wp), Intent(in), Dimension(9) :: rcell
Real(Kind=wp), Intent(In) :: xxx, yyy, zzz, xdc, ydc, zdc
Integer, Intent(In) :: jx, jy, jz, nlx, nly, nlx1s, nlx0e, nly1s, nly0e, nlz1s, nlz0e
Integer, Intent(Out) :: icell
Integer :: ix, iy, iz
Real(Kind=wp) :: x, y, z

!Convert atomic positions from MD cell centered to reduced space coordinates
  x=rcell(1)*xxx+rcell(4)*yyy+rcell(7)*zzz
  y=rcell(2)*xxx+rcell(5)*yyy+rcell(8)*zzz
  z=rcell(3)*xxx+rcell(6)*yyy+rcell(9)*zzz

!  print *, x, y, z, xdc, ydc, zdc, jx, jy, jz

! Get cell coordinates accordingly

     ix = Int(xdc*(x+0.5_wp)) + jx
     iy = Int(ydc*(y+0.5_wp)) + jy
     iz = Int(zdc*(z+0.5_wp)) + jz

! Correction for domain (idnode) only particles (1,natms) but due to
! some tiny numerical inaccuracy kicked into its halo link-cell space
! Put all particles in a bounded link-cell space: lower and upper bounds
! as follows nl_coordinate_0e+1 <= i_coordinate <= nl_coordinate_1s-1 !

     ix = Max( Min( ix , nlx1s-1) , nlx0e+1)
     iy = Max( Min( iy , nly1s-1) , nly0e+1)
     iz = Max( Min( iz , nlz1s-1) , nlz0e+1)

!Convert ix, iy, iz to find the containing cell of the particle

  icell = 1+ ix + ((nlx+2) * (iy + (nly+2)*iz) )


end subroutine domain_particle_to_cell

!!Subroutine just taken straight from link_cells
subroutine halo_particle_to_cell(rcell, xxx, yyy, zzz, xdc, ydc, zdc, jx, jy, jz, nlx0s, nly0s, nlz0s, nlx0e, nly0e, nlz0e, nlx1s, nly1s, nlz1s, nlx1e, nly1e, nlz1e, icell, nlx, nly, nlz )

Implicit None
Real(Kind=wp), Intent(in), Dimension(9) :: rcell
Real(Kind=wp), Intent(In) :: xxx, yyy, zzz, xdc, ydc, zdc
Integer, Intent(In) :: jx, jy, jz, nlx0s, nly0s, nlz0s, nlx0e, nly0e, nlz0e
Integer, Intent(In) :: nlx1s, nly1s, nlz1s, nlx1e, nly1e, nlz1e, nlx, nly, nlz
Integer, Intent(Out) :: icell
Real(Kind=wp) :: x1, y1, z1, dispx, dispy, dispz
Logical :: lx0, lx1, ly0, ly1, lz0, lz1
Real(Kind=wp) :: x, y,  z
Integer :: ix, iy, iz

!Convert atomic positions from MD cell centered to reduced space coordinates
  x=rcell(1)*xxx+rcell(4)*yyy+rcell(7)*zzz
  y=rcell(2)*xxx+rcell(5)*yyy+rcell(8)*zzz
  z=rcell(3)*xxx+rcell(6)*yyy+rcell(9)*zzz

!print *, xxx, yyy, zzz
!print *, x, y, z
 
!Get cell coordinates accordingly
     If (x > -half_plus) Then
        dispx=xdc*(x+0.5_wp)
        ix = Int(dispx) + jx
     Else
        dispx=xdc*Abs(x+0.5_wp)
        ix =-Int(dispx) + jx - 1
     End If
     If (y > -half_plus) Then
        dispy=ydc*(y+0.5_wp)
        iy = Int(dispy) + jy
     Else
        dispy=ydc*Abs(y+0.5_wp)
        iy =-Int(dispy) + jy - 1
     End If
     If( z > -half_plus) Then
        dispz=zdc*(z+0.5_wp)
        iz = Int(dispz) + jz
     Else
        dispz=zdc*Abs(z+0.5_wp)
        iz =-Int(dispz) + jz - 1
     End If
!print *, ix, iy, iz
! Exclude any negatively bound residual halo

     If (ix >= nlx0s .and. iy >= nly0s .and. iz >= nlz0s) Then

! Correction for halo particles (natms+1,nlast) of this domain
! (idnode) but due to some tiny numerical inaccuracy kicked into
! the domain only link-cell space

        lx0=(ix > nlx0e)
        lx1=(ix < nlx1s)
        ly0=(iy > nly0e)
        ly1=(iy < nly1s)
        lz0=(iz > nlz0e)
        lz1=(iz < nlz1s)
        If ( (lx0 .and. lx1) .and. &
             (ly0 .and. ly1) .and. &
             (lz0 .and. lz1) ) Then

! Put the closest to the halo coordinate in the halo

           x1=Abs(x-0.5_wp*Sign(1.0_wp,x))
           y1=Abs(y-0.5_wp*Sign(1.0_wp,y))
           z1=Abs(z-0.5_wp*Sign(1.0_wp,z))
           If      (x1 <= y1 .and. x1 <= z1) Then
              If (x < 0.0_wp) Then
                 ix=nlx0e
              Else
                 ix=nlx1s
              End If
           Else If (y1 <= x1 .and. y1 <= z1) Then
              If (y < 0.0_wp) Then
                 iy=nly0e
              Else
                 iy=nly1s
              End If
           Else
              If (z < 0.0_wp) Then
                 iz=nlz0e
              Else
                 iz=nlz1s
              End If
           End If
        End If

! Check for positively bound residual halo

        lx0=(ix < nlx0s)
        lx1=(ix > nlx1e)
        ly0=(iy < nly0s)
        ly1=(iy > nly1e)
        lz0=(iz < nlz0s)
        lz1=(iz > nlz1e)
        If ( .not. &
             (lx0 .or. lx1 .or. &
              ly0 .or. ly1 .or. &
              lz0 .or. lz1) ) Then

! Hypercube function transformation (counting starts from one
! rather than zero /map_domains/ and 2*nlp more link-cells per
! dimension are accounted /coming from the halo/)

           icell = 1 + ix + (nlx +2)*(iy + (nly+2)*iz)

        Else

! Put possible residual halo in cell=0

           icell = 0

        End If

     Else

! Put possible residual halo in cell=0

        icell = 0

     End If

!print *, icell, ix,nlx1e, iy,nly1e, iz, nlz1e
!print *, "!!!!!!!!!!!!!!!!!!!!!!!!!"

end subroutine halo_particle_to_cell

end module
