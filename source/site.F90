Module site
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global atomic site variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2008
! contrib   - m.a.seaton june 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp,wi
  Use errors_warnings, Only : error
  Implicit None

  Private

  !> Type to store site data
  Type, Public :: site_type
    Private

    !> Maximum size of site arrays
    Integer( Kind = wi ), Public :: max_site

    !> Number of molecule types
    Integer( Kind = wi ), Public :: ntype_mol = 0
    !> Number of atom types
    Integer( Kind = wi ), Public :: ntype_atom = 0
    !> Number of shell types
    Integer( Kind = wi ), Public :: ntype_shell = 0

    !> Molecule names
    Character( Len = 40 ), Allocatable, Public :: mol_name(:)

    !> Site names
    Character( Len = 8 ), Allocatable, Public :: site_name(:)
    !> Unique atom labels
    Character( Len = 8 ), Allocatable, Public :: unique_atom(:)
    !> Unique shell names
    Character( Len = 8 ), Allocatable, Public :: unique_shell(:)

    !> Number of molecules of each type
    Integer( Kind = wi ), Allocatable, Public :: num_mols(:)
    !> NUmber of sites (atoms) in molecules of each type?
    Integer( Kind = wi ), Allocatable, Public :: num_site(:)
    !> Number of frozen molecules of each type
    Integer( Kind = wi ), Allocatable, Public :: num_freeze(:)

    !> Type status
    Integer( Kind = wi ), Allocatable, Public :: type_site(:)
    !> Frozen status
    Integer( Kind = wi ), Allocatable, Public :: freeze_site(:)
    !> Free status?
    Integer( Kind = wi ), Allocatable, Public :: free_site(:)
    !> Weight
    Real( Kind = wp ), Allocatable, Public :: weight_site(:)
    !> Charge
    Real( Kind = wp ), Allocatable, Public :: charge_site(:)
    !> Degrees of freedom
    Real( Kind = wp ), Allocatable, Public :: dof_site(:)

    !> Number of atoms of each type?
    Real( Kind = wp ), Allocatable, Public :: num_type(:)
    !> Number of non-frozen atoms of each type?
    Real( Kind = wp ), Allocatable, Public :: num_type_nf(:)
    !> Density of atoms of each type?
    Real( Kind = wp ), Allocatable, Public :: dens(:)

  Contains
    Private
    Procedure, Public :: init => allocate_site_arrays
    Final :: cleanup
  End Type site_type

Contains

  !> Allocate and initialise site_type arrays
  Subroutine allocate_site_arrays(T,mxtmls,mxatyp)
    Class( site_type ) :: T
    Integer( Kind = wi ), Intent( In    ) :: mxtmls,mxatyp
    Integer, Dimension(1:6) :: fail

    fail = 0

    Allocate (T%molnam(1:mxtmls), stat=fail(1))
    Allocate (T%sitnam(1:T%max_site),T%unqatm(1:T%max_site),T%unqshl(1:T%max_site), stat=fail(2))
    Allocate (T%nummols(1:mxtmls),T%numsit(1:mxtmls),T%numfrz(1:mxtmls), stat=fail(3))
    Allocate (T%typsit(1:T%max_site),T%frzsit(1:T%max_site),T%fresit(1:T%max_site), stat=fail(4))
    Allocate (T%wgtsit(1:T%max_site),T%chgsit(1:T%max_site),T%dofsit(1:T%max_site), stat=fail(5))
    Allocate (T%numtyp(1:mxatyp),T%numtypnf(1:mxatyp),T%dens(1:mxatyp), stat=fail(6))

    If (Any(fail > 0)) Call error(1026)

    T%molnam = ' '
    T%sitnam = ' '
    T%unqatm = ' '
    T%unqshl = ' '

    T%nummols = 0
    T%numsit = 0
    T%numfrz = 0
    T%typsit  = 0
    T%frzsit = 0
    T%fresit = 0

    T%wgtsit = 0.0_wp
    T%chgsit = 0.0_wp
    T%dofsit = 0.0_wp
    T%numtyp = 0.0_wp
    T%numtypnf = 0.0_wp
    T%dens = 0.0_wp
  End Subroutine allocate_site_arrays

  !> Deallocate site_type arrays
  Subroutine cleanup(T)
    Type( site_type ) :: T

    If (Allocated(T%molnam)) Then
      Deallocate(T%molnam)
    End If

    If (Allocated(T%sitnam)) Then
      Deallocate(T%sitnam)
    End If
    If (Allocated(T%unqatm)) Then
      Deallocate(T%unqatm)
    End If
    If (Allocated(T%unqshl)) Then
      Deallocate(T%unqshl)
    End If

    If (Allocated(T%nummols)) Then
      Deallocate(T%nummols)
    End If
    If (Allocated(T%numsit)) Then
      Deallocate(T%numsit)
    End If
    If (Allocated(T%numfrz)) Then
      Deallocate(T%numfrz)
    End If

    If (Allocated(T%typsit)) Then
      Deallocate(T%typsit)
    End If
    If (Allocated(T%frzsit)) Then
      Deallocate(T%frzsit)
    End If
    If (Allocated(T%fresit)) Then
      Deallocate(T%fresit)
    End If

    If (Allocated(T%wgtsit)) Then
      Deallocate(T%wgtsit)
    End If
    If (Allocated(T%chgsit)) Then
      Deallocate(T%chgsit)
    End If
    If (Allocated(T%dofsit)) Then
      Deallocate(T%dofsit)
    End If

    If (Allocated(T%numtyp)) Then
      Deallocate(T%numtyp)
    End If
    If (Allocated(T%numtypnf)) Then
      Deallocate(T%numtypnf)
    End If
    If (Allocated(T%dens)) Then
      Deallocate(T%dens)
    End If
  End Subroutine cleanup
End Module site
