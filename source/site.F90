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

    Allocate (T%mol_name(1:mxtmls), stat=fail(1))
    Allocate (T%site_name(1:T%max_site),T%unique_atom(1:T%max_site),T%unique_shell(1:T%max_site), stat=fail(2))
    Allocate (T%num_mols(1:mxtmls),T%num_site(1:mxtmls),T%num_freeze(1:mxtmls), stat=fail(3))
    Allocate (T%type_site(1:T%max_site),T%freeze_site(1:T%max_site),T%free_site(1:T%max_site), stat=fail(4))
    Allocate (T%weight_site(1:T%max_site),T%charge_site(1:T%max_site),T%dof_site(1:T%max_site), stat=fail(5))
    Allocate (T%num_type(1:mxatyp),T%num_type_nf(1:mxatyp),T%dens(1:mxatyp), stat=fail(6))

    If (Any(fail > 0)) Call error(1026)

    T%mol_name = ' '
    T%site_name = ' '
    T%unique_atom = ' '
    T%unique_shell = ' '

    T%num_mols = 0
    T%num_site = 0
    T%num_freeze = 0
    T%type_site  = 0
    T%freeze_site = 0
    T%free_site = 0

    T%weight_site = 0.0_wp
    T%charge_site = 0.0_wp
    T%dof_site = 0.0_wp
    T%num_type = 0.0_wp
    T%num_type_nf = 0.0_wp
    T%dens = 0.0_wp
  End Subroutine allocate_site_arrays

  !> Deallocate site_type arrays
  Subroutine cleanup(T)
    Type( site_type ) :: T

    If (Allocated(T%mol_name)) Then
      Deallocate(T%mol_name)
    End If

    If (Allocated(T%site_name)) Then
      Deallocate(T%site_name)
    End If
    If (Allocated(T%unique_atom)) Then
      Deallocate(T%unique_atom)
    End If
    If (Allocated(T%unique_shell)) Then
      Deallocate(T%unique_shell)
    End If

    If (Allocated(T%num_mols)) Then
      Deallocate(T%num_mols)
    End If
    If (Allocated(T%num_site)) Then
      Deallocate(T%num_site)
    End If
    If (Allocated(T%num_freeze)) Then
      Deallocate(T%num_freeze)
    End If

    If (Allocated(T%type_site)) Then
      Deallocate(T%type_site)
    End If
    If (Allocated(T%freeze_site)) Then
      Deallocate(T%freeze_site)
    End If
    If (Allocated(T%free_site)) Then
      Deallocate(T%free_site)
    End If

    If (Allocated(T%weight_site)) Then
      Deallocate(T%weight_site)
    End If
    If (Allocated(T%charge_site)) Then
      Deallocate(T%charge_site)
    End If
    If (Allocated(T%dof_site)) Then
      Deallocate(T%dof_site)
    End If

    If (Allocated(T%num_type)) Then
      Deallocate(T%num_type)
    End If
    If (Allocated(T%num_type_nf)) Then
      Deallocate(T%num_type_nf)
    End If
    If (Allocated(T%dens)) Then
      Deallocate(T%dens)
    End If
  End Subroutine cleanup
End Module site
