Module filename
!> Module maintaining access to DL_POLY files
!>
!> Copyright - Daresbury Laboratory
!
!> Author - J. Madge September 2018
!> contrib - a.m.elena October 2018 - use standard integer for units
  Use kinds, Only: wi

  Implicit None

  Private

  !> File data
  Type, Public :: file_type
    Private

    !> Filename
    Character(Len=1024), Public :: filename
    !> Fortran unit number, set with newunit=T%unit_no
    Integer, Public             :: unit_no = -2

  Contains
    Procedure, Public :: init => file_type_init
    Procedure, Public :: rename => file_type_init
    Procedure, Public :: Close => close_file
  End Type file_type

  ! Core file location keys
  !> CONTROL file
  Integer, Parameter, Public :: FILE_CONTROL = 1
  !> OUTPUT file
  Integer, Parameter, Public :: FILE_OUTPUT = 2
  !> CONFIG file
  Integer, Parameter, Public :: FILE_CONFIG = 3
  !> FIELD file
  Integer, Parameter, Public :: FILE_FIELD = 4
  !> STATS file
  Integer, Parameter, Public :: FILE_STATS = 5
  !> HISTORY file
  Integer, Parameter, Public :: FILE_HISTORY = 6
  !> HISTORF file
  Integer, Parameter, Public :: FILE_HISTORF = 7
  !> REVIVE file
  Integer, Parameter, Public :: FILE_REVIVE = 8
  !> REVOLD file
  Integer, Parameter, Public :: FILE_REVOLD = 9
  !> REVCON file
  Integer, Parameter, Public :: FILE_REVCON = 10
  !> CURRENT file
  Integer, Parameter, Public :: FILE_CURRENT = 11
  !> KPOINTS file
  Integer, Parameter, Public :: FILE_KPOINTS = 12
  !> RDF file
  Integer, Parameter, Public :: FILE_RDF = 13
  !> MSD file
  Integer, Parameter, Public :: FILE_MSD = 14

  !> Size of filename array
  Integer(Kind=wi), Parameter, Public :: FILENAME_SIZE = 14

  Public :: default_filenames

Contains

  !> Initialise a file
  Subroutine file_type_init(T, filename)
    Class(file_type)                :: T
    Character(Len=*), Intent(In   ) :: filename

    T%filename = Trim(filename)
  End Subroutine file_type_init

  !> Allocate filename array and set defaults
  Subroutine default_filenames(filenames)
    Type(file_type) :: filenames(FILENAME_SIZE)

    Character(Len=1024), Dimension(FILENAME_SIZE) :: default_names
    Integer(Kind=wi)                              :: file_no

!> Default file names array

    ! Populate default names array
    default_names(FILE_CONTROL) = "CONTROL"
    default_names(FILE_OUTPUT) = "OUTPUT"
    default_names(FILE_CONFIG) = "CONFIG"
    default_names(FILE_FIELD) = "FIELD"
    default_names(FILE_STATS) = "STATIS"
    default_names(FILE_HISTORY) = "HISTORY"
    default_names(FILE_HISTORF) = "HISTORF"
    default_names(FILE_REVIVE) = "REVIVE"
    default_names(FILE_REVOLD) = "REVOLD"
    default_names(FILE_REVCON) = "REVCON"
    default_names(FILE_CURRENT) = "CURRENTS"
    default_names(FILE_KPOINTS) = "KPOINTS"
    default_names(FILE_RDF) = "RDFDAT"
    default_names(FILE_MSD) = "MSDTMP"

    ! Set default filenames
    Do file_no = 1, FILENAME_SIZE
      Call filenames(file_no)%init(default_names(file_no))
    End Do
  End Subroutine default_filenames

  !> close a unit and restore it default value
  Subroutine close_file(T)
    Class(file_type) :: T

    Logical :: is_open

    Inquire (T%unit_no, opened=is_open)
    If (is_open) Then
      Close (T%unit_no)
      T%unit_no = -2
    End If
  End Subroutine close_file

End Module filename
