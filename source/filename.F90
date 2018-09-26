!> Module maintaining access to DL_POLY files
!>
!> Copyright - Daresbury Laboratory
!
!> Author - J. Madge September 2018
Module filename
  Use kinds, only : wi,si
  Implicit None

  Private

  !> File unit number integer kind
  Integer, Parameter, Public :: UNIT_TYPE = si

  !> File data
  Type, Public :: file_type
    Private
    !> Filename
    Character(Len=1024), Public :: filename
    !> Fortran unit number, set with newunit=T%unit_no
    Integer(Kind = UNIT_TYPE), Public :: unit_no = -2
  Contains
    Procedure, Public :: init => file_type_init
    Procedure, Public :: rename => file_type_init
  End Type file_type

  ! Core file location keys
  !> CONTROL file
  Integer(Kind=UNIT_TYPE), Parameter, Public :: FILE_CONTROL = 1
  !> OUTPUT file
  Integer(Kind=UNIT_TYPE), Parameter, Public :: FILE_OUTPUT = 2
  !> CONFIG file
  Integer(Kind=UNIT_TYPE), Parameter, Public :: FILE_CONFIG = 3
  !> FIELD file
  Integer(Kind=UNIT_TYPE), Parameter, Public :: FILE_FIELD = 4
  !> STATS file
  Integer(Kind=UNIT_TYPE), Parameter, Public :: FILE_STATS = 5
  !> HISTORY file
  Integer(Kind=UNIT_TYPE), Parameter, Public :: FILE_HISTORY = 6
  !> HISTORF file
  Integer(Kind=UNIT_TYPE), Parameter, Public :: FILE_HISTORF = 7
  !> REVIVE file
  Integer(Kind=UNIT_TYPE), Parameter, Public :: FILE_REVIVE = 8
  !> REVOLD file
  Integer(Kind=UNIT_TYPE), Parameter, Public :: FILE_REVOLD = 9
  !> REVCON file
  Integer(Kind=UNIT_TYPE), Parameter, Public :: FILE_REVCON = 10

  !> Size of filename array
  Integer(Kind=wi), Parameter :: FILENAME_SIZE = 10

  Public :: default_filenames

Contains

  !> Initialise a file
  Subroutine file_type_init(T, filename)
    Class(file_type) :: T
    Character(Len=*), Intent(In) :: filename

    T%filename = Trim(filename)
  End Subroutine file_type_init

  !> Allocate filename array and set defaults
  Subroutine default_filenames(filenames)
    Type(file_type), Allocatable :: filenames(:)

    !> Default file names array
    Character(Len=1024), Dimension(FILENAME_SIZE) :: default_names

    Integer(Kind=wi) :: file_no

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

    ! Allocate filenames array
    If (Allocated(filenames)) Then
      Deallocate(filenames)
    End If
    Allocate(filenames(FILENAME_SIZE))

    ! Set default filenames
    Do file_no = 1, FILENAME_SIZE
      Call filenames(file_no)%init(default_names(file_no))
    End Do
  End Subroutine default_filenames
End Module filename
