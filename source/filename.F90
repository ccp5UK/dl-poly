!> Module maintaining access to DL_POLY files
!>
!> Copyright - Daresbury Laboratory
!
!> Author - J. Madge September 2018
!> contrib - a.m.elena October 2018 - use standard integer for units
!> contrib - i.Scivetti Aug       2018 - addition of extra files for EVB calculations
Module filename
  Use kinds, Only: wi, STR_FILENAME

  Implicit None

  Private

  !> File data
  Type, Public :: file_type
    Private

    !> Filename
    Character(Len=STR_FILENAME), Public :: filename
    !> Fortran unit number, set with newunit=T%unit_no
    Integer, Public             :: unit_no = -2

  Contains
    Procedure, Public :: init => file_type_init
    Procedure, Public :: rename => file_type_init
    Procedure, Public :: Close => close_file
    Procedure, Public :: is_null => file_type_null
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
  !> TABBND file
  Integer, Parameter, Public :: FILE_TABBND = 15
  !> TABANG file
  Integer, Parameter, Public :: FILE_TABANG = 16
  !> TABDIH file
  Integer, Parameter, Public :: FILE_TABDIH = 17
  !> TABINV file
  Integer, Parameter, Public :: FILE_TABINV = 18
  !> TABVDW file
  Integer, Parameter, Public :: FILE_TABVDW = 19
  !> TABEAM file
  Integer, Parameter, Public :: FILE_TABEAM = 20

  !> SETEVB file
  Integer, Parameter, Public :: FILE_SETEVB = 21
  !> POPEVB file
  Integer, Parameter, Public :: FILE_POPEVB = 22
  !> FIELD2 file
  Integer, Parameter, Public :: FILE_FIELD_2=  23
  !> CONFIG2 file
  Integer, Parameter, Public :: FILE_CONFIG_2 = 24
  !> REVCON2 file
  Integer, Parameter, Public :: FILE_REVCON_2 = 25
  !> FIELD3 file
  Integer, Parameter, Public :: FILE_FIELD_3 =  26
  !> CONFIG3 file
  Integer, Parameter, Public :: FILE_CONFIG_3 = 27
  !> REVCON3 file
  Integer, Parameter, Public :: FILE_REVCON_3 = 28
  !> COR
  Integer, Parameter, Public :: FILE_COR = 29
  !> Size of filename array
  Integer(Kind=wi), Parameter, Public :: FILENAME_SIZE = 29

  Public :: default_filenames

Contains

  !> Initialise a file
  Subroutine file_type_init(T, filename)
    Class(file_type)                :: T
    Character(Len=*), Intent(In   ) :: filename

    T%filename = Trim(filename)
  End Subroutine file_type_init

  !> Determine if file is null
  Function file_type_null(T)
    Class(file_type) :: T
    logical :: file_type_null

    file_type_null = T%filename == 'NONE' .or. &
         T%filename == 'NUL' .or. &
         T%filename == '/dev/null'
  End Function file_type_null

  !> Allocate filename array and set defaults
  Subroutine default_filenames(filenames)
    Type(file_type) :: filenames(FILENAME_SIZE)

    Character(Len=STR_FILENAME), Dimension(FILENAME_SIZE) :: default_names
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
    default_names(FILE_TABBND) = "TABBND"
    default_names(FILE_TABANG) = "TABANG"
    default_names(FILE_TABDIH) = "TABDIH"
    default_names(FILE_TABINV) = "TABINV"
    default_names(FILE_TABVDW) = "TABLE"
    default_names(FILE_TABEAM) = "TABEAM"

    default_names(FILE_SETEVB)   = "SETEVB"
    default_names(FILE_POPEVB)   = "POPEVB"

    default_names(FILE_FIELD_2)  = "FIELD2"
    default_names(FILE_CONFIG_2) = "CONFIG2"
    default_names(FILE_REVCON_2) = "REVCON2"
    default_names(FILE_FIELD_3)  = "FIELD3"
    default_names(FILE_CONFIG_3) = "CONFIG3"
    default_names(FILE_REVCON_3) = "REVCON3"
    default_names(FILE_COR) = "COR"

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
