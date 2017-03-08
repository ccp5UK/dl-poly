Module rdf_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring RDF property variables and arrays
! including USR (umbrella sampling restraint) RDF
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2016
! contrib   - a.b.g.chalk january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer,                        Save :: ncfrdf = 0 , &
                                          ntprdf = 0

  Integer,                        Save :: ncfusr = 0

  Real( Kind = wp ),              Save :: rusr   = 0.0_wp ! USR RDF cutoff

  Integer,           Allocatable, Save :: lstrdf(:)

  Real( Kind = wp ), Allocatable, Save :: rdf(:,:),usr(:)

  Real( Kind = wp ), Allocatable, Save :: block_averages(:,:,:,:)
  Integer, Parameter                   :: num_blocks = 25
  Integer, Save                        :: block_size
  Integer,                        Save :: block_number = 1
  Real( Kind = wp ), Allocatable, Save :: tmp_rdf(:,:,:)
  Logical,                        Save :: tmp_rdf_sync = .FALSE.
  Logical,        Save :: l_block = .FALSE., l_jack = .FALSE.

  Public :: allocate_rdf_arrays, allocate_block_average_array

Contains

  Subroutine allocate_rdf_arrays()

    Use setup_module, Only :mxrdf,mxgrdf,mxgusr

    Implicit None

    Integer, Dimension( 1:3 ) :: fail

    fail = 0

    Allocate (lstrdf(1:mxrdf),       Stat = fail(1))
    Allocate (rdf(1:mxgrdf,1:mxrdf), Stat = fail(2))
    Allocate (usr(1:mxgusr),         Stat = fail(3))

    If (Any(fail > 0)) Call error(1016)

    lstrdf = 0

    rdf = 0.0_wp ; usr = 0.0_wp

  End Subroutine allocate_rdf_arrays

Subroutine allocate_block_average_array(nstrun)

  Use site_module, Only: ntpatm
      Use setup_module, Only : mxrdf,mxgrdf
  Implicit None
  Integer, Intent( In ) :: nstrun
  Integer :: temp1, temp2
  
  Integer, Dimension( 1:2 ) :: fail
  block_size = nstrun/(num_blocks-1)
  if(block_size < 2) then
    block_size = 2
  endif

  temp1 = mxrdf + 16-Mod(mxrdf,16)
  temp2 = mxgrdf + 16-Mod(mxgrdf,16)
  Allocate(block_averages(1:ntpatm,1:ntpatm,1:mxgrdf,1:num_blocks+1), Stat = fail(1))
  Allocate(tmp_rdf( 1:temp2,1:temp1, 1:num_blocks+1 ), Stat = fail(2))

  If (Any(fail > 0)) Call error(1016)
  block_averages = 0.0_wp
  tmp_rdf = 0.0_wp

  End Subroutine allocate_block_average_array


End Module rdf_module
