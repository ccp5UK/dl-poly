#ifdef KIM
#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__
#define TRUEFALSE(TRUTH) merge(1,0,(TRUTH))
#endif

Module kim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global KIM interaction variables and
! arrays
!
! copyright - daresbury laboratory
! author    - r.s.elliott march 2015
! contrib   - h.boateng & i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use, Intrinsic :: iso_c_binding
  Use kinds, Only : wp,wi
  Use neighbours, Only : neighbours_type
  Use errors_warnings, Only : error
#ifdef KIM
  Use KIM_API_F03
  Use domains, Only : map
  Use configuration,  Only : natms,nlast,lsi,lsa,ltg,lsite, &
                             xxx,yyy,zzz,fxx,fyy,fzz
  Use setup,   Only : mxsite,mxatdm,mxbfxp
  Use site,    Only : unqatm,ntpatm,sitnam
  Use comms, Only : comms_type,export_tag,wp_mpi,gsend,gwait,girecv
  Use numerics, Only : local_index
#else  
  Use comms, Only : comms_type
#endif

  Implicit None

  Character( Len = 200 ), Save :: kimim  = ' '      ! KIM IM type for dl_poly
  Real( Kind = wp ),      Save :: rkim = 0.0_wp   ! KIM cutoff for dl_poly
  Integer,                Save :: idhalo(0:2,1:6) ! KIM halo indicator

#ifdef KIM
  Type( c_ptr ),    Save          :: pkim
  Integer( Kind = c_int ), Save, Pointer :: kim_list(:,:)
  Real( Kind = c_double ), Save, Pointer :: kim_Rij(:,:,:)
  Logical,                 Save          :: HalfList, RijNeeded

! variable for use with reverse communication process

  Integer, Parameter                   :: iadd = 4  ! number of scalars per particle to be sent
  Integer, Save                        :: iblock    ! comm buffer half-size
  Real( Kind = wp ), Save, Allocatable :: rev_comm_buffer(:)  ! comm buffer
#endif
  Public :: kim_cutoff,  &
            kim_setup,   &
            kim_forces,  &
            kim_cleanup

#ifdef KIM
  Private :: kim_basic_init, &
             kim_reverse_communication
#endif

Contains

  Subroutine kim_cutoff(num_types,model_types,model_name,cutoff,comm)

!-------------------------------------------------------------------------------
!
! kim_cutoff
!
! This function extracts the cutoff distance for the specified KIM model_name
!
!-------------------------------------------------------------------------------

    Integer( Kind = c_int ), Intent( In    ) :: num_types
    Character( Len = * ),    Intent( In    ) :: model_types(1:num_types)
    Character( Len = * ),    Intent( In    ) :: model_name
    Real( Kind = wp ),       Intent(   Out ) :: cutoff
    Type(comms_type),        Intent( InOut ) :: comm

#ifdef KIM
    Integer( Kind = c_int )          :: ier, idum
    Real( Kind = c_double ), Pointer :: model_cutoff; Type(c_ptr) :: pcut

    ier = KIM_STATUS_OK

    ier = kim_basic_init(num_types, model_types, Trim(model_name))
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_basic_init", ier)
       Call error(0,'kim_basic_init failed in kim_cutoff')
    End If

! Allocate with 1 particle and 1 species, just to get the cutoff

    Call kim_api_allocate(pkim, 1, 1, ier)
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_allocate", ier)
       Call error(0,'kim_api_allocate failed in kim_cutoff')
    End If

    ier = kim_api_model_init(pkim)
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_model_init", ier)
       Call error(0,'kim_api_model_init failed in kim_cutoff')
    End If

    pcut = kim_api_get_data(pkim, "cutoff", ier)
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_getm_data", ier)
       Call error(0,'kim_api_getm_data failed in kim_cutoff')
    End If
    Call c_f_Pointer(pcut, model_cutoff)

    cutoff = model_cutoff

    ier = kim_api_model_destroy(pkim)
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_model_destroy", ier)
       Call error(0,'kim_api_model_destroy failed in kim_cutoff')
    End If

    Call kim_api_free(pkim, ier)
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_free", ier)
       Call error(0,'kim_api_free failed in kim_cutoff')
    End If
    pkim = c_null_ptr
#else
    cutoff = 0.0_wp
    Call kim_message()
#endif
  End Subroutine  kim_cutoff

  Subroutine kim_setup(num_types,model_types,model_name,max_list,comm)

!-------------------------------------------------------------------------------
!
! kim_setup
!
! This function creates the interface via the KIM API to KIM model_name
!
!-------------------------------------------------------------------------------

    Character(Len = *),      Intent( In    ) :: model_name
    Integer( Kind = c_int ), Intent( In    ) :: num_types
    Character( Len = * ),    Intent( In    ) :: model_types(1:num_types)
    Integer( Kind = wi ),    Intent( In    ) :: max_list
    Type(comms_type),        Intent( InOut ) :: comm

#ifdef KIM
    Integer( Kind = c_int ) :: ier,idum
    Integer                 :: fail,i,limit

    Character( Len = KIM_KEY_STRING_LENGTH ) :: ActiveNBC_Method
    Integer( Kind = c_int ), Pointer :: numberOfParticles;           Type( c_ptr ) :: pNOP
    Integer( Kind = c_int ), Pointer :: numberOfSpecies;             Type( c_ptr ) :: pNPT
    Integer( Kind = c_int ), Pointer :: numberContributingParticles; Type( c_ptr ) :: pNCP
    Integer( Kind = c_int ), Pointer :: particleSpecies(:);          Type( c_ptr ) :: pPT

    Character( Len = 256 ) :: message

    ier = KIM_STATUS_OK

    ier = kim_basic_init(num_types, model_types, Trim(model_name))
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_basic_init", ier)
       Call error(0, 'kim_basic_init failed in kim_setup')
    End If

    ier = kim_api_get_NBC_method(pkim, ActiveNBC_Method)
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_get_NBC_method", ier)
       Call error(0, 'kim_api_get_NBC_method failed in kim_setup')
    End If
    If      (Index(ActiveNBC_Method,"NEIGH_RVEC_H") == 1) Then
       HalfList  = .true.
       RijNeeded = .true.
    Else If (Index(ActiveNBC_Method,"NEIGH_RVEC_F") == 1) Then
       HalfList  = .false.
       RijNeeded = .true.
    Else If (Index(ActiveNBC_Method,"NEIGH_PURE_H") == 1) Then
       HalfList  = .true.
       RijNeeded = .false.
    Else If (Index(ActiveNBC_Method,"NEIGH_PURE_F") == 1) Then
       HalfList  = .false.
       RijNeeded = .false.
    End If

    fail=0
    Allocate (kim_list(max_list,mxatdm), Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') 'kim_setup kim_list allocation failure'
       Call error(0,message)
    End If

    If (RijNeeded) Allocate (kim_Rij(3,max_list,mxatdm), Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') 'kim_setup kim_Rij allocation failure'
       Call error(0,message)
    End If

    limit=iadd*mxbfxp

    Allocate (rev_comm_buffer(1:limit), Stat=fail)
    If (fail > 0) Then
       Write(message, '(a)') 'kim_setup rev_comm_buffer allocation failure'
       Call error(0,message)
    End If

! Set buffer limit (half for outgoing data - half for incoming)

    iblock=limit/Merge(2,1,comm%mxnode > 1)

    Call kim_api_allocate(pkim, nlast, 1, ier)
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_allocate", ier)
       Call error(0, 'kim_api_allocate failed in kim_setup')
    End If

    ier = kim_api_model_init(pkim)
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_model_init", ier)
       Call error(0, 'kim_api_model_init failed in kim_setup')
    End If

    Call kim_api_getm_data(pkim, ier,                           &
      "numberOfParticles",           pNOP, 1,                   &
      "numberOfSpecies",             pNPT, 1,                   &
      "numberContributingParticles", pNCP, TRUEFALSE(HalfList), &
      "particleSpecies",             pPT,  1)
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_getm_data", ier)
       Call error(0, 'kim_api_getm_data failed in kim_setup')
    End If
    Call c_f_Pointer(pNOP, numberOfParticles)
    Call c_f_Pointer(pNPT, numberOfSpecies)
    If (HalfList) Call c_f_Pointer(pNCP, numberContributingParticles)
    Call c_f_Pointer(pPT,  particleSpecies, [nlast])

    numberOfParticles = nlast
    numberOfSpecies = ntpatm
    If (HalfList) numberContributingParticles = natms

    Do i=1,nlast
       particleSpecies(i) = kim_api_get_species_code(pkim, &
                                                     Trim(sitnam(lsite(i))), ier)
       If (ier < KIM_STATUS_OK) Then
          idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                      "kim_api_get_species", ier)
          Call error(0, 'kim_api_get_species failed in kim_setup')
       End If
    End Do

    ier = kim_api_set_method(pkim, "get_neigh", 1, c_funloc(get_neigh))
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_set_method", ier)
       Call error(0, 'kim_api_set_method failed in kim_setup')
    End If
#else
    Call kim_message()
#endif
  End Subroutine kim_setup

  Subroutine kim_cleanup(comm)

!-------------------------------------------------------------------------------
!
! kim_cleanup
!
! This function releases resources and cleans up the KIM API interface
!
!-------------------------------------------------------------------------------

    Type(comms_type),        Intent( InOut ) :: comm

#ifdef KIM
    Integer                 :: fail
    Integer( Kind = c_int ) :: ier, idum

    Character( Len = 256 ) :: message

    fail=0
    Deallocate (kim_list, Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') &
            'failure deallocating kim_list in kim_module'
       Call error(0,message)
    End If

    fail=0
    If (RijNeeded) Deallocate (kim_Rij, Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') &
            'failure deallocating kim_Rij in kim_module'
       Call error(0,message)
    End If

    Deallocate (rev_comm_buffer, Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') &
            'failure deallocating rev_comm_buffer in kim_module'
       Call error(0,message)
    End If

    ier = kim_api_model_destroy(pkim)
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_model_destroy", ier)
       Call error(0, 'kim_api_model_destroy failed in kim_cleanup')
    End If

    Call kim_api_free(pkim, ier)
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_free", ier)
       Call error(0, 'kim_api_free failed in kim_cleanup')
    End If
    pkim = c_null_ptr
#else
    Call kim_message()
#endif
  End Subroutine kim_cleanup

  Subroutine kim_forces(engkim,virkim,stress,list,comm)

!-------------------------------------------------------------------------------
!
! kim_forces
!
! This function updates the data needed by the KIM Model and then asks for
! the energy, forces, and virial from the KIM Model, distributes force
! contributions appropriately to all the processors and updates virial
! and pressure values.
!
!-------------------------------------------------------------------------------


    Real( Kind = wp ), Intent( InOut ) :: engkim
    Real( Kind = wp ), Intent( InOut ) :: virkim
    Real( Kind = wp ), Intent( InOut ) :: stress(1:9)
    Integer( Kind = wi), Dimension(-3:,1:) :: list
    Type(comms_type),        Intent( InOut ) :: comm

#ifdef KIM
    Integer( Kind = c_int ) :: i,j,k
    Integer( Kind = c_int ) :: ier,idum

    Real( Kind = c_double ), Pointer :: coordinates(:,:); Type( c_ptr ) :: pCOORD

    Real( Kind = c_double ), Pointer :: energy;           Type( c_ptr ) :: pE
    Real( Kind = c_double ), Pointer :: forces(:,:);      Type( c_ptr ) :: pF
    Real( Kind = c_double ), Pointer :: virial(:);        Type( c_ptr ) :: pV

    Call kim_api_getm_data(pkim, ier,                &
                           "coordinates", pCOORD, 1, &
                           "energy",      pE,     1, &
                           "forces",      pF,     1, &
                           "virial",      pV,     1)
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_getm_data", ier)
       Call error(0, 'kim_api_getm_data failed in kim_forces')
    End If
    Call c_f_Pointer(pCOORD, coordinates, [3,nlast])
    Call c_f_Pointer(pE,     energy)
    Call c_f_Pointer(pF,     forces,      [3,nlast])
    Call c_f_Pointer(pV,     virial,      [6])

    Do i=1,nlast
      coordinates(1,i) = xxx(i)
      coordinates(2,i) = yyy(i)
      coordinates(3,i) = zzz(i)
    End Do

    kim_list = 0
    If (RijNeeded) kim_Rij = 0
    Do i=1,natms
       Do j=1,list(0,i)
          k=list(j,i)
          If (.not.HalfList .or. (i < k)) Then
             kim_list(1,i) = kim_list(1,i) + 1
             kim_list(kim_list(1,i)+1,i) = k
             If (RijNeeded) kim_Rij(:,kim_list(1,i),i) = coordinates(:,k) - coordinates(:,i)
          End If

          If (k <= natms) Then
             If (.not.HalfList .or. (k < i)) Then
                kim_list(1,k) = kim_list(1,k) + 1
                kim_list(kim_list(1,k)+1,k) = i
                If (RijNeeded) kim_Rij(:,kim_list(1,k),k) = coordinates(:,i) - coordinates(:,k)
             End If
          End If
       End Do
    End Do

    ier = kim_api_model_compute(pkim)
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_compute", ier)
       Call error(0, 'kim_api_compute failed in kim_forces')
    End If

    engkim = Real(energy, wp)

    Do i=1,natms
       fxx(i) = fxx(i) + Real(forces(1,i), wp)
       fyy(i) = fyy(i) + Real(forces(2,i), wp)
       fzz(i) = fzz(i) + Real(forces(3,i), wp)
    End Do

! Distribute force contributions on this processors halo particles to their
! respective processors for inclusion in the fxx, fyy, and fzz arrays.

    Call kim_reverse_communication(forces,comm)

    stress(1) = stress(1) - Real(virial(1), wp)
    stress(2) = stress(2) - Real(virial(6), wp)
    stress(3) = stress(3) - Real(virial(5), wp)
    stress(4) = stress(4) - Real(virial(6), wp)
    stress(5) = stress(5) - Real(virial(2), wp)
    stress(6) = stress(6) - Real(virial(4), wp)
    stress(7) = stress(7) - Real(virial(5), wp)
    stress(8) = stress(8) - Real(virial(4), wp)
    stress(9) = stress(9) - Real(virial(3), wp)

    virkim = virkim + Real(virial(1)+virial(2)+virial(3), wp)
#else
    engkim = 0.0_wp
    virkim = 0.0_wp
    stress = 0.0_wp

    Call kim_message()
#endif
  End Subroutine kim_forces

#ifdef KIM
  Subroutine kim_reverse_communication(forces,comm)

!-------------------------------------------------------------------------------
!
! kim_reverse_communication
!
! This function sends force contributions on halo particles to the appropriate
! processors for incorporation in the global force arrays
!
!-------------------------------------------------------------------------------

    Real( Kind = c_double ), Intent( In    ) :: forces(:,:)
    Type( comms_type ), Intent( InOut ) :: comm

    Integer :: i,jdnode,kdnode,j,jj,k,imove,jmove

    Do i=1,6

! jdnode - destination (send to), kdnode - source (receive from)

       If      (i == 1) Then ! Direction -x
          jdnode = map(1)
          kdnode = map(2)
       Else If (i == 2) Then ! Direction +x
          jdnode = map(2)
          kdnode = map(1)
       Else If (i == 3) Then ! Direction -y
          jdnode = map(3)
          kdnode = map(4)
       Else If (i == 4) Then ! Direction +y
          jdnode = map(4)
          kdnode = map(3)
       Else If (i == 5) Then ! Direction -z
          jdnode = map(5)
          kdnode = map(6)
       Else If (i == 6) Then ! Direction +z
          jdnode = map(6)
          kdnode = map(5)
       End If

! pack if need be

       imove=0
       Do j=idhalo(1,i),idhalo(2,i)
          rev_comm_buffer(imove+1)=Real(forces(1,j), wp)
          rev_comm_buffer(imove+2)=Real(forces(2,j), wp)
          rev_comm_buffer(imove+3)=Real(forces(3,j), wp)
          rev_comm_buffer(imove+4)=Real(ltg(j), wp)
          imove=imove+iadd
       End Do

! exchange buffers

       If (comm%mxnode > 1) Then
          jmove=idhalo(0,i)*iadd
          If (jmove > 0) Then
            Call girecv(comm,rev_comm_buffer(iblock+1:iblock+jmove),kdnode,Export_tag)
            Call gsend(comm,rev_comm_buffer(1:imove),jdnode,Export_tag)
            Call gwait(comm)
          End  If
       Else
          jmove=imove
       End If

! unpack if need be

       k=Merge(iblock,0,comm%mxnode > 1)
       Do k=1,jmove/iadd
          jj=local_index(Nint(rev_comm_buffer(j+4)),nlast,lsi,lsa)

          fxx(jj)=fxx(jj)+rev_comm_buffer(j+1)
          fyy(jj)=fyy(jj)+rev_comm_buffer(j+2)
          fzz(jj)=fzz(jj)+rev_comm_buffer(j+3)

          j=j+iadd
       End Do
    End Do

  End Subroutine kim_reverse_communication

  Function get_neigh(pkim,mode,request,atom,numnei,pnei1atom,pRij) Bind(c)

!-------------------------------------------------------------------------------
!
! get_neigh neighbor list access function
!
! This function implements Locator and Iterator mode
!
!-------------------------------------------------------------------------------

    Integer( Kind = c_int )                  :: get_neigh

    Type( c_ptr ),    Intent( In    ) :: pkim
    Integer( Kind = c_int ), Intent( In    ) :: mode
    Integer( Kind = c_int ), Intent( In    ) :: request
    Integer( Kind = c_int ), Intent(   Out ) :: atom
    Integer( Kind = c_int ), Intent(   Out ) :: numnei
    Type( c_ptr ),    Intent(   Out ) :: pnei1atom
    Type( c_ptr ),    Intent(   Out ) :: pRij

    Integer( Kind = c_int ), Save    :: iterVal = 0
    Integer( Kind = c_int )          :: N
    Integer( Kind = c_int )          :: atomToReturn
    Integer( Kind = c_int ), Pointer :: numberOfParticles; Type( c_ptr ) :: pnAtoms
    Integer( Kind = c_int )          :: ier, idum

! unpack number of particles

    pnAtoms = kim_api_get_data(pkim, "numberOfParticles", ier)
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_get_data", ier)
       Call error(0, 'kim_api_get_data failed in get_neigh')
    End If
    Call c_f_Pointer(pnAtoms, numberOfParticles)

! use global access to neighbor list object

    N = natms

! check mode and request

    If (mode == 0) Then ! iterator mode
       If (request == 0) Then ! reset iterator
          iterVal = 0
          get_neigh = KIM_STATUS_NEIGH_ITER_INIT_OK
          Return
       Elseif (request == 1) Then ! increment iterator
          iterVal = iterVal + 1
          If (iterVal > N) Then
             get_neigh = KIM_STATUS_NEIGH_ITER_PAST_END
             Return
          Else
             atomToReturn = iterVal
          End If
       Else
          idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                      "Invalid request in get_neigh", &
                                      KIM_STATUS_NEIGH_INVALID_REQUEST)
          get_neigh = KIM_STATUS_NEIGH_INVALID_REQUEST
          Return
       End If
    Else If (mode == 1) Then ! locator mode
       If (request > N .or. request < 1) Then
          idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                      "Invalid atom ID in get_neigh", &
                                      KIM_STATUS_PARTICLE_INVALID_ID)
          get_neigh = KIM_STATUS_PARTICLE_INVALID_ID
          Return
       Else
          atomToReturn = request
       End If
    Else ! not iterator or locator mode
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "Invalid mode in get_neigh", &
                                   KIM_STATUS_NEIGH_INVALID_MODE)
       get_neigh = KIM_STATUS_NEIGH_INVALID_MODE
       Return
    End If

  ! set the Returned atom

    atom = atomToReturn

    If (atom <= natms) Then
! set the Returned number of neighbors for the returned atom
       numnei = kim_list(1,atom)
! set the location for the Returned neighbor list
       pnei1atom = c_loc(kim_list(2,atom))
    Else
! set the Returned number of neighbors for the returned atom
       numnei = 0
! set the location for the Returned neighbor list
       pnei1atom = c_null_ptr
    End If

! set Pointer to Rij to appropriate value

    If (RijNeeded .and. atom <= natms) Then
       pRij = c_loc(kim_Rij(1,1,atom))
    Else
       pRij = c_null_ptr
    End If

    get_neigh = KIM_STATUS_OK

  End Function get_neigh

  Subroutine Write_KIM_descriptor(NBC_method, num_types, model_types, &
                                  kim_descriptor, ier)

!-------------------------------------------------------------------------------
!
!  Write KIM descriptor file for MiniMol for given NBC and set of
!  support species types
!
!-------------------------------------------------------------------------------

    Character( Len = KIM_KEY_STRING_LENGTH ), Intent( In    ) :: NBC_method(:)
    Integer( Kind = c_int ),                  Intent( In    ) :: num_types
    Character( Len = * ),                     Intent( In    ) :: model_types(num_types)
    Character( Len = 10000 ),                 Intent(   Out ) :: kim_descriptor
    Integer( Kind = c_int ),                  Intent(   Out ) :: ier

    Integer( Kind = c_int)  :: i
    Character( Len = 103  ) :: divider
    Character( Len = 1    ) :: cr
    Character( Len = 1024 ) :: type_line
    Character( Len = 32   ) :: nbcline

! Initialize error flag

    ier = KIM_STATUS_OK

! Define frequently used variables

    cr = Char(10)
    divider = '#######################################################################################################'

! Write Minimol descriptor file into string kim_descriptor

    kim_descriptor = &
     divider                                                                      // cr // &
     '#'                                                                          // cr // &
     '# CDDL HEADER START'                                                        // cr // &
     '#'                                                                          // cr // &
     '# The contents of this file are subject to the terms of the Common Development' // cr // &
     '# and Distribution License Version 1.0 (the "License").'                    // cr // &
     '#'                                                                          // cr // &
     '# You can obtain a copy of the license at'                                  // cr // &
     '# http://www.opensource.org/licenses/CDDL-1.0.  See the License for the'    // cr // &
     '# specific language governing permissions and limitations under the License.' // cr // &
     '#'                                                                          // cr // &
     '# When distributing Covered Code, include this CDDL HEADER in each file and'// cr // &
     '# include the License file in a prominent location with the name LICENSE.CDDL.' // cr // &
     '# If applicable, add the following below this CDDL HEADER, with the fields' // cr // &
     '# enclosed by brackets "[]" replaced with your own identifying information:'// cr // &
     '#'                                                                          // cr // &
     '# Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.' // cr // &
     '#'                                                                          // cr // &
     '# CDDL HEADER END'                                                          // cr // &
     '#'                                                                          // cr // &
                                                                                     cr // &
     '#'                                                                          // cr // &
     '# Copyright (c) 2013--2014, Regents of the University of Minnesota.'        // cr // &
     '# All rights reserved.'                                                     // cr // &
     '#'                                                                          // cr // &
     '# Contributors:'                                                            // cr // &
     '#    Automatically generated by calling Test'                               // cr // &
     '#'                                                                          // cr // &
                                                                                     cr // &
     '#'                                                                          // cr // &
     '# See std/standard.kim for documentation about this file'                   // cr // &
     '#'                                                                          // cr // &
     divider                                                                      // cr // &
                                                                                     cr // &
                                                                                     cr // &
     'KIM_API_Version := 1.6.0'                                                   // cr // &
                                                                                     cr // &
     'Unit_length      := A'                                                      // cr // &
     'Unit_energy      := amu*A^2/(ps)^2'                                         // cr // &
     'Unit_charge      := e'                                                      // cr // &
     'Unit_temperature := K'                                                      // cr // &
     'Unit_time        := ps'                                                     // cr // &
                                                                                     cr // &
                                                                                     cr // &
     divider                                                                      // cr // &
     'PARTICLE_SPECIES:'                                                          // cr // &
     '# Symbol/name               Type                    code'                   // cr

    Do i = 1,num_types
       Write(type_line,'(a,'' '',''spec'',20x,i4)') Trim(model_types(i)),0
       kim_descriptor = Trim(kim_descriptor) // Trim(type_line) // cr
    End Do

    kim_descriptor = Trim(kim_descriptor) // &
                                                                                     cr // &
                                                                                     cr // &
     divider                                                                      // cr // &
     'CONVENTIONS:'                                                               // cr // &
     '# Name                      Type'                                           // cr // &
                                                                                     cr // &
     'OneBasedLists               flag'                                           // cr // &
                                                                                     cr // &
     'Neigh_BothAccess            flag'                                           // cr // &
                                                                                     cr
    Do i=1,Size(NBC_Method,1)
       Write(nbcline,'(a,'' '',''flag'')') Trim(NBC_Method(i))
       kim_descriptor = Trim(kim_descriptor) // nbcline // cr // cr
    End Do

    kim_descriptor = Trim(kim_descriptor) // &
     divider                                                                      // cr // &
     'MODEL_INPUT:'                                                               // cr // &
     '# Name                      Type         Unit       Shape              requirements' // cr // &
     'numberOfParticles           integer      None       []'                     // cr // &
                                                                                     cr // &
     'numberOfSpecies             integer      None       []'                     // cr // &
                                                                                     cr // &
     'particleSpecies             integer      None       [numberOfParticles]'    // cr // &
                                                                                     cr // &
     'coordinates                 double       length     [numberOfParticles,3]'  // cr // &
                                                                                     cr // &
     'get_neigh                   method       None       []'                     // cr // &
                                                                                     cr // &
     'neighObject                 pointer      None       []'                     // cr // &
                                                                                     cr // &
     'numberContributingParticles integer      None       []'                     // cr // &
                                                                                     cr // &
     divider                                                                      // cr // &
     'MODEL_OUTPUT:'                                                              // cr // &
     '# Name                      Type         Unit       Shape              requirements' // cr // &
                                                                                     cr // &
     'destroy                     method       None       []'                     // cr // &
                                                                                     cr // &
     'compute                     method       None       []'                     // cr // &
                                                                                     cr // &
     'cutoff                      double       length     []'                     // cr // &
                                                                                     cr // &
     'energy                      double       energy     []'                     // cr // &
                                                                                     cr // &
     'forces                      double       force      [numberOfParticles,3]'  // cr // &
                                                                                     cr // &
     'virial                      double       energy     [6]'                    // cr // &
                                                                                     cr // &
     divider                                                                            // cr

  End Subroutine Write_KIM_descriptor

  Function kim_basic_init(num_types, model_types, model_name)

!-------------------------------------------------------------------------------
!
! kim_basic_init
!
! This private function is for basic kim setup and is only used by kim_setup
!
!-------------------------------------------------------------------------------

    Integer( Kind = c_int )                  :: kim_basic_init

    Integer( Kind = c_int ), Intent( In    ) :: num_types
    Character( Len = * ),    Intent( In    ) :: model_types(num_types)
    Character( Len = * ),    Intent( In    ) :: model_name

    Integer :: ier, idum
    Character( Len = 10000 ) :: kim_descriptor
    Character( Len = KIM_KEY_STRING_LENGTH ) :: NBC_Method(4)

    Real( Kind = c_double ), Pointer :: model_cutoff; Type( c_ptr ) :: pcut

    NBC_Method = ["NEIGH_RVEC_H", &
                  "NEIGH_PURE_H", &
                  "NEIGH_RVEC_F", &
                  "NEIGH_PURE_F"]

    ier = KIM_STATUS_OK

! Here we assume that all particle types interact via the KIM Model.

    Call Write_KIM_descriptor(NBC_Method, num_types, model_types, &
                              kim_descriptor, ier)
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "Write_KIM_descriptor", ier)
       Call error(0, 'Write_KIM_descriptor failed in kim_basic_init')
    End If

    ier = kim_api_string_init(pkim, Trim(kim_descriptor), Trim(model_name))
    If (ier < KIM_STATUS_OK) Then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_string_init", ier)
       Call error(0, 'kim_api_string_init failed in kim_basic_init')
    End If

    kim_basic_init = ier
  End Function kim_basic_init
#endif

  Subroutine kim_message()
   
#ifndef KIM
    Character( Len = 256 ) :: message
    Write(message,'(1x,a)') "*** warning - kim directive found in FIELD but openKIM not available !!! ***"
    Call error(0,message,.true.)
#endif
  End Subroutine kim_message

End Module kim
