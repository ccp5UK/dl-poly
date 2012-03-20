Subroutine dihedrals_14_check &
           (l_str,ntpmls,nummols,numang,keyang,lstang,numdih,lstdih,prmdih)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to check and resolve conflicting molecular
! forcefield specification for 1-4 interactions
!
! copyright - daresbury laboratory
! author    - w.smith march 1999
! amended   - i.t.todorov march 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : mxtmls,mxtang,mxtdih,mxpdih

  Logical,           Intent( In    ) :: l_str
  Integer,           Intent( In    ) :: ntpmls,nummols(1:mxtmls),numang(1:mxtmls), &
                                        keyang(1:mxtang),lstang(1:3,1:mxtang),     &
                                        numdih(1:mxtmls),lstdih(1:4,1:mxtdih)
  Real( Kind = wp ), Intent( InOut ) :: prmdih(1:mxpdih,1:mxtdih)

  Integer :: kangle,kdihed,itmols,imols,langle,ldihed,mdihed, &
             iang,jang,idih,jdih,kdih,ldih

! Initialise angle and dihedral interaction counters

  kangle=0
  kdihed=0

! loop over molecular types

  Do itmols=1,ntpmls

! loop over molecules in system

     Do imols=1,nummols(itmols)

! check for valence angle on dihedral angle conflicts

        Do langle=1,numang(itmols)

           If (keyang(langle+kangle) > 0) Then

              iang=lstang(1,langle+kangle)
              jang=lstang(3,langle+kangle)

              Do ldihed=1,numdih(itmols)

                 idih=lstdih(1,ldihed+kdihed)
                 jdih=lstdih(4,ldihed+kdihed)

                 If (Min(iang,jang) == Min(idih,jdih) .and. Max(iang,jang) == Max(idih,jdih)) Then

                    prmdih(4,ldihed+kdihed)=0.0_wp
                    prmdih(5,ldihed+kdihed)=0.0_wp
                    If (l_str) Call warning(20,Real(itmols,wp),Real(idih,wp),Real(jdih,wp))

                 End If

              End Do

           End If

        End Do

! check for double dihedral angle conflicts

        Do ldihed=1,numdih(itmols)-1

           idih=lstdih(1,ldihed+kdihed)
           jdih=lstdih(4,ldihed+kdihed)

           Do mdihed=ldihed+1,numdih(itmols)

              kdih=lstdih(1,mdihed+kdihed)
              ldih=lstdih(4,mdihed+kdihed)

              If (Min(kdih,ldih) == Min(idih,jdih) .and. Max(kdih,ldih) == Max(idih,jdih)) Then

                 prmdih(4,mdihed+kdihed)=0.0_wp
                 prmdih(5,mdihed+kdihed)=0.0_wp
                 If (l_str) Call warning(20,Real(itmols,wp),Real(kdih,wp),Real(ldih,wp))

              End If

           End Do

        End Do

    End Do

! Update counters

    kangle=kangle+numang(itmols)
    kdihed=kdihed+numdih(itmols)

  End Do

End Subroutine dihedrals_14_check
