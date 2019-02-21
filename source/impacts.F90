Module impacts
  Use kinds, Only : wp
  Use constants,        Only : eu_ev
  Use comms,               Only : comms_type,gcheck
  Use configuration,       Only : configuration_type
  Use rigid_bodies, Only : rigid_bodies_type
  Use core_shell,   Only : core_shell_type
  Use kinetics,      Only : getvom

  Use numerics, Only : local_index
  Use errors_warnings, Only : error
  Implicit None

  Private
  Type, Public :: impact_type
    Integer          :: imd, tmd
    Real( Kind = wp) :: emd, vmx, vmy, vmz
  End Type 

  Public :: impact

Contains

  Subroutine impact(rigid,cshell,impa,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for setting impact on a particle
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov november 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( impact_type ), Intent( In    ) :: impa
    Type( core_shell_type ), Intent( In    ) :: cshell
    Type( configuration_type ), Intent( InOut ) :: config
    Type( comms_type ) , Intent( InOut ) :: comm
    Logical           :: safe = .true.

    Integer           :: i,j,irgd,jrgd,lrgd,rgdtyp
    Real( Kind = wp ) :: tmp,vom(1:3)


    ! Apply impact to the selected non-frozen, non-massless and non-shell particle

    i=local_index(impa%imd,config%nlast,config%lsi,config%lsa)
    If (i > 0 .and. i <= config%natms) Then
      If (config%lfrzn(i) == 0 .and. config%lfree(i) == 0 .and. All(cshell%listshl(2,1:cshell%ntshl) /= impa%imd)) Then
        tmp=Sqrt(2000.0_wp*impa%emd*eu_ev/config%weight(i)/(impa%vmx**2+impa%vmy**2+impa%vmz**2)) !impa%emd is in keV=1000*eu_ev
        config%vxx(i)=tmp*impa%vmx
        config%vyy(i)=tmp*impa%vmy
        config%vzz(i)=tmp*impa%vmz
      Else
        safe=.false.
      End If
    End If

    Call gcheck(comm,safe)
    If (.not.safe) Call error(610)

    Call config%chvom(.true.) ! Enable COM momentum removal

    ! remove centre of mass motion

    If (rigid%total > 0) Then
      Call getvom(vom,rigid,config,comm)

      Do j=1,config%nfree
        i=config%lstfre(j)

        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
          config%vxx(i) = config%vxx(i) - vom(1)
          config%vyy(i) = config%vyy(i) - vom(2)
          config%vzz(i) = config%vzz(i) - vom(3)
        End If
      End Do

      Do irgd=1,rigid%n_types
        rgdtyp=rigid%list(0,irgd)

        If (rigid%frozen(0,rgdtyp) == 0) Then
          rigid%vxx(irgd) = rigid%vxx(irgd) - vom(1)
          rigid%vyy(irgd) = rigid%vyy(irgd) - vom(2)
          rigid%vzz(irgd) = rigid%vzz(irgd) - vom(3)

          lrgd=rigid%list(-1,irgd)
          Do jrgd=1,lrgd
            i=rigid%index_local(jrgd,irgd) ! local index of particle/site

            If (i <= config%natms) Then
              config%vxx(i) = config%vxx(i) - vom(1)
              config%vyy(i) = config%vyy(i) - vom(2)
              config%vzz(i) = config%vzz(i) - vom(3)
            End If
          End Do
        End If
      End Do
    Else
      Call getvom(vom,config,comm)

      Do i=1,config%natms
        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
          config%vxx(i) = config%vxx(i) - vom(1)
          config%vyy(i) = config%vyy(i) - vom(2)
          config%vzz(i) = config%vzz(i) - vom(3)
        End If
      End Do
    End If

    Call config%chvom() ! default to specification

  End Subroutine impact
End Module impacts
