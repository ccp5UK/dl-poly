Subroutine statistics_connect_set(rcut)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of data between neighbouring
! domains/nodes in order to reconnect some statistical information
! between replayed frames of history
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module, Only : mxnode
  Use setup_module, only : mxatdm,mxstak,zero_plus
  Use domains_module
  Use config_module
  Use statistics_module
  Use msd_module,   Only : l_msd

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: rcut

  Real( Kind = wp ), Save :: cut

  Integer           :: nlx,nly,nlz,i,i0,kk
  Real( Kind = wp ) :: det,celprp(1:10),rcell(1:9),x,y,z, &
                       xdc,ydc,zdc,cwx,cwy,cwz,ecwx,ecwy,ecwz

  If (mxnode > 1) Then

! Define cut

     cut=rcut+1.0e-6_wp

     Call dcell(cell,celprp)
     Call invert(cell,rcell,det)

! calculate link cell dimensions per node

     nlx=Int(celprp(7)/(cut*nprx_r))
     nly=Int(celprp(8)/(cut*npry_r))
     nlz=Int(celprp(9)/(cut*nprz_r))

! Get the total number of link-cells in MD cell per direction

     xdc=Real(nlx*nprx,wp)
     ydc=Real(nly*npry,wp)
     zdc=Real(nlz*nprz,wp)

! link-cell widths in reduced space

     cwx=1.0_wp/xdc
     cwy=1.0_wp/ydc
     cwz=1.0_wp/zdc

! Distance from the - edge of this domain

     ecwx=Nearest( (-0.5_wp+cwx)+Real(idx,wp)*r_nprx , +1.0_wp)+zero_plus
     ecwy=Nearest( (-0.5_wp+cwy)+Real(idy,wp)*r_npry , +1.0_wp)+zero_plus
     ecwz=Nearest( (-0.5_wp+cwz)+Real(idz,wp)*r_nprz , +1.0_wp)+zero_plus

! Distance from the + edge of this domain with a possible
! extension strip for the one linked cell per domain scenario

     cwx=Nearest( (-0.5_wp-cwx)+Real(idx+1,wp)*r_nprx , -1.0_wp)-zero_plus-Merge( cwx*1.0e-10_wp , 0.0_wp , nlx == 1 )
     cwy=Nearest( (-0.5_wp-cwy)+Real(idy+1,wp)*r_npry , -1.0_wp)-zero_plus-Merge( cwy*1.0e-10_wp , 0.0_wp , nly == 1 )
     cwz=Nearest( (-0.5_wp-cwz)+Real(idz+1,wp)*r_nprz , -1.0_wp)-zero_plus-Merge( cwz*1.0e-10_wp , 0.0_wp , nlz == 1 )

     ixyz(1:mxatdm)=0 ! Initialise move (former halo) indicator
     Do i=1,natms
        x=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
        y=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
        z=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)

        If (x <= ecwx) ixyz(i)=ixyz(i)+1
        If (x >=  cwx) ixyz(i)=ixyz(i)+2

        If (y <= ecwy) ixyz(i)=ixyz(i)+10
        If (y >=  cwy) ixyz(i)=ixyz(i)+20

        If (z <= ecwz) ixyz(i)=ixyz(i)+100
        If (z >=  cwz) ixyz(i)=ixyz(i)+200
     End Do

     lsi=0;lsa=0 ! This is a must, unfortunately
     Do i=1,natms
        lsi(i)=i
        lsa(i)=ltg(i)
     End Do
     Call shellsort2(natms,lsi,lsa)

     natms0 = natms
     ltg0(1:natms0) = ltg(1:natms0) !; ltg0(natms0+1: ) = 0
     lsa0(1:natms0) = lsa(1:natms0) !; lsa0(natms0+1: ) = 0
     lsi0(1:natms0) = lsi(1:natms0) !; lsi0(natms0+1: ) = 0

     xin0(1:natms0) = xin(1:natms0) !; xin0(natms0+1: ) = 0 ; xin = 0.0_wp
     yin0(1:natms0) = yin(1:natms0) !; yin0(natms0+1: ) = 0 ; yin = 0.0_wp
     zin0(1:natms0) = zin(1:natms0) !; zin0(natms0+1: ) = 0 ; zin = 0.0_wp

     xto0(1:natms0) = xto(1:natms0) !; xto0(natms0+1: ) = 0
     yto0(1:natms0) = yto(1:natms0) !; yto0(natms0+1: ) = 0
     zto0(1:natms0) = zto(1:natms0) !; zto0(natms0+1: ) = 0

     If (l_msd) Then
        i0=2*natms0
        stpvl00(1:i0)=stpvl0(28:27+i0) !; stpvl00(i0+1: )=0.0_wp
        stpval0(1:i0)=stpval(28:27+i0) !; stpval0(i0+1: )=0.0_wp
        zumval0(1:i0)=zumval(28:27+i0) !; zumval0(i0+1: )=0.0_wp
        ravval0(1:i0)=ravval(28:27+i0) !; ravval0(i0+1: )=0.0_wp
        ssqval0(1:i0)=ssqval(28:27+i0) !; ssqval0(i0+1: )=0.0_wp
        sumval0(1:i0)=sumval(28:27+i0) !; sumval0(i0+1: )=0.0_wp
        Do kk=1,mxstak
           stkval0(kk,1:i0)=stkval(kk,1:i0) !; stkval0(kk,i0+1: )=0.0_wp
        End Do
     End If
  End If

End Subroutine statistics_connect_set
