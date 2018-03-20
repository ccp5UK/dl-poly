Subroutine report_topology               &
           (megatm,megfrz,atmfre,atmfrz, &
           megshl,megcon,megpmf,megrgd,  &
           megtet,megbnd,megang,megdih,meginv)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reporting FIELD topology
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SETUP MODULES

  Use comms_module, Only : idnode
  Use setup_module

! SITE MODULE

  Use site_module

! INTERACTION MODULES

  Use core_shell

  Use constraints
  Use pmf_module

  Use rigid_bodies_module

  Use tethers_module

  Use bonds_module
  Use angles_module
  Use dihedrals
  Use inversions_module

  Implicit None

  Integer, Intent( In    ) :: megatm,megfrz,atmfre,atmfrz, &
                              megshl,megcon,megpmf,megrgd, &
                              megtet,megbnd,megang,megdih,meginv

  Integer :: itmols,nsite,                &
             isite1,isite2,isite3,isite4, &
             iatm1,iatm2,iatm3,iatm4,     &
             ishls,nshels,frzshl,mgfrsh,  &
             mgcon,                       &
             lpmf,ipmf,jpmf,              &
             mgrgd,                       &
             iteth,nteth,frztet,mgfrtt,   &
             ibond,nbonds,frzbnd,mgfrbn,  &
             iang,nangle,frzang,mgfran,   &
             idih,ndihed,frzdih,mgfrdh,   &
             iinv,ninver,frzinv,mgfrin

  Character( Len = 4) :: frzpmf

  nshels = 0
  mgfrsh = 0

  mgcon  = 0

  frzpmf='NONE'

  mgrgd  = 0

  nteth  = 0
  mgfrtt = 0

  nbonds = 0
  mgfrbn = 0

  nangle = 0
  mgfran = 0

  ndihed = 0
  mgfrdh = 0

  ninver = 0
  mgfrin = 0

  nsite  = 0
  Do itmols=1,ntpmls
     frzshl=0
     Do ishls=1,numshl(itmols)
        nshels=nshels+1

        iatm1=lstshl(1,nshels)

        isite1 = nsite + iatm1

        If (frzsit(isite1) == 1) frzshl=frzshl+1
     End Do

     Do lpmf=1,numpmf(itmols) ! numpmf can only be 1 or 0, so the 'Do' loop is used as an 'If' condition
        Do ipmf=1,2
           Do jpmf=1,mxtpmf(ipmf)
              iatm1=lstpmf(jpmf,ipmf)

              isite1 = nsite + iatm1

              If (frzsit(isite1) == 1) frzpmf='ALL'
           End Do
        End Do
     End Do

     frztet=0
     Do iteth=1,numteth(itmols)
        nteth=nteth+1

        iatm1=lsttet(nteth)

        isite1 = nsite + iatm1

        If (frzsit(isite1) == 1) frztet=frztet+1
     End Do

     frzbnd=0
     Do ibond=1,numbonds(itmols)
        nbonds=nbonds+1

        iatm1=lstbnd(1,nbonds)
        iatm2=lstbnd(2,nbonds)

        isite1 = nsite + iatm1
        isite2 = nsite + iatm2

        If (frzsit(isite1)+frzsit(isite2) == 2) frzbnd=frzbnd+1
     End Do

     frzang=0
     Do iang=1,numang(itmols)
        nangle=nangle+1

        iatm1=lstang(1,nangle)
        iatm2=lstang(2,nangle)
        iatm3=lstang(3,nangle)

        isite1 = nsite + iatm1
        isite2 = nsite + iatm2
        isite3 = nsite + iatm3

        If (frzsit(isite1)+frzsit(isite2)+frzsit(isite3) == 3) frzang=frzang+1
     End Do

     frzdih=0
     Do idih=1,numdih(itmols)
        ndihed=ndihed+1

        iatm1=lstdih(1,ndihed)
        iatm2=lstdih(2,ndihed)
        iatm3=lstdih(3,ndihed)
        iatm4=lstdih(4,ndihed)

        isite1 = nsite + iatm1
        isite2 = nsite + iatm2
        isite3 = nsite + iatm3
        isite4 = nsite + iatm4

        If (frzsit(isite1)+frzsit(isite2)+frzsit(isite3)+frzsit(isite4) == 4) frzdih=frzdih+1
     End Do

     frzinv=0
     Do iinv=1,numinv(itmols)
        ninver=ninver+1

        iatm1=lstinv(1,ninver)
        iatm2=lstinv(2,ninver)
        iatm3=lstinv(3,ninver)
        iatm4=lstinv(4,ninver)

        isite1 = nsite + iatm1
        isite2 = nsite + iatm2
        isite3 = nsite + iatm3
        isite4 = nsite + iatm4

        If (frzsit(isite1)+frzsit(isite2)+frzsit(isite3)+frzsit(isite4) == 4) frzinv=frzinv+1
     End Do

     mgcon=mgcon+nummols(itmols)*numcon(itmols)
     mgrgd=mgrgd+nummols(itmols)*numrgd(itmols)

     mgfrsh=mgfrsh+nummols(itmols)*frzshl
     mgfrtt=mgfrtt+nummols(itmols)*frztet
     mgfrbn=mgfrbn+nummols(itmols)*frzbnd
     mgfran=mgfran+nummols(itmols)*frzang
     mgfrdh=mgfrdh+nummols(itmols)*frzdih
     mgfrin=mgfrin+nummols(itmols)*frzinv

     nsite=nsite+numsit(itmols)
  End Do

  If (idnode == 0) Write(nrite,'(/,7(1x,a66,/),4(1x,a29,i11,a8,i11,a7,/),             &
                   & (1x,a29,i11,a15,a4,a7,/),6(1x,a29,i11,a8,i11,a7,/),(1x,a66,/))') &
     "//==============================================================\\",            &
     "||                                                              ||",            &
     "||           SUMMARY  OF  TOPOLOGICAL  DECOMPOSITION            ||",            &
     "||                                                              ||",            &
     "||--------------------------------------------------------------||",            &
     "||  INTERACTION  OR  TYPE  |  GRAND TOTAL | Fully/Partly FROZEN ||",            &
     "||-------------------------+--------------+---------------------||",            &
     "||  all particles/sites    | ", megatm,"  |  F  ",megfrz, "     ||",            &
     "||  free particles         | ", atmfre,"  |  F  ",atmfrz, "     ||",            &
     "||  core-shell units       | ", megshl,"  |  P  ",mgfrsh, "     ||",            &
     "||  constraint bond units  | ", mgcon, "  |  F  ",mgcon-megcon,"     ||",       &
     "||  PMF units              | ", megpmf,"  |  P         ",frzpmf,"     ||",      &
     "||  rigid body units       | ", mgrgd, "  |  F  ",mgrgd-megrgd,"     ||",       &
     "||  tethered atom units    | ", megtet,"  |  F  ",mgfrtt, "     ||",            &
     "||  chemical bond units    | ", megbnd,"  |  F  ",mgfrbn, "     ||",            &
     "||  bond angle units       | ", megang,"  |  F  ",mgfran, "     ||",            &
     "||  dihedral angle units   | ", megdih,"  |  F  ",mgfrdh, "     ||",            &
     "||  inversion angle units  | ", meginv,"  |  F  ",mgfrin, "     ||",            &
     "\\==============================================================//"

End Subroutine report_topology
