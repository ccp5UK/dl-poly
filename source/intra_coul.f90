Subroutine intra_coul(keyfce,rcut,alpha,epsq,chgprd,rrr,rsq,coul,fcoul,safe)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating bond's or 1-4 dihedral
! electrostatics: adjusted by a weighting factor
!
! copyright - daresbury laboratory
! amended   - i.t.todorov march 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : sqrpi,zero_plus

  Implicit None

  Integer,           Intent( In    ) :: keyfce
  Real( Kind = wp ), Intent( In    ) :: chgprd,rcut,alpha,epsq,rrr,rsq
  Real( Kind = wp ), Intent(   Out ) :: coul,fcoul
  Logical,           Intent( InOut ) :: safe

  Logical,           Save :: newjob = .true. , damp
  Real( Kind = wp ), Save :: aa,bb, rfld0,rfld1

  Real( Kind = wp ) :: exp1,tt,erc,fer,b0

  Real( Kind = wp ), Parameter :: aa1 =  0.254829592_wp
  Real( Kind = wp ), Parameter :: aa2 = -0.284496736_wp
  Real( Kind = wp ), Parameter :: aa3 =  1.421413741_wp
  Real( Kind = wp ), Parameter :: aa4 = -1.453152027_wp
  Real( Kind = wp ), Parameter :: aa5 =  1.061405429_wp
  Real( Kind = wp ), Parameter :: pp  =  0.3275911_wp

  If (newjob) Then
     newjob = .false.

! Check for damped force-shifted coulombic and reaction field interactions
! and set force and potential shifting parameters dependingly

     damp=.false.
     If (alpha > zero_plus) Then
        damp=.true.

        exp1= Exp(-(alpha*rcut)**2)
        tt  = 1.0_wp/(1.0_wp+pp*alpha*rcut)

        erc = tt*(aa1+tt*(aa2+tt*(aa3+tt*(aa4+tt*aa5))))*exp1/rcut
        fer = (erc + 2.0_wp*(alpha/sqrpi)*exp1)/rcut**2

        aa  = fer*rcut
        bb  = -(erc + aa*rcut)
     Else If (keyfce == 8) Then
        aa =  1.0_wp/rcut**2
        bb = -2.0_wp/rcut ! = -(1.0_wp/rcut+aa*rcut)
     End If

! set reaction field terms for RFC

     If (keyfce == 10) Then
        b0    = 2.0_wp*(epsq - 1.0_wp)/(2.0_wp*epsq + 1.0_wp)
        rfld0 = b0/rcut**3
        rfld1 = (1.0_wp + 0.5_wp*b0)/rcut
     End If
  End If

! initialise defaults for coulombic energy and force contributions

   coul =0.0_wp
   fcoul=0.0_wp

! Electrostatics by ewald sum = direct coulombic

  If      (keyfce ==  2 .or. keyfce ==  6) Then

     coul = chgprd/rrr
     fcoul= coul/rsq

! distance dependent dielectric

  Else If (keyfce ==  4) Then

     coul = chgprd/rsq
     fcoul= 2.0_wp*coul/rsq

! force shifted coulombic and reaction field

  Else If (keyfce ==  8 .or. keyfce == 10) Then

     If (damp) Then ! calculate damping contributions
        exp1= Exp(-(alpha*rrr)**2)
        tt  = 1.0_wp/(1.0_wp+pp*alpha*rrr)

        erc = tt*(aa1+tt*(aa2+tt*(aa3+tt*(aa4+tt*aa5))))*exp1/rrr
        fer = (erc + 2.0_wp*(alpha/sqrpi)*exp1)/rsq

        coul = chgprd*(erc + aa*rrr + bb)
        fcoul= chgprd*(fer - aa/rrr)
     End If

     If      (keyfce ==  8) Then ! force shifted coulombic
        If (.not.damp) Then ! pure
           coul = chgprd*(1.0_wp/rrr + aa*rrr+ bb)
           fcoul= chgprd*(1.0_wp/rsq - aa)/rrr
        Else                ! damped
           coul = coul
           fcoul= fcoul
        End If
     Else If (keyfce == 10) Then ! reaction field
        If (.not.damp) Then ! pure
           coul = chgprd*(1.0_wp/rrr + 0.5_wp*rfld0*rsq - rfld1)
           fcoul= chgprd*(1.0_wp/rsq/rrr - rfld0)
        Else                ! damped
           coul = coul  + chgprd*(0.5_wp*rfld0*rsq - rfld1)
           fcoul= fcoul + chgprd*(- rfld0)
        End If
     End If

  Else

     safe = .false.

  End If

End Subroutine intra_coul
