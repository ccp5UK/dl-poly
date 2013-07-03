Subroutine metal_table_derivatives(ityp,buffer,v2d,vvv)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to calculate numerical derivatives of tabulated
! EAM metal potentials
!
! copyright - daresbury laboratory
! author    - w.smith march 2006
! amended   - i.t.todorov june 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : zero_plus,mxgrid

  Implicit None

  Integer,           Intent( In    ) :: ityp,v2d
  Real( Kind = wp ), Intent( In    ) :: buffer(1:mxgrid)
  Real( Kind = wp ), Intent( InOut ) :: vvv(1:mxgrid,1:v2d,1:2)

  Integer           :: i,v_end,i_start,i_end
  Real( Kind = wp ) :: delmet,aa0,aa1,aa2,aa3,aa4, &
                       d1y,d2y,d3y,d4y,f0,f1,f2,f3,f4

! interpolation parameters

  vvv(1,ityp,2)=buffer(1)
  vvv(2,ityp,2)=buffer(2)
  vvv(3,ityp,2)=buffer(3)
  vvv(4,ityp,2)=buffer(4)

! construct interpolation table

  delmet=buffer(4)
  v_end=Nint(buffer(1))
  i_start=5    +2
  i_end  =v_end-2
  Do i=i_start,i_end
     aa0=buffer(i)
     If (Abs(aa0) <= zero_plus) Then
        f0=0.0_wp
        f1=0.0_wp
        f2=0.0_wp
        f3=0.0_wp
        f4=0.0_wp
     Else
        f0=buffer(i-2)/aa0
        f1=buffer(i-1)/aa0
        f2=1.0_wp
        f3=buffer(i+1)/aa0
        f4=buffer(i+2)/aa0
     End If

! calculate numerical differences for 5-point interpolation

     d1y=(f1-f0)
     d2y=(f2-f1)-(f1-f0)
     d3y=(f3-f0)+3.0_wp*(f1-f2)
     d4y=(f4-f3)+3.0_wp*(f2-f3)+3.0_wp*(f2-f1)+(f0-f1)

! calculate polynomial coefficients

     aa0=aa0/delmet
     aa4=d4y/24.0_wp
     aa3=(d3y+12.0_wp*aa4)/6.0_wp
     aa2=(d2y+6.0_wp*aa3-14.0_wp*aa4)/2.0_wp
     aa1=d1y+3.0_wp*aa2-7.0_wp*aa3+15.0_wp*aa4

! calculate derivatives

     vvv(i,ityp,2)=aa1*aa0

! derivatives at extremes of range

     If      (i == i_start) Then
        vvv(i_start-2,ityp,2)=(aa1-4.0_wp*aa2+12.0_wp*aa3-32.0_wp*aa4)*aa0
        vvv(i_start-1,ityp,2)=(aa1-2.0_wp*aa2+3.0_wp*aa3-4.0_wp*aa4)*aa0
     Else If (i == i_end  ) Then
        vvv(i_end  +1,ityp,2)=(aa1+2.0_wp*aa2+3.0_wp*aa3+4.0_wp*aa4)*aa0
        vvv(i_end  +2,ityp,2)=(aa1+4.0_wp*aa2+12.0_wp*aa3+32.0_wp*aa4)*aa0
     End If
  End Do

! set derivatives to zero beyond end point of function

  Do i=v_end+3,mxgrid
     vvv(i,ityp,2)=0.0_wp
  End Do

End Subroutine metal_table_derivatives
