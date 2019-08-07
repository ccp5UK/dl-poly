!> Module containing magic numbers
!> the usage of these numbers comes from nothing else
!> than the whimsical arbitrariness of the author  
!>
!> Copyright - Daresbury Laboratory
!
!> Author - A.M. Elena August 2019

Module magic_numbers
  Use kinds, Only : wp
  Implicit none
  Private

  Real(Kind=wp), Public :: m1 = 1.0e-10_wp ! some small number
  Real(Kind=wp), Public :: m2 = 1.0e-6_wp ! some small number but not as small as m1
  Real(Kind=wp), Public :: m3 = 1.7_wp ! 1.7
  Real(Kind=wp), Public :: m4 = 0.45_wp ! 
  Integer, Public :: m5 = 1005
  Real(Kind=wp), Public :: m6 = 0.05_wp  
  Real(Kind=wp), Public :: m7 = 0.005_wp  
  Real(Kind=wp), Public :: m8 = 0.02_wp  
  Real(Kind=wp), Public :: m9 = 0.95_wp  
  Real(Kind=wp), Public :: m10 = 0.65_wp  
  Real(Kind=wp), Public :: m11 = 0.35_wp  
  Real(Kind=wp), Public :: m12 = 0.5_wp  
  Real(Kind=wp), Public :: m13 = 2.5_wp  
  Real(Kind=wp), Public :: m14 = 1.25_wp  
  Real(Kind=wp), Public :: m15 = 8.0_wp  
  Real(Kind=wp), Public :: m16 = 7.0_wp  
  Real(Kind=wp), Public :: m17 = 8.0_wp  
  Real(Kind=wp), Public :: m18 = 0.2_wp  
End Module magic_numbers
