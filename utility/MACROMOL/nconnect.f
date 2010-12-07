      subroutine nconnect(icn,name)

c***********************************************************************
c
c     DL_POLY utility to assign number of neighbours for an atom type
c     entries are for the AMBER and GROMOS  force fields
c
c     author t.forester feb 1996
c     copyright daresbury laboratory 1996
c
c     itt
c     2010-10-30 17:20:53
c     1.3
c     Exp
c
c***********************************************************************

      character*4 name

      if(name(1:1).eq."C") then
         
        icn = 4
        if(name(2:2).le."F") then 
          if(name.eq."C   ") icn = 3   
          if(name.eq."C*  ") icn = 3    
          if(name.eq."C2  ") icn = 2   
          if(name.eq."C3  ") icn = 1
          if(name.eq."CA  ") icn = 3
          if(name.eq."CB  ") icn = 3
          if(name.eq."CC  ") icn = 3
          if(name.eq."CD  ") icn = 2
          if(name.eq."CE  ") icn = 2
          if(name.eq."CF  ") icn = 2

        elseif(name(2:2).le."M") then 

          if(name.eq."CG  ") icn = 2
          if(name.eq."CH  ") icn = 3
          if(name.eq."CH1 ") icn = 3
          if(name.eq."CH2 ") icn = 2
          if(name.eq."CH3 ") icn = 1
          if(name.eq."CI  ") icn = 2
          if(name.eq."CJ  ") icn = 2
          if(name.eq."CK  ") icn = 3
          if(name.eq."CM  ") icn = 3

        else

          if(name.eq."CN  ") icn = 3
          if(name.eq."CP  ") icn = 2
          if(name.eq."CQ  ") icn = 3
          if(name.eq."CR  ") icn = 3
          if(name.eq."CR51") icn = 2
          if(name.eq."CR61") icn = 2
          if(name.eq."CS1 ") icn = 3
          if(name.eq."CS2 ") icn = 2
          if(name.eq."CT  ") icn = 4
          if(name.eq."CV  ") icn = 3
          if(name.eq."CW  ") icn = 3

        endif
          
      elseif (name(1:1).eq."N") then
          
        if(name.eq."N   ") icn = 3
        if(name.eq."N*  ") icn = 3
        if(name.eq."N2  ") icn = 3
        if(name.eq."N3  ") icn = 4
        if(name.eq."NA  ") icn = 3     
        if(name.eq."NB  ") icn = 2
        if(name.eq."NC  ") icn = 2
        if(name.eq."NE  ") icn = 3
        if(name.eq."NL  ") icn = 4
        if(name.eq."NP  ") icn = 3
        if(name.eq."NR5 ") icn = 2
        if(name.eq."NR5*") icn = 3
        if(name.eq."NR6 ") icn = 2
        if(name.eq."NR6*") icn = 3
        if(name.eq."NT  ") icn = 3
        if(name.eq."NZ  ") icn = 3
          
      elseif (name(1:1).eq."O") then
          
        if(name.eq."O   ") icn = 1
        if(name.eq."O2  ") icn = 1
        if(name.eq."OA  ") icn = 2
        if(name.eq."OH  ") icn = 2
        if(name.eq."OM  ") icn = 1
        if(name.eq."OS  ") icn = 2
        if(name.eq."OW  ") icn = 2
          
      elseif (name(1:1).eq."H") then
          
        icn = 1
        if(name(1:2).eq."HW") icn =2
          
      elseif (name(1:1).eq."S") then
          
        icn = 4
        if(name(2:2).eq.'g') icn = 2 ! gromos sulphur
          
      elseif (name(1:1).eq."P") then
          
        icn = 4
          
      elseif (name(1:2).eq."LP") then
          
        icn = 1
        
      else

         icn = 0
       
       endif
       return 
       end

