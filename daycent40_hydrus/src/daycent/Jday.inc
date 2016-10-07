c************************************************************************
c  jday.inc
c  For use with the daily version of century (in getwth & simsom)
c
c  Melannie D. Hartman
c
c************************************************************************

      common /jdays/ timstep, usexdrvrs, leapyr, dysimo, frstdy, lstdy

      real timstep
      integer usexdrvrs
      logical leapyr

      integer dysimo(12)
      integer frstdy(12)
      integer lstdy(12)

      integer k
      integer idysimo(12)
      data (idysimo(k), k=1,12) /31,28,31,30,31,30,31,31,30,31,30,31/

      integer ilstdy(12)
      data (ilstdy(k), k=1,12) 
     &                 /31,59,90,120,151,181,212,243,273,304,334,365/

      integer ifrstdy(12)
      data (ifrstdy(k), k=1,12) 
     &                 /1,32,60,91,121,152,182,213,244,274,305,335/

      save /jdays/
