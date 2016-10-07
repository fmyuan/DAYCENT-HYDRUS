
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine phshift(time)

      implicit none
      include 'param.inc'

c ... Argument declarations
      real time

c ... Compute the shift in pH.

c ... Function declarations
      real      ramp
      external  ramp

c ... Calculate the change in pH

c ... Calculate the new pH value
      if (time .le. phtm(1)) then
c ..... use pH value as read from <site>.100 file
        ph = phstart

      else if (time .ge. phtm(2)) then
c ..... we have reached the final pH value
        ph = phend

      else
c ..... ramping
        ph = ramp(time,phtm(1),phstart,phtm(2),phend)
      endif

      return
      end
