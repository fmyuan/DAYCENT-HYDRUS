
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... RAMP.F

      real function ramp(x, x1, y1, x2, y2)

      implicit none

c ... Argument declarations
      real      x, x1, y1, x2, y2

c ... This function models a "ramp":
c ...             /-----
c ...            /
c ...      -----/

c ... Function declarations
      real      line
      external  line

      if (x .lt. x1) then 
        ramp = y1
      else if (x .ge. x2) then
        ramp = y2
      else
        ramp = line(x, x1, y1, x2, y2)
c ..... where line = (y2 - y1) / (x2 - x1) * (x - x2) + y2
      endif

      return
      end
