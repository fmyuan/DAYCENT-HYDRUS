
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      real function line(x, x1, y1, x2, y2)

      implicit none

c ... Argument declarations
      real x, x1, y1, x2, y2

c ... This function is the generic equation of a line from
c ... two points.
c ... Given 2 known points and a new X, calculate Y.
c ... slope = (y2 - y1) / (x2 - x1)
c ... y = slope * (x - x2) + y2
c ...
c ... CALLED FROM:  co2eff.f
c ...               froota.f
c ...               grem.f
c ...               prelim.f
c ...               ramp.f

      line = (y2 - y1) / (x2 - x1) * (x - x2) + y2

      return
      end
