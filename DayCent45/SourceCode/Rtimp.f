
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... RTIMP.F

      real function rtimp(riint, rictrl, rootc)

      implicit none

c ... Argument declarations
      real      riint, rictrl, rootc

c ... This function calculates and returns a value between 0-1 which is the impact
c ... of root biomass on available nutrients.  It is used in the calculation of
c ... total plant production in RESTRP.

c ... Called From:     GROWTH
c ...                  CROPDYNC
c ...                  TREEDYNC
c ...                  TREEGROW

c ... Check added to handle underflow potential in exp intrinsic
      if ((rictrl * rootc * 2.5) .gt. 33) then
        rtimp = 1.0
      else
        rtimp = (1.0 - riint * exp(-rictrl * rootc * 2.5))
      endif

      return
      end
