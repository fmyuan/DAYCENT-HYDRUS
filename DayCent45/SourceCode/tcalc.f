
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      real function tcalc(stemp, teff)

      implicit none

c ... Argument declarations
      real stemp, teff(4)

c ... This function computes the effect of temperature on
c ... decomposition.  It is an exponential function.  Older
c ... versions of Century used a density function.
c ... Created 10/95 - rm
c ...
c ... The temperature effect is now being computed using an
c ... arctangent curve.  CAK - 03/16/01
c ...
c ... Called From:  calcdefac
c ...
c ... Variables
c ...   STEMP:    soil temperature
c ...   TEFF(1):  "x" location of inflection point
c ...   TEFF(2):  "y" location of inflection point
c ...   TEFF(3):  step size (distance from the maximum point to the minimum point)
c ...   TEFF(4):  slope of line at inflection point

c ... Allow C code to call Fortran function
      !MS$ATTRIBUTES C :: tcalc

c ... Function declarations
      real      catanf
      external  catanf

c ... Local variables
      real      normalizer

c      tcalc = teff(1) + teff(2) * exp(teff(3) * stemp)

c ... The normalizer is the value of the numerator at 30 deg C

      normalizer = catanf(30.0, teff(1), teff(2), teff(3), teff(4))

      tcalc = max(0.01,
     &            catanf(stemp, teff(1), teff(2), teff(3), teff(4)) /
     &            normalizer)

      return
      end
