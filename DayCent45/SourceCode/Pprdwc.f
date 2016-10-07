
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      real function pprdwc(wc,x,pprpts)

      implicit none

c ... Argument declarations
      real     wc, x, pprpts(3)

c ... This funtion returns a value for potential plant production
c ... due to water content.  Basically you have an equation of a
c ... line with a moveable y-intercept depending on the soil type.
c ... The value passed in for x is ((avh2o(1) + prcurr(month) + irract)/pet)

c ... pprpts(1):  The minimum ratio of available water to pet which
c ...             would completely limit production assuming wc=0.
c ... pprpts(2):  The effect of wc on the intercept, allows the
c ...             user to increase the value of the intercept and
c ...             thereby increase the slope of the line.
c ... pprpts(3):  The lowest ratio of available water to pet at which
c ...             there is no restriction on production.

c ... Local variables
      real     intcpt, slope

c ... The equation for the y-intercept (intcpt) is A+B*WC.  A and B
c ... determine the effect of soil texture on plant production based
c ... on the soil water holding capacity.

c ... Old way:
c      intcpt = 0.0 + 1.0 * wc
cc ... The second point in the equation is (.8,1.0)
c      slope = (1.0-0.0)/(.8-intcpt)
c      pprdwc = 1.0+slope*(x-.8)

c ... New way:
      intcpt = pprpts(1) + (pprpts(2) * wc)
      slope = 1.0 / (pprpts(3) - intcpt)
      pprdwc = 1.0 + slope * (x - pprpts(3))

      if (pprdwc .gt. 1.0) then
        pprdwc = 1.0
      elseif (pprdwc .lt. 0.01) then
        pprdwc = 0.01
      endif

      return
      end
