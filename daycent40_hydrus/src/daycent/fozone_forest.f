c
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... FOZONE.F

c ... ******************************************************************
      real function fozone_prd(x)
c
      implicit none
c
c ... Argument declarations
      real    x   ! W126 (unit: ppm-h)
c
c ...
c ... function to estimate ozone reduction factor to production
c ...
c ... This routine is based on a review by J.M. Pye. 1988. Impact of Ozone
c ... on the Growth and Yield of Trees: A Review. J. of Environmental
c ... Quality, 17: 347 - 360.

c ... 09/30/2012 by Fengming Yuan

c ... empirical parameters for the regression
      real a, b, c

c ... multiple linear regression
      a = 1.6e-5
      b = -0.00014
      c = 3.0

      fozone_prd = 0.
      if (x .gt. 0.) then
        fozone_prd = (a*x*x+b*x+c)/100.0
      endif

      if (fozone_prd .gt. 1.0) fozone_prd = 1.0
      if (fozone_prd .lt. 0.) fozone_prd = 0.

      return
      end



c ... ******************************************************************

      real function fozone_dfoliage(x)

      implicit none

c ... Argument declarations
      real    x   ! daily concentration (unit: ppb per hour, 24-h basis)

c ...
c ... function to estimate ozone reduction factor to production
c ...
c ... This routine is based on limited study across SBM ozone gradient
c ... by Grulke and Balduman 1999. Deciduous conifers: high N deposition
c ... and O3 exposure effects on growth and biomass allocation.
c ... Water, Air and Soil Pollution 116: 235-248.

c ... It's assumed that Chlorotic mottle percentage as defoliage rate

c ... 09/30/2012 by Fengming Yuan

c ... empirical parameters for the regression
      real a, b, c,d

c ... multiple linear regression
      a = 0.0011
      b = -0.17
      c = 8.7
      d = -1400.0

      fozone_dfoliage = 0.
      if (x .gt. 0.) then
        fozone_dfoliage = (a*x*x*x+b*x*x+c*x+d)/100.0
      endif

c ... 50% is the max. defoliage from the datasets
      if (fozone_dfoliage .gt. 0.50) fozone_dfoliage = 0.50
      if (fozone_dfoliage .lt. 0.) fozone_dfoliage = 0.

      return
      end

