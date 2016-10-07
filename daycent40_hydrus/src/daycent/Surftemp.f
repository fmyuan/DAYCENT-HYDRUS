
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine surftemp(aglivb, sfclit, stdead, elitst, pmxtmp, 
     &                    pmntmp, pmxbio, tmax, tmin, tmxs, tmns,
     &                    stemp)

      implicit none

c ... Argument declarations
      real aglivb, sfclit, stdead, elitst, pmxtmp, pmntmp
      real pmxbio, tmax, tmin
      real tmxs, tmns, stemp

c ... Compute the minimum, maximum, and average soil surface temperatures

c ... Inputs:
c ...   aglivb - above ground live biomass (g/m^2)
c ...   elitst - effect of litter on soil temperature relative
c ...            to live and standing dead biomass (fix.100)
c ...   pmntmp - effect of biomass on minimum surface temperature (fix.100)
c ...   pmxbio - maximum biomass for soil temperature calculations (fix.100)
c ...   pmxtmp - effect of biomass on maximum surface temperature (fix.100)
c ...   sfclit - surface litter biomass (g/m^2)
c ...   stdead - above ground standing dead biomass (g/m^2)
c ...   tmax   - maximum air temperature (deg C)
c ...   tmin   - minimum air temperature (deg C)

c ... Outputs:
c ...   stemp - average soil surface temperature
c ...   tmns  - minimum soil surface temperature
c ...   tmxs  - maximum soil surface temperature

c ... Local Variables
      real bio

      bio = aglivb + stdead + elitst * sfclit
      bio = min(bio,pmxbio)

c ... Maximum temperature
      tmxs = tmax+(25.4/(1.+18.*exp(-.20*tmax)))*
     &       (exp(pmxtmp*bio)-.13)

c ... Minimum temperature
      tmns = tmin+pmntmp*bio-1.78

c ... Average surface temperature
c ... Note: Soil temperature used to calculate potential production does not
c ...       take into account the effect of snow (AKM)
c      stemp = (tmxs+tmns)/2.
c ... This code captures the diurnal temperature variation
      stemp = 0.41*tmxs + 0.59*tmns

      return
      end
