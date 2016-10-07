
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine surftemp(elitst, pmxtmp, pmntmp, pmxbio, tmax, tmin,
     &                    tmxs, tmns, stemp, jdaywk)

      implicit none
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'site.inc'
      include 'zztim.inc'

c ... Argument declarations
      real elitst, pmxtmp, pmntmp
      real pmxbio, tmax, tmin
      real tmxs, tmns, stemp
      integer jdaywk

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
c ...   woodb  - wood biomass, fine branch + large wood (g/m^2)

c ... Outputs:
c ...   stemp - average soil surface temperature
c ...   tmns  - minimum soil surface temperature
c ...   tmxs  - maximum soil surface temperature

c ... Fortran to C prototype
      INTERFACE
        SUBROUTINE daylen(jdaywk, sitlat, daylength)
          !MS$ATTRIBUTES ALIAS:'_daylen' :: daylen
          INTEGER jdaywk
          REAL    sitlat
          REAL    daylength
        END SUBROUTINE daylen
      END INTERFACE

c ... Local Variables
      real bio, woodbio
      real daylength, tmns_mlt, tmxs_mlt
      real tmns_leaf, tmns_wood, tmxs_leaf, tmxs_wood

c ... Leaf biomass
      bio = aglivb + stdead + elitst * sfclit
      bio = min(bio, pmxbio)

c ... Wood biomass
      woodbio = min(woodb, 5000.0)

c ... Maximum temperature with leaf shading
      tmxs_leaf = tmax+(25.4/(1.+18.*exp(-.20*tmax)))*
     &            (exp(pmxtmp*bio)-.13)
c ... Minimum temperature with leaf shading
      tmns_leaf = tmin+pmntmp*bio-1.78

c ... Maximum temperature with wood shading
      tmxs_wood = tmax+(25.4/(1.+18.*exp(-.20*tmax)))*
     &            (exp(pmxtmp*0.1*woodbio)-.13)
c ... Minimum temperature with wood shading
      tmns_wood = tmin+pmntmp*0.1*woodbio-1.78

      tmxs = min(tmxs_leaf, tmxs_wood)
      tmns = max(tmns_leaf, tmns_wood)

c ... Average surface temperature
c ... Note: Soil temperature used to calculate potential production does not
c ...       take into account the effect of snow (AKM)
c      stemp = (tmxs+tmns)/2.
c ... This code captures the diurnal temperature variation
c      stemp = 0.41*tmxs + 0.59*tmns
c ... Use day length to compute the values for the multipliers on minimum
c ... and maximum soil surface temperature, cak - 04/25/03
c ... Compute daylength
      call daylen(jdaywk, sitlat, daylength)
      if (daylength .lt. 12.0) then
        tmns_mlt = ((12.0 - daylength) * 3.0 + 12.0) / 24.0
      else
        tmns_mlt = ((12.0 - daylength) * 1.2 + 12.0) / 24.0
      endif
      tmns_mlt = min(0.95, tmns_mlt)
      tmns_mlt = max(0.05, tmns_mlt)
      tmxs_mlt = 1.0 - tmns_mlt
      stemp = tmxs_mlt*tmxs + tmns_mlt*tmns

c ... Added code to allow warming of the soil surface temperature,
c ... cak - 07/02/03
      if (stsys .gt. 0 .and. time .ge. ststart) then
        stemp = stemp + stamt
      endif

      return
      end
