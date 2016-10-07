
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine cropDynC(agprod, rtsh)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'potent.inc'

c ... Argument declarations
      real agprod, rtsh

c ... Compute carbon allocation fractions for aboveground and belowground
c ... plant parts.
c ...   agprod - estimated aboveground plant production
c ...   rtsh   - root/shoot ratio

c ... Function declarations
      real     froota, rtimp
      external froota, rtimp

c ... Local variables
      integer  iel, lyr
      real     availm(MAXIEL), demand, eavail(MAXIEL), fracrc, maxNfix,
     &         rimpct, totale

c ... Estimate the fraction of carbon going to the roots
      fracrc = (frtc(1)+frtc(2))/2.0

c ... Estimate total production
      tgprod = agprod / (1 - fracrc)

c ... Determine nutrients available to plants for growth.
      do 20 iel = 1, nelem
        availm(iel) = 0.0
        do 30 lyr = 1, nlayer
          availm(iel) = availm(iel) + minerl(lyr, iel)
30      continue
20    continue

c ... Calculate impact of root biomass on available nutrients
      rimpct = rtimp(riint, rictrl, bglivc)

c ... Calculate soil available nutrients, based on a maximum fraction
c ... (favail) and the impact of root biomass (rimpct), adding storage.
      do 45 iel = 1, nelem
        eavail(iel) = (availm(iel) * favail(iel) * rimpct) + 
     &                 crpstg(iel)
45    continue
        
c ... Compute the minimum and maximum C/E ratios
      call fltce(nelem, aglivc, co2cce)

c ... Estimate the demand
      do 50  iel = 1, nelem
c ..... Initialize fixation to 0
        maxNfix = 0.0
c ..... N FIXATION
        if (iel .eq. N) then
          maxNfix = snfxmx(CRPSYS) * (tgprod / 2.5)
        endif
c ..... DEMAND based on the maximum E/C ratio.
        demand = 0.0
        demand = demand + ((agprod / 2.5) *
     &                    (1.0 / cercrp(IMIN,ABOVE,iel)))
        demand = demand + ((tgprod - agprod) / 2.5) *
     &                    (1.0 / cercrp(IMIN,BELOW,iel))
        totale = eavail(iel) + maxNfix

c ..... New calculation -mdh 5/10/01
        crop_a2drat(iel) = min(1.0, totale / demand)
        crop_a2drat(iel) = max(0.0, crop_a2drat(iel))
50    continue

c ... New way of calculating fracrc, the fraction of root carbon. -mdh 3/21/00
c ... crop_a2drat(iel) = ratio of available mineral to mineral demand

      fracrc = froota(crop_a2drat,h2ogef(1),CRPSYS)
      rtsh = fracrc/(1 - fracrc)

      return
      end