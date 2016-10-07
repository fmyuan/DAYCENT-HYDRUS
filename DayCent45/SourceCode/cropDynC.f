
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine cropDynC(rtsh, fracrc)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'potent.inc'

c ... Argument declarations
      real fracrc, rtsh

c ... Compute carbon allocation fractions for aboveground and belowground
c ... plant parts.
c ...   agprod - estimated aboveground plant production
c ...   rtsh   - root/shoot ratio

c ... Function declarations
      real     froota, rtimp
      external froota, rtimp

c ... Local variables
      integer  iel, lyr
      real     agprod, availm(MAXIEL), demand, eavail(MAXIEL), maxNfix,
     &         rimpct, totale

c ... Estimate the fraction of carbon going to the roots
      if (frtcindx .eq. 0) then
c ..... Use Great Plains equation for root to shoot ratio
        rtsh = (bgppa + grwprc*bgppb)/(agppa + grwprc*agppb)
        fracrc = 1.0 / (1.0/rtsh + 1.0)
      elseif (frtcindx .eq. 1 .or. frtcindx .eq. 3) then
c ..... A perennial plant (grass)
        fracrc = (cfrtcw(1)+cfrtcw(2)+cfrtcn(1)+cfrtcn(2))/4.0
      elseif (frtcindx .eq. 2 .or. frtcindx .ge. 4) then
c ..... An annual plant (crop)
        fracrc = (frtc(1)+frtc(2))/2.0
      endif

c ... Estimate total production
c      tgprod = agprod / (1 - fracrc)
c ... Estimate aboveground production, cak - 08/22/03
      agprod = tgprod * (1 - fracrc)

c ... Determine nutrients available to plants for growth.
      do 20 iel = 1, nelem
        availm(iel) = 0.0
c ..... Nutrients available to grasses/crops are in the top claypg layers,
c ..... cak 01/29/03
c        do 30 lyr = 1, nlayer
        do 30 lyr = 1, claypg
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
