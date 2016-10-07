
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... SOMDEC.F

      subroutine somdec(amovdly, dtm, newminrl)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'

c ... Argument declarations
      real      amovdly(MAXLYR)
      real      dtm
      real*8    newminrl

c ... Soil Organic Matter Decomposition                written by vek, 04/91
c ... Decompose SOM1 (surface and soil), SOM2, and SOM3.
c ... defac = decomposition factor based on water and
c ...         temperature computed in prelim and in cycle
c ... dtm   = time step (dt/ntspm)

c ... Fortran to C prototype
      INTERFACE
        SUBROUTINE flow(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow' :: flow
          REAL from
          REAL to
          REAL when
          REAL howmuch
        END SUBROUTINE flow
      END INTERFACE

c ... Function declarations
      real      bgdrat
      logical   candec
      real      catanf
      external  catanf

c ... Local variables
      integer   iel
      real      accum(ISOS), cf3co2, cfs1s2, cfs1s3, cfs2s1, cfs2s3, 
     &          cfs3s1, cfsfs2, cleach, co2los, linten, mnrflo, orgflow, 
     &          pheff, radds1, rceof1, rceto1(MAXIEL), rceto2(MAXIEL), 
     &          rceto3(MAXIEL), rcetob, tcflow
      real      somc(2)

c ... *******************************************************************
c ... Initialize ACCUM even though it is not really used
      accum(UNLABL) = 0.0
      accum(LABELD) = 0.0
      somc(1) = 0.0
      somc(2) = 0.0

c ... Surface SOM1 decomposes to SOM2 with CO2 loss
      if (som1c(SRFC) .gt. 1.e-07) then

c ..... Determine C/E ratios for flows to som2
        do 10 iel=1,nelem
          radds1 = rad1p(1,iel) + rad1p(2,iel) *
     &             ((som1c(1)/som1e(1,iel))-pcemic(2,iel))
          rceto2(iel) = som1c(SRFC)/som1e(SRFC,iel) + radds1
          rceto2(iel) = max(rceto2(iel), rad1p(3,iel))
10      continue

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

c ..... Compute total C flow out of surface microbes.
c ..... Add pH effect on decomposition to calculation, cak - 08/02/02
c        tcflow = som1c(SRFC) * defac * dec3(SRFC) * dtm
        tcflow = som1c(SRFC) * defac * dec3(SRFC) * dtm * pheff
c ..... where
c .....   som1c(SRFC) =  unlabeled and labeled C in surface microbes
c .....                  (som1ci(SRFC,1)+som1ci(SRFC,2))
c .....   dec3(SRFC)  =  intrinsic decomposition rate of
c .....                  surface microbes

c ..... If decomposition can occur, schedule flows associated with respiration 
c ..... and decomposition 
        if (candec(nelem,aminrl,som1c(SRFC),som1e,2,SRFC,rceto2)) then

c ....... CO2 loss - Compute and schedule respiration flows.
          co2los = tcflow * p1co2(SRFC)
c ....... where 
c .......   p1co2(SRFC)  = set to p1co2a(SRFC) in prelim.
c .......   p1co2a(SRFC) = a fixed parameter;
c .......                  intercept parameter which controls flow from soil
c .......                  organic matter with fast turnover to CO2 (fraction
c .......                  of carbon lost to CO2 when there is no sand in the
c .......                  soil)
c ....... Changed csrsnk to s11c2 (10/92)

          call respir(co2los,2,SRFC,som1c,som1ci,s11c2,resp,
     &                som1e,minerl,gromin,s1mnr,newminrl)

c ....... Decompose Surface SOM1 to SOM2

c ....... cfsfs2 is C Flow from SurFace som1 to Som2
          cfsfs2 = tcflow - co2los

c ....... Partition and schedule C flows by isotope
          call csched(cfsfs2,som1ci(SRFC,LABELD),som1c(SRFC),
     &                som1ci(SRFC,UNLABL),som2ci(UNLABL),
     &                som1ci(SRFC,LABELD),som2ci(LABELD),
     &                1.0,accum)

c ....... Compute and schedule N, P, and S flows.

c ....... Update mineralization accumulators.
          do 20 iel=1,nelem
            call esched(cfsfs2,som1c(SRFC),rceto2(iel),
     &                  som1e(SRFC,iel),som2e(iel),
     &                  minerl(SRFC,iel),mnrflo)
            call mnracc(mnrflo,gromin(iel),s1mnr(SRFC,iel))
c ......... newminrl should be updated only for nitrogen
c ......... akm via cak 07/31/01
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
20        continue
        endif
      endif

c ... End of SOM1 (surface layer) Decomposition

c ... *******************************************************************

c ... Soil SOM1 decomposes to SOM2 and SOM3 with CO2 loss and possible leaching 
c ... of organics.
      if (som1c(SOIL) .gt. 1.e-07) then

c ..... Determine C/E ratios for flows to som2
        do 30 iel=1,nelem
          rceto2(iel) = bgdrat(aminrl,varat2,iel)
30      continue

c ..... Compute total C flow out of soil microbes.

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 4.8, 0.5, 1.14, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

c ..... Added impact of soil anaerobic conditions -rm 12/91
c ..... Add pH effect on decomposition to calculation, cak - 08/02/02
c        tcflow = som1c(SOIL) * defac *dec3(SOIL) *cltfac(1) * eftext *
c     &           anerb * dtm
        tcflow = som1c(SOIL) * defac *dec3(SOIL) *cltfac(1) * eftext *
     &           anerb * dtm * pheff
c ..... where
c .....   som1c(SOIL) = unlabeled and labeled C in soil microbes
c .....                 (som1ci(SOIL,1)+som1ci(SOIL,2))
c .....   dec3(SOIL)  = intrinsic decomposition rate of
c .....                 soil microbes
c .....   cltfac(1)   = cultivation factor for som1
c .....                 (set in cycle)
c .....   eftext      = effect of soil texture on the soil microbe
c .....                 decomposition rate (computed in prelim)

c ..... If soil som1 can decompose to som2, it will also go to som3.
c ..... If it can't go to som2, it can't decompose at all.

c ..... If decomposition can occur,
        if (candec(nelem,aminrl,som1c(SOIL),som1e,2,SOIL,rceto2)) then

c ....... CO2 Loss - Compute and schedule respiration flows
          co2los = tcflow * p1co2(SOIL)
c ....... where 
c .......   p1co2(SOIL) is computed in prelim as a function of fixed parameters 
c .......   p1co2a(SOIL) and p1co2b(SOIL) and soil texture.
c ....... Changed csrsnk to s21c2 (10/92)
          call respir(co2los,2,SOIL,som1c,som1ci,s21c2,resp,
     &                som1e,minerl,gromin,s1mnr,newminrl)

c ....... Decompose Soil SOM1 to SOM3
c ....... The fraction of tcflow that goes to SOM3 is a function of 
c ....... clay content (fps1s3 is computed in prelim).
          cfs1s3 = tcflow * fps1s3 * (1.0 + animpt * 
     &             (1.0 -anerb))

c ....... Partition and schedule C flows by isotope
          call csched(cfs1s3,som1ci(SOIL,LABELD),som1c(SOIL),
     &                som1ci(SOIL,UNLABL),som3ci(UNLABL),
     &                som1ci(SOIL,LABELD),som3ci(LABELD),
     &                1.0,accum)

c ....... Compute and schedule N, P, and S flows and update mineralization 
c ....... accumulators.
          do 40 iel=1,nelem
            rceto3(iel) = bgdrat(aminrl,varat3,iel)
            call esched(cfs1s3,som1c(SOIL),rceto3(iel),
     &                  som1e(SOIL,iel),som3e(iel),
     &                  minerl(SRFC,iel),mnrflo)
            call mnracc(mnrflo,gromin(iel),s1mnr(SOIL,iel))
c ......... newminrl should be updated only for nitrogen
c ......... akm via cak 07/31/01
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
40        continue

c ....... Leaching of Organics
c ....... This only occurs when the water flow out of water layer 2
c ....... exceeds a critical value.  Use the same C/N, C/P, and C/S
c ....... ratios as for the flow to SOM3.

c ....... Removed organic leaching sink and replaced it with the
c ....... stream flows. -rm 2/92
          if(amovdly(2) .gt. 0.0) then
            linten = min(1.0-(omlech(3)-amovdly(2))/
     &                   omlech(3), 1.0)
            cleach = tcflow * orglch * linten
c ......... Partition and schedule C flows by isotope
            call csched(cleach,som1ci(SOIL,LABELD),som1c(SOIL),
     &                  som1ci(SOIL,UNLABL),strm5u,
     &                  som1ci(SOIL,LABELD),strm5l,
     &                  1.0,accum)

c ......... Compute and schedule N, P, and S flows and update mineralization 
c ......... accumulators.
            do 50 iel=1,nelem

c ........... Need to use the ratio for som1 for organic leaching     -rm 3/92
              rceof1 = som1c(SOIL) / som1e(SOIL,iel)
              orgflow = cleach / rceof1
              call flow(som1e(SOIL,iel),stream(iel+5),time,
     &                  orgflow)
50          continue

c ....... No leaching this time step
          else
            cleach = 0.
          endif

c ....... Decompose Soil SOM1 to SOM2.
c ....... SOM2 gets what's left of tcflow.
          cfs1s2 = tcflow - co2los - cfs1s3 - cleach

c ....... Partition and schedule C flows by isotope
          call csched(cfs1s2,som1ci(SOIL,LABELD),som1c(SOIL),
     &                som1ci(SOIL,UNLABL),som2ci(UNLABL),
     &                som1ci(SOIL,LABELD),som2ci(LABELD),
     &                1.0,accum)

c ....... Compute and schedule N, P, and S flows and update mineralization 
c ....... accumulators.
          do 60 iel=1,nelem
            call esched(cfs1s2,som1c(SOIL),rceto2(iel),
     &                  som1e(SOIL,iel),som2e(iel),
     &                  minerl(SRFC,iel),mnrflo)
            call mnracc(mnrflo,gromin(iel),s1mnr(SOIL,iel))
c ......... newminrl should be updated only for nitrogen
c ......... akm via cak 07/31/01
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
60        continue
        endif
      endif

c ... End of Soil SOM1 decomposition

c ... *****************************************************************

c ... SOM2 decomposes to soil SOM1 and SOM3 with CO2 loss
      if (som2c .gt. 1.e-07) then

c ..... Determine C/E ratios for flows to SOM1
        do 70 iel=1,nelem
          rceto1(iel) = bgdrat(aminrl,varat1,iel)
70      continue

c ..... Compute total C flow out of SOM2C

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 3.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

c ..... Added impact of soil anaerobic conditions -rm 12/91
c ..... Add pH effect on decomposition to calculation, cak - 08/02/02
c        tcflow = som2c * defac * dec5 * cltfac(2) * anerb * dtm
        tcflow = som2c * defac * dec5 * cltfac(2) * anerb * dtm * pheff
c ..... where
c .....   dec5      = intrinsic decomposition rate of som2
c .....   cltfac(2) = cultivation factor for som2 (set in cycle)

c ..... If som2 can decompose to som1, it will also go to som3.
c ..... If it can't go to som1, it can't decompose at all.

c ..... If decomposition can occur,
        if (candec(nelem,aminrl,som2c,som2e,1,1,rceto1)) then

c ....... CO2 loss - Compute and schedule respiration flows
          co2los = tcflow * p2co2

c ....... Changed csrsnk to s2c2 (10/92)
          somc(1) = som2c
          call respir(co2los,1,SRFC,somc,som2ci,s2c2,resp,
     &                som2e,minerl,gromin,s2mnr,newminrl)

c ....... Rearranged the order of calculation.  There is a calculated
c ....... impact of clay on the decomposition of som2 to som3.  -rm 12/91

c ....... Decompose SOM2 to SOM3, SOM3 gets what's left of tcflow.
c ....... Added impact of soil anaerobic conditions -rm 12/91
          cfs2s3 = tcflow * fps2s3 * (1.0 + animpt * (1.0 - anerb))

c ....... Partition and schedule C flows by isotope
          call csched(cfs2s3,som2ci(LABELD),som2c,
     &                som2ci(UNLABL),som3ci(UNLABL),
     &                som2ci(LABELD),som3ci(LABELD),
     &                1.0,accum)

c ....... Compute and schedule N, P, and S flows and update mineralization 
c ....... accumulators.
          do 80 iel=1,nelem
            rcetob = bgdrat(aminrl,varat3,iel)
            call esched(cfs2s3,som2c,rcetob,som2e(iel),
     &                  som3e(iel),minerl(SRFC,iel),mnrflo)
            call mnracc(mnrflo,gromin(iel),s2mnr(iel))
c ......... newminrl should be updated only for nitrogen
c ......... akm via cak 07/31/01
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
80        continue

c ....... Decompose SOM2 to SOM1

c ....... Added impact of soil anaerobic conditions -rm 12/91
          cfs2s1 = tcflow  - co2los - cfs2s3

c ....... Partition and schedule C flows by isotope
          call csched(cfs2s1,som2ci(LABELD),som2c,
     &                som2ci(UNLABL),som1ci(SOIL,UNLABL),
     &                som2ci(LABELD),som1ci(SOIL,LABELD),
     &                1.0,accum)

c ....... Compute and schedule N, P, and S flows and update mineralization 
c ....... accumulators.
          do 90 iel=1,nelem
            call esched(cfs2s1,som2c,rceto1(iel),som2e(iel),
     &                  som1e(SOIL,iel),minerl(SRFC,iel),mnrflo)
            call mnracc(mnrflo,gromin(iel),s2mnr(iel))
c ......... newminrl should be updated only for nitrogen
c ......... akm via cak 07/31/01
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
90        continue
        endif
      endif

c ... End of SOM2 Decompositon

c ... **********************************************************************

c ... SOM3 decomposes to soil SOM1 with CO2 loss.
      if (som3c .gt. 1.e-07) then

c ..... Determine C/E ratios for flows to SOM1.
        do 100 iel=1,nelem
          rceto1(iel) = bgdrat(aminrl,varat1,iel)
100     continue

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 3.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

c ..... Compute total C flow out of SOM3C
c ..... Add pH effect on decomposition to calculation, cak - 08/02/02
c        tcflow = som3c * defac * dec4 * cltfac(3) * anerb * dtm
        tcflow = som3c * defac * dec4 * cltfac(3) * anerb * dtm * pheff
c ..... where
c .....   dec4      = intrinsic decomposition rate of som2
c .....   cltfac(3) = cultivation factor for som3 (set in cycle)

c ..... If decomposition can occur,
        if (candec(nelem,aminrl,som3c,som3e,1,1,rceto1)) then

c ....... CO2 loss - Compute and schedule respiration flows.
          cf3co2 = tcflow * p3co2 * anerb

c ....... Changed csrsnk to s3c2 (10/92)
          somc(1) = som3c
          call respir(cf3co2,1,SRFC,somc,som3ci,s3c2,resp,
     &                som3e,minerl,gromin,s3mnr,newminrl)

c ....... Decompose SOM3 to soil SOM1
          cfs3s1 = tcflow - cf3co2

c ....... Partition and schedule C flows by isotope
          call csched(cfs3s1,som3ci(LABELD),som3c,
     &                som3ci(UNLABL),som1ci(SOIL,UNLABL),
     &                som3ci(LABELD),som1ci(SOIL,LABELD),
     &                1.0,accum)

c ....... Compute and schedule N, P, and S flows and update mineralization 
c ....... accumulators.
          do 110 iel=1,nelem
            call esched(cfs3s1,som3c,rceto1(iel),som3e(iel),
     &                  som1e(SOIL,iel),minerl(SRFC,iel),mnrflo)
            call mnracc(mnrflo,gromin(iel),s3mnr(iel))
c ......... newminrl should be updated only for nitrogen
c ......... akm via cak 07/31/01
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
110       continue
        endif
      endif

c ... End of SOM3 Decomposition

      return
      end
