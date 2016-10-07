
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... LITDEC.F

      subroutine litdec(dtm, newminrl)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'

c ... Argument declarations
      real      dtm
      real*8    newminrl

c ... Litter Decomposition
c ... Decompose structural and metabolic material for surface and soil.
c ... written by vek 04/91

c ... Function declarations
      logical  candec
      real     agdrat, bgdrat, catanf
      external agdrat, bgdrat, catanf

c ... Local variables
      integer  iel, lyr
      real     accum(ISOS), biocnv, cfmes1, co2los, mnrflo, pheff,
     &         rceto1(3), tcflow

c ... Factor to convert C to biomass is 2.5 for everything but wood.
      parameter (biocnv = 2.5)

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

c ... Surface STRUCTURAL Material
      if(strucc(SRFC) .gt. 1.e-07) then

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

c ..... Compute total C flow out of structural in layer SRFC
c ..... Add pH effect on decomposition to calculation, cak - 08/02/02
c        tcflow = min(strucc(SRFC),strmax(SRFC)) * defac * dec1(SRFC) *
c     &               exp(-pligst(SRFC)*strlig(SRFC)) * dtm
        tcflow = min(strucc(SRFC),strmax(SRFC)) * defac * dec1(SRFC) *
     &               exp(-pligst(SRFC)*strlig(SRFC)) * dtm * pheff
c ..... where
c .....   tcflow       = grams C
c .....   strucc(SRFC) = the current value for total structural C
c .....                  (grams C) in layer SRFC
c .....   strmax(SRFC) = the maximum amount of structural C in 
c .....                  layer SRFC that will decompose (grams C).
c .....   defac        = the decomposition factor based on water and 
c .....                  temperature computed in prelim and in cycle
c .....   dec1(SRFC)   = the intrinsic decomposition rate of 
c .....                  structural C (a fixed parameter)
c .....   pligst       = a fixed parameter that represents the effect 
c .....                  of lignin-to-structural-ratio on structural 
c .....                  decomposition
c .....   strlig(SRFC) = the lignin content of surface structural 
c .....                  residue (grams lignin/grams biomass)
c .....   dtm          = the time step in years

c ..... Decompose structural into som1 and som2 with CO2 loss.
c ..... Changed csrsnk to st1c2 (10/92)

        call declig(aminrl,strlig(SRFC),SRFC,nelem,2,ps1co2,rnewas,
     &              rsplig,tcflow,strucc,st1c2,strcis,struce,
     &              gromin,minerl,strmnr,resp,som1ci,som1e,som2ci,
     &              som2e,newminrl)

      endif

c ... Soil STRUCTURAL Material
      if (strucc(SOIL) .gt. 1.e-07) then

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

c ..... Compute total C flow out of structural in layer SOIL
c ..... Added impact of soil anerobic conditions -rm 12/91
c ..... Add pH effect on decomposition to calculation, cak - 08/02/02
c        tcflow = min(strucc(SOIL),strmax(SOIL)) * defac * dec1(SOIL) *
c     &               exp(-pligst(SOIL)*strlig(SOIL)) * cltfac(4) * 
c     &               anerb * dtm
        tcflow = min(strucc(SOIL),strmax(SOIL)) * defac * dec1(SOIL) *
     &               exp(-pligst(SOIL)*strlig(SOIL)) * cltfac(4) * 
     &               anerb * dtm * pheff
c ..... where
c .....   tcflow       = grams C
c .....   strucc(SOIL) = the current value for total structural C
c .....                  (grams C) in layer SOIL
c .....   strmax(SOIL) = the maximum amount of structural C in 
c .....                  layer SOIL that will decompose (grams C).
c .....   defac        = the decomposition factor based on water and 
c .....                  temperature computed in prelim and in cycle
c .....   dec1(SOIL)   = the intrinsic decomposition rate of 
c .....                  structural C (a fixed parameter)
c .....   pligst       = a fixed parameter that represents the effect 
c .....                  of lignin-to-structural-ratio on structural 
c .....                  decomposition
c .....   strlig(SOIL) = the lignin content of soil structural 
c .....                  residue (grams lignin/grams biomass)
c .....   cltfac(4)    = the cultivation factor for soil
c .....                  structural material (set in cycle)
c .....   dtm          = the time step in years (= dt/ntspm)

c ..... Decompose structural into som1 and som2 with CO2 loss.
c ..... Changed csrsnk to st2c2 (10/92)

        call declig(aminrl,strlig(SOIL),SOIL,nelem,2,ps1co2,rnewbs,
     &              rsplig,tcflow,strucc,st2c2,strcis,struce,
     &              gromin,minerl,strmnr,resp,som1ci,som1e,som2ci,
     &              som2e,newminrl)

      endif
c ... End of Structural Decomposition

c ... METABOLIC Material

c ... Process each layer
      do 30 lyr = SRFC, SOIL

        if (metabc(lyr) .gt. 1.e-07) then

c ....... Determine C/E ratios for flows to SOM1
          do 10 iel=1,nelem

c ......... Compute ratios for surface metabolic residue
            if (lyr .eq. SRFC) then
              rceto1(iel) = agdrat(metabe(lyr,iel),metabc(lyr),
     &                             biocnv,iel,pcemic,cemicb)

c ......... Compute ratios for soil metabolic residue
            else
              rceto1(iel) = bgdrat(aminrl,varat1,iel)
            endif
10        continue

c ....... Compute pH effect on decomposition
          pheff = catanf(ph, 4.8, 0.5, 1.14, 0.7)
          pheff = min(pheff, 1.0)
          pheff = max(pheff, 0.0)

c ....... Compute total C flow out of metabolic in layer lyr
c ....... Add pH effect on decomposition to calculation, cak - 08/02/02
c          tcflow = metabc(lyr) * defac * dec2(lyr) * dtm
          tcflow = metabc(lyr) * defac * dec2(lyr) * dtm * pheff

c ....... Added impact of soil anerobic conditions -rm 12/91
          if (lyr .eq. SOIL) then
            tcflow = tcflow * anerb
          endif
c ....... where:
c .......   tcflow      = grams C
c .......   metabc(lyr) = the current value for total metabolic C
c .......                 (grams C) in layer lyr
c .......   defac       = the decomposition factor
c .......   dec2(lyr)   = the intrinsic decomposition rate of 
c .......                 metabolic C
c .......   dtm         = the time step in years

c ....... Make sure metab does not go negative.
          if (tcflow .gt. metabc(lyr)) then
            tcflow = metabc(lyr)
          endif

c ....... If decomposition can occur,
          if (candec(nelem,aminrl,metabc(lyr),metabe,2,lyr,rceto1))
     &        then

c ......... CO2 loss
            co2los = tcflow * pmco2(lyr)

c ......... Changed csrsnk to mt1c2, mt2c2 (10/92)
            if (lyr .eq. SRFC) then
              call respir(co2los,2,lyr,metabc,metcis,mt1c2,
     &                    resp,metabe,minerl,gromin,metmnr,newminrl)
            else
              call respir(co2los,2,lyr,metabc,metcis,mt2c2,
     &                    resp,metabe,minerl,gromin,metmnr,newminrl)
            endif

c ......... Decompose metabolic into som1
            cfmes1 = tcflow - co2los

c ......... Partition and schedule C flows by isotope
            call csched (cfmes1,metcis(lyr,LABELD),metabc(lyr),
     &                   metcis(lyr,UNLABL),som1ci(lyr,UNLABL),
     &                   metcis(lyr,LABELD),som1ci(lyr,LABELD),
     &                   1.0,accum)

c ......... Compute and schedule N, P, and S flows and update
c ......... mineralization accumulators.
            do 20 iel = 1, nelem
             
              call esched(cfmes1,metabc(lyr),rceto1(iel),
     &                    metabe(lyr,iel),som1e(lyr,iel),
     &                    minerl(SRFC,iel),mnrflo)
              call mnracc(mnrflo,gromin(iel),metmnr(lyr,iel))
c ........... newminrl should be updated only for nitrogen
c ........... akm via cak 07/31/01
              if (iel .eq. N) then
                newminrl = newminrl + mnrflo
              endif
20          continue
          endif
        endif

c ..... End of Metabolic Decomposition

c ... Next layer
30    continue

      return
      end
