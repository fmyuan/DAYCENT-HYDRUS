
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine erosn(psloss,bulkd,edepth,enrich,
     &                 lhzci,lhze,nelem)

      implicit none
      include 'const.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'

c ... Argument declarations
      integer   nelem
      real      psloss, bulkd, edepth, enrich, lhzci(3,2), lhze(3,3)

c ... Soil is removed from the system.

c ... psloss is the rate of soil loss (kg/m**2/m) or per week -mdh 1/97
c ... bulkd is the bulk density of the soil (kg/liter)
c ... edepth is the depth of the soil (m)
c ... enrich is the enrichment factor for SOM losses
c ... scloss is the total carbon loss from soil organic
c ...   matter and below ground litter for this month or week
c ... sclosa is the accumulator for scloss
c ... lhzci(pool,iso) is the carbon in the lower horizon pools
c ... lhze(pool,iel) is the organic N,P,S in the lower horizon pools
c ... lhzcac is the accumulator for carbon inputs from the
c ...   lower horizon pools
c ... lhzeac(iel) is the accumulator for N,P,S inputs from the
c ...   lower horizon pools

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

c ... Local variables
      integer   iel, iso
      real      flost, input

c ... Compute fraction of soil which is lost.
      flost = psloss/(1000.*bulkd*edepth) * enrich

c ... Total carbon loss
      scloss = somtc * flost
      sclosa = sclosa + scloss

c ... Soil losses from below ground som1, som2, som3, below ground
c ... metabolic and structural pools
c ... Argument after nelem tells how many layers of som1, som2, or
c ... som3 are being modeled.  Only the soil layer is used in soilos.
c ... vek 08/91
      call soilos(time,nelem,2,flost,som1c(SOIL),som1ci,csrsnk,som1e,
     &            esrsnk)
      call soilos(time,nelem,1,flost,som2c,som2ci,csrsnk,som2e,
     &            esrsnk)
      call soilos(time,nelem,1,flost,som3c,som3ci,csrsnk,som3e,
     &            esrsnk)
      call soilos(time,nelem,2,flost,metabc(SOIL),metcis,csrsnk,
     &            metabe,esrsnk)
      call soilos(time,nelem,2,flost,strucc(SOIL),strcis,csrsnk,
     &            struce,esrsnk)

c ... This section commented out. It is assumed that losses from mineral
c ... pools are replaced by an equivalent amount from the next layer
c ... Calculate soil losses from mineral pools based on edepth
c      real eloss, sum
c      flost = psloss/(1000.*bulkd*edepth)
c      do 10 iel = 1, nelem
c        eloss = parent(iel)*flost
c        call flow(parent(iel),esrsnk(iel),time,eloss)
c        eloss = secndy(iel)*flost
c        call flow(secndy(iel),esrsnk(iel),time,eloss)
c        if (iel .eq. P) then
c          eloss = occlud*flost
c          call flow(occlud,esrsnk(iel),time,eloss)
c        endif
c10    continue
c
c ... Calculate soil losses from mineral pools based on adep(1)
c ... Compute fraction of soil which is lost.
c      flost = psloss/(1000.*bulkd*adep1)
c
c      do 20 iel = 1, nelem
c        eloss = minerl(SRFC,iel)*flost
c        call flow(minerl(SRFC,iel),esrsnk(iel),time,eloss)
c20    continue

c ... Calculate input of organic matter from next soil horizon
c ... using an equivalent depth of soil to that eroded
      do 30 iso = UNLABL, LABELD 
        input = flost*lhzci(1,iso)
        lhzci(1,iso) = lhzci(1,iso) - input
        lhzcac = lhzcac + input
        call flow(csrsnk(iso),som1ci(SOIL,iso),time,input)
        input = flost*lhzci(2,iso)
        lhzci(2,iso) = lhzci(2,iso) - input
        lhzcac = lhzcac + input
        call flow(csrsnk(iso),som2ci(iso),time,input)
        input = flost*lhzci(3,iso)
        lhzci(3,iso) = lhzci(3,iso) - input
        lhzcac = lhzcac + input
        call flow(csrsnk(iso),som3ci(iso),time,input)
30    continue

      do 40 iel = 1, nelem
        input = flost*lhze(1,iel)
        lhze(1,iel) = lhze(1,iel) - input
        lhzeac(iel) = lhzeac(iel) + input
        call flow(esrsnk(iel),som1e(SOIL,iel),time,input)
        input = flost*lhze(2,iel)
        lhze(2,iel) = lhze(2,iel) - input
        lhzeac(iel) = lhzeac(iel) + input
        call flow(esrsnk(iel),som2e(iel),time,input)
        input = flost*lhze(3,iel)
        lhze(3,iel) = lhze(3,iel) - input
        lhzeac(iel) = lhzeac(iel) + input
        call flow(esrsnk(iel),som3e(iel),time,input)
40    continue

      return
      end
