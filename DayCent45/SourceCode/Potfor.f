
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine potfor(month,irractwk,pptwk,petwk,tfrac,tavemth,
     &                  jdaywk)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'dynam.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'site.inc'

c ... Argument declarations
      integer month, jdaywk
      real    irractwk, pptwk, petwk, tfrac, tavemth

c ... Compute monthly potential production for forest
c ...
c ... Outputs
c ...   pforc gC/m^2/time
c ...   h2ogef(2)

c ... RESPPT is an array of respiration values for the forest system
c ... production components.

c ... Added savanna model, pcropc now pforc
c ...                      tgprod now tfprod (BO)

c ... Function declarations
      real     gpdf, frespr, laprod, pprdwc, shwave
      external gpdf, frespr, laprod, pprdwc, shwave

c ... Local variables
      real     fcmax, frlive, lai, potprd, tmoist

c ... Estimate potential production based on temp & h2o

      if (stemp .gt. 0.0) then

c ..... Calculate temperature effect on growth.
        potprd = gpdf(stemp, ppdf(1,2), ppdf(2,2), ppdf(3,2), 
     &                ppdf(4,2))

c ..... Added to match version 3.0 -lh 4/93
        potprd = potprd * .8

c ..... Calculate moisture effect on growth -
c ..... Value for potential plant production is now calculated from the
c ..... equation of a line whose intercept changes depending on water
c ..... content based on soil type.  The function PPRDWC contains the
c ..... equation for the line.  -rm 9/91
 
        tmoist = pptwk + irractwk
 
        if (petwk .ge. .01) then
c ....... Allow trees access to water in the profile -mdh 11/96.
c          pptprd = (avh2o(1) + tmoist) / petwk
c          h2ogef(2) = (avh2o(2) + tmoist) / petwk
c ....... The potential plant growth was not responding as expected to
c ....... water stress.  The suspicion is that the original potential
c ....... production equation was optimized for a monthly time step.  The
c ....... code modification is designed to work for weekly time step
c ....... cak - 11/12/01, change suggested by Steve Del Grosso
c          h2ogef(2) = (avh2o(2) + (tmoist/tfrac))/(petwk/tfrac)
c ....... Calculate potential growth based on the relative water content
c ....... of the wettest soil layer, cak - 12/06/04
          h2ogef(2) = 1.0/(1.0 + 30.0 * exp(-9.0 * twstress))
        else
          h2ogef(2) = 0.01
        endif

c        h2ogef(2) = pprdwc(wc,h2ogef(2),pprpts)

c ..... For large wood, calculate the percentage which is live (sapwood)
        frlive = sapk / (sapk + rlwodc)
 
c ..... Calculate LAI and then use in function for LAI reducer on
c ..... production.  -rm 5/91
c ..... prdx(2) is maximum GROSS production in g biomass/m^2/mo
c ..... Calculate theoretical maximum for LAI based on large wood biomass,
c ..... cak - 07/24/02
c ... Include fine branch carbon as part of the woody component in the
c ... LAI calculation, cak - 10/20/2006
c        call lacalc(lai, rleavc, rlwodc, btolai, maxlai, klai)
c        call lacalc(lai, rlwodc, maxlai, klai)
        call lacalc(lai, fbrchc, rlwodc, maxlai, klai)
c ..... No longer using the maximum gross production parameter to calculate
c ..... potential production, cak - 07/01/02
c        tfprod = prdx(2)*tfrac*potprd*h2ogef(2)*laprod(lai,laitop)*
c     &           co2cpr(FORSYS)

c ..... Calculate monthly maintenance respiration
c ..... Maintenance respiration is now being calculated in the treegrow
c ..... subroutine
c        resppt(LEAF) = frespr(tavewk,rleave(N))*tfrac
c        resppt(FROOT) = frespr(stemp,froote(N))*tfrac
c        resppt(FBRCH) = frespr(tavewk,fbrche(N))*tfrac
c        resppt(LWOOD) = frespr(tavewk,frlive*rlwode(N))*tfrac
c        resppt(CROOT) = frespr(stemp,frlive*croote(N))*tfrac

c        sumrsp = 0.0
c        do 10 ii = 1, FPARTS
c          sumrsp = sumrsp + resppt(ii)
c10      continue

c ..... Use 2.0 to convert from biomass to carbon in forest system
c ..... Mike Ryan & Dennis Ojima
c ..... Added calculation of fcmax  -rm  11/91
c ..... prdx(3) is maximum NET production in g C/m^2/mo

c ..... Since we are no longer using the maximum gross production equation the
c ..... pforc calculation now now uses only the maximum net production value,
c ..... cak - 07/01/02
c        pforc = (tfprod / 2.0) - sumrsp
c ..... Add call to new subroutine for calculating photo period effect
c ..... on growth, in the fall, when the day length is decreasing,
c ..... growth will slow down.  The definition of PRDX(2) has been
c ..... changed, it now represents the coefficient for calculating the
c ..... potential monthly forest production as a function of solar
c ..... radiation outside of the atmosphere.  The PRDX(3) parameter has
c ..... been removed from the tree.100 file.  cak - 08/26/02
c        fcmax = prdx(3) * tfrac * potprd * h2ogef(2) *
c     &          laprod(lai,laitop) * co2cpr(FORSYS)
        fcmax = shwave(month,sitlat,jdaywk) * prdx(2) * tfrac *
     &          potprd * h2ogef(2) * laprod(lai,laitop) *
     &          co2cpr(FORSYS)
        pforc = fcmax

c ..... Compute carbon allocation fractions for each tree part
c ..... mdh 5/11/01
c ..... Call moved from treegrow subroutine, cak - 07/01/02
        if (pforc .gt. 0.01) then
          call treeDynC(pforc, tavemth, tree_cfrac)
        else
          pforc = 0.0
        endif

      else
        pforc = 0.
      endif

      return
      end
