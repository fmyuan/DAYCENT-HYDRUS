
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine mthacc(agdefacsum, bgdefacsum)

      implicit none
      include 'const.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

c ... Argument declarations
      real      agdefacsum, bgdefacsum

c ... Reset monthly accumulators.

c ... Local variables
      integer iel, ii, ilyr, iso

c ... Initialize monthly accumulators
      do 20 ii = 1, 8
        stream(ii) = 0
20    continue
      pet = 0
      evap = 0
      tran = 0
      pttr = 0
      rain = 0
      agdefacsum = 0.0
      bgdefacsum = 0.0
      irract = 0.0
      runoff = 0.0
      do 102 ilyr = 1,nlayer
        amov(ilyr) = 0
102   continue

c ... Initialize monthly co2 accumlators (10/92)
      do 25 iso = 1, 2
        st1c2(iso) = 0.0
        st2c2(iso) = 0.0
        mt1c2(iso) = 0.0
        mt2c2(iso) = 0.0
        s11c2(iso) = 0.0
        s21c2(iso) = 0.0
        s2c2(iso)  = 0.0
        s3c2(iso)  = 0.0
        wd1c2(iso) = 0.0
        wd2c2(iso) = 0.0
        wd3c2(iso) = 0.0
25    continue

c ... Initialize monthly accumulator for volatilization of N during
c ... harvest, senescence, and return from grazing animal waste,
c ... cak 01/02
      volpl = 0.0

c ... Initialize monthly accumulator for symbiotic N fixation to track
c ... fixation for both grasses and trees as necessary, cak - 10/15/02
      nfix = 0.0

c ... Initialize monthly C production, cak - 11/20/03
      cprodc = 0.0
      cprodf = 0.0
      do 30 iel = 1, MAXIEL
        eprodc(iel) = 0.0
        eprodf(iel) = 0.0
30    continue

c ... Initialize monthly accumulator for soil surface temperature,
c ... cak - 11/20/03
      stempmth = 0.0

c ... Initialize monthly accumulator for maintenance respiration
c ... cak - 05/14/04
      mrspflow(CRPSYS) = 0.0
      mrspflow(FORSYS) = 0.0
      sumrsp = 0.0

c ... Initialize monthly trace gas accumulator output variables,
c ... cak - 05/14/04
      N2O_month = 0.0
      NO_month = 0.0
      N2_month = 0.0
      ch4_month = 0.0
      nit_amt_month = 0.0
      pptmonth = 0.0

      return
      end
