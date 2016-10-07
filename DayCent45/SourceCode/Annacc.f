
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine annacc

      implicit none
      include 'const.inc'
      include 'monprd.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'param.inc'

c ... Reset annual accumulators.
c ... NOTE: The annet annual accumulator is reset in eachyr as it is
c ...       before used in the calculation for non-symbiotic soil N
c ...       fixation being reset

c ... Local variables
      integer iel, ii

c ... Initialize annual removal accumulators
      prcann = 0.0
      petann = 0.0
      nfixac = 0.0
      cgrain = 0.0
      cgracc = 0.0
      snfxac(CRPSYS) = 0.0
      snfxac(FORSYS) = 0.0
      crmvst = 0.0
      accrst = 0.0
      shrema = 0.0
      shrmai(UNLABL) = 0.0
      shrmai(LABELD) = 0.0
      sdrema = 0.0
      sdrmai(UNLABL) = 0.0
      sdrmai(LABELD) = 0.0
      creta = 0.0
      resp(UNLABL) = 0.0
      resp(LABELD) = 0.0
      do 10 iel = 1, nelem
        ereta(iel) = 0.0
        shrmae(iel) = 0.0
        sdrmae(iel) = 0.0
        egracc(iel) = 0.0
c ..... Initialize mineralization accumulators
        tnetmn(iel) = 0.0
        sumnrs(iel) = 0.0
        soilnm(iel) = 0.0
10    continue

c ... Initialize annual C production
      cproda = 0.0

c ... Initialize cinputs
      cinput = 0.0

c ... Reset minimum total non-living C, an annual value
      totc = 1000000
      
c ... Initialize co2 accumulators (10/92)
      ast1c2 = 0.0
      ast2c2 = 0.0
      amt1c2 = 0.0
      amt2c2 = 0.0
      as11c2 = 0.0
      as21c2 = 0.0
      as2c2 = 0.0
      as3c2 = 0.0

c ... Initialize the annual maintenance respiration accumulators
      mrspann(CRPSYS) = 0.0
      mrspann(FORSYS) = 0.0

c ... Initialize stream accumulators (04/03)
      do 20 ii = 1, 8
        strmac(ii) = 0.0
20    continue

c ... Initialize fertilizer addition accumulator (04/12/04)
      do 30 ii = 1, 3
        fertac(ii) = 0.0
30    continue

c ... Initialize annual accumulators for N volatilization (04/14/04)
      voleac = 0.0
      volgac = 0.0
      volpac = 0.0

c ... Initialize E return from grazing accumulators (04/14/04)
      do 40 ii = 1, 3
        tgzrte(ii) = 0.0
40    continue

c ... Initialize accumulators for N deposition and non-symbiotic soil N
c ... fixation
      wdfxas = 0.0
      wdfxaa = 0.0
      wdfxa = 0.0

c ... Reset accumulators for yearly trace gas output, cak - 09/23/02
      N2O_year = 0.0
      NO_year = 0.0
      N2_year = 0.0
      CH4_year = 0.0
      nit_amt_year = 0.0
      annppt = 0.0

c ... Initialize OMAD addition accumulators, cak - 07/13/2006
      omadac = 0.0
      do 50 ii = 1, 3
        omadae(ii) = 0.0
50    continue

      return
      end
