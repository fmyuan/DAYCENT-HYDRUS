
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

      subroutine annacc

      implicit none
      include 'const.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'param.inc'

c ... Reset annual accumulators.

c ... Local variables
      integer iel

c ... Initialize annual removal accumulators
      prcann = 0.0
      petann = 0.0
      nfixac = 0.0
      cgrain = 0.0
      snfxac(CRPSYS) = 0.0
      snfxac(FORSYS) = 0.0
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

c ..... Initialize mineralization accumulators
        tnetmn(iel) = 0.0
        sumnrs(iel) = 0.0
        soilnm(iel) = 0.0
10    continue

c ... Initialize annual C production
      cproda = 0.0
      cprodc = 0.0
      cprodf = 0.0

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

      return
      end
