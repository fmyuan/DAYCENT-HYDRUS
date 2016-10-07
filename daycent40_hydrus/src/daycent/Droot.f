
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine droot(pltlig,tfrac)

      implicit none
      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'

c ... Argument declarations
      real      pltlig(2)
      real      tfrac

c ... Simulate death of roots for the month.

c ... Local variables
      integer   iel
      real      accum(ISOS), fr14, recres(MAXIEL), rdeath, rtdh,
     &          mrspstgFracRemoved, mrspstgLoss, mRespStorage

c ... Death of roots

c ... Mod. added (stemp .lt. 2.0) conditional since roots don't really
c ... die in winter because it's too cold for physiological activity. 
c ... Also added rtdh term for drought conditions.
c ... -rm  9-12-90
      if ((bglivc .le. 0.) .or. (stemp .lt. rtdtmp)) then
        goto 999
      endif

c ... This is representing the death of fine surface roots.  They depend
c ... on moisture in the top soil layers, so use avh2o(1).
      rtdh = 1.0 - avh2o(1)/(deck5+avh2o(1))
      rdeath = rdr * tfrac * rtdh
      if (rdeath .gt. 0.95) then
        rdeath = 0.95
      endif
      rdeath = rdeath * bglivc
      do 10 iel = 1, nelem
        recres(iel) = bglive(iel)/bglivc
10    continue
      fr14 = bglcis(LABELD)/bglivc
c ... A fraction of the maintenance respiration storage pool is
c ... removed in proportion to the amount of live roots removed,
c ... compute before bglcis is updated, mdh - 09/01
      mrspstgFracRemoved = rdeath / bglivc
      call partit(rdeath,recres,2,bglcis,bglive,pltlig(BELOW),fr14)
c ... When live roots die, a proportional amount of the maintenance
c ... respiration storage carbon pool flows to the same destination
c ... as the roots, mdh - 09/01
      mRespStorage = mrspstg(CRPSYS,UNLABL) + mrspstg(CRPSYS,LABELD)
      mrspstgLoss = mrspstgFracRemoved * mRespStorage
      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0
      if (mrspstg(CRPSYS,UNLABL) .lt. 0.0) then
        write(*,*) 'Error in droot, mrspstg(CRPSYS,UNLABL) < 0.0'
        write(*,*) 'mrspstg(CRPSYS,UNLABL) = ',mrspstg(CRPSYS,UNLABL)
        STOP
      endif
      if (mrspstg(CRPSYS,LABELD) .lt. 0.0) then
        write(*,*) 'Error in droot, mrspstg(CRPSYS,UNLABL) < 0.0'
        write(*,*) 'mrspstg(CRPSYS,LABELD) = ',mrspstg(CRPSYS,LABELD)
        STOP
      endif
      call csched(mrspstgLoss, mrspstg(CRPSYS,LABELD), mRespStorage,
     &            mrspstg(CRPSYS,UNLABL), metcis(SOIL,UNLABL),
     &            mrspstg(CRPSYS,LABELD), metcis(SOIL,LABELD),
     &            1.0, accum)
c ... Need this?
c      cinput = cinput + accum(LABELD) + accum(UNLABL)

999   continue

      return
      end
