
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine harvst(month,pltlig)

      implicit none
      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'

c ... Argument declarations
      integer   month
      real      pltlig(2)

c ... Harvest the crop

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE flow(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow' :: flow
          REAL from
          REAL to
          REAL when
          REAL howmuch
        END SUBROUTINE flow

        SUBROUTINE flowup(time)
          !MS$ATTRIBUTES ALIAS:'_flowup' :: flowup
          REAL time
        END SUBROUTINE flowup

        SUBROUTINE flowup_double(time)
          !MS$ATTRIBUTES ALIAS:'_flowup_double' :: flowup_double
          REAL time
        END SUBROUTINE flowup_double

        SUBROUTINE flowup_double_in(time)
          !MS$ATTRIBUTES ALIAS:'_flowup_double_in' :: flowup_double_in
          REAL time
        END SUBROUTINE flowup_double_in

        SUBROUTINE flowup_double_out(time)
          !MS$ATTRIBUTES ALIAS:'_flowup_double_out' :: flowup_double_out
          REAL time
        END SUBROUTINE flowup_double_out

      END INTERFACE

c ... Local variables
      integer   iel, mm, hmonth(2)
      real      accum(ISOS), addsdc, addsde, bgd, cisbgd(ISOS), 
     &          cstraw, ctubes, etubes, fr14, recres(MAXIEL), 
     &          resid, sumpttr, sumtran, harv_volpl,
     &          mrspstgFracRemoved, mrspstgLoss, mRespStorage

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

c ... Check that there is material to harvest
      if (aglivc .le. 0.001) then
        cgrain = 0.0
        crmvst = 0.0
        do 5 iel = 1, nelem
          egrain(iel) = 0.0
          ermvst(iel) = 0.0
5       continue
        goto 999
      endif

c ... Carbon

c ... Grain
c ... (Alister's new way of calculating:)
      if (flghrv .eq. 1) then
        sumtran = 0
        sumpttr = 0
        hmonth(1) = month - himon(1)
        hmonth(2) = month - himon(2)
        if (hmonth(1) .lt. 1) then
          hmonth(1) = hmonth(1) + MONTHS
        endif
        if (hmonth(2) .lt. 1)  then
          hmonth(2) = hmonth(2) + MONTHS
        endif
        if (hmonth(2) .ge. hmonth(1)) then
          do 10 mm = hmonth(1), hmonth(2)
            sumtran = sumtran + htran(mm)
            sumpttr = sumpttr + hpttr(mm)
10        continue
        else
          do 15 mm = hmonth(1), MONTHS
            sumtran = sumtran + htran(mm)
            sumpttr = sumpttr + hpttr(mm)
15        continue
          do 16 mm = 1, hmonth(2)
            sumtran = sumtran + htran(mm)
            sumpttr = sumpttr + hpttr(mm)
16        continue
        endif
        hi = himax * (1.0 - hiwsf * (1.0 - (sumtran / sumpttr)))
      else
        hi = 0.0
      endif
      cgrain = hi * aglivc * (1.0 - aglrem)
      call csched(cgrain,aglcis(LABELD),aglivc,
     &            aglcis(UNLABL),csrsnk(UNLABL),
     &            aglcis(LABELD),csrsnk(LABELD),
     &            1.0,cisgra)

c ... Straw
      cstraw = aglivc * (1.0 - aglrem) - cgrain

c ... Straw removal
      crmvst = rmvstr * cstraw
      accrst = accrst + crmvst

      call csched(crmvst,aglcis(LABELD),aglivc,
     &            aglcis(UNLABL),csrsnk(UNLABL),
     &            aglcis(LABELD),csrsnk(LABELD),
     &            1.0,accum)

c ... Some straw will remain as standing dead
      addsdc = remwsd * (cstraw-crmvst)
      call csched(addsdc,aglcis(LABELD),aglivc,
     &            aglcis(UNLABL),stdcis(UNLABL),
     &            aglcis(LABELD),stdcis(LABELD),
     &            1.0,accum)

c ... Other elements
      do 20 iel = 1, nelem

c ..... Grain
        if (flghrv .eq. 1) then
          egrain(iel) = efrgrn(iel) * aglive(iel) * (1.0 - aglrem) *
     &                  sqrt(hi/himax)
          call flow(aglive(iel),esrsnk(iel),time,egrain(iel))
        else
          egrain(iel) = 0.0
        endif

c ..... Volatilization of N from plants
c ..... Use the local variable harv_volpl so that volatilization that
c ..... occurs at harvest and senescence can both be tracked, see dshoot
c ..... cak - 01/02
        if (iel .eq. N) then
          harv_volpl = vlossp * aglive(iel)
          call flow(aglive(iel),esrsnk(iel),time,harv_volpl)
          volpl = volpl + harv_volpl
          volpla = volpla + harv_volpl

c ....... N/C ratio in straw
          recres(iel) = ((aglive(iel) - harv_volpl) * (1.0 - aglrem) - 
     &                    egrain(iel)) / cstraw
        else
c ....... P/C, or S/C ratio in straw
          recres(iel) = (aglive(iel) * (1.0 - aglrem) - 
     &                   egrain(iel)) / cstraw
        endif

c ..... Straw removal
        ermvst(iel) = crmvst * recres(iel)
        call flow(aglive(iel),esrsnk(iel),time,ermvst(iel))

c ..... Some straw remains as standing dead
        addsde = addsdc * recres(iel)
        call flow(aglive(iel),stdede(iel),time,addsde)
20    continue

c ... Partition c, n, p, and s in remaining straw into top layer
c ... of structural and metabolic
      resid = cstraw - crmvst - addsdc
      fr14 = aglcis(LABELD)/aglivc
      call partit(resid,recres,1,aglcis,aglive,pltlig(ABOVE),fr14)

c ... Below ground removal (root harvest) -lh 8/91

      ctubes = hibg * bglivc * (1.0 - bglrem)
      cgrain = cgrain + ctubes
c ... A fraction of the maintenance respiration storage pool is
c ... removed in proportion to the amount of live roots removed,
c ... compute before bglcis is updated, mdh - 09/01
      mrspstgFracRemoved = ctubes / bglivc
      call csched(ctubes,bglcis(LABELD),bglivc,
     &            bglcis(UNLABL), csrsnk(UNLABL),
     &            bglcis(LABELD), csrsnk(LABELD),
     &            1.0,cisgra)
c ... When live roots die, a proportional amount of the maintenance
c ... respiration storage carbon pool flows to the same destination
c ... as the roots, mdh - 09/01
      mRespStorage = mrspstg(CRPSYS,UNLABL) + mrspstg(CRPSYS,LABELD)
      mrspstgLoss = mrspstgFracRemoved * mRespStorage
      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0
      if (mrspstg(CRPSYS,UNLABL) .lt. 0.0) then
        write(*,*) 'Error in droot, mrspstg(CRPSYS,UNLABL) < 0.0'
        STOP
      endif
      if (mrspstg(CRPSYS,LABELD) .lt. 0.0) then
        write(*,*) 'Error in droot, mrspstg(CRPSYS,UNLABL) < 0.0'
        STOP
      endif
      call csched(mrspstgLoss, mrspstg(CRPSYS,LABELD), mRespStorage,
     &            mrspstg(CRPSYS,UNLABL), csrsnk(UNLABL),
     &            mrspstg(CRPSYS,LABELD), csrsnk(LABELD),
     &            1.0, accum)
c ... Need this?
c      cinput = cinput + accum(LABELD) + accum(UNLABL)
      if (bglivc .gt. 0.0001) then
        cisbgd(LABELD) = (bglivc*(1.0 - bglrem)-ctubes) *
     &                   (bglcis(LABELD)/bglivc)
      else
        cisbgd(LABELD) = 0.0
      endif
      cisbgd(UNLABL) = (bglivc * (1.0 - bglrem) - ctubes) -
     &                 cisbgd(LABELD)

      do 30 iel = 1, nelem
        etubes = hibg * (1.0 - bglrem) * bglive(iel)
        call flow(bglive(iel), esrsnk(iel), time, etubes)
        egrain(iel) = egrain(iel) + etubes
        recres(iel) = bglive(iel) / bglivc

c ..... Remove from crop storage as well
        etubes = hibg * (1.0 - bglrem) * crpstg(iel)
        call flow(bglive(iel), esrsnk(iel), time, etubes)

30    continue

c ... Calculation of accumulator for grain production
      cgracc = cgracc + cgrain
      do 40 iel = 1, nelem
        egracc(iel) = egracc(iel)+egrain(iel)
40    continue

c ... Partition c, n, p, and s in remaining roots into bottom layer of
c ... structural and metabolic
      bgd = cisbgd(LABELD) + cisbgd(UNLABL)
      call partit(bgd,recres,2,bglcis,bglive,pltlig(BELOW),fr14)

c ... Update state variables and accumulators.
      call flowup(time)
      call flowup_double(time)
      call flowup_double_in(time)
      call flowup_double_out(time)
      call sumcar

c ... Check status of values to make sure that everything
c ... has been reset correctly

      if (aglcis(UNLABL)+aglcis(LABELD) .lt. 1.e-05) then
        aglcis(UNLABL) = 0.0
        aglcis(LABELD) = 0.0
      endif
      if (aglcis(UNLABL)+aglcis(LABELD) .eq. 0) then
        do 50 iel = 1, MAXIEL
          aglive(iel) = 0.0
50      continue
      endif
      
      if (bglcis(UNLABL)+bglcis(LABELD) .lt. 1.e-05) then
        bglcis(UNLABL) = 0.0
        bglcis(LABELD) = 0.0
      endif  
      if (bglcis(UNLABL)+bglcis(LABELD) .eq. 0) then
        do 60 iel = 1, MAXIEL
          bglive(iel) = 0.0
60      continue
      endif

      if (stdcis(UNLABL)+stdcis(LABELD) .lt. 1.e-05)then
        stdcis(UNLABL) = 0.0
        stdcis(LABELD) = 0.0
      endif  
      if (stdcis(UNLABL)+stdcis(LABELD) .eq. 0) then
        do 70 iel = 1, MAXIEL
          stdede(iel) = 0
70      continue
      endif
          
      do 80 iel = 1, MAXIEL
        if (aglive(iel) .lt. 1.e-05) then
          aglive(iel) = 0.0
        endif
        if (bglive(iel) .lt. 1.e-05) then
          bglive(iel) = 0.0
        endif
        if (stdede(iel) .lt. 1.e-05) then
          stdede(iel) = 0.0
        endif
80    continue

999   continue

      return
      end
