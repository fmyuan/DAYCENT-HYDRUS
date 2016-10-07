
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine treegrow (tfrac,tavewk,startofmo)

      implicit none
      include 'const.inc'
      include 'dynam.inc'
      include 'fertil.inc'
      include 'isovar.inc'
      include 'monprd.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'zztim.inc'

c ... Argument declarations
      real    tfrac, tavewk
      logical startofmo

c ... Simulate forest production associated with tree growth.
c ... This function was a part of TREES.F.
c ...   tfrac  - fraction of month over which current production event
c ...            occurs (0-1)
c ...   tavewk - mean air temperature over production period (deg C)

c ... Fortran to C prototype
      INTERFACE
        SUBROUTINE cmpnfrac(lyr, ammonium, nitrate, minerl,
     &                      frac_nh4, frac_no3)
          !MS$ATTRIBUTES ALIAS:'_cmpnfrac' :: cmpnfrac
          INTEGER lyr
          REAL*8  ammonium
          REAL*8  nitrate(*)
          REAL    minerl(*)
          REAL*8  frac_nh4
          REAL*8  frac_no3
        END SUBROUTINE cmpnfrac

        SUBROUTINE flow(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow' :: flow
          REAL from
          REAL to
          REAL when
          REAL howmuch
        END SUBROUTINE flow

        SUBROUTINE update_npool(clyr, amt, frac_nh4, frac_no3, 
     &             ammonium, nitrate, subname)
          !MS$ATTRIBUTES ALIAS:'_update_npool' :: update_npool
          INTEGER   clyr
          REAL      amt
          REAL*8    frac_nh4
          REAL*8    frac_no3
          REAL*8    ammonium
          REAL*8    nitrate(*)
          CHARACTER subname*10
        END SUBROUTINE update_npool

      END INTERFACE

c ... Function declarations
      real      leafa, rtimp
      external  leafa, rtimp

c ... Local variables
      integer   iel, ipart, iptr, lyr
      real      accum(ISOS), availm(MAXIEL), amt
      real      calcup, cprodfLeft
      real      euf(FPARTS), fsol
      real      mrspTempEffect
      real      namt
      real      remCfrac, rimpct
      real      sum_cfrac, toler, totCup
!      real      uptake(4,MAXIEL)
      real      uptake(5,MAXIEL)
      real      treeNfix
      real*8    frac_nh4, frac_no3
      character subname*10

      if (tfrac .lt. 0.0 .or. tfrac .gt. 1.0) then
        write(*,*) 'Error in treegrow, tfrac = ', tfrac
        STOP
      endif
      if (tavewk .le. -999.0) then
        write(*,*) 'Error in treegrow, tavewk = ', tavewk
        STOP
      endif

      subname = 'treegrow  '
      toler = 1.0E-30

c ... Initialize monthly production, maintenance respiration, and
c ... uptake variables
      accum(UNLABL) = 0.0
      accum(LABELD) = 0.0
      mrspwkflow(FORSYS) = 0.0
      if (startofmo) then
        mrspflow(FORSYS) = 0.0
      endif
      do 10 ipart=1, FPARTS
        mfprd(ipart) = 0.0
10    continue
      do 15 iel = 1, nelem
        uptake(ESTOR,iel) = 0.0
        uptake(ESOIL,iel) = 0.0
        uptake(ENFIX,iel) = 0.0
        uptake(EFERT,iel) = 0.0
        uptake(EFOLI,iel) = 0.0
15    continue

c ... Previous flowup, in trees, should have updated mineral pools.
c ... Determine nutrients available to plants for growth.
      do 20 iel = 1, nelem
        availm(iel) = 0.0
        do 30 lyr = 1, nlayer
          availm(iel) = availm(iel) + minerl(lyr, iel)
30      continue
20    continue

c ... Determine old or new forest
c ... iptr points to new forest carbon allocation fractions (iptr = 1) or
c ... mature forest carbon allocation fractions (iptr = 2) for each of the
c ... tree parts; leaves, fine roots, fine branches, large wood, and
c ... coarse roots.  Switch from new forest allocation fractions to old
c ... forest allocation fractions at time = swold

C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
      if ((time-time0) .le. swold) then             ! time0: time when trees start to grow/regrow
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
!      if (time .le. swold) then
c ..... Use juvenile forest C allocation fractions
        iptr = 1
      else
c ..... Use mature forest C allocation fractions
        iptr = 2
      endif

c ... If growth can occur
      if (forgrw .eq. 1 .and. pforc .gt. 0.0) then

c ..... Calculate actual production values
c ..... Calculate impact of root biomass on available nutrients
        rimpct = rtimp(riint, rictrl, frootc)
c ..... Determine actual production values, restricting the C/E ratios
c ..... When calling restrp we are only looking at allocation to fine roots
c ..... and leaves, cak - 07/02/02
        call restrp(elimit, nelem, availm, ccefor, 2, tree_cfrac,
     &              pforc, rimpct, forstg, snfxmx(FORSYS), cprodf,
     &              eprodf, uptake, tree_a2drat, treeNfix, relyld)
      else
        cprodf = 0.0
      endif

c ... If growth occurs...
      if (cprodf .gt. 0.) then

c ..... Compute carbon allocation fractions for each tree part
c ..... Calculate how much of the carbon the roots use
        cprodfLeft = cprodf - (cprodf * tree_cfrac(FROOT))

c ..... Calculate how much of the carbon the leaves use, we allocate leaves
c ..... up to a optimal LAI
        tree_cfrac(LEAF) = leafa(rleavc, rlwodc, cprodfLeft, cprodf)
        remCfrac = 1.0 - tree_cfrac(FROOT) - tree_cfrac(LEAF)

c ..... If we have leftover carbon allocate it to the woody plant parts
c ..... using a weighted average
        if (remCfrac .lt. 1.0E-05) then
c ..... for FBRCH, LWOOD, and CROOT ...
          do 60 ipart = FBRCH, CROOT
            tree_cfrac(ipart) = 0.0
60        continue
        else
c ....... for FBRCH, LWOOD, and CROOT ...
          totCup = 0.0
          do 70 ipart = FBRCH, CROOT
            tree_cfrac(ipart) = fcfrac(ipart, iptr)
            totCup = totCup + tree_cfrac(ipart)
70        continue
          if (totCup .gt. 0.0) then
c ....... for FBRCH, LWOOD, and CROOT ...
            do 80 ipart = FBRCH, CROOT
              tree_cfrac(ipart) = tree_cfrac(ipart) / totCup * remCfrac
80          continue
          else
            write(*,*) 'Error in treegrow'
            write(*,*) 'fcfrac(FBRCH)+fcfrac(LWOOD)+fcfrac(CROOT) <= 0'
            write(*,*) 'should be > 0 for tree'
            STOP
          endif
        endif

c ... Error checking
        sum_cfrac = 0.0
        do 90 ipart = 1, FPARTS
          if (tree_cfrac(ipart) .lt. 0.0) then
            write(*,*) 'Error in treegrow, tree_cfrac(ipart) < 0'
            STOP
          else if (tree_cfrac(ipart) .gt. 1.0) then
            write(*,*) 'Error in treegrow, tree_cfrac(ipart) > 1'
            STOP
          else
            sum_cfrac = sum_cfrac + tree_cfrac(ipart)
          endif
90      continue
        if (abs(1.0 - sum_cfrac) .gt. 0.000001) then
          write(*,*) "Error in tree carbon allocation fractions!"
          STOP
        endif

c ..... Recalculate actual production values with updated C-allocation
c ..... fractions, restricting the C/E ratios -mdh 5/11/01
c ..... Calculate impact of root biomass on available nutrients
        rimpct = rtimp(riint, rictrl, frootc)
c ..... Determine actual production values, restricting the C/E ratios
        call restrp(elimit, nelem, availm, ccefor, FPARTS, tree_cfrac,
     &              pforc, rimpct, forstg, snfxmx(FORSYS), cprodf,
     &              eprodf, uptake, tree_a2drat, treeNfix, relyld)

c ..... Calculations for symbiotic N fixation accumulators moved from
c ..... nutrlm subroutine, cak - 10/17/02
c ..... Compute N fixation which actually occurs and add to the
c ..... N fixation accumulator.
        nfix = nfix + treeNfix
        snfxac(FORSYS) = snfxac(FORSYS) + treeNfix
c ..... Add computation for nfixac -mdh 1/16/02
        nfixac = nfixac + treeNfix

c ..... Update accumulators for N, P, and S uptake
        do 85 iel = 1, nelem
          eupacc(iel) = eupacc(iel) + eprodf(iel)
          do 75 ipart = 1, FPARTS
            eupprt(ipart,iel) = eupprt(ipart,iel) + eup(ipart,iel)
75        continue
85      continue

c ..... C/N ratio for production
        if (eprodf(N) .eq. 0.0) then
          write(*,*) 'Error in treegrowth, eprodf(N) = 0.0)'
          STOP
        endif
        tcnpro = cprodf/eprodf(N)

c ..... Calculate production for each tree part
c ..... New variable MFPRD added for gridded output - 6/96 rm
c ..... Added maintenance respiration (mrspflow) calculation. -mdh 2/99
c ..... Added mrspwkflow variable for weekly timestep, cak - 08/07/02
        do 95 ipart = 1, FPARTS
          mfprd(ipart) = tree_cfrac(ipart) * cprodf
          mrspwkflow(FORSYS) = mrspwkflow(FORSYS) +
     &                         mfprd(ipart) * kmrsp(FORSYS)
          mrspflow(FORSYS) = mrspflow(FORSYS) +
     &                       mfprd(ipart) * kmrsp(FORSYS)
95      continue

c ..... Growth of leaves; split into labeled & unlabeled parts
        call csched(mfprd(LEAF),cisotf,1.0,
     &              csrsnk(UNLABL),rlvcis(UNLABL),
     &              csrsnk(LABELD),rlvcis(LABELD),
     &              1.0,alvcis)

c ..... Growth of fine roots; split into labeled & unlabeled parts
        call csched(mfprd(FROOT),cisotf,1.0,
     &              csrsnk(UNLABL),frtcis(UNLABL),
     &              csrsnk(LABELD),frtcis(LABELD),
     &              1.0,afrcis)

c ..... Growth of fine branches; split into labeled & unlabeled parts
        call csched(mfprd(FBRCH),cisotf,1.0,
     &              csrsnk(UNLABL),fbrcis(UNLABL),
     &              csrsnk(LABELD),fbrcis(LABELD),
     &              1.0,afbcis)

c ..... Growth of large wood; split into labeled & unlabeled parts
        call csched(mfprd(LWOOD),cisotf,1.0,
     &              csrsnk(UNLABL),rlwcis(UNLABL),
     &              csrsnk(LABELD),rlwcis(LABELD),
     &              1.0,alwcis)

c ..... Growth of coarse roots; split into labeled & unlabeled parts
        call csched(mfprd(CROOT),cisotf,1.0,
     &              csrsnk(UNLABL),crtcis(UNLABL),
     &              csrsnk(LABELD),crtcis(LABELD),
     &              1.0,acrcis)

c ..... Added maintenance respiration (mrspflow) flow. -mdh 2/99
c ..... Maintenance respiration flows are added to the maintenance
c ..... respiration storage pool.
c ..... Use the weekly flow, cak - 08/07/02
        call csched(mrspwkflow(FORSYS),cisotf,1.0,
     &              csrsnk(UNLABL),mrspstg(FORSYS,UNLABL),
     &              csrsnk(LABELD),mrspstg(FORSYS,LABELD),
     &              1.0,accum)

c ..... Maintenance respiration flux calculation added, mdh - 9/4/01
        mrspTempEffect = 0.1 * exp(0.07 * tavewk)
c ..... Bound maintenance respiration temperature effect between 0.0 and 1.0,
c ..... cak - 09/16/02
        mrspTempEffect = min(1.0, mrspTempEffect)
        mrspTempEffect = max(0.0, mrspTempEffect)
        fmrspflux(LEAF) = fkmrspmx(LEAF) * mrspTempEffect * rleavc *
     &                    tfrac
        fmrspflux(FROOT) = fkmrspmx(FROOT) * mrspTempEffect * frootc *
     &                     tfrac
        fmrspflux(FBRCH) = fkmrspmx(FBRCH) * mrspTempEffect * fbrchc *
     &                     tfrac
        fmrspflux(LWOOD) = fkmrspmx(LWOOD) * mrspTempEffect * rlwodc *
     &                     tfrac
        fmrspflux(CROOT) = fkmrspmx(CROOT) * mrspTempEffect * crootc *
     &                     tfrac
        do 150 ipart = 1, FPARTS
          mrspann(FORSYS) = mrspann(FORSYS) + fmrspflux(ipart)
c ....... Monthly accumulator
          sumrsp = sumrsp + fmrspflux(ipart)
          call csched(fmrspflux(ipart),cisotf,1.0,
     &                mrspstg(FORSYS,UNLABL),csrsnk(UNLABL),
     &                mrspstg(FORSYS,LABELD),csrsnk(LABELD),
     &                1.0,accum)
150     continue

c ..... Actual Uptake
        do 110 iel = 1, nelem
          if (eprodf(iel) .eq. 0.0) then
            write(*,*) 'Divide by zero in treegrow, eprodf(iel) = 0'
            STOP
          endif
          euf(LEAF) = eup(LEAF,iel) / eprodf(iel)
          euf(FROOT) = eup(FROOT,iel) / eprodf(iel)
          euf(FBRCH) = eup(FBRCH,iel) / eprodf(iel)
          euf(LWOOD) = eup(LWOOD,iel) / eprodf(iel)
          euf(CROOT) = eup(CROOT,iel) / eprodf(iel)

c ....... Take up nutrients from internal storage pool
          amt = uptake(ESTOR,iel) * euf(LEAF)
          call flow(forstg(iel),rleave(iel),time,amt)
          amt = uptake(ESTOR,iel) * euf(FROOT)
          call flow(forstg(iel),froote(iel),time,amt)
          amt = uptake(ESTOR,iel) * euf(FBRCH)
          call flow(forstg(iel),fbrche(iel),time,amt)
          amt = uptake(ESTOR,iel) * euf(LWOOD)
          call flow(forstg(iel),rlwode(iel),time,amt)
          amt = uptake(ESTOR,iel) * euf(CROOT)
          call flow(forstg(iel),croote(iel),time,amt)

c ....... Take up nutrients from soil
          do 100 lyr = 1, nlayer
            if (minerl(lyr,iel) .gt. toler) then
              fsol = 1.0
c ........... The fsol calculation for P is not needed here to compute
c ........... the weighted average, cak - 04/05/02
c              if (iel .eq. P) then
c                fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
c              endif
              call cmpnfrac(lyr,ammonium,nitrate,minerl,
     &                      frac_nh4,frac_no3)
              calcup = uptake(ESOIL,iel)*minerl(lyr,iel)*
     &                 fsol/availm(iel)
              amt = calcup * euf(LEAF)
              if (iel .eq. N) then
                namt = -1.0*amt
                call update_npool(lyr, namt, frac_nh4, frac_no3, 
     &                            ammonium, nitrate, subname)
              endif
              call flow(minerl(lyr,iel),rleave(iel),time,amt)
              amt = calcup * euf(FROOT)
              if (iel .eq. N) then
                namt = -1.0*amt
                call update_npool(lyr, namt, frac_nh4, frac_no3, 
     &                            ammonium, nitrate, subname)
              endif
              call flow(minerl(lyr,iel),froote(iel),time,amt)
              amt = calcup * euf(FBRCH)
              if (iel .eq. N) then
                namt = -1.0*amt
                call update_npool(lyr, namt, frac_nh4, frac_no3, 
     &                            ammonium, nitrate, subname)
              endif
              call flow(minerl(lyr,iel),fbrche(iel),time,amt)
              amt = calcup * euf(LWOOD)
              if (iel .eq. N) then
                namt = -1.0*amt
                call update_npool(lyr, namt, frac_nh4, frac_no3, 
     &                            ammonium, nitrate, subname)
              endif
              call flow(minerl(lyr,iel),rlwode(iel),time,amt)
              amt = calcup * euf(CROOT)
              if (iel .eq. N) then
                namt = -1.0*amt
                call update_npool(lyr, namt, frac_nh4, frac_no3, 
     &                            ammonium, nitrate, subname)
              endif
              call flow(minerl(lyr,iel),croote(iel),time,amt)
            endif
100       continue

c ....... Take up nutrients from nitrogen fixation
          if (iel .eq. N .and. treeNfix .gt. 0) then
            amt = uptake(ENFIX,iel) * euf(LEAF)
            call flow(esrsnk(iel),rleave(iel),time,amt)
            amt = uptake(ENFIX,iel) * euf(FROOT)
            call flow(esrsnk(iel),froote(iel),time,amt)
            amt = uptake(ENFIX,iel) * euf(FBRCH)
            call flow(esrsnk(iel),fbrche(iel),time,amt)
            amt = uptake(ENFIX,iel) * euf(LWOOD)
            call flow(esrsnk(iel),rlwode(iel),time,amt)
            amt = uptake(ENFIX,iel) * euf(CROOT)
            call flow(esrsnk(iel),croote(iel),time,amt)
          endif
 
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
          if (uptake(EFOLI,iel) .gt. 0) then
!            amt = uptake(EFOLI,iel)
            amt = uptake(EFOLI,iel) * euf(LEAF)
            call flow(esrsnk(iel),rleave(iel),time,amt)
            amt = uptake(EFOLI,iel) * euf(FROOT)
            call flow(esrsnk(iel),froote(iel),time,amt)
            amt = uptake(EFOLI,iel) * euf(FBRCH)
            call flow(esrsnk(iel),fbrche(iel),time,amt)
            amt = uptake(EFOLI,iel) * euf(LWOOD)
            call flow(esrsnk(iel),rlwode(iel),time,amt)
            amt = uptake(EFOLI,iel) * euf(CROOT)
            call flow(esrsnk(iel),croote(iel),time,amt)
          endif
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------

c ....... Take up nutrients from fertilizer
          if (aufert .ne. 0) then
            if (uptake(EFERT,iel) .gt. 0.) then

c ........... Automatic fertilizer added to plant pools
              amt = uptake(EFERT,iel) * euf(LEAF)
              call flow(esrsnk(iel),rleave(iel),time,amt)
              amt = uptake(EFERT,iel) * euf(FROOT)
              call flow(esrsnk(iel),froote(iel),time,amt)
              amt = uptake(EFERT,iel) * euf(FBRCH)
              call flow(esrsnk(iel),fbrche(iel),time,amt)
              amt = uptake(EFERT,iel) * euf(LWOOD)
              call flow(esrsnk(iel),rlwode(iel),time,amt)
              amt = uptake(EFERT,iel) * euf(CROOT)
              call flow(esrsnk(iel),croote(iel),time,amt)

c ........... Automatic fertilizer added to mineral pool
              frac_nh4 = 0.5
              frac_no3 = 0.5
              if (favail(iel) .eq. 0.0) then
                write(*,*) 'Error in treegrow, favail(iel) = 0'
                STOP
              endif
              amt = uptake(EFERT,iel) * (1.0/favail(iel) - 1.0)
              fertot(iel) = fertot(iel) + uptake(EFERT,iel) + amt
              if (iel .eq. N) then
                call update_npool(lyr, amt, frac_nh4, frac_no3, 
     &                            ammonium, nitrate, subname)
              endif
              call flow(esrsnk(iel),minerl(SRFC,iel),time,amt)
            endif
          endif
110     continue
 
c ... Else there is no production this month
      else  
        cprodf = 0.0
        do 140 iel = 1, MAXIEL
          eprodf(iel) = 0.0
          do 130 ipart = 1, FPARTS
            eup(ipart,iel) = 0.0
130       continue
140     continue
      endif

      return
      end
