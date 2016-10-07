
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... WDEATH

      subroutine wdeath (tavewk, bgwfunc, tfrac)

      implicit none
      include 'const.inc'
      include 'forrem.inc'
      include 'isovar.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'pheno.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'site.inc'
      include 'timvar.inc'
      include 'zztim.inc'

c ... Argument declarations
      real tavewk, bgwfunc, tfrac

c ... Death of leaves, fine branches, large wood, fine roots, and coarse roots.

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
      integer iel, drpwks
      logical drpdlv
      real    accum(ISOS), ctodie, etodie, fr14, recres(MAXIEL), tostore
      data    drpwks /0/

c ... Saved variables
      save drpdlv
      save drpwks

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

      if ((time - strtyr) .le. 0.00001) drpdlv = .FALSE.

c ... Death of leaves
c ... NOTE:  WOODDR(1)   - the death rate in fall for deciduous forests
c ...        LEAFDR(MTH) - the monthly death rate for leaves in every
c ...                      case except for fall in deciduous forests.
      if (rleavc .gt. 0.0001) then
        if (decid .ge. 1) then

c ....... Deciduous forest
c ....... If the daylight hours are increasing - it must be spring
          if ((hrsinc) .and. (tavewk .gt. tmplff)) drpdlv = .FALSE.

c ....... If daylight hours are decreasing and the temperature is low
c ....... enough drop leaves for fall season.
c ....... Add check for number of daylight hours to conditional for
c ....... determining if leaf drop should occur, cak - 06/30/03
c ....... If leaf drop has not occurred by the time the winter solstice
c ....... is reached force leaf drop to occur, cak - 10/28/04
          if (decid .eq. 1) then
            if (((tavewk .lt. tmplff) .and. (.not. drpdlv) .and.
     &           (.not. hrsinc) .and. (dayhrs .lt. 12.0)) .or.
     &          ((.not.drpdlv) .and. (sitlat .ge. 0) .and.
     &           (month .eq. 12)) .or.
     &          ((.not.drpdlv) .and. (sitlat .lt. 0) .and.
     &           (month .eq. 6))) then
c ........... Allow leaf drop to occur over a 4 week period, dropping
c ........... all of the remaining leaves in the last week.
              if (drpwks .lt. 4) then
                ctodie = rleavc * wooddr(LEAF)*tfrac
                drpwks = drpwks + 1
                decidgrow = .FALSE.
               else
                ctodie = rleavc * wooddr(LEAF)
                drpdlv = .TRUE.
                drpwks = 0
              endif
            else
              ctodie = rleavc * leafdr(month)*tfrac
            endif
          elseif (decid .eq. 2) then
c ......... Drought deciduous forest
c ......... Compute death for drought deciduous forests
            ctodie = rleavc * (1. - bgwfunc) * wooddr(LEAF)*tfrac
          endif
        else
c ....... Continuous forest
c ....... Use leaf death rate multiplier from EACHYR
          ctodie = rleavc * leafdr(month)*tfrac * ldrmlt
        endif

c ..... Compute E/C ratios
        do 10 iel = 1, nelem
          recres(iel) = rleave(iel) / rleavc

c ....... Compute flow to retranslocation storage
          tostore = recres(iel) * ctodie * forrtf(iel)
          call flow(rleave(iel), forstg(iel), time, tostore)

c ....... Decrease E/C by the amount that is retranslocated
          recres(iel) = recres(iel) * (1 - forrtf(iel))
10      continue

c ..... If evntyp is greater than 1 the leaves go to the source/sink
c ..... rather than to the litter, cak - 02/07/2006
        if (evntyp .lt. 2) then
          fr14 = rlvcis(LABELD) / rleavc
          call partit(ctodie, recres, 1, rlvcis, rleave, wdlig(LEAF),
     &                fr14)
        else
          call csched(ctodie, cisotf, 1.0,
     &                rlvcis(UNLABL), csrsnk(UNLABL),
     &                rlvcis(LABELD), csrsnk(LABELD),
     &                1.0, accum)
          do 15 iel = 1, nelem
            etodie = ctodie * (rleave(iel) / rleavc)
            if (iel .eq. N) then
              call flow(rleave(iel), esrsnk(iel), time,
     &                  etodie - tostore)
            else
              call flow(rleave(iel), esrsnk(iel), time, etodie)
            endif
15        continue
        endif
      endif

c ... Death of fine roots
      if (frootc .gt. 0.0) then
        ctodie = frootc * wooddr(FROOT)*tfrac
        do 20 iel = 1, nelem
          recres(iel) = froote(iel) / frootc
20      continue
        fr14 = frtcis(LABELD) / frootc
        call partit(ctodie, recres, 2, frtcis, froote, wdlig(FROOT),
     &              fr14)
      endif

c ... Fine Branches, Large Wood, and Coarse Roots go to the dead wood
c ... compartments: WOOD1, WOOD2, WOOD3

c ... Death of fine branches
      if (fbrchc .gt. 0.0) then
        ctodie = fbrchc * wooddr(FBRCH)*tfrac
        call csched(ctodie, cisotf, 1.0,
     &              fbrcis(UNLABL), wd1cis(UNLABL),
     &              fbrcis(LABELD), wd1cis(LABELD),
     &              1.0, accum)

        do 30 iel = 1, nelem
          etodie = ctodie * (fbrche(iel) / fbrchc)
          call flow(fbrche(iel), wood1e(iel), time, etodie)
30      continue
      endif

c ... Death of large wood
      if (rlwodc .gt. 0.0) then
        ctodie = rlwodc * wooddr(LWOOD)*tfrac
        call csched(ctodie, cisotf, 1.0,
     &              rlwcis(UNLABL), wd2cis(UNLABL),
     &              rlwcis(LABELD), wd2cis(LABELD),
     &              1.0, accum)

        do 40 iel = 1, nelem
          etodie = ctodie * (rlwode(iel) / rlwodc)
          call flow(rlwode(iel), wood2e(iel), time, etodie)
40      continue
      endif

c ... Death of coarse roots
      if (crootc .gt. 0.0) then
        ctodie = crootc * wooddr(CROOT)*tfrac
        call csched(ctodie, cisotf, 1.0,
     &              crtcis(UNLABL), wd3cis(UNLABL),
     &              crtcis(LABELD), wd3cis(LABELD),
     &              1.0, accum)

        do 50 iel = 1, nelem
          etodie = ctodie * (croote(iel) / crootc)
          call flow(croote(iel), wood3e(iel), time, etodie)
50      continue
      endif

      return
      end
